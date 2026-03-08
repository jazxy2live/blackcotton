#!/usr/bin/env python3
"""
robustness_analyzer.py — Uncertainty-aware candidate ranking
=============================================================

Moves the project from single-point predictions to robust design selection by:
  1) Loading top deterministic candidates from tradeoff optimization
  2) Running Monte Carlo perturbations over biological/model uncertainty
  3) Estimating success probability and failure risks
  4) Ranking candidates by robust score (not just nominal best case)

Usage:
    python -m src.robustness_analyzer
"""

import copy
import json
from pathlib import Path
from typing import Any

import numpy as np

from src.expression_model import active_promoter_keys, promoter_config_for_gene
from src.failure_risk_model import resolved_failure_risks
from src.tradeoff_optimizer import evaluate_candidate_proxy, load_params

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"

DEFAULT_THRESHOLDS = {
    "max_color_L": 25.0,
    "min_strength_g_tex": 28.0,
    "min_yield_index": 0.85,
    "min_temporal_gap_days": 0.0,
    "max_toxicity_pre_cellulose": 0.09,
}

DEFAULT_NOISE = {
    "timing_sigma_dpa": 0.75,     # absolute DPA perturbation
    "strength_cv": 0.08,          # multiplicative coefficient of variation
    "competition_cv": 0.12,
    "efficiency_cv": 0.10,
    "retention_sigma": 0.05,      # additive perturbation (bounded to [0, 0.9])
    "hill_cv": 0.08,
    "leakage_cv": 0.12,
    "copper_sigma": 0.04,
    "activation_sigma": 0.05,
    "ros_cv": 0.10,
    "silencing_sigma": 0.03,
    "event_cv_sigma": 0.04,
}

DEFAULT_CORRELATED_PROFILE = {
    "name": "construct_bundle_v1",
    "shared_expression_drop": {
        "probability": 0.22,
        "sigma": 0.75,
    },
    "timing_overlap": {
        "probability": 0.18,
        "sigma": 0.70,
    },
    "storage_slip": {
        "probability": 0.14,
        "sigma": 0.85,
    },
}


def _bounded(value: float, lo: float, hi: float) -> float:
    return float(np.clip(value, lo, hi))


def _lognormal_factor(rng: np.random.Generator, cv: float) -> float:
    """Sample multiplicative uncertainty with mean ~1.0."""
    if cv <= 0:
        return 1.0
    sigma2 = np.log(1.0 + cv * cv)
    sigma = float(np.sqrt(sigma2))
    mu = -0.5 * sigma2
    return float(np.exp(rng.normal(mu, sigma)))


def correlated_profile_config(profile: str | dict[str, Any] | None) -> dict[str, Any] | None:
    if profile is None or profile == "" or profile == "none":
        return None
    if isinstance(profile, dict):
        return copy.deepcopy(profile)
    key = str(profile)
    profiles = {
        DEFAULT_CORRELATED_PROFILE["name"]: DEFAULT_CORRELATED_PROFILE,
        "construct_bundle": DEFAULT_CORRELATED_PROFILE,
        "anti_silencing_stress": DEFAULT_CORRELATED_PROFILE,
    }
    if key not in profiles:
        raise ValueError(f"Unknown correlated profile: {profile}")
    return copy.deepcopy(profiles[key])


def _shock_severity(
    rng: np.random.Generator,
    probability: float,
    sigma: float,
) -> float:
    if probability <= 0.0 or sigma <= 0.0:
        return 0.0
    if float(rng.random()) > float(probability):
        return 0.0
    return float(np.clip(abs(rng.normal(0.0, sigma)), 0.0, 3.0))


def sample_candidate_parameters(
    candidate: dict[str, Any],
    rng: np.random.Generator,
    noise: dict[str, float],
) -> dict[str, float]:
    """Perturb candidate parameters to model biological/process variability."""
    mat_activation = _bounded(
        float(candidate["mat_activation_dpa"]) + rng.normal(0.0, noise["timing_sigma_dpa"]),
        28.0,
        45.0,
    )
    scw_activation = _bounded(
        float(candidate["scw_activation_dpa"]) + rng.normal(0.0, noise["timing_sigma_dpa"]),
        22.0,
        40.0,
    )

    return {
        "mat_activation_dpa": mat_activation,
        "scw_activation_dpa": scw_activation,
        "mat_strength": _bounded(
            float(candidate["mat_strength"]) * _lognormal_factor(rng, noise["strength_cv"]),
            0.4,
            2.0,
        ),
        "scw_strength": _bounded(
            float(candidate["scw_strength"]) * _lognormal_factor(rng, noise["strength_cv"]),
            0.4,
            2.0,
        ),
        "k_competition": _bounded(
            float(candidate["k_competition"]) * _lognormal_factor(rng, noise["competition_cv"]),
            0.01,
            1.0,
        ),
        "melanin_efficiency": _bounded(
            float(candidate["melanin_efficiency"]) * _lognormal_factor(rng, noise["efficiency_cv"]),
            0.4,
            2.0,
        ),
        "late_retention_factor": _bounded(
            float(candidate["late_retention_factor"]) + rng.normal(0.0, noise["retention_sigma"]),
            0.0,
            0.9,
        ),
    }


def jitter_promoters(
    params: dict[str, Any],
    rng: np.random.Generator,
    noise: dict[str, float],
) -> dict[str, Any]:
    """Perturb promoter shape/leakage to test timing robustness."""
    perturbed = copy.deepcopy(params)
    for key in active_promoter_keys(perturbed):
        promoter = perturbed["promoters"][key]
        promoter["hill_coefficient"] = max(
            1.0,
            float(promoter["hill_coefficient"]) * _lognormal_factor(rng, noise["hill_cv"]),
        )
        promoter["leakage_fraction"] = _bounded(
            float(promoter["leakage_fraction"]) * _lognormal_factor(rng, noise["leakage_cv"]),
            0.0,
            0.25,
        )
    return perturbed


def _failure_risk_block(params: dict[str, Any]) -> dict[str, float]:
    return resolved_failure_risks(params)


def jitter_failure_risks(
    params: dict[str, Any],
    rng: np.random.Generator,
    noise: dict[str, float],
) -> dict[str, Any]:
    perturbed = copy.deepcopy(params)
    base = _failure_risk_block(params)
    fr = perturbed.setdefault("failure_risks", {})

    fr["copper_loading_fraction"] = _bounded(
        base["copper_loading_fraction"] + rng.normal(0.0, float(noise.get("copper_sigma", 0.0))),
        0.0,
        1.0,
    )
    fr["tyrosinase_activation_fraction"] = _bounded(
        base["tyrosinase_activation_fraction"] + rng.normal(0.0, float(noise.get("activation_sigma", 0.0))),
        0.0,
        1.0,
    )
    fr["ros_buffer_capacity"] = _bounded(
        base["ros_buffer_capacity"] * _lognormal_factor(rng, float(noise.get("ros_cv", 0.0))),
        0.1,
        1.5,
    )
    fr["silencing_probability"] = _bounded(
        base["silencing_probability"] + rng.normal(0.0, float(noise.get("silencing_sigma", 0.0))),
        0.0,
        0.80,
    )
    fr["event_expression_cv"] = _bounded(
        base["event_expression_cv"] + rng.normal(0.0, float(noise.get("event_cv_sigma", 0.0))),
        0.0,
        0.80,
    )
    fr["quinone_stress_threshold"] = max(base["quinone_stress_threshold"], 1e-6)
    return perturbed


def apply_correlated_failure_profile(
    sampled: dict[str, float],
    perturbed_params: dict[str, Any],
    rng: np.random.Generator,
    correlated_profile: str | dict[str, Any] | None,
) -> tuple[dict[str, float], dict[str, Any], dict[str, float]]:
    profile = correlated_profile_config(correlated_profile)
    if not profile:
        return sampled, perturbed_params, {}

    out_sampled = dict(sampled)
    out_params = copy.deepcopy(perturbed_params)
    fr = out_params.setdefault("failure_risks", {})
    comp = out_params.setdefault("melanin_pathway", {}).setdefault("compartmentalization", {})
    prom_mel = out_params["promoters"]["pGhMat1"]
    scw_promoters = [
        out_params["promoters"][key]
        for key in active_promoter_keys(out_params)
        if key != "pGhMat1"
    ]
    meta = {"profile_name": str(profile["name"])}

    shared_sev = _shock_severity(
        rng,
        float(profile.get("shared_expression_drop", {}).get("probability", 0.0)),
        float(profile.get("shared_expression_drop", {}).get("sigma", 0.0)),
    )
    meta["shared_expression_drop"] = float(shared_sev)
    if shared_sev > 0.0:
        out_sampled["mat_strength"] = _bounded(
            out_sampled["mat_strength"] * max(0.50, 1.0 - 0.12 * shared_sev),
            0.4,
            2.0,
        )
        out_sampled["scw_strength"] = _bounded(
            out_sampled["scw_strength"] * max(0.48, 1.0 - 0.16 * shared_sev),
            0.4,
            2.0,
        )
        out_sampled["melanin_efficiency"] = _bounded(
            out_sampled["melanin_efficiency"] * max(0.45, 1.0 - 0.14 * shared_sev),
            0.4,
            2.0,
        )
        out_sampled["late_retention_factor"] = _bounded(
            out_sampled["late_retention_factor"] - 0.05 * shared_sev,
            0.0,
            0.9,
        )
        fr["copper_loading_fraction"] = _bounded(
            float(fr.get("copper_loading_fraction", 1.0)) - 0.05 * shared_sev,
            0.0,
            1.0,
        )
        fr["tyrosinase_activation_fraction"] = _bounded(
            float(fr.get("tyrosinase_activation_fraction", 1.0)) - 0.05 * shared_sev,
            0.0,
            1.0,
        )
        fr["silencing_probability"] = _bounded(
            float(fr.get("silencing_probability", 0.0)) + 0.08 * shared_sev,
            0.0,
            0.80,
        )
        fr["event_expression_cv"] = _bounded(
            float(fr.get("event_expression_cv", 0.0)) + 0.09 * shared_sev,
            0.0,
            0.80,
        )
        fr["ros_buffer_capacity"] = _bounded(
            float(fr.get("ros_buffer_capacity", 1.0)) * max(0.75, 1.0 - 0.08 * shared_sev),
            0.1,
            1.5,
        )
        prom_mel["leakage_fraction"] = _bounded(
            float(prom_mel["leakage_fraction"]) + 0.010 * shared_sev,
            0.0,
            0.25,
        )
        for promoter in scw_promoters:
            promoter["leakage_fraction"] = _bounded(
                float(promoter["leakage_fraction"]) + 0.012 * shared_sev,
                0.0,
                0.25,
            )

    timing_sev = _shock_severity(
        rng,
        float(profile.get("timing_overlap", {}).get("probability", 0.0)),
        float(profile.get("timing_overlap", {}).get("sigma", 0.0)),
    )
    meta["timing_overlap"] = float(timing_sev)
    if timing_sev > 0.0:
        out_sampled["mat_activation_dpa"] = _bounded(
            out_sampled["mat_activation_dpa"] - 0.70 * timing_sev,
            28.0,
            45.0,
        )
        out_sampled["scw_activation_dpa"] = _bounded(
            out_sampled["scw_activation_dpa"] + 0.35 * timing_sev,
            22.0,
            40.0,
        )
        out_sampled["late_retention_factor"] = _bounded(
            out_sampled["late_retention_factor"] - 0.04 * timing_sev,
            0.0,
            0.9,
        )
        prom_mel["hill_coefficient"] = max(
            1.0,
            float(prom_mel["hill_coefficient"]) * max(0.72, 1.0 - 0.10 * timing_sev),
        )
        prom_mel["leakage_fraction"] = _bounded(
            float(prom_mel["leakage_fraction"]) + 0.018 * timing_sev,
            0.0,
            0.25,
        )
        for promoter in scw_promoters:
            promoter["hill_coefficient"] = max(
                1.0,
                float(promoter["hill_coefficient"]) * max(0.80, 1.0 - 0.06 * timing_sev),
            )

    storage_sev = _shock_severity(
        rng,
        float(profile.get("storage_slip", {}).get("probability", 0.0)),
        float(profile.get("storage_slip", {}).get("sigma", 0.0)),
    )
    meta["storage_slip"] = float(storage_sev)
    if storage_sev > 0.0:
        if bool(out_params.get("construct", {}).get("use_vacuolar_transit_peptides", False)):
            comp["vacuolar_sequestration_fraction"] = _bounded(
                float(comp.get("vacuolar_sequestration_fraction", 0.0)) - 0.15 * storage_sev,
                0.0,
                0.95,
            )
            comp["cytosolic_quinone_leak_fraction"] = _bounded(
                float(comp.get("cytosolic_quinone_leak_fraction", 1.0)) + 0.12 * storage_sev,
                0.05,
                1.0,
            )
        fr["ros_buffer_capacity"] = _bounded(
            float(fr.get("ros_buffer_capacity", 1.0)) * max(0.65, 1.0 - 0.10 * storage_sev),
            0.1,
            1.5,
        )
        out_sampled["melanin_efficiency"] = _bounded(
            out_sampled["melanin_efficiency"] * max(0.60, 1.0 - 0.06 * storage_sev),
            0.4,
            2.0,
        )
        out_sampled["late_retention_factor"] = _bounded(
            out_sampled["late_retention_factor"] - 0.03 * storage_sev,
            0.0,
            0.9,
        )

    meta["total_severity"] = float(shared_sev + timing_sev + storage_sev)
    return out_sampled, out_params, meta


def is_success(trial: dict[str, Any], thresholds: dict[str, float]) -> bool:
    tox_thr = float(thresholds.get("max_toxicity_pre_cellulose", 1e9))
    tox_val = float(trial.get("toxicity_pre_cellulose", 0.0))
    return (
        float(trial["color_L"]) <= float(thresholds["max_color_L"])
        and float(trial["strength_g_tex"]) >= float(thresholds["min_strength_g_tex"])
        and float(trial["yield_index"]) >= float(thresholds["min_yield_index"])
        and float(trial["temporal_gap_days"]) >= float(thresholds["min_temporal_gap_days"])
        and tox_val <= tox_thr
    )


def _percentiles(values: np.ndarray) -> tuple[float, float, float]:
    p10, p50, p90 = np.percentile(values, [10, 50, 90])
    return float(p10), float(p50), float(p90)


def summarize_trials(trials: list[dict[str, Any]], thresholds: dict[str, float]) -> dict[str, float]:
    if not trials:
        return {
            "n_trials": 0,
            "success_rate": 0.0,
            "robust_score": 0.0,
            "risk_temporal_overlap": 1.0,
            "risk_strength_failure": 1.0,
            "risk_darkness_failure": 1.0,
            "risk_yield_failure": 1.0,
            "risk_toxicity_failure": 1.0,
            "fragility_index": 1.0,
        }

    arr_color = np.array([float(t["color_L"]) for t in trials], dtype=float)
    arr_strength = np.array([float(t["strength_g_tex"]) for t in trials], dtype=float)
    arr_yield = np.array([float(t["yield_index"]) for t in trials], dtype=float)
    arr_gap = np.array([float(t["temporal_gap_days"]) for t in trials], dtype=float)
    arr_toxicity = np.array([float(t.get("toxicity_pre_cellulose", 0.0)) for t in trials], dtype=float)
    arr_comp = np.array([float(t.get("composite_score", 0.0)) for t in trials], dtype=float)

    success_flags = np.array([1.0 if is_success(t, thresholds) else 0.0 for t in trials], dtype=float)
    success_rate = float(success_flags.mean())

    comp_p10, comp_p50, comp_p90 = _percentiles(arr_comp)
    color_p10, color_p50, color_p90 = _percentiles(arr_color)
    strength_p10, strength_p50, strength_p90 = _percentiles(arr_strength)
    yield_p10, yield_p50, yield_p90 = _percentiles(arr_yield)
    gap_p10, gap_p50, gap_p90 = _percentiles(arr_gap)
    tox_p10, tox_p50, tox_p90 = _percentiles(arr_toxicity)

    risk_temporal_overlap = float(np.mean(arr_gap < thresholds["min_temporal_gap_days"]))
    risk_strength_failure = float(np.mean(arr_strength < thresholds["min_strength_g_tex"]))
    risk_darkness_failure = float(np.mean(arr_color > thresholds["max_color_L"]))
    risk_yield_failure = float(np.mean(arr_yield < thresholds["min_yield_index"]))
    tox_thr = float(thresholds.get("max_toxicity_pre_cellulose", 1e9))
    risk_toxicity_failure = float(np.mean(arr_toxicity > tox_thr))

    fragility_index = float(
        np.clip(
            0.30 * risk_darkness_failure
            + 0.25 * risk_temporal_overlap
            + 0.25 * risk_toxicity_failure
            + 0.10 * risk_strength_failure
            + 0.10 * risk_yield_failure,
            0.0,
            1.0,
        )
    )

    # Weighted objective: mostly "chance of meeting all thresholds", partially
    # resilient quality among pass/fail boundary regions, while discouraging
    # high-tail fragility even if average performance is strong.
    robust_score = (
        0.68 * success_rate
        + 0.27 * float(np.clip(comp_p50, 0.0, 1.0))
        + 0.05 * (1.0 - fragility_index)
    )

    return {
        "n_trials": int(len(trials)),
        "success_rate": success_rate,
        "robust_score": float(robust_score),
        "risk_temporal_overlap": risk_temporal_overlap,
        "risk_strength_failure": risk_strength_failure,
        "risk_darkness_failure": risk_darkness_failure,
        "risk_yield_failure": risk_yield_failure,
        "risk_toxicity_failure": risk_toxicity_failure,
        "fragility_index": fragility_index,
        "p10_color_L": color_p10,
        "p50_color_L": color_p50,
        "p90_color_L": color_p90,
        "p10_strength_g_tex": strength_p10,
        "p50_strength_g_tex": strength_p50,
        "p90_strength_g_tex": strength_p90,
        "p10_yield_index": yield_p10,
        "p50_yield_index": yield_p50,
        "p90_yield_index": yield_p90,
        "p10_temporal_gap_days": gap_p10,
        "p50_temporal_gap_days": gap_p50,
        "p90_temporal_gap_days": gap_p90,
        "p10_toxicity_pre_cellulose": tox_p10,
        "p50_toxicity_pre_cellulose": tox_p50,
        "p90_toxicity_pre_cellulose": tox_p90,
        "p10_composite_score": comp_p10,
        "p50_composite_score": comp_p50,
        "p90_composite_score": comp_p90,
    }


def run_robustness_analysis(
    params: dict[str, Any],
    candidates: list[dict[str, Any]],
    n_trials: int = 120,
    seed: int = 42,
    thresholds: dict[str, float] | None = None,
    noise: dict[str, float] | None = None,
    correlated_profile: str | dict[str, Any] | None = None,
    collect_trial_records: bool = False,
) -> list[dict[str, Any]]:
    thresholds = thresholds or DEFAULT_THRESHOLDS
    noise = noise or DEFAULT_NOISE
    rng = np.random.default_rng(seed)

    robust_results: list[dict[str, Any]] = []
    for idx, candidate in enumerate(candidates, 1):
        summary, trial_records = evaluate_candidate_trials(
            params=params,
            candidate=candidate,
            n_trials=n_trials,
            rng=rng,
            thresholds=thresholds,
            noise=noise,
            correlated_profile=correlated_profile,
            candidate_rank=idx,
            collect_trial_records=collect_trial_records,
        )
        result = {
            "input_rank": idx,
            "base_candidate": {
                "color_L": float(candidate["color_L"]),
                "strength_g_tex": float(candidate["strength_g_tex"]),
                "yield_index": float(candidate["yield_index"]),
                "temporal_gap_days": float(candidate["temporal_gap_days"]),
                "toxicity_pre_cellulose": float(candidate.get("toxicity_pre_cellulose", 0.0)),
                "mat_activation_dpa": float(candidate["mat_activation_dpa"]),
                "scw_activation_dpa": float(candidate["scw_activation_dpa"]),
                "mat_strength": float(candidate["mat_strength"]),
                "scw_strength": float(candidate["scw_strength"]),
                "k_competition": float(candidate["k_competition"]),
                "melanin_efficiency": float(candidate["melanin_efficiency"]),
                "late_retention_factor": float(candidate["late_retention_factor"]),
            },
            **summary,
        }
        if collect_trial_records:
            result["trial_records"] = trial_records
        robust_results.append(result)

    robust_results.sort(
        key=lambda r: (
            r["robust_score"],
            r["success_rate"],
            1.0 - float(r.get("fragility_index", 1.0)),
            r["p50_composite_score"],
            -r["p50_color_L"],  # lower L* preferred
            r["p50_strength_g_tex"],
        ),
        reverse=True,
    )
    return robust_results


def evaluate_candidate_trials(
    params: dict[str, Any],
    candidate: dict[str, Any],
    n_trials: int,
    rng: np.random.Generator,
    thresholds: dict[str, float],
    noise: dict[str, float],
    correlated_profile: str | dict[str, Any] | None,
    candidate_rank: int,
    collect_trial_records: bool = False,
) -> tuple[dict[str, float], list[dict[str, Any]]]:
    """Run Monte Carlo trials for a single candidate."""
    trials: list[dict[str, Any]] = []
    trial_records: list[dict[str, Any]] = []

    for trial_idx in range(n_trials):
        sampled = sample_candidate_parameters(candidate, rng, noise)
        perturbed_params = jitter_promoters(params, rng, noise)
        perturbed_params = jitter_failure_risks(perturbed_params, rng, noise)
        sampled, perturbed_params, correlated_meta = apply_correlated_failure_profile(
            sampled=sampled,
            perturbed_params=perturbed_params,
            rng=rng,
            correlated_profile=correlated_profile,
        )
        sampled_risk = _failure_risk_block(perturbed_params)
        trial = evaluate_candidate_proxy(
            params=perturbed_params,
            mat_activation_dpa=sampled["mat_activation_dpa"],
            scw_activation_dpa=sampled["scw_activation_dpa"],
            mat_strength=sampled["mat_strength"],
            scw_strength=sampled["scw_strength"],
            k_competition=sampled["k_competition"],
            melanin_efficiency=sampled["melanin_efficiency"],
            late_retention_factor=sampled["late_retention_factor"],
        )
        trials.append(trial)

        if collect_trial_records:
            success = is_success(trial, thresholds)
            tyrp1_promoter = promoter_config_for_gene(perturbed_params, "TYRP1")
            dct_promoter = promoter_config_for_gene(perturbed_params, "DCT")
            trial_records.append(
                {
                    "candidate_input_rank": int(candidate_rank),
                    "trial_index": int(trial_idx),
                    "sampled_mat_activation_dpa": float(sampled["mat_activation_dpa"]),
                    "sampled_scw_activation_dpa": float(sampled["scw_activation_dpa"]),
                    "sampled_mat_strength": float(sampled["mat_strength"]),
                    "sampled_scw_strength": float(sampled["scw_strength"]),
                    "sampled_k_competition": float(sampled["k_competition"]),
                    "sampled_melanin_efficiency": float(sampled["melanin_efficiency"]),
                    "sampled_late_retention_factor": float(sampled["late_retention_factor"]),
                    "sampled_hill_melA": float(perturbed_params["promoters"]["pGhMat1"]["hill_coefficient"]),
                    "sampled_hill_scw": float(perturbed_params["promoters"]["pGhSCW_late"]["hill_coefficient"]),
                    "sampled_hill_tyrp1": float(tyrp1_promoter["hill_coefficient"]),
                    "sampled_hill_dct": float(dct_promoter["hill_coefficient"]),
                    "sampled_leak_melA": float(perturbed_params["promoters"]["pGhMat1"]["leakage_fraction"]),
                    "sampled_leak_scw": float(perturbed_params["promoters"]["pGhSCW_late"]["leakage_fraction"]),
                    "sampled_leak_tyrp1": float(tyrp1_promoter["leakage_fraction"]),
                    "sampled_leak_dct": float(dct_promoter["leakage_fraction"]),
                    "sampled_activation_tyrp1_dpa": float(tyrp1_promoter["activation_dpa"]),
                    "sampled_activation_dct_dpa": float(dct_promoter["activation_dpa"]),
                    "sampled_strength_tyrp1": float(tyrp1_promoter["strength_relative"]),
                    "sampled_strength_dct": float(dct_promoter["strength_relative"]),
                    "sampled_copper_loading_fraction": float(sampled_risk["copper_loading_fraction"]),
                    "sampled_tyrosinase_activation_fraction": float(sampled_risk["tyrosinase_activation_fraction"]),
                    "sampled_ros_buffer_capacity": float(sampled_risk["ros_buffer_capacity"]),
                    "sampled_silencing_probability": float(sampled_risk["silencing_probability"]),
                    "sampled_event_expression_cv": float(sampled_risk["event_expression_cv"]),
                    "sampled_vacuolar_sequestration_fraction": float(
                        perturbed_params.get("melanin_pathway", {})
                        .get("compartmentalization", {})
                        .get("vacuolar_sequestration_fraction", 0.0)
                    ),
                    "sampled_cytosolic_quinone_leak_fraction": float(
                        perturbed_params.get("melanin_pathway", {})
                        .get("compartmentalization", {})
                        .get("cytosolic_quinone_leak_fraction", 1.0)
                    ),
                    "correlated_profile_name": str(correlated_meta.get("profile_name", "")),
                    "correlated_shared_expression_drop": float(correlated_meta.get("shared_expression_drop", 0.0)),
                    "correlated_timing_overlap": float(correlated_meta.get("timing_overlap", 0.0)),
                    "correlated_storage_slip": float(correlated_meta.get("storage_slip", 0.0)),
                    "correlated_total_severity": float(correlated_meta.get("total_severity", 0.0)),
                    "color_L": float(trial["color_L"]),
                    "strength_g_tex": float(trial["strength_g_tex"]),
                    "yield_index": float(trial["yield_index"]),
                    "temporal_gap_days": float(trial["temporal_gap_days"]),
                    "toxicity_pre_cellulose": float(trial.get("toxicity_pre_cellulose", 0.0)),
                    "toxicity_ratio": float(trial.get("toxicity_ratio", 0.0)),
                    "toxicity_gate_pass": bool(trial.get("toxicity_gate_pass", True)),
                    "composite_score": float(trial["composite_score"]),
                    "success": bool(success),
                }
            )

    summary = summarize_trials(trials, thresholds)
    return summary, trial_records


def save_outputs(
    robust_results: list[dict[str, Any]],
    n_trials: int,
    thresholds: dict[str, float],
    noise: dict[str, float],
    seed: int,
) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    summary = {
        "n_input_candidates": len(robust_results),
        "n_trials_per_candidate": int(n_trials),
        "seed": int(seed),
        "thresholds": thresholds,
        "noise_model": noise,
        "best_input_rank": int(robust_results[0]["input_rank"]) if robust_results else None,
        "best_robust_score": float(robust_results[0]["robust_score"]) if robust_results else 0.0,
        "best_success_rate": float(robust_results[0]["success_rate"]) if robust_results else 0.0,
    }

    with open(RESULTS_DIR / "robustness_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    with open(RESULTS_DIR / "robustness_top_candidates.json", "w") as f:
        json.dump(robust_results, f, indent=2)


def print_report(robust_results: list[dict[str, Any]]) -> None:
    print("\n" + "=" * 124)
    print("  BLACKCOTTON ROBUSTNESS ANALYSIS")
    print("=" * 124)
    print("\n  Ranking by robust score (uncertainty-aware), not only deterministic best case.")
    print("\n  " + "-" * 122)
    print("  #  in_rank  robust  success  fragility  risk_gap  risk_dark  risk_tox  p50_L*  p50_Str  p50_Yield")
    print("  " + "-" * 122)
    for i, row in enumerate(robust_results[:12], 1):
        print(
            f"  {i:>2}    {row['input_rank']:>2}    {row['robust_score']:>6.3f}   "
            f"{row['success_rate']:>6.3f}    {row.get('fragility_index', 0.0):>7.3f}   "
            f"{row['risk_temporal_overlap']:>7.3f}   {row['risk_darkness_failure']:>8.3f}   "
            f"{row.get('risk_toxicity_failure', 0.0):>8.3f}   "
            f"{row['p50_color_L']:>6.2f}  {row['p50_strength_g_tex']:>7.2f}   {row['p50_yield_index']:>8.3f}"
        )
    print("  " + "-" * 122)
    print("\n  Saved:")
    print("    - results/robustness_summary.json")
    print("    - results/robustness_top_candidates.json")
    print("\n" + "=" * 124)


if __name__ == "__main__":
    print("\n🧬 BlackCotton Robustness Analyzer")
    print("=" * 50)

    top_path = RESULTS_DIR / "top_candidates.json"
    if not top_path.exists():
        raise FileNotFoundError(
            "Missing results/top_candidates.json. Run `python -m src.tradeoff_optimizer` first."
        )

    with open(top_path) as f:
        top_candidates = json.load(f)

    params = load_params()
    n_trials = 120
    seed = 42
    thresholds = DEFAULT_THRESHOLDS
    noise = DEFAULT_NOISE

    robust_results = run_robustness_analysis(
        params=params,
        candidates=top_candidates[:12],
        n_trials=n_trials,
        seed=seed,
        thresholds=thresholds,
        noise=noise,
    )
    save_outputs(robust_results, n_trials=n_trials, thresholds=thresholds, noise=noise, seed=seed)
    print_report(robust_results)
    print("\n✅ Robustness analysis complete!")
