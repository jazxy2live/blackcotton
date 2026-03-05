#!/usr/bin/env python3
"""
tradeoff_optimizer.py — Multi-objective black-cotton optimization
==================================================================

Two-stage strategy:
  1) Coarse sweep with an explicit hard constraint: temporal_gap_days >= 0.
  2) ODE refinement (expression + melanin pathway) for top coarse candidates.

This version also introduces a late-only pigment retention mechanism:
  - Pigment load drives structural penalties.
  - Late-retained pigment boosts visible darkness without extra cellulose penalty.
"""

import copy
import io
import json
from contextlib import redirect_stdout
from itertools import product
from pathlib import Path
from typing import Any

import numpy as np

from src.config_loader import load_config
from src.expression_model import (
    cellulose_accumulation,
    promoter_activity,
    run_simulation,
)
from src.fiber_model import model_engineered_black
from src.melanin_pathway import run_melanin_simulation

BASE_DIR = Path(__file__).resolve().parent.parent
CONFIG_DIR = BASE_DIR / "config"
RESULTS_DIR = BASE_DIR / "results"

# Calibration constants to avoid over-saturation in optical darkness mapping.
PROXY_LOAD_SCALE = 0.82
REFINED_LOAD_SCALE = 0.88
MAX_STRUCTURAL_LOAD = 0.90
MAX_VISIBLE_PIGMENT = 0.98
RETENTION_GAIN = 1.00
DEFAULT_FAILURE_RISKS = {
    "copper_loading_fraction": 1.0,
    "tyrosinase_activation_fraction": 1.0,
    "quinone_stress_threshold": 0.08,
    "ros_buffer_capacity": 1.0,
    "silencing_probability": 0.0,
    "event_expression_cv": 0.0,
}


def load_params() -> dict:
    return load_config()


def sigmoid(x: np.ndarray, midpoint: float, slope: float) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-slope * (x - midpoint)))


def failure_risk_params(params: dict[str, Any]) -> dict[str, float]:
    cfg = dict(DEFAULT_FAILURE_RISKS)
    cfg.update(params.get("failure_risks", {}))
    return {
        "copper_loading_fraction": float(np.clip(cfg["copper_loading_fraction"], 0.0, 1.0)),
        "tyrosinase_activation_fraction": float(np.clip(cfg["tyrosinase_activation_fraction"], 0.0, 1.0)),
        "quinone_stress_threshold": float(max(cfg["quinone_stress_threshold"], 1e-6)),
        "ros_buffer_capacity": float(np.clip(cfg["ros_buffer_capacity"], 0.0, 1.5)),
        "silencing_probability": float(np.clip(cfg["silencing_probability"], 0.0, 0.80)),
        "event_expression_cv": float(np.clip(cfg["event_expression_cv"], 0.0, 0.80)),
    }


def effective_melanin_efficiency(
    params: dict[str, Any],
    melanin_efficiency: float,
) -> tuple[float, dict[str, float]]:
    risk = failure_risk_params(params)
    cofactor_scale = float(risk["copper_loading_fraction"] * risk["tyrosinase_activation_fraction"])
    silencing_scale = float(max(0.05, 1.0 - risk["silencing_probability"]))
    event_scale = float(np.exp(-0.5 * (risk["event_expression_cv"] ** 2)))
    ros_scale = float(np.clip(0.75 + 0.30 * risk["ros_buffer_capacity"], 0.55, 1.10))

    total_scale = cofactor_scale * silencing_scale * event_scale * ros_scale
    adjusted = float(np.clip(float(melanin_efficiency) * total_scale, 0.0, 2.2))
    return adjusted, {
        "cofactor_scale": float(cofactor_scale),
        "silencing_scale": float(silencing_scale),
        "event_scale": float(event_scale),
        "ros_scale": float(ros_scale),
        "total_scale": float(total_scale),
    }


def toxicity_threshold(params: dict[str, Any]) -> float:
    risk = failure_risk_params(params)
    # More ROS buffering increases tolerated intermediate load.
    return float(max(1e-6, risk["quinone_stress_threshold"] * (0.70 + 0.60 * risk["ros_buffer_capacity"])))


def compartmentalization_profile(params: dict[str, Any]) -> dict[str, float]:
    """
    Translate compartment-targeting settings into a cytosolic toxicity scale.

    A higher sequestration fraction and lower leak fraction reduce quinone burden
    experienced by the cytosol, where pre-cellulose toxicity is most damaging.
    """
    mp = params.get("melanin_pathway", {})
    comp = mp.get("compartmentalization", {})
    construct_cfg = params.get("construct", {})
    has_transit_targeting = bool(construct_cfg.get("use_vacuolar_transit_peptides", False))

    # Transit targeting is required for sequestration assumptions to apply.
    if has_transit_targeting:
        sequestration = float(np.clip(comp.get("vacuolar_sequestration_fraction", 0.0), 0.0, 0.95))
        leak = float(np.clip(comp.get("cytosolic_quinone_leak_fraction", 1.0), 0.05, 1.0))
    else:
        sequestration = 0.0
        leak = 1.0
    cytosolic_scale = float(np.clip((1.0 - sequestration) * leak, 0.02, 1.0))
    return {
        "use_vacuolar_transit_peptides": bool(has_transit_targeting),
        "vacuolar_sequestration_fraction": sequestration,
        "cytosolic_quinone_leak_fraction": leak,
        "cytosolic_quinone_scale": cytosolic_scale,
    }


def effective_pigment_levels(
    pigment_load: float,
    late_fraction: float,
    late_retention_factor: float,
) -> tuple[float, float]:
    """
    Convert structural pigment load into visible pigment density.

    Structural penalties use `structural_load`, while optical darkness uses
    `visible_pigment` with a late-retention bonus that cannot grow unbounded.
    """
    structural_load = float(np.clip(pigment_load, 0.0, MAX_STRUCTURAL_LOAD))
    late_fraction = float(np.clip(late_fraction, 0.0, 1.0))
    late_retention_factor = max(float(late_retention_factor), 0.0)

    retention_bonus = late_retention_factor * late_fraction * (1.0 - structural_load) * RETENTION_GAIN
    visible_pigment = float(np.clip(structural_load + retention_bonus, 0.0, MAX_VISIBLE_PIGMENT))
    return structural_load, visible_pigment


def compute_scores(fiber, temporal_gap_days: float, params: dict) -> tuple[float, float, float, float]:
    fd = params["fiber_development"]
    rel_len = fiber.uhml_mm / float(fd["target_length_mm"])
    rel_strength = fiber.strength_g_tex / float(fd["target_strength_g_tex"])
    rel_uniformity = fiber.uniformity_pct / float(fd["target_uniformity"])
    yield_index = float(np.clip(rel_len * rel_strength * (fiber.cellulose_pct / 95.0), 0.0, 1.2))

    darkness_score = float(np.clip((82.0 - fiber.color_L) / (82.0 - 5.0), 0.0, 1.0))
    quality_score = float(np.clip((rel_len + rel_strength + rel_uniformity) / 3.0, 0.0, 1.2))

    temporal_penalty = min(max(-temporal_gap_days, 0.0) / 10.0, 1.0)
    composite_score = (
        0.48 * darkness_score +
        0.30 * min(quality_score, 1.0) +
        0.22 * min(yield_index, 1.0) -
        0.20 * temporal_penalty
    )
    return yield_index, darkness_score, quality_score, composite_score


def compute_proxy_signals(
    params: dict,
    mat_activation_dpa: float,
    scw_activation_dpa: float,
    mat_strength: float,
    scw_strength: float,
    melanin_efficiency: float,
) -> tuple[float, float, float, float]:
    """
    Returns:
        overlap_index
        pigment_load_proxy (0-1)
        temporal_gap_days
        late_fraction_proxy (0-1)
    """
    total_days = params["fiber_development"]["total_development_days"]
    t_dpa = np.linspace(0.0, float(total_days), int(total_days * 8) + 1)  # 3-hour resolution

    prom = params["promoters"]
    mat_cfg = prom["pGhMat1"]
    scw_cfg = prom["pGhSCW_late"]

    mat_peak_delta = max(float(mat_cfg["peak_dpa"]) - float(mat_cfg["activation_dpa"]), 1.0)
    scw_peak_delta = max(float(scw_cfg["peak_dpa"]) - float(scw_cfg["activation_dpa"]), 1.0)
    mat_peak = mat_activation_dpa + mat_peak_delta
    scw_peak = scw_activation_dpa + scw_peak_delta

    p_melA = np.array([
        promoter_activity(
            t,
            mat_activation_dpa,
            mat_peak,
            float(mat_cfg["hill_coefficient"]),
            float(mat_cfg["leakage_fraction"]),
            mat_strength,
        )
        for t in t_dpa
    ])
    p_scw = np.array([
        promoter_activity(
            t,
            scw_activation_dpa,
            scw_peak,
            float(scw_cfg["hill_coefficient"]),
            float(scw_cfg["leakage_fraction"]),
            scw_strength,
        )
        for t in t_dpa
    ])

    max_strength = max(mat_strength, scw_strength, 1.0)
    synergy = np.sqrt(np.clip(p_melA * p_scw, 0.0, None))
    melanin_signal = np.clip((0.55 * p_melA + 0.25 * p_scw + 0.20 * synergy) / max_strength, 0.0, 1.0)

    cellulose = np.array([cellulose_accumulation(t, params) for t in t_dpa])
    cellulose_max = float(cellulose.max()) if cellulose.size else 0.0
    if cellulose_max <= 0:
        return 0.0, 0.0, -999.0, 0.0

    cellulose_norm = cellulose / cellulose_max
    cellulose_rate = np.clip(np.gradient(cellulose_norm, t_dpa), 0.0, None)
    if float(cellulose_rate.sum()) <= 0:
        cellulose_weight = np.full_like(cellulose_rate, 1.0 / len(cellulose_rate))
    else:
        cellulose_weight = cellulose_rate / cellulose_rate.sum()

    overlap_index = float(np.sum(cellulose_weight * melanin_signal))

    late_window = sigmoid(t_dpa, midpoint=37.0, slope=1.1)
    late_window = late_window / float(late_window.sum())

    early_window = cellulose_weight
    early_activity = float(np.sum(melanin_signal * early_window))
    late_activity = float(np.sum(melanin_signal * late_window))

    pigment_load = float(np.clip(melanin_efficiency * max(late_activity - 0.55 * early_activity, 0.0), 0.0, 1.0))
    late_fraction_proxy = float(np.clip(late_activity / (late_activity + early_activity + 1e-9), 0.0, 1.0))

    cellulose_90_idx = int(np.argmax(cellulose_norm >= 0.9))
    melanin_50_idx = int(np.argmax(melanin_signal >= 0.5))
    temporal_gap_days = float(t_dpa[melanin_50_idx] - t_dpa[cellulose_90_idx])

    return overlap_index, pigment_load, temporal_gap_days, late_fraction_proxy


def evaluate_candidate_proxy(
    params: dict,
    mat_activation_dpa: float,
    scw_activation_dpa: float,
    mat_strength: float,
    scw_strength: float,
    k_competition: float,
    melanin_efficiency: float,
    late_retention_factor: float,
) -> dict:
    effective_eff, eff_meta = effective_melanin_efficiency(params, melanin_efficiency)
    comp_profile = compartmentalization_profile(params)
    overlap_index, pigment_load, temporal_gap_days, late_fraction_proxy = compute_proxy_signals(
        params,
        mat_activation_dpa,
        scw_activation_dpa,
        mat_strength,
        scw_strength,
        effective_eff,
    )

    scaled_load = pigment_load * PROXY_LOAD_SCALE
    structural_load, effective_color_pigment = effective_pigment_levels(
        scaled_load,
        late_fraction_proxy,
        late_retention_factor,
    )

    fiber = model_engineered_black(
        params,
        melanin=effective_color_pigment,
        overlap=overlap_index,
        k_competition=k_competition,
        melanin_load=structural_load,
    )
    yield_index, darkness_score, quality_score, composite_score = compute_scores(
        fiber,
        temporal_gap_days,
        params,
    )
    tox_threshold = toxicity_threshold(params)
    toxicity_pre_cellulose_proxy = float(
        overlap_index
        * scaled_load
        * (1.0 + 0.5 * max(effective_eff, 0.0))
        * comp_profile["cytosolic_quinone_scale"]
    )
    toxicity_gate_pass = bool(toxicity_pre_cellulose_proxy <= tox_threshold)
    toxicity_ratio = float(toxicity_pre_cellulose_proxy / max(tox_threshold, 1e-9))
    toxicity_penalty = 0.12 * max(toxicity_ratio - 1.0, 0.0)
    composite_score = float(composite_score - toxicity_penalty)

    return {
        "model_stage": "proxy",
        "mat_activation_dpa": float(mat_activation_dpa),
        "scw_activation_dpa": float(scw_activation_dpa),
        "mat_strength": float(mat_strength),
        "scw_strength": float(scw_strength),
        "k_competition": float(k_competition),
        "melanin_efficiency": float(melanin_efficiency),
        "effective_melanin_efficiency": float(effective_eff),
        "efficiency_scale_total": float(eff_meta["total_scale"]),
        "efficiency_scale_cofactor": float(eff_meta["cofactor_scale"]),
        "efficiency_scale_silencing": float(eff_meta["silencing_scale"]),
        "efficiency_scale_event": float(eff_meta["event_scale"]),
        "efficiency_scale_ros": float(eff_meta["ros_scale"]),
        "late_retention_factor": float(late_retention_factor),
        "overlap_index": float(overlap_index),
        "temporal_gap_days": float(temporal_gap_days),
        "pigment_load": float(structural_load),
        "pigment_density": float(effective_color_pigment),
        "late_fraction": float(late_fraction_proxy),
        "uhml_mm": float(fiber.uhml_mm),
        "strength_g_tex": float(fiber.strength_g_tex),
        "uniformity_pct": float(fiber.uniformity_pct),
        "micronaire": float(fiber.micronaire),
        "color_L": float(fiber.color_L),
        "cellulose_pct": float(fiber.cellulose_pct),
        "yield_index": float(yield_index),
        "darkness_score": float(darkness_score),
        "quality_score": float(quality_score),
        "composite_score": float(composite_score),
        "toxicity_pre_cellulose": float(toxicity_pre_cellulose_proxy),
        "toxicity_threshold": float(tox_threshold),
        "toxicity_ratio": float(toxicity_ratio),
        "toxicity_gate_pass": bool(toxicity_gate_pass),
        "use_vacuolar_transit_peptides": bool(comp_profile["use_vacuolar_transit_peptides"]),
        "vacuolar_sequestration_fraction": float(comp_profile["vacuolar_sequestration_fraction"]),
        "cytosolic_quinone_leak_fraction": float(comp_profile["cytosolic_quinone_leak_fraction"]),
        "cytosolic_quinone_scale": float(comp_profile["cytosolic_quinone_scale"]),
        "grade": fiber.grade(),
    }


def build_candidate_params(base_params: dict, candidate: dict) -> dict:
    params = copy.deepcopy(base_params)
    prom = params["promoters"]

    base_mat = base_params["promoters"]["pGhMat1"]
    base_scw = base_params["promoters"]["pGhSCW_late"]
    total_days = float(params["fiber_development"]["total_development_days"])

    mat_peak_delta = max(float(base_mat["peak_dpa"]) - float(base_mat["activation_dpa"]), 1.0)
    scw_peak_delta = max(float(base_scw["peak_dpa"]) - float(base_scw["activation_dpa"]), 1.0)

    prom["pGhMat1"]["activation_dpa"] = float(candidate["mat_activation_dpa"])
    prom["pGhMat1"]["peak_dpa"] = float(min(candidate["mat_activation_dpa"] + mat_peak_delta, total_days))
    prom["pGhMat1"]["strength_relative"] = float(candidate["mat_strength"])

    prom["pGhSCW_late"]["activation_dpa"] = float(candidate["scw_activation_dpa"])
    prom["pGhSCW_late"]["peak_dpa"] = float(min(candidate["scw_activation_dpa"] + scw_peak_delta, total_days))
    prom["pGhSCW_late"]["strength_relative"] = float(candidate["scw_strength"])
    return params


def run_odes_silent(func, *args, **kwargs):
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        return func(*args, **kwargs)


def simulate_baseline_melanin(params: dict) -> float:
    expr = run_odes_silent(run_simulation, params)
    mel = run_odes_silent(run_melanin_simulation, params, expr)
    return max(float(mel["melanin"][-1]), 1e-9)


def simulate_candidate_profile(
    base_params: dict,
    candidate: dict,
    baseline_melanin: float,
    cache: dict,
) -> dict:
    key = (
        float(candidate["mat_activation_dpa"]),
        float(candidate["scw_activation_dpa"]),
        float(candidate["mat_strength"]),
        float(candidate["scw_strength"]),
    )
    if key in cache:
        return cache[key]

    params = build_candidate_params(base_params, candidate)
    comp_profile = compartmentalization_profile(params)
    expr = run_odes_silent(run_simulation, params)
    mel = run_odes_silent(run_melanin_simulation, params, expr)

    t_expr = expr["t_dpa"]
    cellulose = expr["cellulose"]
    melA = expr["protein_melA"]

    cellulose_max = float(cellulose.max()) if cellulose.size else 0.0
    if cellulose_max <= 0:
        profile = {
            "overlap_index": 0.0,
            "temporal_gap_days": -999.0,
            "final_melanin_mM": 0.0,
            "pigment_from_ode": 0.0,
            "late_fraction": 0.0,
            "toxicity_pre_cellulose": 0.0,
        }
        cache[key] = profile
        return profile

    cellulose_norm = cellulose / cellulose_max
    cellulose_90_idx = int(np.argmax(cellulose_norm >= 0.9))
    cellulose_90_dpa = float(t_expr[cellulose_90_idx])

    melA_max = float(melA.max()) if melA.size else 0.0
    if melA_max <= 0:
        melA_half_dpa = 0.0
        melA_norm = np.zeros_like(melA)
    else:
        melA_half_idx = int(np.argmax(melA >= 0.5 * melA_max))
        melA_half_dpa = float(t_expr[melA_half_idx])
        melA_norm = melA / melA_max

    temporal_gap_days = melA_half_dpa - cellulose_90_dpa

    cellulose_rate = np.clip(np.gradient(cellulose_norm, t_expr), 0.0, None)
    if float(cellulose_rate.sum()) <= 0:
        cellulose_weight = np.full_like(cellulose_rate, 1.0 / len(cellulose_rate))
    else:
        cellulose_weight = cellulose_rate / cellulose_rate.sum()
    overlap_index = float(np.sum(cellulose_weight * melA_norm))

    # Expression-timing efficiency drives how much melanin contributes to
    # visible pigment without over-penalizing cellulose structure.
    expr_total = float(np.trapezoid(melA_norm, t_expr))
    expr_late = float(np.trapezoid(melA_norm[cellulose_90_idx:], t_expr[cellulose_90_idx:])) if expr_total > 0 else 0.0
    expression_late_fraction = float(expr_late / expr_total) if expr_total > 0 else 0.0

    final_melanin = float(mel["melanin"][-1])
    timing_efficiency = np.clip(
        0.35 + 0.90 * expression_late_fraction - 0.85 * overlap_index,
        0.20,
        1.05,
    )
    pigment_from_ode = float(np.clip(
        (final_melanin / baseline_melanin) * timing_efficiency * 0.75,
        0.0,
        1.0,
    ))

    t_mel = mel["t_dpa"]
    melanin_arr = mel["melanin"]
    idx_after_90 = int(np.argmax(t_mel >= cellulose_90_dpa))
    late_gain = max(float(melanin_arr[-1] - melanin_arr[idx_after_90]), 0.0)
    pathway_late_fraction = float(late_gain / melanin_arr[-1]) if float(melanin_arr[-1]) > 0 else 0.0
    late_fraction = max(expression_late_fraction, pathway_late_fraction)

    tox_slice = slice(0, idx_after_90 + 1)
    toxicity_pre_cellulose_raw = float(max(
        np.max(mel["dopaquinone"][tox_slice]),
        np.max(mel["dopachrome"][tox_slice]),
    ))
    toxicity_pre_cellulose = float(toxicity_pre_cellulose_raw * comp_profile["cytosolic_quinone_scale"])

    profile = {
        "overlap_index": float(overlap_index),
        "temporal_gap_days": float(temporal_gap_days),
        "final_melanin_mM": float(final_melanin),
        "pigment_from_ode": float(pigment_from_ode),
        "late_fraction": float(late_fraction),
        "expression_late_fraction": float(expression_late_fraction),
        "pathway_late_fraction": float(pathway_late_fraction),
        "toxicity_pre_cellulose_raw": float(toxicity_pre_cellulose_raw),
        "toxicity_pre_cellulose": float(toxicity_pre_cellulose),
        "use_vacuolar_transit_peptides": bool(comp_profile["use_vacuolar_transit_peptides"]),
        "vacuolar_sequestration_fraction": float(comp_profile["vacuolar_sequestration_fraction"]),
        "cytosolic_quinone_leak_fraction": float(comp_profile["cytosolic_quinone_leak_fraction"]),
        "cytosolic_quinone_scale": float(comp_profile["cytosolic_quinone_scale"]),
    }
    cache[key] = profile
    return profile


def refine_candidates_with_odes(
    params: dict,
    coarse_candidates: list[dict],
    max_refine: int = 80,
) -> list[dict]:
    if not coarse_candidates:
        return []

    baseline_melanin = simulate_baseline_melanin(params)
    tox_threshold = toxicity_threshold(params)
    cache = {}
    refined = []

    for candidate in coarse_candidates[:max_refine]:
        profile = simulate_candidate_profile(params, candidate, baseline_melanin, cache)

        # Hard constraint: only keep timing-safe solutions.
        if profile["temporal_gap_days"] < 0:
            continue
        # Fail-fast filter for early toxic intermediate accumulation.
        if profile["toxicity_pre_cellulose"] > tox_threshold:
            continue

        effective_eff, eff_meta = effective_melanin_efficiency(params, candidate["melanin_efficiency"])
        scaled_load = profile["pigment_from_ode"] * effective_eff * REFINED_LOAD_SCALE
        pigment_load, effective_color_pigment = effective_pigment_levels(
            scaled_load,
            profile["late_fraction"],
            candidate["late_retention_factor"],
        )

        fiber = model_engineered_black(
            params,
            melanin=effective_color_pigment,
            overlap=profile["overlap_index"],
            k_competition=candidate["k_competition"],
            melanin_load=pigment_load,
        )
        yield_index, darkness_score, quality_score, composite_score = compute_scores(
            fiber,
            profile["temporal_gap_days"],
            params,
        )
        tox_ratio = float(profile["toxicity_pre_cellulose"] / max(tox_threshold, 1e-9))
        toxicity_penalty = 0.12 * max(tox_ratio - 1.0, 0.0)
        composite_score = float(composite_score - toxicity_penalty)

        refined.append({
            **candidate,
            "model_stage": "ode_refined",
            "overlap_index": float(profile["overlap_index"]),
            "temporal_gap_days": float(profile["temporal_gap_days"]),
            "final_melanin_mM": float(profile["final_melanin_mM"]),
            "late_fraction": float(profile["late_fraction"]),
            "effective_melanin_efficiency": float(effective_eff),
            "efficiency_scale_total": float(eff_meta["total_scale"]),
            "efficiency_scale_cofactor": float(eff_meta["cofactor_scale"]),
            "efficiency_scale_silencing": float(eff_meta["silencing_scale"]),
            "efficiency_scale_event": float(eff_meta["event_scale"]),
            "efficiency_scale_ros": float(eff_meta["ros_scale"]),
            "toxicity_pre_cellulose_raw": float(profile.get("toxicity_pre_cellulose_raw", profile["toxicity_pre_cellulose"])),
            "toxicity_pre_cellulose": float(profile["toxicity_pre_cellulose"]),
            "toxicity_threshold": float(tox_threshold),
            "toxicity_ratio": float(tox_ratio),
            "toxicity_gate_pass": True,
            "use_vacuolar_transit_peptides": bool(profile.get("use_vacuolar_transit_peptides", False)),
            "vacuolar_sequestration_fraction": float(profile.get("vacuolar_sequestration_fraction", 0.0)),
            "cytosolic_quinone_leak_fraction": float(profile.get("cytosolic_quinone_leak_fraction", 1.0)),
            "cytosolic_quinone_scale": float(profile.get("cytosolic_quinone_scale", 1.0)),
            "pigment_load": float(pigment_load),
            "pigment_density": float(effective_color_pigment),
            "uhml_mm": float(fiber.uhml_mm),
            "strength_g_tex": float(fiber.strength_g_tex),
            "uniformity_pct": float(fiber.uniformity_pct),
            "micronaire": float(fiber.micronaire),
            "color_L": float(fiber.color_L),
            "cellulose_pct": float(fiber.cellulose_pct),
            "yield_index": float(yield_index),
            "darkness_score": float(darkness_score),
            "quality_score": float(quality_score),
            "composite_score": float(composite_score),
            "grade": fiber.grade(),
        })

    return refined


def dominates(a: dict, b: dict) -> bool:
    not_worse = (
        a["color_L"] <= b["color_L"] and
        a["strength_g_tex"] >= b["strength_g_tex"] and
        a["yield_index"] >= b["yield_index"]
    )
    strictly_better = (
        a["color_L"] < b["color_L"] or
        a["strength_g_tex"] > b["strength_g_tex"] or
        a["yield_index"] > b["yield_index"]
    )
    return not_worse and strictly_better


def pareto_front(candidates: list[dict]) -> list[dict]:
    front = []
    for i, candidate in enumerate(candidates):
        dominated = False
        for j, other in enumerate(candidates):
            if i == j:
                continue
            if dominates(other, candidate):
                dominated = True
                break
        if not dominated:
            front.append(candidate)
    deduped = []
    seen = set()
    for c in front:
        key = (
            round(c["color_L"], 4),
            round(c["strength_g_tex"], 4),
            round(c["yield_index"], 4),
            round(c["temporal_gap_days"], 4),
        )
        if key in seen:
            continue
        seen.add(key)
        deduped.append(c)
    return deduped


def select_top_candidates(candidates: list[dict], n: int = 12) -> list[dict]:
    viable = [
        c for c in candidates
        if (
            c["color_L"] <= 32.0
            and c["strength_g_tex"] >= 27.8
            and c["yield_index"] >= 0.82
            and bool(c.get("toxicity_gate_pass", True))
        )
    ]
    pool = viable if viable else candidates
    ranked = sorted(pool, key=lambda c: c["composite_score"], reverse=True)

    # Keep top list diverse across objective space.
    selected = []
    seen = set()
    for c in ranked:
        key = (
            round(c["color_L"], 2),
            round(c["strength_g_tex"], 2),
            round(c["yield_index"], 3),
            round(c["temporal_gap_days"], 2),
        )
        if key in seen:
            continue
        seen.add(key)
        selected.append(c)
        if len(selected) >= n:
            break
    return selected


def select_seed_candidates(candidates: list[dict], n: int = 80) -> list[dict]:
    if not candidates:
        return []

    ordered = sorted(candidates, key=lambda c: c["composite_score"], reverse=True)
    selected = []
    seen = set()

    def add_candidate(c):
        key = (
            c["mat_activation_dpa"],
            c["scw_activation_dpa"],
            c["mat_strength"],
            c["scw_strength"],
            c["k_competition"],
            c["melanin_efficiency"],
            c["late_retention_factor"],
        )
        if key in seen:
            return
        selected.append(c)
        seen.add(key)

    # Best overall.
    for c in ordered[: max(n // 3, 1)]:
        add_candidate(c)
        if len(selected) >= n:
            return selected[:n]

    # Force coverage across melanin efficiency, retention, and competition regimes.
    for field, quota_divisor in [
        ("melanin_efficiency", 10),
        ("late_retention_factor", 12),
        ("k_competition", 12),
    ]:
        values = sorted({c[field] for c in ordered})
        per_value_quota = max(n // max(len(values), 1) // quota_divisor + 1, 2)
        for value in values:
            subset = [c for c in ordered if c[field] == value]
            for c in subset[:per_value_quota]:
                add_candidate(c)
                if len(selected) >= n:
                    return selected[:n]

    # Add darkest safe candidates if capacity remains.
    for c in sorted(ordered, key=lambda c: c["color_L"]):
        add_candidate(c)
        if len(selected) >= n:
            return selected[:n]

    return selected[:n]


def run_sweep(params: dict) -> list[dict]:
    mat_activation_values = np.arange(32.0, 40.1, 2.0)  # 32-40 DPA
    scw_activation_values = np.arange(26.0, 34.1, 2.0)  # 26-34 DPA
    mat_strength_values = [0.8, 1.0, 1.2, 1.4]
    scw_strength_values = [0.6, 0.8, 1.0, 1.2]
    k_competition_values = [0.10, 0.20, 0.30, 0.40]
    melanin_eff_values = [0.85, 1.00, 1.15, 1.30]
    late_retention_values = [0.0, 0.15, 0.30, 0.45]

    candidates = []
    for mat_act, scw_act, mat_str, scw_str, k_comp, mel_eff, late_ret in product(
        mat_activation_values,
        scw_activation_values,
        mat_strength_values,
        scw_strength_values,
        k_competition_values,
        melanin_eff_values,
        late_retention_values,
    ):
        candidate = evaluate_candidate_proxy(
            params=params,
            mat_activation_dpa=float(mat_act),
            scw_activation_dpa=float(scw_act),
            mat_strength=float(mat_str),
            scw_strength=float(scw_str),
            k_competition=float(k_comp),
            melanin_efficiency=float(mel_eff),
            late_retention_factor=float(late_ret),
        )

        # Hard constraint in coarse stage: no overlap-first regimes.
        if candidate["temporal_gap_days"] < 0.0:
            continue
        # Hard constraint in coarse stage: reject pre-cellulose toxicity tails.
        if not bool(candidate.get("toxicity_gate_pass", True)):
            continue

        candidates.append(candidate)

    return candidates


def save_outputs(
    coarse_candidates: list[dict],
    refined_candidates: list[dict],
    front: list[dict],
    top: list[dict],
) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    with open(RESULTS_DIR / "optimization_sweep.json", "w") as f:
        json.dump(coarse_candidates, f, indent=2)
    with open(RESULTS_DIR / "optimization_refined.json", "w") as f:
        json.dump(refined_candidates, f, indent=2)
    with open(RESULTS_DIR / "pareto_front.json", "w") as f:
        json.dump(front, f, indent=2)
    with open(RESULTS_DIR / "top_candidates.json", "w") as f:
        json.dump(top, f, indent=2)

    summary = {
        "constraint": "temporal_gap_days >= 0 AND toxicity_pre_cellulose <= threshold (hard)",
        "num_candidates_coarse": len(coarse_candidates),
        "num_candidates_refined": len(refined_candidates),
        "num_pareto_refined": len(front),
        "num_top": len(top),
        "best_composite_score": max((c["composite_score"] for c in refined_candidates), default=0.0),
        "best_darkness_L": min((c["color_L"] for c in refined_candidates), default=82.0),
        "best_strength_g_tex": max((c["strength_g_tex"] for c in refined_candidates), default=0.0),
        "best_yield_index": max((c["yield_index"] for c in refined_candidates), default=0.0),
    }
    with open(RESULTS_DIR / "optimization_summary.json", "w") as f:
        json.dump(summary, f, indent=2)


def print_candidate_table(title: str, candidates: list[dict]) -> None:
    print(f"\n  {title}")
    print("  " + "-" * 106)
    print(
        f"  {'#':>2} {'L*':>6} {'Str':>6} {'Yield':>7} {'Gap(d)':>7} {'ovlp':>6} "
        f"{'late%':>6} {'ret':>5}  {'mat/scw act':<13} {'mat/scw str':<13} {'k':>4}"
    )
    print("  " + "-" * 106)
    for i, c in enumerate(candidates, 1):
        print(
            f"  {i:>2} {c['color_L']:>6.1f} {c['strength_g_tex']:>6.1f} {c['yield_index']:>7.3f} "
            f"{c['temporal_gap_days']:>7.2f} {c['overlap_index']:>6.3f} "
            f"{100.0 * c.get('late_fraction', 0.0):>6.1f} {c.get('late_retention_factor', 0.0):>5.2f}  "
            f"{c['mat_activation_dpa']:.0f}/{c['scw_activation_dpa']:.0f} DPA   "
            f"{c['mat_strength']:.1f}/{c['scw_strength']:.1f}      {c['k_competition']:>4.2f}"
        )
    print("  " + "-" * 106)


def print_report(
    coarse_candidates: list[dict],
    refined_candidates: list[dict],
    front: list[dict],
    top: list[dict],
) -> None:
    print("\n" + "=" * 110)
    print("  BLACKCOTTON TRADEOFF OPTIMIZER RESULTS")
    print("=" * 110)
    print("\n  Hard constraints enforced: temporal_gap_days >= 0 and toxicity_pre_cellulose <= threshold")
    print(f"  Coarse candidates (safe region): {len(coarse_candidates)}")
    print(f"  ODE-refined candidates:          {len(refined_candidates)}")
    print(f"  Pareto-optimal refined:          {len(front)}")
    print(f"  Top candidates reported:         {len(top)}")
    print_candidate_table("TOP CANDIDATES (ODE-Refined, Safe-Timing Only)", top)
    print("\n  Saved:")
    print("    - results/optimization_sweep.json")
    print("    - results/optimization_refined.json")
    print("    - results/pareto_front.json")
    print("    - results/top_candidates.json")
    print("    - results/optimization_summary.json")
    print("\n" + "=" * 110)


if __name__ == "__main__":
    print("\n🧬 BlackCotton Tradeoff Optimizer")
    print("=" * 50)

    params = load_params()

    coarse_candidates = run_sweep(params)
    coarse_seed = select_seed_candidates(coarse_candidates, n=120)
    refined_candidates = refine_candidates_with_odes(
        params=params,
        coarse_candidates=coarse_seed,
        max_refine=120,
    )

    front = pareto_front(refined_candidates)
    top = select_top_candidates(refined_candidates, n=12)

    save_outputs(coarse_candidates, refined_candidates, front, top)
    print_report(coarse_candidates, refined_candidates, front, top)

    print("\n✅ Optimization complete!")
