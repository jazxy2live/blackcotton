#!/usr/bin/env python3
"""
anti_gibberish_audit.py — Ablation + holdout sanity audit for BlackCotton
==========================================================================

Purpose:
  1) Ablation test: isolate whether key model upgrades actually drive robustness
  2) Holdout stress test: evaluate top candidates on scenarios not used in tuning
  3) Confidence report: separate strong evidence from unresolved assumptions

Usage:
  python -m src.anti_gibberish_audit
"""

import argparse
import copy
import json
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
from typing import Any

import numpy as np

from src.adversarial_robustness_suite import scale_noise, strict_thresholds
from src.construct_designer import (
    GeneticElement,
    SEQ_DIR,
    maybe_apply_transit_peptides,
    resolve_melA_cds,
    resolve_strict_cds,
)
from src.robustness_analyzer import DEFAULT_NOISE, DEFAULT_THRESHOLDS, run_robustness_analysis
from src.tradeoff_optimizer import evaluate_candidate_proxy, load_params
from src.worst_case_hardening_sprint import pick_top_unique_robust_rows, select_seed_candidates, strict_margin_score

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"

AUDIT_GRID = {
    "mat_activation_dpa": [37.0, 38.0, 39.0, 40.0, 41.0, 42.0],
    "scw_activation_dpa": [32.0, 33.0, 34.0, 35.0],
    "mat_strength": [1.1, 1.3],
    "scw_strength": [1.0, 1.2, 1.3],
    "k_competition": [0.03, 0.07, 0.12],
    "melanin_efficiency": [1.5, 1.6, 1.7],
    "late_retention_factor": [0.65, 0.75, 0.85],
}

VARIANTS = [
    {
        "id": "all_off",
        "label": "All Off (control)",
        "traffic": False,
        "storage": False,
        "transit": False,
    },
    {
        "id": "transit_only",
        "label": "Transit Only",
        "traffic": False,
        "storage": False,
        "transit": True,
    },
    {
        "id": "traffic_only",
        "label": "Traffic Only",
        "traffic": True,
        "storage": False,
        "transit": False,
    },
    {
        "id": "storage_only",
        "label": "Storage Only (Transit+Compartment)",
        "traffic": False,
        "storage": True,
        "transit": True,
    },
    {
        "id": "all_on",
        "label": "All On",
        "traffic": True,
        "storage": True,
        "transit": True,
    },
]


def strict_high_noise_scenario() -> dict[str, Any]:
    return {
        "name": "strict_high_noise",
        "thresholds": strict_thresholds(dict(DEFAULT_THRESHOLDS)),
        "noise": scale_noise(dict(DEFAULT_NOISE), factor=1.8),
    }


def holdout_scenarios() -> list[dict[str, Any]]:
    base_strict = strict_thresholds(dict(DEFAULT_THRESHOLDS))

    holdout_a_thr = dict(base_strict)
    holdout_a_thr["max_color_L"] = 21.5
    holdout_a_thr["min_temporal_gap_days"] = 0.75
    holdout_a_thr["max_toxicity_pre_cellulose"] = 0.075
    holdout_a_noise = scale_noise(dict(DEFAULT_NOISE), factor=2.2)
    holdout_a_noise["timing_sigma_dpa"] = max(float(holdout_a_noise["timing_sigma_dpa"]), 1.40)

    holdout_b_thr = dict(base_strict)
    holdout_b_thr["max_color_L"] = 21.0
    holdout_b_thr["min_strength_g_tex"] = 28.6
    holdout_b_thr["min_yield_index"] = 0.89
    holdout_b_thr["min_temporal_gap_days"] = 0.75
    holdout_b_thr["max_toxicity_pre_cellulose"] = 0.080
    holdout_b_noise = scale_noise(dict(DEFAULT_NOISE), factor=1.9)
    holdout_b_noise["competition_cv"] = float(holdout_b_noise["competition_cv"]) + 0.06
    holdout_b_noise["efficiency_cv"] = float(holdout_b_noise["efficiency_cv"]) + 0.05

    return [
        {
            "name": "holdout_ultra_noise_gate",
            "description": "Unseen: stricter darkness/toxicity gate + heavier timing/process noise",
            "thresholds": holdout_a_thr,
            "noise": holdout_a_noise,
        },
        {
            "name": "holdout_dark_yield_clamp",
            "description": "Unseen: deeper-black + higher strength/yield under elevated uncertainty",
            "thresholds": holdout_b_thr,
            "noise": holdout_b_noise,
        },
    ]


def _clone_gene(gene: GeneticElement) -> GeneticElement:
    return GeneticElement(
        name=str(gene.name),
        element_type=str(gene.element_type),
        sequence=str(gene.sequence),
        description=str(gene.description),
    )


def build_gene_templates() -> dict[str, GeneticElement]:
    melA_seq, _ = resolve_melA_cds()
    tyrp1_seq, _ = resolve_strict_cds(
        real_path=SEQ_DIR / "TYRP1_CDS_real.fasta",
        fallback_path=SEQ_DIR / "GhTYRP1_optimized.fasta",
    )
    dct_seq, _ = resolve_strict_cds(
        real_path=SEQ_DIR / "DCT_CDS_real.fasta",
        fallback_path=SEQ_DIR / "GhDCT_optimized.fasta",
    )
    return {
        "melA": GeneticElement(
            name="melA",
            element_type="gene",
            sequence=melA_seq,
            description="Tyrosinase",
        ),
        "TYRP1": GeneticElement(
            name="TYRP1",
            element_type="gene",
            sequence=tyrp1_seq,
            description="TYRP1",
        ),
        "DCT": GeneticElement(
            name="DCT",
            element_type="gene",
            sequence=dct_seq,
            description="DCT",
        ),
    }


def traffic_off_block() -> dict[str, float]:
    return {
        "melA_entry_scale": 1.0,
        "dopaquinone_drain_scale": 1.0,
        "dopachrome_drain_scale": 1.0,
        "dct_flux_scale": 1.0,
        "tyrp1_flux_scale": 1.0,
        "eumelanin_push_scale": 1.0,
    }


def apply_variant(base_params: dict[str, Any], variant: dict[str, Any]) -> dict[str, Any]:
    params = copy.deepcopy(base_params)
    mp = params.setdefault("melanin_pathway", {})
    construct = params.setdefault("construct", {})

    base_tc = copy.deepcopy(base_params.get("melanin_pathway", {}).get("traffic_control", {}))
    base_comp = copy.deepcopy(base_params.get("melanin_pathway", {}).get("compartmentalization", {}))

    if bool(variant["traffic"]):
        mp["traffic_control"] = base_tc if base_tc else traffic_off_block()
    else:
        mp["traffic_control"] = traffic_off_block()

    if bool(variant["storage"]):
        mp["compartmentalization"] = {
            "vacuolar_sequestration_fraction": float(base_comp.get("vacuolar_sequestration_fraction", 0.70)),
            "cytosolic_quinone_leak_fraction": float(base_comp.get("cytosolic_quinone_leak_fraction", 0.35)),
        }
    else:
        mp["compartmentalization"] = {
            "vacuolar_sequestration_fraction": 0.0,
            "cytosolic_quinone_leak_fraction": 1.0,
        }

    construct["use_vacuolar_transit_peptides"] = bool(variant["transit"])
    return params


def build_deterministic_pool_from_grid(
    params: dict[str, Any],
    strict_thresholds_cfg: dict[str, float],
    min_temporal_gap_days: float,
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for values in product(
        AUDIT_GRID["mat_activation_dpa"],
        AUDIT_GRID["scw_activation_dpa"],
        AUDIT_GRID["mat_strength"],
        AUDIT_GRID["scw_strength"],
        AUDIT_GRID["k_competition"],
        AUDIT_GRID["melanin_efficiency"],
        AUDIT_GRID["late_retention_factor"],
    ):
        row = evaluate_candidate_proxy(
            params=params,
            mat_activation_dpa=float(values[0]),
            scw_activation_dpa=float(values[1]),
            mat_strength=float(values[2]),
            scw_strength=float(values[3]),
            k_competition=float(values[4]),
            melanin_efficiency=float(values[5]),
            late_retention_factor=float(values[6]),
        )
        if float(row["temporal_gap_days"]) < float(min_temporal_gap_days):
            continue
        if not bool(row.get("toxicity_gate_pass", True)):
            continue
        margin_score, strict_pass = strict_margin_score(row, strict_thresholds_cfg)
        row["_strict_margin_score"] = float(margin_score)
        row["_strict_pass"] = bool(strict_pass)
        out.append(row)
    return out


def transit_status_for_variant(
    params: dict[str, Any],
    templates: dict[str, GeneticElement],
) -> dict[str, Any]:
    try:
        melA, tyrp1, dct = maybe_apply_transit_peptides(
            params,
            _clone_gene(templates["melA"]),
            _clone_gene(templates["TYRP1"]),
            _clone_gene(templates["DCT"]),
        )
        names = [melA.name, tyrp1.name, dct.name]
        return {
            "enabled": bool(params.get("construct", {}).get("use_vacuolar_transit_peptides", False)),
            "fused": all(str(name).endswith("_vTP") for name in names),
            "gene_names": names,
            "gene_lengths": [int(melA.length), int(tyrp1.length), int(dct.length)],
            "error": "",
        }
    except Exception as exc:
        return {
            "enabled": bool(params.get("construct", {}).get("use_vacuolar_transit_peptides", False)),
            "fused": False,
            "gene_names": [],
            "gene_lengths": [],
            "error": str(exc),
        }


def evaluate_variant(
    params: dict[str, Any],
    variant: dict[str, Any],
    strict_cfg: dict[str, Any],
    holdouts: list[dict[str, Any]],
    templates: dict[str, GeneticElement],
    n_trials_strict: int,
    n_trials_holdout: int,
    seed_candidates_limit: int,
    min_temporal_gap_days: float,
    seed: int,
) -> dict[str, Any]:
    deterministic_pool = build_deterministic_pool_from_grid(
        params=params,
        strict_thresholds_cfg=dict(strict_cfg["thresholds"]),
        min_temporal_gap_days=float(min_temporal_gap_days),
    )
    seeds = select_seed_candidates(deterministic_pool, limit=int(seed_candidates_limit))

    strict_rows = run_robustness_analysis(
        params=params,
        candidates=seeds,
        n_trials=int(n_trials_strict),
        seed=int(seed),
        thresholds=dict(strict_cfg["thresholds"]),
        noise=dict(strict_cfg["noise"]),
        collect_trial_records=False,
    )
    strict_top = strict_rows[0] if strict_rows else {}

    top3_rows = pick_top_unique_robust_rows(
        strict_rows,
        n=3,
        min_success_rate=0.20,
        max_fragility_index=0.40,
    )
    if not top3_rows:
        top3_rows = strict_rows[:3]
    top3_candidates = [dict(row["base_candidate"]) for row in top3_rows]

    holdout_rows: list[dict[str, Any]] = []
    for i, holdout in enumerate(holdouts):
        rows = run_robustness_analysis(
            params=params,
            candidates=top3_candidates,
            n_trials=int(n_trials_holdout),
            seed=int(seed) + (i + 1) * 101,
            thresholds=dict(holdout["thresholds"]),
            noise=dict(holdout["noise"]),
            collect_trial_records=False,
        )
        top = rows[0] if rows else {}
        success_values = [float(r.get("success_rate", 0.0)) for r in rows]
        holdout_rows.append(
            {
                "name": str(holdout["name"]),
                "description": str(holdout["description"]),
                "top_success_rate": float(top.get("success_rate", 0.0)) if top else 0.0,
                "top_robust_score": float(top.get("robust_score", 0.0)) if top else 0.0,
                "mean_success_rate": float(np.mean(success_values)) if success_values else 0.0,
                "min_success_rate": float(min(success_values)) if success_values else 0.0,
            }
        )

    holdout_top_success = [float(r["top_success_rate"]) for r in holdout_rows]
    holdout_min_top_success = min(holdout_top_success) if holdout_top_success else 0.0
    holdout_mean_top_success = float(np.mean(holdout_top_success)) if holdout_top_success else 0.0

    construct_status = transit_status_for_variant(params, templates)

    top_base = dict(strict_top.get("base_candidate", {})) if strict_top else {}
    return {
        "variant_id": str(variant["id"]),
        "label": str(variant["label"]),
        "switches": {
            "traffic": bool(variant["traffic"]),
            "storage": bool(variant["storage"]),
            "transit": bool(variant["transit"]),
        },
        "strict_pool_size": int(len(deterministic_pool)),
        "strict_seed_candidates": int(len(seeds)),
        "strict_top_success_rate": float(strict_top.get("success_rate", 0.0)) if strict_top else 0.0,
        "strict_top_robust_score": float(strict_top.get("robust_score", 0.0)) if strict_top else 0.0,
        "strict_top_fragility": float(strict_top.get("fragility_index", 1.0)) if strict_top else 1.0,
        "strict_top_candidate": top_base,
        "strict_top3": [
            {
                "success_rate": float(r.get("success_rate", 0.0)),
                "robust_score": float(r.get("robust_score", 0.0)),
                "fragility_index": float(r.get("fragility_index", 0.0)),
                "base_candidate": dict(r.get("base_candidate", {})),
            }
            for r in top3_rows
        ],
        "holdout": holdout_rows,
        "holdout_min_top_success_rate": float(holdout_min_top_success),
        "holdout_mean_top_success_rate": float(holdout_mean_top_success),
        "construct_transit_status": construct_status,
    }


def confidence_summary(variant_rows: list[dict[str, Any]]) -> dict[str, Any]:
    by_id = {str(row["variant_id"]): row for row in variant_rows}

    base = by_id.get("all_off", {})
    transit_only = by_id.get("transit_only", {})
    traffic_only = by_id.get("traffic_only", {})
    storage_only = by_id.get("storage_only", {})
    all_on = by_id.get("all_on", {})

    strict_gain = float(all_on.get("strict_top_success_rate", 0.0)) - float(base.get("strict_top_success_rate", 0.0))
    holdout_gain = float(all_on.get("holdout_min_top_success_rate", 0.0)) - float(
        base.get("holdout_min_top_success_rate", 0.0)
    )
    traffic_gain = float(traffic_only.get("strict_top_success_rate", 0.0)) - float(base.get("strict_top_success_rate", 0.0))
    storage_gain = float(storage_only.get("strict_top_success_rate", 0.0)) - float(
        transit_only.get("strict_top_success_rate", 0.0)
    )
    transit_gain = float(transit_only.get("strict_top_success_rate", 0.0)) - float(base.get("strict_top_success_rate", 0.0))

    rating = "LOW"
    if strict_gain >= 0.20 and holdout_gain >= 0.10:
        rating = "HIGH"
    elif strict_gain >= 0.10 and holdout_gain >= 0.05:
        rating = "MEDIUM"

    strong_points: list[str] = []
    caveats: list[str] = []

    if strict_gain > 0.0:
        strong_points.append(
            f"Combined upgrades improve strict-high-noise top success by {strict_gain:+.3f} versus control."
        )
    if holdout_gain > 0.0:
        strong_points.append(
            f"Combined upgrades improve worst holdout top success by {holdout_gain:+.3f} versus control."
        )
    if traffic_gain > 0.0:
        strong_points.append(f"Traffic tuning alone contributes positive gain ({traffic_gain:+.3f}).")
    if storage_gain > 0.0:
        strong_points.append(f"Storage targeting (with transit) contributes positive gain ({storage_gain:+.3f}).")

    if transit_gain <= 0.01:
        caveats.append(
            "Transit-only change has near-zero direct robustness impact unless paired with compartment/toxicity assumptions."
        )
    if float(all_on.get("strict_top_success_rate", 0.0)) < 0.12:
        caveats.append("Combined setup did not clear the 0.12 strict-high-noise minimum objective.")
    if float(all_on.get("holdout_min_top_success_rate", 0.0)) < 0.20:
        caveats.append("Holdout worst-case success is still low; model remains fragile under unseen stress.")

    return {
        "rating": rating,
        "strict_success_gain_all_on_vs_control": float(strict_gain),
        "holdout_worst_top_success_gain_all_on_vs_control": float(holdout_gain),
        "marginal_gain_traffic_only_vs_control": float(traffic_gain),
        "marginal_gain_storage_only_vs_transit_only": float(storage_gain),
        "marginal_gain_transit_only_vs_control": float(transit_gain),
        "strong_points": strong_points,
        "caveats": caveats,
    }


def write_outputs(
    out_prefix: str,
    n_trials_strict: int,
    n_trials_holdout: int,
    min_temporal_gap_days: float,
    seed_candidates_limit: int,
    seed: int,
    strict_cfg: dict[str, Any],
    holdouts: list[dict[str, Any]],
    variants: list[dict[str, Any]],
    summary: dict[str, Any],
) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    variants_path = RESULTS_DIR / f"{out_prefix}_variants.json"
    summary_path = RESULTS_DIR / f"{out_prefix}_summary.json"
    report_path = RESULTS_DIR / f"{out_prefix}_report.md"

    variants_path.write_text(json.dumps(variants, indent=2) + "\n")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    lines: list[str] = []
    lines.append("# Anti-Gibberish Audit")
    lines.append("")
    lines.append("Focused audit: ablation + unseen holdouts + confidence scoring.")
    lines.append("")
    lines.append("## Run Setup")
    lines.append("")
    lines.append(f"- Generated: `{summary['generated_at']}`")
    lines.append(f"- Strict trials/candidate: `{n_trials_strict}`")
    lines.append(f"- Holdout trials/candidate: `{n_trials_holdout}`")
    lines.append(f"- Seed candidates per variant: `{seed_candidates_limit}`")
    lines.append(f"- Min temporal gap filter: `{min_temporal_gap_days:.2f}` days")
    lines.append(f"- RNG seed: `{seed}`")
    lines.append(f"- Strict scenario: `{strict_cfg['name']}`")
    lines.append(f"- Holdouts: `{len(holdouts)}`")
    lines.append("")

    lines.append("## Ablation Results")
    lines.append("")
    lines.append("| Variant | Traffic | Storage | Transit | Strict Top Success | Strict Top Robust | Holdout Worst Top Success |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for row in variants:
        sw = row["switches"]
        lines.append(
            f"| {row['label']} | {int(bool(sw['traffic']))} | {int(bool(sw['storage']))} | {int(bool(sw['transit']))} | "
            f"{row['strict_top_success_rate']:.3f} | {row['strict_top_robust_score']:.3f} | {row['holdout_min_top_success_rate']:.3f} |"
        )
    lines.append("")

    lines.append("## Holdout Leaderboard (Top Success by Variant)")
    lines.append("")
    for holdout in holdouts:
        lines.append(f"### {holdout['name']}")
        lines.append("")
        lines.append("| Variant | Top Success | Top Robust |")
        lines.append("|---|---:|---:|")
        sorted_rows = sorted(
            variants,
            key=lambda r: next(
                (
                    float(h["top_success_rate"])
                    for h in r["holdout"]
                    if str(h["name"]) == str(holdout["name"])
                ),
                0.0,
            ),
            reverse=True,
        )
        for row in sorted_rows:
            h = next((x for x in row["holdout"] if str(x["name"]) == str(holdout["name"])), {})
            lines.append(
                f"| {row['label']} | {float(h.get('top_success_rate', 0.0)):.3f} | "
                f"{float(h.get('top_robust_score', 0.0)):.3f} |"
            )
        lines.append("")

    lines.append("## Confidence")
    lines.append("")
    lines.append(f"- Confidence rating: `{summary['confidence']['rating']}`")
    lines.append(
        "- Combined strict gain vs control: "
        f"`{summary['confidence']['strict_success_gain_all_on_vs_control']:+.3f}`"
    )
    lines.append(
        "- Combined holdout worst-top gain vs control: "
        f"`{summary['confidence']['holdout_worst_top_success_gain_all_on_vs_control']:+.3f}`"
    )
    lines.append(
        "- Marginal traffic gain: "
        f"`{summary['confidence']['marginal_gain_traffic_only_vs_control']:+.3f}`"
    )
    lines.append(
        "- Marginal storage gain (vs transit-only): "
        f"`{summary['confidence']['marginal_gain_storage_only_vs_transit_only']:+.3f}`"
    )
    lines.append(
        "- Marginal transit-only gain: "
        f"`{summary['confidence']['marginal_gain_transit_only_vs_control']:+.3f}`"
    )
    lines.append("")

    lines.append("### Strong Evidence")
    lines.append("")
    if summary["confidence"]["strong_points"]:
        for point in summary["confidence"]["strong_points"]:
            lines.append(f"- {point}")
    else:
        lines.append("- No strong evidence criteria were met.")
    lines.append("")

    lines.append("### Caveats")
    lines.append("")
    if summary["confidence"]["caveats"]:
        for caveat in summary["confidence"]["caveats"]:
            lines.append(f"- {caveat}")
    else:
        lines.append("- No major caveats flagged by this audit run.")
    lines.append("")

    report_path.write_text("\n".join(lines) + "\n")


def run(
    n_trials_strict: int,
    n_trials_holdout: int,
    seed_candidates_limit: int,
    min_temporal_gap_days: float,
    seed: int,
    out_prefix: str,
) -> None:
    base_params = load_params()
    strict_cfg = strict_high_noise_scenario()
    holdouts = holdout_scenarios()
    templates = build_gene_templates()

    results: list[dict[str, Any]] = []
    for idx, variant in enumerate(VARIANTS):
        params_variant = apply_variant(base_params, variant)
        row = evaluate_variant(
            params=params_variant,
            variant=variant,
            strict_cfg=strict_cfg,
            holdouts=holdouts,
            templates=templates,
            n_trials_strict=int(n_trials_strict),
            n_trials_holdout=int(n_trials_holdout),
            seed_candidates_limit=int(seed_candidates_limit),
            min_temporal_gap_days=float(min_temporal_gap_days),
            seed=int(seed) + idx * 701,
        )
        results.append(row)

    conf = confidence_summary(results)
    summary = {
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "out_prefix": out_prefix,
        "n_variants": len(results),
        "variants": [r["variant_id"] for r in results],
        "strict_scenario": strict_cfg,
        "holdout_scenarios": holdouts,
        "confidence": conf,
    }

    write_outputs(
        out_prefix=out_prefix,
        n_trials_strict=n_trials_strict,
        n_trials_holdout=n_trials_holdout,
        min_temporal_gap_days=min_temporal_gap_days,
        seed_candidates_limit=seed_candidates_limit,
        seed=seed,
        strict_cfg=strict_cfg,
        holdouts=holdouts,
        variants=results,
        summary=summary,
    )

    best_row = max(results, key=lambda r: float(r.get("strict_top_success_rate", 0.0))) if results else {}
    print("\n🧪 Anti-Gibberish Audit")
    print("=" * 34)
    print(f"Variants evaluated: {len(results)}")
    print(
        "Best strict-high-noise top success: "
        f"{best_row.get('label', 'n/a')} = {float(best_row.get('strict_top_success_rate', 0.0)):.3f}"
    )
    print(f"Confidence rating: {conf['rating']}")
    print("Saved:")
    print(f"  - results/{out_prefix}_variants.json")
    print(f"  - results/{out_prefix}_summary.json")
    print(f"  - results/{out_prefix}_report.md")
    print("=" * 34)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ablation + holdout anti-gibberish audit.")
    parser.add_argument("--n-trials-strict", type=int, default=120, help="Robustness trials per candidate in strict scenario.")
    parser.add_argument("--n-trials-holdout", type=int, default=140, help="Robustness trials per candidate in holdout scenarios.")
    parser.add_argument("--seed-candidates", type=int, default=120, help="Number of seed candidates per variant.")
    parser.add_argument("--min-temporal-gap-days", type=float, default=0.5, help="Deterministic minimum temporal gap filter.")
    parser.add_argument("--seed", type=int, default=20260304, help="RNG seed base.")
    parser.add_argument(
        "--out-prefix",
        default="anti_gibberish_audit_2026_03_04",
        help="Output prefix under results/.",
    )
    args = parser.parse_args()

    run(
        n_trials_strict=max(int(args.n_trials_strict), 60),
        n_trials_holdout=max(int(args.n_trials_holdout), 80),
        seed_candidates_limit=max(int(args.seed_candidates), 40),
        min_temporal_gap_days=float(args.min_temporal_gap_days),
        seed=int(args.seed),
        out_prefix=str(args.out_prefix),
    )
