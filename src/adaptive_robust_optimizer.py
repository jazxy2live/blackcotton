#!/usr/bin/env python3
"""
adaptive_robust_optimizer.py — Sensitivity-guided robust optimization
======================================================================

Uses sensitivity insights to retune the search space and increase robust success:
  - Focuses on ranges that improve darkness while preserving safe timing/quality
  - Runs uncertainty-aware scoring on targeted candidates
  - Exports a lab-priority shortlist

Usage:
    python -m src.adaptive_robust_optimizer
"""

import json
from itertools import product
from pathlib import Path
from typing import Any

from src.robustness_analyzer import (
    DEFAULT_NOISE,
    DEFAULT_THRESHOLDS,
    run_robustness_analysis,
)
from src.tradeoff_optimizer import evaluate_candidate_proxy, load_params

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"


def targeted_parameter_grid() -> dict[str, list[float]]:
    """
    Sensitivity-guided ranges:
      - push melanin efficiency and retention upward (darkness bottleneck)
      - keep competition moderate to protect strength/yield
      - keep timing in safe-late window
    """
    return {
        "mat_activation_dpa": [38.0, 40.0, 42.0],
        "scw_activation_dpa": [30.0, 32.0, 34.0],
        "mat_strength": [0.9, 1.1, 1.3],
        "scw_strength": [0.9, 1.1, 1.3],
        "k_competition": [0.05, 0.10, 0.15, 0.20],
        "melanin_efficiency": [1.30, 1.45, 1.60],
        "late_retention_factor": [0.45, 0.60, 0.75],
    }


def run_targeted_sweep(params: dict[str, Any]) -> list[dict[str, Any]]:
    grid = targeted_parameter_grid()
    candidates: list[dict[str, Any]] = []

    for (
        mat_act,
        scw_act,
        mat_str,
        scw_str,
        k_comp,
        mel_eff,
        late_ret,
    ) in product(
        grid["mat_activation_dpa"],
        grid["scw_activation_dpa"],
        grid["mat_strength"],
        grid["scw_strength"],
        grid["k_competition"],
        grid["melanin_efficiency"],
        grid["late_retention_factor"],
    ):
        row = evaluate_candidate_proxy(
            params=params,
            mat_activation_dpa=mat_act,
            scw_activation_dpa=scw_act,
            mat_strength=mat_str,
            scw_strength=scw_str,
            k_competition=k_comp,
            melanin_efficiency=mel_eff,
            late_retention_factor=late_ret,
        )
        if float(row["temporal_gap_days"]) < 0.0:
            continue
        candidates.append(row)

    # Keep best deterministic subset for robust analysis.
    candidates.sort(key=lambda c: c["composite_score"], reverse=True)
    return candidates[:260]


def pick_shortlist(
    robust_results: list[dict[str, Any]],
    n: int = 12,
    min_success: float = 0.15,
) -> list[dict[str, Any]]:
    pool = [r for r in robust_results if float(r["success_rate"]) >= float(min_success)]
    if not pool:
        pool = robust_results
    return pool[:n]


def save_outputs(
    robust_all: list[dict[str, Any]],
    shortlist: list[dict[str, Any]],
    summary: dict[str, Any],
) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    stripped = []
    trial_records = []
    for row in robust_all:
        item = dict(row)
        trial_records.extend(item.pop("trial_records", []))
        stripped.append(item)
    shortlist_stripped = []
    for row in shortlist:
        item = dict(row)
        item.pop("trial_records", None)
        shortlist_stripped.append(item)

    with open(RESULTS_DIR / "adaptive_robust_candidates.json", "w") as f:
        json.dump(stripped, f, indent=2)
    with open(RESULTS_DIR / "adaptive_robust_top_candidates.json", "w") as f:
        json.dump(shortlist_stripped, f, indent=2)
    with open(RESULTS_DIR / "adaptive_robust_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    with open(RESULTS_DIR / "adaptive_robust_trial_records.jsonl", "w") as f:
        for rec in trial_records:
            f.write(json.dumps(rec) + "\n")

    lines = []
    lines.append("# Adaptive Robust Optimization Report")
    lines.append("")
    lines.append("Targeted robust optimization after sensitivity-guided retuning.")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Deterministic targeted seeds: `{summary['num_targeted_seed_candidates']}`")
    lines.append(f"- Trials per seed: `{summary['n_trials_per_seed']}`")
    lines.append(f"- Best robust score: `{summary['best_robust_score']:.3f}`")
    lines.append(f"- Best success rate: `{summary['best_success_rate']:.3f}`")
    lines.append("")
    lines.append("## Top 12")
    lines.append("")
    lines.append("| Rank | Success | Robust | p50 L* | p50 Str | p50 Yield | Risk Gap | Risk Dark | Params |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for i, row in enumerate(shortlist, 1):
        b = row["base_candidate"]
        lines.append(
            f"| {i} | {row['success_rate']:.3f} | {row['robust_score']:.3f} | "
            f"{row['p50_color_L']:.2f} | {row['p50_strength_g_tex']:.2f} | {row['p50_yield_index']:.3f} | "
            f"{row['risk_temporal_overlap']:.3f} | {row['risk_darkness_failure']:.3f} | "
            f"{b['mat_activation_dpa']:.0f}/{b['scw_activation_dpa']:.0f}, "
            f"k={b['k_competition']:.2f}, eff={b['melanin_efficiency']:.2f}, ret={b['late_retention_factor']:.2f} |"
        )
    lines.append("")
    with open(RESULTS_DIR / "adaptive_robust_report.md", "w") as f:
        f.write("\n".join(lines) + "\n")


def print_report(shortlist: list[dict[str, Any]], summary: dict[str, Any]) -> None:
    print("\n" + "=" * 112)
    print("  BLACKCOTTON ADAPTIVE ROBUST OPTIMIZER")
    print("=" * 112)
    print(f"\n  Targeted deterministic seeds: {summary['num_targeted_seed_candidates']}")
    print(f"  Trials per seed:              {summary['n_trials_per_seed']}")
    print(f"  Best robust score:            {summary['best_robust_score']:.3f}")
    print(f"  Best success rate:            {summary['best_success_rate']:.3f}")
    print("\n  TOP ADAPTIVE ROBUST CANDIDATES")
    print("  " + "-" * 110)
    print("  #  succ   robust   p50 L*   p50 Str   p50 Yield   risk_gap  risk_dark  mat/scw  k   eff  ret")
    print("  " + "-" * 110)
    for i, row in enumerate(shortlist, 1):
        b = row["base_candidate"]
        print(
            f"  {i:>2} {row['success_rate']:>6.3f} {row['robust_score']:>7.3f} "
            f"{row['p50_color_L']:>7.2f} {row['p50_strength_g_tex']:>9.2f} {row['p50_yield_index']:>10.3f} "
            f"{row['risk_temporal_overlap']:>9.3f} {row['risk_darkness_failure']:>10.3f} "
            f"{b['mat_activation_dpa']:.0f}/{b['scw_activation_dpa']:.0f} "
            f"{b['k_competition']:.2f} {b['melanin_efficiency']:.2f} {b['late_retention_factor']:.2f}"
        )
    print("  " + "-" * 110)
    print("\n  Saved:")
    print("    - results/adaptive_robust_candidates.json")
    print("    - results/adaptive_robust_top_candidates.json")
    print("    - results/adaptive_robust_summary.json")
    print("    - results/adaptive_robust_trial_records.jsonl")
    print("    - results/adaptive_robust_report.md")
    print("\n" + "=" * 112)


if __name__ == "__main__":
    print("\n🧬 BlackCotton Adaptive Robust Optimizer")
    print("=" * 50)

    params = load_params()
    targeted_candidates = run_targeted_sweep(params)

    thresholds = dict(DEFAULT_THRESHOLDS)
    noise = dict(DEFAULT_NOISE)
    n_trials = 120
    seed = 123

    robust_all = run_robustness_analysis(
        params=params,
        candidates=targeted_candidates,
        n_trials=n_trials,
        seed=seed,
        thresholds=thresholds,
        noise=noise,
        collect_trial_records=True,
    )
    shortlist = pick_shortlist(robust_all, n=12, min_success=0.15)

    summary = {
        "num_targeted_seed_candidates": len(targeted_candidates),
        "n_trials_per_seed": n_trials,
        "seed": seed,
        "thresholds": thresholds,
        "noise_model": noise,
        "num_shortlisted": len(shortlist),
        "best_robust_score": float(shortlist[0]["robust_score"]) if shortlist else 0.0,
        "best_success_rate": float(shortlist[0]["success_rate"]) if shortlist else 0.0,
    }

    save_outputs(robust_all=robust_all, shortlist=shortlist, summary=summary)
    print_report(shortlist=shortlist, summary=summary)
    print("\n✅ Adaptive robust optimization complete!")
