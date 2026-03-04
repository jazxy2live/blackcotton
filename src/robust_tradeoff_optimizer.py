#!/usr/bin/env python3
"""
robust_tradeoff_optimizer.py — Robust candidate search under uncertainty
========================================================================

Builds on deterministic tradeoff optimization and robustness analysis:
  1) Generate safe coarse candidates (temporal_gap_days >= 0)
  2) Select a diverse seed subset
  3) Run Monte Carlo uncertainty analysis on the seed subset
  4) Rank by robust score / success rate and export robust shortlist

Usage:
    python -m src.robust_tradeoff_optimizer
"""

import json
from pathlib import Path
from typing import Any

from src.robustness_analyzer import (
    DEFAULT_NOISE,
    DEFAULT_THRESHOLDS,
    run_robustness_analysis,
)
from src.tradeoff_optimizer import (
    load_params,
    run_sweep,
    select_seed_candidates,
)

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"


def select_robust_top(
    robust_results: list[dict[str, Any]],
    n: int = 12,
    min_success_rate: float = 0.05,
) -> list[dict[str, Any]]:
    """
    Choose top robust candidates while keeping diversity.
    """
    pool = [r for r in robust_results if float(r["success_rate"]) >= float(min_success_rate)]
    if not pool:
        pool = robust_results

    selected: list[dict[str, Any]] = []
    seen = set()
    for row in pool:
        base = row["base_candidate"]
        key = (
            round(float(base["color_L"]), 2),
            round(float(base["strength_g_tex"]), 2),
            round(float(base["yield_index"]), 3),
            round(float(base["temporal_gap_days"]), 2),
            round(float(base["mat_activation_dpa"]), 2),
            round(float(base["scw_activation_dpa"]), 2),
            round(float(base["k_competition"]), 2),
        )
        if key in seen:
            continue
        seen.add(key)
        selected.append(row)
        if len(selected) >= n:
            break
    return selected


def save_outputs(
    robust_all: list[dict[str, Any]],
    robust_top: list[dict[str, Any]],
    summary: dict[str, Any],
) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    stripped_candidates = []
    trial_records = []
    for row in robust_all:
        item = dict(row)
        trials = item.pop("trial_records", [])
        if trials:
            trial_records.extend(trials)
        stripped_candidates.append(item)

    with open(RESULTS_DIR / "robust_optimization_candidates.json", "w") as f:
        json.dump(stripped_candidates, f, indent=2)
    with open(RESULTS_DIR / "robust_top_candidates.json", "w") as f:
        json.dump(robust_top, f, indent=2)
    with open(RESULTS_DIR / "robust_optimization_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    with open(RESULTS_DIR / "robust_trial_records.jsonl", "w") as f:
        for record in trial_records:
            f.write(json.dumps(record) + "\n")
    write_markdown_report(robust_top=robust_top, summary=summary)


def write_markdown_report(
    robust_top: list[dict[str, Any]],
    summary: dict[str, Any],
) -> None:
    lines = []
    lines.append("# Robust Optimization Report")
    lines.append("")
    lines.append("Uncertainty-aware ranking of black cotton candidates.")
    lines.append("")
    lines.append("## Run Summary")
    lines.append("")
    lines.append(f"- Coarse safe candidates: `{summary['num_coarse_candidates']}`")
    lines.append(f"- Robust seeds evaluated: `{summary['num_seed_candidates']}`")
    lines.append(f"- Trials per seed: `{summary['n_trials_per_seed']}`")
    lines.append(f"- Best robust score: `{summary['best_robust_score']:.3f}`")
    lines.append(f"- Best success rate: `{summary['best_success_rate']:.3f}`")
    lines.append("")
    lines.append("## Top 12 (Robust)")
    lines.append("")
    lines.append("| Rank | Success | Robust Score | p50 L* | p50 Strength | p50 Yield | Risk Gap | Risk Darkness | Params (mat/scw, k, eff, ret) |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for i, row in enumerate(robust_top, 1):
        base = row["base_candidate"]
        lines.append(
            f"| {i} | {row['success_rate']:.3f} | {row['robust_score']:.3f} | "
            f"{row['p50_color_L']:.2f} | {row['p50_strength_g_tex']:.2f} | {row['p50_yield_index']:.3f} | "
            f"{row['risk_temporal_overlap']:.3f} | {row['risk_darkness_failure']:.3f} | "
            f"{base['mat_activation_dpa']:.0f}/{base['scw_activation_dpa']:.0f}, "
            f"{base['k_competition']:.2f}, {base['melanin_efficiency']:.2f}, {base['late_retention_factor']:.2f} |"
        )
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("- Robust winners are less dark in median conditions but survive uncertainty better.")
    lines.append("- Main failure mode is still darkness variability; strength/yield are comparatively stable.")
    lines.append("- Next optimization cycles should focus on reducing darkness-risk without creating timing overlap.")
    lines.append("")
    with open(RESULTS_DIR / "robust_optimization_report.md", "w") as f:
        f.write("\n".join(lines) + "\n")


def print_report(
    robust_top: list[dict[str, Any]],
    summary: dict[str, Any],
) -> None:
    print("\n" + "=" * 112)
    print("  BLACKCOTTON ROBUST TRADEOFF OPTIMIZER")
    print("=" * 112)
    print(f"\n  Coarse safe candidates:      {summary['num_coarse_candidates']}")
    print(f"  Robust seeds evaluated:      {summary['num_seed_candidates']}")
    print(f"  Trials per seed:             {summary['n_trials_per_seed']}")
    print(f"  Minimum success for shortlist: {summary['min_success_rate_for_shortlist']:.3f}")
    print(f"  Best robust score:           {summary['best_robust_score']:.3f}")
    print(f"  Best success rate:           {summary['best_success_rate']:.3f}")
    print("\n  ROBUST TOP CANDIDATES")
    print("  " + "-" * 110)
    print("  #  succ   robust   p50 L*   p50 Str   p50 Yield   risk_gap  risk_dark  mat/scw  k   eff  ret")
    print("  " + "-" * 110)
    for i, row in enumerate(robust_top, 1):
        base = row["base_candidate"]
        print(
            f"  {i:>2} {row['success_rate']:>6.3f} {row['robust_score']:>7.3f} "
            f"{row['p50_color_L']:>7.2f} {row['p50_strength_g_tex']:>9.2f} {row['p50_yield_index']:>10.3f} "
            f"{row['risk_temporal_overlap']:>9.3f} {row['risk_darkness_failure']:>10.3f} "
            f"{base['mat_activation_dpa']:.0f}/{base['scw_activation_dpa']:.0f} "
            f"{base['k_competition']:.2f} {base['melanin_efficiency']:.2f} {base['late_retention_factor']:.2f}"
        )
    print("  " + "-" * 110)
    print("\n  Saved:")
    print("    - results/robust_optimization_candidates.json")
    print("    - results/robust_top_candidates.json")
    print("    - results/robust_optimization_summary.json")
    print("    - results/robust_trial_records.jsonl")
    print("    - results/robust_optimization_report.md")
    print("\n" + "=" * 112)


if __name__ == "__main__":
    print("\n🧬 BlackCotton Robust Tradeoff Optimizer")
    print("=" * 50)

    params = load_params()

    # Stage 1: deterministic safe candidate generation.
    coarse_candidates = run_sweep(params)

    # Stage 2: robust evaluation on a diverse subset.
    seed_candidates = select_seed_candidates(coarse_candidates, n=160)

    thresholds = dict(DEFAULT_THRESHOLDS)
    noise = dict(DEFAULT_NOISE)
    n_trials = 100
    seed = 42
    min_success_rate_for_shortlist = 0.08

    robust_all = run_robustness_analysis(
        params=params,
        candidates=seed_candidates,
        n_trials=n_trials,
        seed=seed,
        thresholds=thresholds,
        noise=noise,
        collect_trial_records=True,
    )
    robust_top = select_robust_top(
        robust_all,
        n=12,
        min_success_rate=min_success_rate_for_shortlist,
    )

    summary = {
        "num_coarse_candidates": len(coarse_candidates),
        "num_seed_candidates": len(seed_candidates),
        "n_trials_per_seed": n_trials,
        "seed": seed,
        "thresholds": thresholds,
        "noise_model": noise,
        "min_success_rate_for_shortlist": min_success_rate_for_shortlist,
        "num_shortlisted": len(robust_top),
        "num_trial_records": sum(len(r.get("trial_records", [])) for r in robust_all),
        "best_robust_score": float(robust_top[0]["robust_score"]) if robust_top else 0.0,
        "best_success_rate": float(robust_top[0]["success_rate"]) if robust_top else 0.0,
        "best_shortlisted_input_rank": int(robust_top[0]["input_rank"]) if robust_top else None,
    }

    save_outputs(robust_all=robust_all, robust_top=robust_top, summary=summary)
    print_report(robust_top=robust_top, summary=summary)
    print("\n✅ Robust tradeoff optimization complete!")
