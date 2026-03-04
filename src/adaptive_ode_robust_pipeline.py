#!/usr/bin/env python3
"""
adaptive_ode_robust_pipeline.py — Finalize lab shortlist with ODE + robustness
===============================================================================

Pipeline:
  1) Load adaptive robust candidates (proxy stage shortlist)
  2) ODE-refine them with full expression + melanin dynamics
  3) Re-run robustness on ODE-refined candidates
  4) Freeze a final lab Top 3

Usage:
    python -m src.adaptive_ode_robust_pipeline
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
    pareto_front,
    refine_candidates_with_odes,
    select_top_candidates,
)

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"


def normalize_candidate_row(row: dict[str, Any]) -> dict[str, float]:
    """Normalize candidate row from either robust or base-candidate format."""
    base = row.get("base_candidate", row)
    return {
        "mat_activation_dpa": float(base["mat_activation_dpa"]),
        "scw_activation_dpa": float(base["scw_activation_dpa"]),
        "mat_strength": float(base["mat_strength"]),
        "scw_strength": float(base["scw_strength"]),
        "k_competition": float(base["k_competition"]),
        "melanin_efficiency": float(base["melanin_efficiency"]),
        "late_retention_factor": float(base["late_retention_factor"]),
        "color_L": float(base.get("color_L", 82.0)),
        "strength_g_tex": float(base.get("strength_g_tex", 0.0)),
        "yield_index": float(base.get("yield_index", 0.0)),
        "temporal_gap_days": float(base.get("temporal_gap_days", 0.0)),
        "composite_score": float(row.get("robust_score", row.get("composite_score", 0.0))),
    }


def dedupe_param_candidates(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    seen = set()
    out = []
    for c in candidates:
        key = (
            round(float(c["mat_activation_dpa"]), 4),
            round(float(c["scw_activation_dpa"]), 4),
            round(float(c["mat_strength"]), 4),
            round(float(c["scw_strength"]), 4),
            round(float(c["k_competition"]), 4),
            round(float(c["melanin_efficiency"]), 4),
            round(float(c["late_retention_factor"]), 4),
        )
        if key in seen:
            continue
        seen.add(key)
        out.append(c)
    return out


def select_final_top_candidates(
    robust_rows: list[dict[str, Any]],
    n: int = 3,
    min_success_rate: float = 0.20,
) -> list[dict[str, Any]]:
    pool = [r for r in robust_rows if float(r["success_rate"]) >= min_success_rate]
    if not pool:
        pool = robust_rows

    selected = []
    seen = set()
    for row in pool:
        b = row["base_candidate"]
        key = (
            round(float(b["color_L"]), 2),
            round(float(b["strength_g_tex"]), 2),
            round(float(b["yield_index"]), 3),
            round(float(b["temporal_gap_days"]), 2),
        )
        if key in seen:
            continue
        seen.add(key)
        selected.append(row)
        if len(selected) >= n:
            break
    return selected


def save_outputs(
    ode_refined: list[dict[str, Any]],
    ode_pareto: list[dict[str, Any]],
    ode_top: list[dict[str, Any]],
    robust_rows: list[dict[str, Any]],
    final_top3: list[dict[str, Any]],
    summary: dict[str, Any],
) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    stripped_robust = []
    trial_records = []
    for row in robust_rows:
        item = dict(row)
        trial_records.extend(item.pop("trial_records", []))
        stripped_robust.append(item)
    final_top3_stripped = []
    for row in final_top3:
        item = dict(row)
        item.pop("trial_records", None)
        final_top3_stripped.append(item)

    with open(RESULTS_DIR / "adaptive_ode_refined_candidates.json", "w") as f:
        json.dump(ode_refined, f, indent=2)
    with open(RESULTS_DIR / "adaptive_ode_pareto.json", "w") as f:
        json.dump(ode_pareto, f, indent=2)
    with open(RESULTS_DIR / "adaptive_ode_top_candidates.json", "w") as f:
        json.dump(ode_top, f, indent=2)
    with open(RESULTS_DIR / "adaptive_ode_robust_candidates.json", "w") as f:
        json.dump(stripped_robust, f, indent=2)
    with open(RESULTS_DIR / "adaptive_ode_robust_top3.json", "w") as f:
        json.dump(final_top3_stripped, f, indent=2)
    with open(RESULTS_DIR / "final_lab_top3.json", "w") as f:
        json.dump(final_top3_stripped, f, indent=2)
    with open(RESULTS_DIR / "adaptive_ode_robust_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    with open(RESULTS_DIR / "adaptive_ode_robust_trial_records.jsonl", "w") as f:
        for record in trial_records:
            f.write(json.dumps(record) + "\n")

    lines = []
    lines.append("# Adaptive ODE Robust Final Report")
    lines.append("")
    lines.append("Finalization pipeline: adaptive robust shortlist -> ODE refinement -> robustness -> lab Top 3.")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Input adaptive candidates: `{summary['num_input_adaptive_candidates']}`")
    lines.append(f"- ODE refined candidates: `{summary['num_ode_refined_candidates']}`")
    lines.append(f"- ODE pareto candidates: `{summary['num_ode_pareto']}`")
    lines.append(f"- ODE top candidates for robustness: `{summary['num_ode_top_for_robustness']}`")
    lines.append(f"- Robust trials per candidate: `{summary['n_trials_per_candidate']}`")
    lines.append(f"- Best robust success rate: `{summary['best_success_rate']:.3f}`")
    lines.append(f"- Best robust score: `{summary['best_robust_score']:.3f}`")
    lines.append("")
    lines.append("## Final Lab Top 3")
    lines.append("")
    lines.append("| Rank | Success | Robust | p50 L* | p50 Str | p50 Yield | Risk Gap | Risk Dark | Params |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for i, row in enumerate(final_top3, 1):
        b = row["base_candidate"]
        lines.append(
            f"| {i} | {row['success_rate']:.3f} | {row['robust_score']:.3f} | "
            f"{row['p50_color_L']:.2f} | {row['p50_strength_g_tex']:.2f} | {row['p50_yield_index']:.3f} | "
            f"{row['risk_temporal_overlap']:.3f} | {row['risk_darkness_failure']:.3f} | "
            f"{b['mat_activation_dpa']:.0f}/{b['scw_activation_dpa']:.0f}, "
            f"k={b['k_competition']:.2f}, eff={b['melanin_efficiency']:.2f}, ret={b['late_retention_factor']:.2f} |"
        )
    lines.append("")

    with open(RESULTS_DIR / "adaptive_ode_robust_report.md", "w") as f:
        f.write("\n".join(lines) + "\n")


def print_report(summary: dict[str, Any], final_top3: list[dict[str, Any]]) -> None:
    print("\n" + "=" * 116)
    print("  ADAPTIVE ODE ROBUST FINALIZATION")
    print("=" * 116)
    print(f"\n  Input adaptive candidates:     {summary['num_input_adaptive_candidates']}")
    print(f"  ODE refined candidates:        {summary['num_ode_refined_candidates']}")
    print(f"  ODE pareto candidates:         {summary['num_ode_pareto']}")
    print(f"  ODE top for robustness:        {summary['num_ode_top_for_robustness']}")
    print(f"  Robust trials per candidate:   {summary['n_trials_per_candidate']}")
    print(f"  Best robust success rate:      {summary['best_success_rate']:.3f}")
    print(f"  Best robust score:             {summary['best_robust_score']:.3f}")
    print("\n  FINAL LAB TOP 3")
    print("  " + "-" * 114)
    print("  #  succ   robust   p50 L*   p50 Str   p50 Yield   risk_gap  risk_dark  mat/scw  k   eff  ret")
    print("  " + "-" * 114)
    for i, row in enumerate(final_top3, 1):
        b = row["base_candidate"]
        print(
            f"  {i:>2} {row['success_rate']:>6.3f} {row['robust_score']:>7.3f} "
            f"{row['p50_color_L']:>7.2f} {row['p50_strength_g_tex']:>9.2f} {row['p50_yield_index']:>10.3f} "
            f"{row['risk_temporal_overlap']:>9.3f} {row['risk_darkness_failure']:>10.3f} "
            f"{b['mat_activation_dpa']:.0f}/{b['scw_activation_dpa']:.0f} "
            f"{b['k_competition']:.2f} {b['melanin_efficiency']:.2f} {b['late_retention_factor']:.2f}"
        )
    print("  " + "-" * 114)
    print("\n  Saved:")
    print("    - results/adaptive_ode_refined_candidates.json")
    print("    - results/adaptive_ode_pareto.json")
    print("    - results/adaptive_ode_top_candidates.json")
    print("    - results/adaptive_ode_robust_candidates.json")
    print("    - results/adaptive_ode_robust_top3.json")
    print("    - results/final_lab_top3.json")
    print("    - results/adaptive_ode_robust_summary.json")
    print("    - results/adaptive_ode_robust_trial_records.jsonl")
    print("    - results/adaptive_ode_robust_report.md")
    print("\n" + "=" * 116)


if __name__ == "__main__":
    print("\n🧬 BlackCotton Adaptive ODE Robust Pipeline")
    print("=" * 50)

    input_path = RESULTS_DIR / "adaptive_robust_top_candidates.json"
    if not input_path.exists():
        raise FileNotFoundError(
            "Missing results/adaptive_robust_top_candidates.json. "
            "Run `python -m src.adaptive_robust_optimizer` first."
        )

    with open(input_path) as f:
        rows = json.load(f)

    input_candidates = dedupe_param_candidates([normalize_candidate_row(r) for r in rows])
    params = load_params()

    ode_refined = refine_candidates_with_odes(
        params=params,
        coarse_candidates=input_candidates,
        max_refine=len(input_candidates),
    )
    ode_pareto = pareto_front(ode_refined)
    ode_top = select_top_candidates(ode_refined, n=12)

    robust_rows = run_robustness_analysis(
        params=params,
        candidates=ode_top,
        n_trials=120,
        seed=777,
        thresholds=dict(DEFAULT_THRESHOLDS),
        noise=dict(DEFAULT_NOISE),
        collect_trial_records=True,
    )
    final_top3 = select_final_top_candidates(robust_rows, n=3, min_success_rate=0.20)

    summary = {
        "num_input_adaptive_candidates": len(input_candidates),
        "num_ode_refined_candidates": len(ode_refined),
        "num_ode_pareto": len(ode_pareto),
        "num_ode_top_for_robustness": len(ode_top),
        "n_trials_per_candidate": 120,
        "seed": 777,
        "thresholds": dict(DEFAULT_THRESHOLDS),
        "noise_model": dict(DEFAULT_NOISE),
        "num_robust_candidates": len(robust_rows),
        "num_final_top3": len(final_top3),
        "best_success_rate": float(final_top3[0]["success_rate"]) if final_top3 else 0.0,
        "best_robust_score": float(final_top3[0]["robust_score"]) if final_top3 else 0.0,
    }

    save_outputs(
        ode_refined=ode_refined,
        ode_pareto=ode_pareto,
        ode_top=ode_top,
        robust_rows=robust_rows,
        final_top3=final_top3,
        summary=summary,
    )
    print_report(summary=summary, final_top3=final_top3)
    print("\n✅ Adaptive ODE robust finalization complete!")
