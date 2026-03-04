#!/usr/bin/env python3
"""
sensitivity_analyzer.py — Parameter sensitivity from robust trial records
=========================================================================

Analyzes Monte Carlo trial records to quantify which parameters drive:
  - overall success/failure
  - darkness failures
  - temporal overlap failures

Usage:
    python -m src.sensitivity_analyzer
"""

import json
import argparse
from pathlib import Path
from typing import Any

import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"

FEATURE_COLUMNS = [
    "sampled_mat_activation_dpa",
    "sampled_scw_activation_dpa",
    "sampled_mat_strength",
    "sampled_scw_strength",
    "sampled_k_competition",
    "sampled_melanin_efficiency",
    "sampled_late_retention_factor",
    "sampled_hill_melA",
    "sampled_hill_scw",
    "sampled_leak_melA",
    "sampled_leak_scw",
    "sampled_copper_loading_fraction",
    "sampled_tyrosinase_activation_fraction",
    "sampled_ros_buffer_capacity",
    "sampled_silencing_probability",
    "sampled_event_expression_cv",
]


def _safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    if x.size == 0 or y.size == 0:
        return 0.0
    if np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def load_trial_records(path: Path) -> list[dict[str, Any]]:
    records = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            records.append(json.loads(line))
    return records


def recommend_direction(
    corr_success: float,
    corr_dark_fail: float,
    corr_gap_fail: float,
) -> str:
    # Prioritize reducing failure channels over raw success correlation.
    if abs(corr_dark_fail) >= abs(corr_gap_fail) and abs(corr_dark_fail) >= 0.05:
        return "decrease" if corr_dark_fail > 0 else "increase"
    if abs(corr_gap_fail) >= 0.05:
        return "decrease" if corr_gap_fail > 0 else "increase"
    if abs(corr_success) >= 0.03:
        return "increase" if corr_success > 0 else "decrease"
    return "neutral"


def rank_feature_impacts(
    records: list[dict[str, Any]],
    thresholds: dict[str, float],
) -> list[dict[str, Any]]:
    if not records:
        return []

    color = np.array([float(r["color_L"]) for r in records], dtype=float)
    strength = np.array([float(r["strength_g_tex"]) for r in records], dtype=float)
    yld = np.array([float(r["yield_index"]) for r in records], dtype=float)
    gap = np.array([float(r["temporal_gap_days"]) for r in records], dtype=float)
    success = np.array([1.0 if bool(r["success"]) else 0.0 for r in records], dtype=float)

    dark_fail = (color > float(thresholds["max_color_L"])).astype(float)
    gap_fail = (gap < float(thresholds["min_temporal_gap_days"])).astype(float)

    ranked = []
    for col in FEATURE_COLUMNS:
        x = np.array([float(r.get(col, 0.0)) for r in records], dtype=float)
        corr_success = _safe_corr(x, success)
        corr_dark_fail = _safe_corr(x, dark_fail)
        corr_gap_fail = _safe_corr(x, gap_fail)
        corr_color = _safe_corr(x, color)
        corr_strength = _safe_corr(x, strength)
        corr_yield = _safe_corr(x, yld)

        # Weight failure channels highest; they are the gating risks.
        impact_score = (
            1.4 * abs(corr_dark_fail)
            + 1.2 * abs(corr_gap_fail)
            + 1.0 * abs(corr_success)
            + 0.4 * abs(corr_strength)
            + 0.4 * abs(corr_yield)
        )

        ranked.append(
            {
                "feature": col,
                "impact_score": float(impact_score),
                "corr_success": float(corr_success),
                "corr_darkness_failure": float(corr_dark_fail),
                "corr_temporal_overlap_failure": float(corr_gap_fail),
                "corr_color_L": float(corr_color),
                "corr_strength_g_tex": float(corr_strength),
                "corr_yield_index": float(corr_yield),
                "recommended_direction": recommend_direction(corr_success, corr_dark_fail, corr_gap_fail),
            }
        )

    ranked.sort(key=lambda row: row["impact_score"], reverse=True)
    return ranked


def candidate_failure_table(
    records: list[dict[str, Any]],
    thresholds: dict[str, float],
) -> list[dict[str, Any]]:
    if not records:
        return []

    grouped: dict[int, list[dict[str, Any]]] = {}
    for r in records:
        grouped.setdefault(int(r["candidate_input_rank"]), []).append(r)

    table = []
    for rank, rows in grouped.items():
        arr_color = np.array([float(r["color_L"]) for r in rows], dtype=float)
        arr_strength = np.array([float(r["strength_g_tex"]) for r in rows], dtype=float)
        arr_yield = np.array([float(r["yield_index"]) for r in rows], dtype=float)
        arr_gap = np.array([float(r["temporal_gap_days"]) for r in rows], dtype=float)
        arr_success = np.array([1.0 if bool(r["success"]) else 0.0 for r in rows], dtype=float)

        table.append(
            {
                "candidate_input_rank": int(rank),
                "n_trials": int(len(rows)),
                "success_rate": float(arr_success.mean()),
                "risk_darkness_failure": float(np.mean(arr_color > float(thresholds["max_color_L"]))),
                "risk_temporal_overlap": float(np.mean(arr_gap < float(thresholds["min_temporal_gap_days"]))),
                "risk_strength_failure": float(np.mean(arr_strength < float(thresholds["min_strength_g_tex"]))),
                "risk_yield_failure": float(np.mean(arr_yield < float(thresholds["min_yield_index"]))),
                "p50_color_L": float(np.percentile(arr_color, 50)),
                "p50_strength_g_tex": float(np.percentile(arr_strength, 50)),
                "p50_yield_index": float(np.percentile(arr_yield, 50)),
                "p50_temporal_gap_days": float(np.percentile(arr_gap, 50)),
            }
        )

    table.sort(key=lambda row: row["success_rate"], reverse=True)
    return table


def save_outputs(
    summary: dict[str, Any],
    ranked_features: list[dict[str, Any]],
    candidate_table: list[dict[str, Any]],
    out_prefix: str = "sensitivity",
) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_DIR / f"{out_prefix}_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    with open(RESULTS_DIR / f"{out_prefix}_ranked_features.json", "w") as f:
        json.dump(ranked_features, f, indent=2)
    with open(RESULTS_DIR / f"{out_prefix}_candidate_failures.json", "w") as f:
        json.dump(candidate_table, f, indent=2)

    lines = []
    lines.append("# Sensitivity Analysis Report")
    lines.append("")
    lines.append("Global sensitivity derived from robust trial records.")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Trial records analyzed: `{summary['n_trial_records']}`")
    lines.append(f"- Candidate ranks analyzed: `{summary['n_candidate_ranks']}`")
    lines.append(f"- Overall success rate: `{summary['overall_success_rate']:.3f}`")
    lines.append(f"- Main failure mode: `{summary['primary_failure_mode']}`")
    lines.append("")
    lines.append("## Top Feature Impacts")
    lines.append("")
    lines.append("| Rank | Feature | Impact | Corr Success | Corr Dark Fail | Corr Gap Fail | Recommended Direction |")
    lines.append("|---:|---|---:|---:|---:|---:|---|")
    for i, row in enumerate(ranked_features[:10], 1):
        lines.append(
            f"| {i} | {row['feature']} | {row['impact_score']:.3f} | {row['corr_success']:.3f} | "
            f"{row['corr_darkness_failure']:.3f} | {row['corr_temporal_overlap_failure']:.3f} | {row['recommended_direction']} |"
        )
    lines.append("")
    lines.append("## Actionable Tuning")
    lines.append("")
    for row in ranked_features[:5]:
        lines.append(
            f"- `{row['feature']}`: `{row['recommended_direction']}` "
            f"(impact {row['impact_score']:.3f}, dark-fail corr {row['corr_darkness_failure']:.3f}, "
            f"gap-fail corr {row['corr_temporal_overlap_failure']:.3f})"
        )
    lines.append("")
    with open(RESULTS_DIR / f"{out_prefix}_report.md", "w") as f:
        f.write("\n".join(lines) + "\n")


def print_report(summary: dict[str, Any], ranked_features: list[dict[str, Any]], out_prefix: str = "sensitivity") -> None:
    print("\n" + "=" * 104)
    print("  BLACKCOTTON SENSITIVITY ANALYSIS")
    print("=" * 104)
    print(f"\n  Trial records:         {summary['n_trial_records']}")
    print(f"  Candidate ranks:       {summary['n_candidate_ranks']}")
    print(f"  Overall success rate:  {summary['overall_success_rate']:.3f}")
    print(f"  Primary failure mode:  {summary['primary_failure_mode']}")
    print("\n  TOP FEATURE IMPACTS")
    print("  " + "-" * 102)
    print("  #  feature                         impact   corr_succ   corr_dark   corr_gap   direction")
    print("  " + "-" * 102)
    for i, row in enumerate(ranked_features[:12], 1):
        print(
            f"  {i:>2} {row['feature']:<30} {row['impact_score']:>7.3f} "
            f"{row['corr_success']:>10.3f} {row['corr_darkness_failure']:>11.3f} "
            f"{row['corr_temporal_overlap_failure']:>10.3f} {row['recommended_direction']:>10}"
        )
    print("  " + "-" * 102)
    print("\n  Saved:")
    print(f"    - results/{out_prefix}_summary.json")
    print(f"    - results/{out_prefix}_ranked_features.json")
    print(f"    - results/{out_prefix}_candidate_failures.json")
    print(f"    - results/{out_prefix}_report.md")
    print("\n" + "=" * 104)


if __name__ == "__main__":
    print("\n🧬 BlackCotton Sensitivity Analyzer")
    print("=" * 50)

    parser = argparse.ArgumentParser(description="Sensitivity analysis from robust trial records.")
    parser.add_argument(
        "--summary",
        default=str(RESULTS_DIR / "robust_optimization_summary.json"),
        help="Path to summary JSON containing thresholds.",
    )
    parser.add_argument(
        "--trials",
        default=str(RESULTS_DIR / "robust_trial_records.jsonl"),
        help="Path to JSONL trial records.",
    )
    parser.add_argument(
        "--out-prefix",
        default="sensitivity",
        help="Output filename prefix under results/.",
    )
    args = parser.parse_args()

    summary_path = Path(args.summary)
    trials_path = Path(args.trials)
    if not summary_path.exists() or not trials_path.exists():
        raise FileNotFoundError(
            "Missing robust optimization artifacts. "
            "Run `python -m src.robust_tradeoff_optimizer` first."
        )

    with open(summary_path) as f:
        robust_summary = json.load(f)
    thresholds = robust_summary["thresholds"]

    records = load_trial_records(trials_path)
    ranked_features = rank_feature_impacts(records, thresholds)
    candidate_table = candidate_failure_table(records, thresholds)

    dark_fails = np.mean([float(r["color_L"]) > float(thresholds["max_color_L"]) for r in records]) if records else 0.0
    gap_fails = np.mean([float(r["temporal_gap_days"]) < float(thresholds["min_temporal_gap_days"]) for r in records]) if records else 0.0
    str_fails = np.mean([float(r["strength_g_tex"]) < float(thresholds["min_strength_g_tex"]) for r in records]) if records else 0.0
    yld_fails = np.mean([float(r["yield_index"]) < float(thresholds["min_yield_index"]) for r in records]) if records else 0.0
    failure_modes = {
        "darkness_failure": float(dark_fails),
        "temporal_overlap_failure": float(gap_fails),
        "strength_failure": float(str_fails),
        "yield_failure": float(yld_fails),
    }
    primary_failure_mode = max(failure_modes.items(), key=lambda kv: kv[1])[0]

    summary = {
        "n_trial_records": len(records),
        "n_candidate_ranks": len({int(r["candidate_input_rank"]) for r in records}),
        "overall_success_rate": float(np.mean([1.0 if r["success"] else 0.0 for r in records])) if records else 0.0,
        "failure_mode_rates": failure_modes,
        "primary_failure_mode": primary_failure_mode,
    }

    save_outputs(
        summary=summary,
        ranked_features=ranked_features,
        candidate_table=candidate_table,
        out_prefix=args.out_prefix,
    )
    print_report(summary=summary, ranked_features=ranked_features, out_prefix=args.out_prefix)
    print("\n✅ Sensitivity analysis complete!")
