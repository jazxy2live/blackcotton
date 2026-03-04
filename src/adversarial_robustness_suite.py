#!/usr/bin/env python3
"""
adversarial_robustness_suite.py — Red-team stress tests for final candidates
=============================================================================

Runs the frozen Top 3 through stricter threshold gates and harsher uncertainty
scenarios to estimate worst-case robustness before wet-lab execution.

Usage:
    python -m src.adversarial_robustness_suite
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from src.robustness_analyzer import (
    DEFAULT_NOISE,
    DEFAULT_THRESHOLDS,
    run_robustness_analysis,
)
from src.tradeoff_optimizer import load_params

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"


def load_json(path: Path, default: Any) -> Any:
    if not path.exists():
        return default
    with open(path) as f:
        return json.load(f)


def scale_noise(noise: dict[str, float], factor: float) -> dict[str, float]:
    scaled = dict(noise)
    for key, value in noise.items():
        scaled[key] = float(value) * float(factor)
    # Keep additive retention jitter within reasonable bounds.
    if "retention_sigma" in scaled:
        scaled["retention_sigma"] = min(0.20, scaled["retention_sigma"])
    return scaled


def strict_thresholds(base: dict[str, float]) -> dict[str, float]:
    out = dict(base)
    out["max_color_L"] = 22.0
    out["min_strength_g_tex"] = 28.5
    out["min_yield_index"] = 0.88
    out["min_temporal_gap_days"] = 0.5
    if "max_toxicity_pre_cellulose" in out:
        out["max_toxicity_pre_cellulose"] = 0.08
    return out


def normalize_top_candidate_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    candidates = []
    for i, row in enumerate(rows, 1):
        base = row.get("base_candidate", row)
        candidates.append(
            {
                "candidate_arm": f"CAND_{i}",
                "input_rank": int(row.get("input_rank", i)),
                "color_L": float(base.get("color_L", row.get("p50_color_L", 0.0))),
                "strength_g_tex": float(base.get("strength_g_tex", row.get("p50_strength_g_tex", 0.0))),
                "yield_index": float(base.get("yield_index", row.get("p50_yield_index", 0.0))),
                "temporal_gap_days": float(base.get("temporal_gap_days", row.get("p50_temporal_gap_days", 0.0))),
                "toxicity_pre_cellulose": float(
                    base.get("toxicity_pre_cellulose", row.get("p50_toxicity_pre_cellulose", 0.0))
                ),
                "mat_activation_dpa": float(base.get("mat_activation_dpa")),
                "scw_activation_dpa": float(base.get("scw_activation_dpa")),
                "mat_strength": float(base.get("mat_strength")),
                "scw_strength": float(base.get("scw_strength")),
                "k_competition": float(base.get("k_competition")),
                "melanin_efficiency": float(base.get("melanin_efficiency")),
                "late_retention_factor": float(base.get("late_retention_factor")),
                "composite_score": float(row.get("robust_score", row.get("composite_score", 0.0))),
            }
        )
    return candidates


def scenario_configs() -> list[dict[str, Any]]:
    base_noise = dict(DEFAULT_NOISE)
    base_thr = dict(DEFAULT_THRESHOLDS)
    strict_thr = strict_thresholds(base_thr)
    return [
        {
            "name": "baseline",
            "thresholds": base_thr,
            "noise": base_noise,
            "description": "Current production thresholds and uncertainty model.",
        },
        {
            "name": "strict_gate",
            "thresholds": strict_thr,
            "noise": base_noise,
            "description": "Stricter pass/fail gates, same uncertainty model.",
        },
        {
            "name": "high_noise",
            "thresholds": base_thr,
            "noise": scale_noise(base_noise, factor=1.8),
            "description": "Current gates with heavier biological/process uncertainty.",
        },
        {
            "name": "strict_high_noise",
            "thresholds": strict_thr,
            "noise": scale_noise(base_noise, factor=1.8),
            "description": "Strict gates plus heavy uncertainty (worst-case in-silico).",
        },
    ]


def scenario_results_rows(
    scenario_name: str,
    rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    out = []
    for row in rows:
        out.append(
            {
                "scenario": scenario_name,
                "candidate_input_rank": int(row["input_rank"]),
                "success_rate": float(row["success_rate"]),
                "robust_score": float(row["robust_score"]),
                "risk_temporal_overlap": float(row["risk_temporal_overlap"]),
                "risk_darkness_failure": float(row["risk_darkness_failure"]),
                "risk_strength_failure": float(row["risk_strength_failure"]),
                "risk_yield_failure": float(row["risk_yield_failure"]),
                "risk_toxicity_failure": float(row.get("risk_toxicity_failure", 0.0)),
                "fragility_index": float(row.get("fragility_index", 0.0)),
                "p50_color_L": float(row["p50_color_L"]),
                "p50_strength_g_tex": float(row["p50_strength_g_tex"]),
                "p50_yield_index": float(row["p50_yield_index"]),
                "p50_temporal_gap_days": float(row["p50_temporal_gap_days"]),
                "p50_toxicity_pre_cellulose": float(row.get("p50_toxicity_pre_cellulose", 0.0)),
            }
        )
    return out


def aggregate_worst_case(per_scenario_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_candidate: dict[int, list[dict[str, Any]]] = {}
    for row in per_scenario_rows:
        by_candidate.setdefault(int(row["candidate_input_rank"]), []).append(row)

    summary = []
    for rank, rows in sorted(by_candidate.items()):
        success_vals = [float(r["success_rate"]) for r in rows]
        robust_vals = [float(r["robust_score"]) for r in rows]
        worst_row = min(rows, key=lambda r: (float(r["success_rate"]), float(r["robust_score"])))
        summary.append(
            {
                "candidate_input_rank": int(rank),
                "num_scenarios": int(len(rows)),
                "min_success_rate": float(min(success_vals)),
                "mean_success_rate": float(sum(success_vals) / len(success_vals)),
                "min_robust_score": float(min(robust_vals)),
                "mean_robust_score": float(sum(robust_vals) / len(robust_vals)),
                "worst_scenario": str(worst_row["scenario"]),
                "worst_scenario_success_rate": float(worst_row["success_rate"]),
                "worst_scenario_robust_score": float(worst_row["robust_score"]),
            }
        )
    summary.sort(key=lambda x: (x["min_success_rate"], x["min_robust_score"]), reverse=True)
    return summary


def write_report(
    run_summary: dict[str, Any],
    scenario_summary: list[dict[str, Any]],
    candidate_summary: list[dict[str, Any]],
    out_prefix: str,
) -> None:
    lines: list[str] = []
    lines.append("# Adversarial Robustness Report")
    lines.append("")
    lines.append("Red-team robustness check across strict thresholds and high-noise scenarios.")
    lines.append("")
    lines.append("## Run Summary")
    lines.append("")
    lines.append(f"- Generated: `{run_summary['generated_at']}`")
    lines.append(f"- Candidates tested: `{run_summary['num_candidates']}`")
    lines.append(f"- Scenarios: `{run_summary['num_scenarios']}`")
    lines.append(f"- Trials per scenario/candidate: `{run_summary['n_trials_per_candidate']}`")
    lines.append(f"- Worst-case best candidate input rank: `{run_summary['worst_case_best_candidate_input_rank']}`")
    lines.append("")
    lines.append("## Scenario Leaderboard")
    lines.append("")
    lines.append("| Scenario | Top Candidate Rank | Top Success | Top Robust |")
    lines.append("|---|---:|---:|---:|")
    for row in scenario_summary:
        lines.append(
            f"| {row['scenario']} | {row['top_candidate_input_rank']} | "
            f"{row['top_success_rate']:.3f} | {row['top_robust_score']:.3f} |"
        )
    lines.append("")
    lines.append("## Worst-Case Candidate Ranking")
    lines.append("")
    lines.append("| Candidate Rank | Min Success | Mean Success | Min Robust | Worst Scenario |")
    lines.append("|---:|---:|---:|---:|---|")
    for row in candidate_summary:
        lines.append(
            f"| {row['candidate_input_rank']} | {row['min_success_rate']:.3f} | "
            f"{row['mean_success_rate']:.3f} | {row['min_robust_score']:.3f} | {row['worst_scenario']} |"
        )
    lines.append("")

    (RESULTS_DIR / f"{out_prefix}_report.md").write_text("\n".join(lines) + "\n")


def run(n_trials: int, seed: int, out_prefix: str) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    top_rows = load_json(RESULTS_DIR / "final_lab_top3.json", [])
    if not top_rows:
        raise FileNotFoundError(
            "Missing results/final_lab_top3.json. Run `python -m src.adaptive_ode_robust_pipeline` first."
        )
    candidates = normalize_top_candidate_rows(top_rows)
    params = load_params()

    all_rows: list[dict[str, Any]] = []
    scenario_summary: list[dict[str, Any]] = []
    scenario_meta = scenario_configs()

    for i, scenario in enumerate(scenario_meta):
        robust_rows = run_robustness_analysis(
            params=params,
            candidates=candidates,
            n_trials=n_trials,
            seed=seed + i * 101,
            thresholds=dict(scenario["thresholds"]),
            noise=dict(scenario["noise"]),
            collect_trial_records=False,
        )
        scenario_rows = scenario_results_rows(scenario["name"], robust_rows)
        all_rows.extend(scenario_rows)

        top = robust_rows[0] if robust_rows else {}
        scenario_summary.append(
            {
                "scenario": scenario["name"],
                "description": scenario["description"],
                "thresholds": scenario["thresholds"],
                "noise": scenario["noise"],
                "top_candidate_input_rank": int(top.get("input_rank", 0)) if top else 0,
                "top_success_rate": float(top.get("success_rate", 0.0)) if top else 0.0,
                "top_robust_score": float(top.get("robust_score", 0.0)) if top else 0.0,
            }
        )

    candidate_summary = aggregate_worst_case(all_rows)
    best_worst_case = candidate_summary[0] if candidate_summary else {}

    run_summary = {
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "num_candidates": len(candidates),
        "num_scenarios": len(scenario_meta),
        "n_trials_per_candidate": int(n_trials),
        "seed": int(seed),
        "worst_case_best_candidate_input_rank": int(best_worst_case.get("candidate_input_rank", 0)) if best_worst_case else 0,
    }

    with open(RESULTS_DIR / f"{out_prefix}_summary.json", "w") as f:
        json.dump(run_summary, f, indent=2)
    with open(RESULTS_DIR / f"{out_prefix}_scenarios.json", "w") as f:
        json.dump(scenario_summary, f, indent=2)
    with open(RESULTS_DIR / f"{out_prefix}_rows.json", "w") as f:
        json.dump(all_rows, f, indent=2)
    with open(RESULTS_DIR / f"{out_prefix}_candidate_summary.json", "w") as f:
        json.dump(candidate_summary, f, indent=2)
    write_report(run_summary, scenario_summary, candidate_summary, out_prefix=out_prefix)

    print("\n🛡️ BlackCotton Adversarial Robustness Suite")
    print("=" * 49)
    print(f"Candidates tested: {len(candidates)}")
    print(f"Scenarios: {len(scenario_meta)}")
    print(f"Trials per scenario/candidate: {n_trials}")
    if best_worst_case:
        print(
            f"Worst-case best candidate rank: {best_worst_case['candidate_input_rank']} "
            f"(min success={best_worst_case['min_success_rate']:.3f}, "
            f"min robust={best_worst_case['min_robust_score']:.3f})"
        )
    print("Saved:")
    print(f"  - results/{out_prefix}_summary.json")
    print(f"  - results/{out_prefix}_scenarios.json")
    print(f"  - results/{out_prefix}_rows.json")
    print(f"  - results/{out_prefix}_candidate_summary.json")
    print(f"  - results/{out_prefix}_report.md")
    print("=" * 49)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run adversarial robustness stress tests.")
    parser.add_argument("--n-trials", type=int, default=240, help="Trials per scenario/candidate.")
    parser.add_argument("--seed", type=int, default=20260302, help="RNG seed.")
    parser.add_argument("--out-prefix", default="adversarial_robustness", help="Output prefix under results/.")
    args = parser.parse_args()

    run(
        n_trials=max(50, int(args.n_trials)),
        seed=int(args.seed),
        out_prefix=str(args.out_prefix),
    )
