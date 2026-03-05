#!/usr/bin/env python3
"""
calibration_impact_analyzer.py — End-to-end baseline vs calibrated comparison
=============================================================================

Compares two parameter configurations (baseline and calibrated) through:
  1) Expression timing snapshot
  2) Constrained optimization (coarse -> ODE refinement)
  3) Robustness ranking under default uncertainty
  4) Worst-case robustness under strict + high-noise stress

Usage:
    python -m src.calibration_impact_analyzer
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

from src.adversarial_robustness_suite import scale_noise, strict_thresholds
from src.config_loader import resolve_config_path
from src.expression_model import run_simulation
from src.robustness_analyzer import DEFAULT_NOISE, DEFAULT_THRESHOLDS, run_robustness_analysis
from src.tradeoff_optimizer import (
    pareto_front,
    refine_candidates_with_odes,
    run_sweep,
    select_seed_candidates,
    select_top_candidates,
)
from src.transcriptome_calibrator import expression_snapshot

BASE_DIR = Path(__file__).resolve().parent.parent
CONFIG_DIR = BASE_DIR / "config"
RESULTS_DIR = BASE_DIR / "results"


def load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Missing config file: {path}")
    with open(path) as f:
        return yaml.safe_load(f)


def config_eval(
    label: str,
    params: dict[str, Any],
    max_seed: int,
    max_refine: int,
    top_for_robust: int,
    n_trials: int,
    seed: int,
) -> dict[str, Any]:
    expr = run_simulation(params)
    expr_snap = expression_snapshot(expr)

    coarse = run_sweep(params)
    seed_pool = select_seed_candidates(coarse, n=max_seed)
    refined = refine_candidates_with_odes(params=params, coarse_candidates=seed_pool, max_refine=max_refine)
    front = pareto_front(refined)
    top = select_top_candidates(refined, n=max(12, top_for_robust))
    robust_top_pool = top[:top_for_robust]

    robust_default = run_robustness_analysis(
        params=params,
        candidates=robust_top_pool,
        n_trials=n_trials,
        seed=seed,
        thresholds=dict(DEFAULT_THRESHOLDS),
        noise=dict(DEFAULT_NOISE),
        collect_trial_records=False,
    )
    robust_worst = run_robustness_analysis(
        params=params,
        candidates=robust_top_pool,
        n_trials=n_trials,
        seed=seed + 202,
        thresholds=strict_thresholds(dict(DEFAULT_THRESHOLDS)),
        noise=scale_noise(dict(DEFAULT_NOISE), 1.8),
        collect_trial_records=False,
    )

    det_best = top[0] if top else {}
    robust_best = robust_default[0] if robust_default else {}
    worst_best = robust_worst[0] if robust_worst else {}

    return {
        "label": label,
        "n_coarse": int(len(coarse)),
        "n_seed": int(len(seed_pool)),
        "n_refined": int(len(refined)),
        "n_pareto": int(len(front)),
        "n_top": int(len(top)),
        "expression_snapshot": expr_snap,
        "deterministic_best": {
            "color_L": float(det_best.get("color_L", 82.0)),
            "strength_g_tex": float(det_best.get("strength_g_tex", 0.0)),
            "yield_index": float(det_best.get("yield_index", 0.0)),
            "temporal_gap_days": float(det_best.get("temporal_gap_days", 0.0)),
            "composite_score": float(det_best.get("composite_score", 0.0)),
        },
        "robust_default_best": {
            "success_rate": float(robust_best.get("success_rate", 0.0)),
            "robust_score": float(robust_best.get("robust_score", 0.0)),
            "p50_color_L": float(robust_best.get("p50_color_L", 82.0)),
            "p50_strength_g_tex": float(robust_best.get("p50_strength_g_tex", 0.0)),
            "p50_yield_index": float(robust_best.get("p50_yield_index", 0.0)),
            "p50_temporal_gap_days": float(robust_best.get("p50_temporal_gap_days", 0.0)),
        },
        "robust_worst_case_best": {
            "success_rate": float(worst_best.get("success_rate", 0.0)),
            "robust_score": float(worst_best.get("robust_score", 0.0)),
            "p50_color_L": float(worst_best.get("p50_color_L", 82.0)),
            "p50_strength_g_tex": float(worst_best.get("p50_strength_g_tex", 0.0)),
            "p50_yield_index": float(worst_best.get("p50_yield_index", 0.0)),
            "p50_temporal_gap_days": float(worst_best.get("p50_temporal_gap_days", 0.0)),
        },
    }


def delta(calibrated: float, baseline: float) -> float:
    return float(calibrated - baseline)


def pick_winner(baseline: dict[str, Any], calibrated: dict[str, Any]) -> dict[str, Any]:
    b = baseline["robust_worst_case_best"]
    c = calibrated["robust_worst_case_best"]
    b_key = (float(b["success_rate"]), float(b["robust_score"]))
    c_key = (float(c["success_rate"]), float(c["robust_score"]))
    if c_key > b_key:
        return {
            "winner": "calibrated",
            "reason": "better worst-case success/robust score",
        }
    if c_key < b_key:
        return {
            "winner": "baseline",
            "reason": "better worst-case success/robust score",
        }
    return {
        "winner": "tie",
        "reason": "equal worst-case success and robust score",
    }


def write_report(
    summary: dict[str, Any],
    out_prefix: str,
) -> None:
    b = summary["baseline"]
    c = summary["calibrated"]
    d = summary["delta_cal_minus_base"]
    decision = summary["decision"]

    lines: list[str] = []
    lines.append("# Calibration Impact Report")
    lines.append("")
    lines.append("End-to-end comparison of baseline vs calibrated parameters.")
    lines.append("")
    lines.append("## Decision")
    lines.append("")
    lines.append(f"- Winner: `{decision['winner']}`")
    lines.append(f"- Reason: `{decision['reason']}`")
    lines.append("")
    lines.append("## Optimization Footprint")
    lines.append("")
    lines.append("| Metric | Baseline | Calibrated | Delta (cal-base) |")
    lines.append("|---|---:|---:|---:|")
    for key in ("n_coarse", "n_seed", "n_refined", "n_pareto", "n_top"):
        lines.append(f"| {key} | {b[key]} | {c[key]} | {d[key]:.3f} |")
    lines.append("")
    lines.append("## Robustness Comparison")
    lines.append("")
    lines.append("| Metric | Baseline | Calibrated | Delta (cal-base) |")
    lines.append("|---|---:|---:|---:|")
    for key in ("success_rate", "robust_score", "p50_color_L", "p50_strength_g_tex", "p50_yield_index", "p50_temporal_gap_days"):
        lines.append(
            f"| default.{key} | {b['robust_default_best'][key]:.3f} | {c['robust_default_best'][key]:.3f} | "
            f"{d['robust_default_best'][key]:+.3f} |"
        )
    for key in ("success_rate", "robust_score", "p50_color_L", "p50_strength_g_tex", "p50_yield_index", "p50_temporal_gap_days"):
        lines.append(
            f"| worst.{key} | {b['robust_worst_case_best'][key]:.3f} | {c['robust_worst_case_best'][key]:.3f} | "
            f"{d['robust_worst_case_best'][key]:+.3f} |"
        )
    lines.append("")
    lines.append("## Expression Timing Snapshot")
    lines.append("")
    for key in ("cellulose_90pct_dpa", "melA_halfmax_dpa", "temporal_gap_days"):
        lines.append(
            f"- `{key}` baseline `{b['expression_snapshot'][key]:.3f}` vs calibrated `{c['expression_snapshot'][key]:.3f}` "
            f"(delta `{d['expression_snapshot'][key]:+.3f}`)"
        )
    lines.append("")

    (RESULTS_DIR / f"{out_prefix}_report.md").write_text("\n".join(lines) + "\n")


def run(
    baseline_config: Path,
    calibrated_config: Path,
    max_seed: int,
    max_refine: int,
    top_for_robust: int,
    n_trials: int,
    seed: int,
    out_prefix: str,
) -> dict[str, Any]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    baseline_params = load_yaml(baseline_config)
    calibrated_params = load_yaml(calibrated_config)

    baseline = config_eval(
        label="baseline",
        params=baseline_params,
        max_seed=max_seed,
        max_refine=max_refine,
        top_for_robust=top_for_robust,
        n_trials=n_trials,
        seed=seed,
    )
    calibrated = config_eval(
        label="calibrated",
        params=calibrated_params,
        max_seed=max_seed,
        max_refine=max_refine,
        top_for_robust=top_for_robust,
        n_trials=n_trials,
        seed=seed + 999,
    )

    delta_obj = {
        "n_coarse": delta(calibrated["n_coarse"], baseline["n_coarse"]),
        "n_seed": delta(calibrated["n_seed"], baseline["n_seed"]),
        "n_refined": delta(calibrated["n_refined"], baseline["n_refined"]),
        "n_pareto": delta(calibrated["n_pareto"], baseline["n_pareto"]),
        "n_top": delta(calibrated["n_top"], baseline["n_top"]),
        "expression_snapshot": {
            k: delta(calibrated["expression_snapshot"][k], baseline["expression_snapshot"][k])
            for k in baseline["expression_snapshot"].keys()
        },
        "robust_default_best": {
            k: delta(calibrated["robust_default_best"][k], baseline["robust_default_best"][k])
            for k in baseline["robust_default_best"].keys()
        },
        "robust_worst_case_best": {
            k: delta(calibrated["robust_worst_case_best"][k], baseline["robust_worst_case_best"][k])
            for k in baseline["robust_worst_case_best"].keys()
        },
    }
    decision = pick_winner(baseline=baseline, calibrated=calibrated)

    summary = {
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "baseline_config": str(baseline_config.relative_to(BASE_DIR)),
        "calibrated_config": str(calibrated_config.relative_to(BASE_DIR)),
        "settings": {
            "max_seed": int(max_seed),
            "max_refine": int(max_refine),
            "top_for_robust": int(top_for_robust),
            "n_trials": int(n_trials),
            "seed": int(seed),
        },
        "baseline": baseline,
        "calibrated": calibrated,
        "delta_cal_minus_base": delta_obj,
        "decision": decision,
    }

    with open(RESULTS_DIR / f"{out_prefix}_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    with open(RESULTS_DIR / f"{out_prefix}_baseline.json", "w") as f:
        json.dump(baseline, f, indent=2)
    with open(RESULTS_DIR / f"{out_prefix}_calibrated.json", "w") as f:
        json.dump(calibrated, f, indent=2)

    write_report(summary, out_prefix=out_prefix)

    print("\n📊 BlackCotton Calibration Impact Analyzer")
    print("=" * 46)
    print(f"Baseline config:   {summary['baseline_config']}")
    print(f"Calibrated config: {summary['calibrated_config']}")
    print(f"Winner:            {decision['winner']} ({decision['reason']})")
    print(
        f"Worst-case success (base -> cal): "
        f"{baseline['robust_worst_case_best']['success_rate']:.3f} -> "
        f"{calibrated['robust_worst_case_best']['success_rate']:.3f}"
    )
    print(
        f"Worst-case robust score (base -> cal): "
        f"{baseline['robust_worst_case_best']['robust_score']:.3f} -> "
        f"{calibrated['robust_worst_case_best']['robust_score']:.3f}"
    )
    print("Saved:")
    print(f"  - results/{out_prefix}_summary.json")
    print(f"  - results/{out_prefix}_baseline.json")
    print(f"  - results/{out_prefix}_calibrated.json")
    print(f"  - results/{out_prefix}_report.md")
    print("=" * 46)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare baseline vs calibrated parameter configs.")
    parser.add_argument(
        "--baseline-config",
        default=str(resolve_config_path()),
        help="Baseline YAML config path.",
    )
    parser.add_argument(
        "--calibrated-config",
        default=str(CONFIG_DIR / "parameters_calibrated.yaml"),
        help="Calibrated YAML config path.",
    )
    parser.add_argument("--max-seed", type=int, default=96, help="Seed candidates to retain from coarse sweep.")
    parser.add_argument("--max-refine", type=int, default=96, help="Candidates to ODE-refine.")
    parser.add_argument("--top-for-robust", type=int, default=8, help="Top refined candidates for robustness analysis.")
    parser.add_argument("--n-trials", type=int, default=120, help="Monte Carlo trials per candidate.")
    parser.add_argument("--seed", type=int, default=20260302, help="Random seed base.")
    parser.add_argument("--out-prefix", default="calibration_impact", help="Output prefix under results/.")
    args = parser.parse_args()

    run(
        baseline_config=Path(args.baseline_config),
        calibrated_config=Path(args.calibrated_config),
        max_seed=max(20, int(args.max_seed)),
        max_refine=max(20, int(args.max_refine)),
        top_for_robust=max(3, int(args.top_for_robust)),
        n_trials=max(40, int(args.n_trials)),
        seed=int(args.seed),
        out_prefix=str(args.out_prefix),
    )
