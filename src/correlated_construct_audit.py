#!/usr/bin/env python3
"""
correlated_construct_audit.py — Correlated-failure audit for construct stability.
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from src.adversarial_robustness_suite import scale_noise, strict_thresholds
from src.config_loader import load_config, resolve_config_path
from src.failure_risk_model import stability_design_metadata
from src.robustness_analyzer import (
    DEFAULT_NOISE,
    DEFAULT_THRESHOLDS,
    run_robustness_analysis,
)
from src.worst_case_hardening_sprint import (
    build_deterministic_pool,
    pick_top_unique_robust_rows,
    select_seed_candidates,
)

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"
DEFAULT_VARIANT_PATH = BASE_DIR / "config" / "variants" / "variant_anti_silencing_construct.yaml"


def _rel(path: Path) -> str:
    return str(path.resolve().relative_to(BASE_DIR))


def run_arm(
    label: str,
    config_path: Path,
    n_trials: int,
    seed: int,
    seed_candidates_limit: int,
    min_temporal_gap_days: float,
    correlated_profile: str,
) -> dict[str, Any]:
    params = load_config(config_path)
    thresholds = strict_thresholds(dict(DEFAULT_THRESHOLDS))
    noise = scale_noise(dict(DEFAULT_NOISE), factor=1.8)

    deterministic_pool = build_deterministic_pool(
        params=params,
        strict_thresholds=thresholds,
        min_temporal_gap_days=float(min_temporal_gap_days),
    )
    seed_candidates = select_seed_candidates(deterministic_pool, limit=int(seed_candidates_limit))

    independent_rows = run_robustness_analysis(
        params=params,
        candidates=seed_candidates,
        n_trials=int(n_trials),
        seed=int(seed),
        thresholds=thresholds,
        noise=noise,
        collect_trial_records=False,
    )
    correlated_rows = run_robustness_analysis(
        params=params,
        candidates=seed_candidates,
        n_trials=int(n_trials),
        seed=int(seed) + 101,
        thresholds=thresholds,
        noise=noise,
        correlated_profile=correlated_profile,
        collect_trial_records=False,
    )

    top3 = pick_top_unique_robust_rows(
        correlated_rows,
        n=3,
        min_success_rate=0.20,
        max_fragility_index=0.40,
    )
    if not top3:
        top3 = correlated_rows[:3]

    independent_top = independent_rows[0] if independent_rows else {}
    correlated_top = correlated_rows[0] if correlated_rows else {}
    stability_meta = stability_design_metadata(params)

    return {
        "label": label,
        "config_path": _rel(config_path),
        "strict_thresholds": thresholds,
        "noise": noise,
        "correlated_profile": correlated_profile,
        "stability_design": stability_meta["design"],
        "stability_silencing_scale": float(stability_meta["silencing_scale"]),
        "stability_event_cv_scale": float(stability_meta["event_expression_cv_scale"]),
        "deterministic_pool_size": int(len(deterministic_pool)),
        "seed_candidates": int(len(seed_candidates)),
        "independent_top_success_rate": float(independent_top.get("success_rate", 0.0)),
        "independent_top_robust_score": float(independent_top.get("robust_score", 0.0)),
        "correlated_top_success_rate": float(correlated_top.get("success_rate", 0.0)),
        "correlated_top_robust_score": float(correlated_top.get("robust_score", 0.0)),
        "correlated_fragility_index": float(correlated_top.get("fragility_index", 1.0)),
        "correlated_penalty_success_rate": float(correlated_top.get("success_rate", 0.0))
        - float(independent_top.get("success_rate", 0.0)),
        "correlated_penalty_robust_score": float(correlated_top.get("robust_score", 0.0))
        - float(independent_top.get("robust_score", 0.0)),
        "correlated_top3": top3,
    }


def write_report(path: Path, summary: dict[str, Any]) -> None:
    baseline = summary["baseline"]
    variant = summary["variant"]

    lines: list[str] = []
    lines.append("# Correlated Construct Audit")
    lines.append("")
    lines.append(
        f"Objective: compare {baseline['label']} against {variant['label']} under strict "
        "high-noise stress with correlated failure bundles."
    )
    lines.append("")
    lines.append("## Headline")
    lines.append("")
    lines.append(
        f"- Baseline correlated top success: `{baseline['correlated_top_success_rate']:.3f}`"
    )
    lines.append(
        f"- Variant correlated top success: `{variant['correlated_top_success_rate']:.3f}`"
    )
    lines.append(
        f"- Delta: `{summary['delta_correlated_top_success_rate']:+.3f}`"
    )
    lines.append(
        f"- Correlated profile: `{summary['correlated_profile']}`"
    )
    lines.append("")
    lines.append("## Comparison")
    lines.append("")
    lines.append("| Arm | Independent Top Success | Correlated Top Success | Correlated Top Robust | Fragility |")
    lines.append("|---|---:|---:|---:|---:|")
    for row in (baseline, variant):
        lines.append(
            f"| {row['label']} | {row['independent_top_success_rate']:.3f} | "
            f"{row['correlated_top_success_rate']:.3f} | {row['correlated_top_robust_score']:.3f} | "
            f"{row['correlated_fragility_index']:.3f} |"
        )
    lines.append("")
    lines.append(f"## {variant['label']} Top 3")
    lines.append("")
    lines.append("| Rank | Success | Robust | Fragility | p50 L* | p50 Strength | p50 Yield | Params |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---|")
    for i, row in enumerate(variant["correlated_top3"], 1):
        base = row["base_candidate"]
        lines.append(
            f"| {i} | {row['success_rate']:.3f} | {row['robust_score']:.3f} | "
            f"{row.get('fragility_index', 0.0):.3f} | {row['p50_color_L']:.2f} | "
            f"{row['p50_strength_g_tex']:.2f} | {row['p50_yield_index']:.3f} | "
            f"{base['mat_activation_dpa']:.0f}/{base['scw_activation_dpa']:.0f}, "
            f"k={base['k_competition']:.2f}, eff={base['melanin_efficiency']:.2f}, "
            f"ret={base['late_retention_factor']:.2f} |"
        )
    lines.append("")
    lines.append("## Construct Stability Assumptions")
    lines.append("")
    lines.append(
        f"- Variant silencing scale: `{variant['stability_silencing_scale']:.3f}`"
    )
    lines.append(
        f"- Variant event-CV scale: `{variant['stability_event_cv_scale']:.3f}`"
    )
    lines.append("")

    path.write_text("\n".join(lines) + "\n")


def run(
    baseline_config: Path,
    variant_config: Path,
    n_trials: int,
    seed: int,
    seed_candidates_limit: int,
    min_temporal_gap_days: float,
    correlated_profile: str,
    out_prefix: str,
    baseline_label: str,
    variant_label: str,
) -> dict[str, Any]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    baseline = run_arm(
        label=baseline_label,
        config_path=baseline_config,
        n_trials=n_trials,
        seed=seed,
        seed_candidates_limit=seed_candidates_limit,
        min_temporal_gap_days=min_temporal_gap_days,
        correlated_profile=correlated_profile,
    )
    variant = run_arm(
        label=variant_label,
        config_path=variant_config,
        n_trials=n_trials,
        seed=seed + 1000,
        seed_candidates_limit=seed_candidates_limit,
        min_temporal_gap_days=min_temporal_gap_days,
        correlated_profile=correlated_profile,
    )

    summary = {
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "baseline": baseline,
        "variant": variant,
        "correlated_profile": correlated_profile,
        "n_trials": int(n_trials),
        "seed_candidates_limit": int(seed_candidates_limit),
        "min_temporal_gap_days": float(min_temporal_gap_days),
        "delta_correlated_top_success_rate": float(variant["correlated_top_success_rate"] - baseline["correlated_top_success_rate"]),
        "delta_correlated_top_robust_score": float(variant["correlated_top_robust_score"] - baseline["correlated_top_robust_score"]),
    }

    summary_path = RESULTS_DIR / f"{out_prefix}_summary.json"
    report_path = RESULTS_DIR / f"{out_prefix}_report.md"
    top3_path = RESULTS_DIR / f"{out_prefix}_top3.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    top3_path.write_text(
        json.dumps(
            {
                "baseline_top3": baseline["correlated_top3"],
                "variant_top3": variant["correlated_top3"],
            },
            indent=2,
        )
        + "\n"
    )
    write_report(report_path, summary)
    return {
        "summary_path": summary_path,
        "report_path": report_path,
        "top3_path": top3_path,
        "summary": summary,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run correlated-failure construct audit.")
    parser.add_argument(
        "--baseline-config",
        default=str(resolve_config_path()),
        help="Baseline config path.",
    )
    parser.add_argument(
        "--variant-config",
        default=str(DEFAULT_VARIANT_PATH),
        help="Variant config path.",
    )
    parser.add_argument("--baseline-label", default="Production", help="Label for the baseline arm.")
    parser.add_argument("--variant-label", default="Anti-Silencing", help="Label for the variant arm.")
    parser.add_argument("--n-trials", type=int, default=140, help="Trials per candidate.")
    parser.add_argument("--seed", type=int, default=20260307, help="RNG seed.")
    parser.add_argument("--seed-candidates", type=int, default=220, help="Number of deterministic seeds.")
    parser.add_argument("--min-temporal-gap-days", type=float, default=0.5, help="Deterministic minimum gap filter.")
    parser.add_argument(
        "--correlated-profile",
        default="construct_bundle_v1",
        help="Correlated failure profile name.",
    )
    parser.add_argument(
        "--out-prefix",
        default="correlated_construct_audit_2026_03_07",
        help="Output prefix under results/.",
    )
    args = parser.parse_args()

    output = run(
        baseline_config=Path(args.baseline_config),
        variant_config=Path(args.variant_config),
        n_trials=max(int(args.n_trials), 60),
        seed=int(args.seed),
        seed_candidates_limit=max(int(args.seed_candidates), 40),
        min_temporal_gap_days=float(args.min_temporal_gap_days),
        correlated_profile=str(args.correlated_profile),
        out_prefix=str(args.out_prefix),
        baseline_label=str(args.baseline_label),
        variant_label=str(args.variant_label),
    )

    summary = output["summary"]
    print("\nCorrelated construct audit complete")
    print("=" * 36)
    print(f"Baseline correlated top success: {summary['baseline']['correlated_top_success_rate']:.3f}")
    print(f"Variant correlated top success: {summary['variant']['correlated_top_success_rate']:.3f}")
    print(f"Delta: {summary['delta_correlated_top_success_rate']:+.3f}")
    print(f"Saved: {output['summary_path'].relative_to(BASE_DIR)}")
    print(f"Saved: {output['report_path'].relative_to(BASE_DIR)}")
    print(f"Saved: {output['top3_path'].relative_to(BASE_DIR)}")
