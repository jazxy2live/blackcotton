#!/usr/bin/env python3
"""
lab_validation_planner.py — Build a practical lab package from final Top 3
===========================================================================

Generates a concrete first-cycle experimental package from the frozen robust
shortlist so wet-lab validation can start with minimal ambiguity.

Usage:
    python -m src.lab_validation_planner
"""

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"


def load_json(path: Path, default: Any) -> Any:
    if not path.exists():
        return default
    with open(path) as f:
        return json.load(f)


def as_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def build_candidate_rows(
    final_top3: list[dict[str, Any]],
    replicates: int,
    run_id: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for rank, row in enumerate(final_top3, 1):
        base = row.get("base_candidate", {})
        rows.append(
            {
                "run_id": run_id,
                "arm_type": "candidate",
                "arm_id": f"CAND_{rank}",
                "priority_rank": rank,
                "source_input_rank": int(row.get("input_rank", rank)),
                "replicates": replicates,
                "mat_activation_dpa": as_float(base.get("mat_activation_dpa")),
                "scw_activation_dpa": as_float(base.get("scw_activation_dpa")),
                "mat_strength": as_float(base.get("mat_strength")),
                "scw_strength": as_float(base.get("scw_strength")),
                "k_competition": as_float(base.get("k_competition")),
                "melanin_efficiency": as_float(base.get("melanin_efficiency")),
                "late_retention_factor": as_float(base.get("late_retention_factor")),
                "pred_success_rate": as_float(row.get("success_rate")),
                "pred_robust_score": as_float(row.get("robust_score")),
                "pred_p50_color_L": as_float(row.get("p50_color_L")),
                "pred_p50_strength_g_tex": as_float(row.get("p50_strength_g_tex")),
                "pred_p50_yield_index": as_float(row.get("p50_yield_index")),
                "pred_p50_temporal_gap_days": as_float(row.get("p50_temporal_gap_days")),
                "pred_p50_toxicity_pre_cellulose": as_float(row.get("p50_toxicity_pre_cellulose")),
                "notes": "Top robust candidate from adaptive ODE pipeline",
            }
        )
    return rows


def build_control_rows(replicates: int, run_id: str) -> list[dict[str, Any]]:
    return [
        {
            "run_id": run_id,
            "arm_type": "control",
            "arm_id": "CTRL_WHITE_WT",
            "priority_rank": 0,
            "source_input_rank": 0,
            "replicates": replicates,
            "mat_activation_dpa": "",
            "scw_activation_dpa": "",
            "mat_strength": "",
            "scw_strength": "",
            "k_competition": "",
            "melanin_efficiency": "",
            "late_retention_factor": "",
            "pred_success_rate": "",
            "pred_robust_score": "",
            "pred_p50_color_L": "",
            "pred_p50_strength_g_tex": "",
            "pred_p50_yield_index": "",
            "pred_p50_temporal_gap_days": "",
            "notes": "Wild-type white cotton control",
        },
        {
            "run_id": run_id,
            "arm_type": "control",
            "arm_id": "CTRL_DYED_BLACK",
            "priority_rank": 0,
            "source_input_rank": 0,
            "replicates": replicates,
            "mat_activation_dpa": "",
            "scw_activation_dpa": "",
            "mat_strength": "",
            "scw_strength": "",
            "k_competition": "",
            "melanin_efficiency": "",
            "late_retention_factor": "",
            "pred_success_rate": "",
            "pred_robust_score": "",
            "pred_p50_color_L": "",
            "pred_p50_strength_g_tex": "",
            "pred_p50_yield_index": "",
            "pred_p50_temporal_gap_days": "",
            "notes": "Commercially dyed black cotton benchmark",
        },
    ]


def evaluate_gate(metrics: dict[str, float], thresholds: dict[str, float]) -> bool:
    tox_thr = thresholds.get("max_toxicity_pre_cellulose")
    tox_ok = True
    if tox_thr is not None:
        tox_ok = as_float(metrics.get("toxicity_pre_cellulose"), default=0.0) <= as_float(
            tox_thr, default=1e9
        )
    return (
        as_float(metrics.get("color_L"), default=1e9) <= as_float(thresholds.get("max_color_L"), default=25.0)
        and as_float(metrics.get("strength_g_tex"), default=-1e9) >= as_float(thresholds.get("min_strength_g_tex"), default=28.0)
        and as_float(metrics.get("yield_index"), default=-1e9) >= as_float(thresholds.get("min_yield_index"), default=0.85)
        and as_float(metrics.get("temporal_gap_days"), default=-1e9) >= as_float(thresholds.get("min_temporal_gap_days"), default=0.0)
        and tox_ok
    )


def build_measurement_template(arms: list[dict[str, Any]], replicates: int, run_id: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for arm in arms:
        for rep in range(1, replicates + 1):
            rows.append(
                {
                    "run_id": run_id,
                    "arm_id": arm["arm_id"],
                    "arm_type": arm["arm_type"],
                    "replicate_id": rep,
                    "sample_id": f"{arm['arm_id']}_R{rep:02d}",
                    "stage_dpa": "",
                    "assay_batch": "",
                    "melA_expr_norm": "",
                    "TYRP1_expr_norm": "",
                    "DCT_expr_norm": "",
                    "cellulose_marker_expr_norm": "",
                    "melanin_content_mg_g": "",
                    "color_L": "",
                    "strength_g_tex": "",
                    "yield_index": "",
                    "temporal_gap_days_est": "",
                    "wash_cycles": "",
                    "postwash_color_L": "",
                    "pass_timing": "",
                    "pass_darkness": "",
                    "pass_strength": "",
                    "pass_yield": "",
                    "pass_overall": "",
                    "operator": "",
                    "timestamp_iso": "",
                    "notes": "",
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_pass_fail_rules(path: Path, thresholds: dict[str, Any], run_id: str, final_top3: list[dict[str, Any]]) -> None:
    rules = {
        "run_id": run_id,
        "phase1_gate": {
            "description": "Minimum gate to continue candidate in this project.",
            "max_color_L": as_float(thresholds.get("max_color_L"), 25.0),
            "min_strength_g_tex": as_float(thresholds.get("min_strength_g_tex"), 28.0),
            "min_yield_index": as_float(thresholds.get("min_yield_index"), 0.85),
            "min_temporal_gap_days": as_float(thresholds.get("min_temporal_gap_days"), 0.0),
        },
        "phase2_stretch_gate": {
            "description": "Stricter target for partner demo readiness.",
            "max_color_L": 22.0,
            "min_strength_g_tex": 28.5,
            "min_yield_index": 0.88,
            "min_temporal_gap_days": 0.5,
        },
        "candidate_expectations": [
            {
                "arm_id": f"CAND_{i}",
                "pred_success_rate": as_float(row.get("success_rate")),
                "pred_robust_score": as_float(row.get("robust_score")),
                "pred_p50_color_L": as_float(row.get("p50_color_L")),
                "pred_p50_strength_g_tex": as_float(row.get("p50_strength_g_tex")),
                "pred_p50_yield_index": as_float(row.get("p50_yield_index")),
                "pred_p50_temporal_gap_days": as_float(row.get("p50_temporal_gap_days")),
            }
            for i, row in enumerate(final_top3, 1)
        ],
    }
    with open(path, "w") as f:
        json.dump(rules, f, indent=2)


def write_plan_markdown(
    path: Path,
    run_id: str,
    replicates: int,
    thresholds: dict[str, Any],
    candidate_rows: list[dict[str, Any]],
) -> None:
    lines: list[str] = []
    lines.append("# Lab Validation Plan — Final Top 3")
    lines.append("")
    lines.append(f"- Run ID: `{run_id}`")
    lines.append(f"- Generated (UTC): `{datetime.now(timezone.utc).replace(microsecond=0).isoformat()}`")
    lines.append(f"- Candidate replicates per arm: `{replicates}`")
    lines.append("")
    lines.append("## Objective")
    lines.append("")
    lines.append("Validate whether the simulated Top 3 preserve quality while reaching dark fiber under timing-safe conditions.")
    lines.append("")
    lines.append("## Phase 1 Pass/Fail Gate")
    lines.append("")
    lines.append(f"- `color_L <= {as_float(thresholds.get('max_color_L'), 25.0):.2f}`")
    lines.append(f"- `strength_g_tex >= {as_float(thresholds.get('min_strength_g_tex'), 28.0):.2f}`")
    lines.append(f"- `yield_index >= {as_float(thresholds.get('min_yield_index'), 0.85):.3f}`")
    lines.append(f"- `temporal_gap_days >= {as_float(thresholds.get('min_temporal_gap_days'), 0.0):.2f}`")
    if "max_toxicity_pre_cellulose" in thresholds:
        lines.append(
            f"- `toxicity_pre_cellulose <= {as_float(thresholds.get('max_toxicity_pre_cellulose'), 0.0):.3f}`"
        )
    lines.append("")
    lines.append("## Candidate Arms")
    lines.append("")
    lines.append("| Arm | Pred Success | Pred Robust | Pred p50 L* | Pred p50 Strength | Pred p50 Yield | Timing (mat/scw) | Key Params |")
    lines.append("|---|---:|---:|---:|---:|---:|---|---|")
    for row in candidate_rows:
        lines.append(
            f"| {row['arm_id']} | {as_float(row['pred_success_rate']):.3f} | {as_float(row['pred_robust_score']):.3f} | "
            f"{as_float(row['pred_p50_color_L']):.2f} | {as_float(row['pred_p50_strength_g_tex']):.2f} | "
            f"{as_float(row['pred_p50_yield_index']):.3f} | "
            f"{as_float(row['mat_activation_dpa']):.0f}/{as_float(row['scw_activation_dpa']):.0f} DPA | "
            f"k={as_float(row['k_competition']):.2f}, eff={as_float(row['melanin_efficiency']):.2f}, ret={as_float(row['late_retention_factor']):.2f} |"
        )
    lines.append("")
    lines.append("## Required Measurements")
    lines.append("")
    lines.append("- Expression timing: melA, TYRP1, DCT, cellulose markers across DPA stages")
    lines.append("- Pigment endpoint: melanin content and color L*")
    lines.append("- Quality endpoint: strength (g/tex), yield index")
    lines.append("- Durability endpoint: color L* shift after wash cycles")
    lines.append("")
    lines.append("## Generated Files")
    lines.append("")
    lines.append("- `results/lab_validation_matrix.csv`")
    lines.append("- `results/lab_measurements_template.csv`")
    lines.append("- `results/lab_pass_fail_rules.json`")
    lines.append("")

    path.write_text("\n".join(lines) + "\n")


def run(replicates: int, run_id: str) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    final_top3 = load_json(RESULTS_DIR / "final_lab_top3.json", [])
    if not final_top3:
        raise FileNotFoundError(
            "Missing results/final_lab_top3.json. Run `python -m src.adaptive_ode_robust_pipeline` first."
        )
    robust_summary = load_json(RESULTS_DIR / "adaptive_ode_robust_summary.json", {})
    thresholds = robust_summary.get("thresholds", {})

    candidate_rows = build_candidate_rows(final_top3=final_top3, replicates=replicates, run_id=run_id)
    control_rows = build_control_rows(replicates=replicates, run_id=run_id)
    all_arms = candidate_rows + control_rows
    template_rows = build_measurement_template(arms=all_arms, replicates=replicates, run_id=run_id)

    write_csv(RESULTS_DIR / "lab_validation_matrix.csv", all_arms)
    write_csv(RESULTS_DIR / "lab_measurements_template.csv", template_rows)
    write_pass_fail_rules(
        path=RESULTS_DIR / "lab_pass_fail_rules.json",
        thresholds=thresholds,
        run_id=run_id,
        final_top3=final_top3,
    )
    write_plan_markdown(
        path=RESULTS_DIR / "lab_validation_plan.md",
        run_id=run_id,
        replicates=replicates,
        thresholds=thresholds,
        candidate_rows=candidate_rows,
    )

    phase1_ok = evaluate_gate(
        {
            "color_L": as_float(final_top3[0].get("p50_color_L")),
            "strength_g_tex": as_float(final_top3[0].get("p50_strength_g_tex")),
            "yield_index": as_float(final_top3[0].get("p50_yield_index")),
            "temporal_gap_days": as_float(final_top3[0].get("p50_temporal_gap_days")),
        },
        thresholds=thresholds,
    )

    print("\n🧪 BlackCotton Lab Validation Planner")
    print("=" * 44)
    print(f"Run ID: {run_id}")
    print(f"Replicates per arm: {replicates}")
    print(f"Arms generated: {len(all_arms)} (candidates={len(candidate_rows)}, controls={len(control_rows)})")
    print(f"Measurement template rows: {len(template_rows)}")
    print(f"Lead candidate passes phase-1 gate by p50 prediction: {phase1_ok}")
    print("Saved:")
    print("  - results/lab_validation_plan.md")
    print("  - results/lab_validation_matrix.csv")
    print("  - results/lab_measurements_template.csv")
    print("  - results/lab_pass_fail_rules.json")
    print("=" * 44)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create first-cycle lab validation package from Top 3.")
    parser.add_argument(
        "--replicates",
        type=int,
        default=8,
        help="Biological replicates per arm (default: 8).",
    )
    parser.add_argument(
        "--run-id",
        default=f"BC-LAB-{datetime.now(timezone.utc).strftime('%Y-%m-%d')}",
        help="Run identifier used in generated files.",
    )
    args = parser.parse_args()
    run(replicates=max(1, int(args.replicates)), run_id=str(args.run_id))
