#!/usr/bin/env python3
"""
plant_vs_bacterial_decision_report.py — Consolidate plant-path vs bacterial fallback results.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"


def load_json(path: Path) -> Any:
    with open(path) as f:
        return json.load(f)


def maybe_load_json(path: Path) -> Any | None:
    if not path.exists():
        return None
    return load_json(path)


def round_float(value: float | None, digits: int = 3) -> float | None:
    if value is None:
        return None
    return round(float(value), int(digits))


def build_payload(hardening_prefix: str | None) -> dict[str, Any]:
    optimization_summary = load_json(RESULTS_DIR / "optimization_summary.json")
    top_candidates = load_json(RESULTS_DIR / "top_candidates.json")
    chemical_validation = load_json(RESULTS_DIR / "chemical_gate_validation_2026_03_07.json")
    correlated_summary = load_json(RESULTS_DIR / "adversarial_robustness_correlated_2026_03_07_candidate_summary.json")
    final_lab_top3 = load_json(RESULTS_DIR / "final_lab_top3.json")
    bacterial_focus = load_json(RESULTS_DIR / "bacterial_alternative_single_strain_focus_2026_03_07_summary.json")

    hardening_summary = None
    hardening_top3 = None
    if hardening_prefix:
        hardening_summary = maybe_load_json(RESULTS_DIR / f"{hardening_prefix}_summary.json")
        hardening_top3 = maybe_load_json(RESULTS_DIR / f"{hardening_prefix}_top3.json")

    best_correlated = correlated_summary[0]
    reference_final_lab = final_lab_top3[1] if len(final_lab_top3) > 1 else final_lab_top3[0]
    best_opt_strength = max(float(row["strength_g_tex"]) for row in top_candidates)
    best_opt_yield = max(float(row["yield_index"]) for row in top_candidates)
    best_opt_ros = min(float(row["chemical_peak_ros_ratio"]) for row in top_candidates)

    bacterial_arch = bacterial_focus["architectures"][0]
    bacterial_rec = bacterial_arch["recommended_design"]
    bacterial_target_binder = bacterial_arch["required_binding_factor_at_target"]
    bacterial_best_strength = max(float(row["wet_matrix_strength"]) for row in bacterial_arch["top_safe_designs"])

    recommendation = {
        "primary_path": "plant_chemical",
        "secondary_path": "bacterial_single_strain_binder_engineering",
        "deprioritized_path": "csc_structural_metamaterial",
        "reason": (
            "Plant path still clears cotton-like strength with an active chemical gate, "
            "while the bacterial fallback remains materially below cotton spinning strength "
            "even after binder rescue."
        ),
    }

    payload: dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "recommendation": recommendation,
        "plant": {
            "reference_candidate": {
                "mat_activation_dpa": float(reference_final_lab["base_candidate"]["mat_activation_dpa"]),
                "scw_activation_dpa": float(reference_final_lab["base_candidate"]["scw_activation_dpa"]),
                "k_competition": float(reference_final_lab["base_candidate"]["k_competition"]),
                "melanin_efficiency": float(reference_final_lab["base_candidate"]["melanin_efficiency"]),
                "late_retention_factor": float(reference_final_lab["base_candidate"]["late_retention_factor"]),
                "success_rate": float(reference_final_lab["success_rate"]),
                "robust_score": float(reference_final_lab["robust_score"]),
                "strength_g_tex": float(reference_final_lab["base_candidate"]["strength_g_tex"]),
                "yield_index": float(reference_final_lab["base_candidate"]["yield_index"]),
            },
            "chemical_gate": {
                "all_refined_candidates_safe": int(optimization_summary["num_candidates_refined"])
                == int(optimization_summary["num_candidates_refined_chemical_safe"]),
                "num_refined_candidates": int(optimization_summary["num_candidates_refined"]),
                "num_refined_candidates_chemical_safe": int(optimization_summary["num_candidates_refined_chemical_safe"]),
                "best_chemical_peak_ros_ratio": float(optimization_summary["best_chemical_peak_ros_ratio"]),
                "reference_ros_ratio": float(chemical_validation["safe_current"]["peak_ros_ratio"]),
                "reference_kill_before_opacity_90": bool(chemical_validation["safe_current"]["kill_before_opacity_90"]),
                "unsafe_melA_multiplier": 2.25,
                "unsafe_kill_opacity_fraction": float(chemical_validation["unsafe_melA_2_25x"]["kill_opacity_fraction"]),
            },
            "worst_case_correlated": {
                "best_min_success_rate": float(best_correlated["min_success_rate"]),
                "best_mean_success_rate": float(best_correlated["mean_success_rate"]),
                "worst_scenario": str(best_correlated["worst_scenario"]),
                "worst_scenario_success_rate": float(best_correlated["worst_scenario_success_rate"]),
            },
            "current_search_frontier": {
                "best_strength_g_tex": best_opt_strength,
                "best_yield_index": best_opt_yield,
                "best_ros_ratio_seen": best_opt_ros,
            },
        },
        "bacterial": {
            "architecture": str(bacterial_arch["label"]),
            "target_ratio": float(bacterial_arch["target_ratio"]),
            "required_binding_factor_at_target": float(bacterial_target_binder),
            "recommended_design": {
                "ratio": float(bacterial_rec["melanin_to_cellulose_ratio"]),
                "binding_protein_factor": float(bacterial_rec["binding_protein_factor"]),
                "hbond_retention": float(bacterial_rec["hbond_retention"]),
                "wet_matrix_strength": float(bacterial_rec["wet_matrix_strength"]),
            },
            "best_safe_strength_in_focus_window": float(bacterial_best_strength),
        },
    }

    if hardening_summary is not None:
        top1 = (hardening_top3 or [None])[0]
        payload["plant"]["hardening"] = {
            "prefix": hardening_prefix,
            "baseline_strict_high_noise_top_success_rate": float(
                hardening_summary["baseline_strict_high_noise_top_success_rate"]
            ),
            "hardened_strict_high_noise_top_success_rate": float(
                hardening_summary["hardened_strict_high_noise_top_success_rate"]
            ),
            "delta_strict_high_noise_top_success_rate": float(
                hardening_summary["delta_strict_high_noise_top_success_rate"]
            ),
            "baseline_strict_high_noise_top_robust_score": float(
                hardening_summary["baseline_strict_high_noise_top_robust_score"]
            ),
            "hardened_strict_high_noise_top_robust_score": float(
                hardening_summary["hardened_strict_high_noise_top_robust_score"]
            ),
            "delta_strict_high_noise_top_robust_score": float(
                hardening_summary["delta_strict_high_noise_top_robust_score"]
            ),
            "seed_candidates_before_chemical_gate": int(hardening_summary["seed_candidates_before_chemical_gate"]),
            "seed_candidates_after_chemical_gate": int(hardening_summary["seed_candidates_after_chemical_gate"]),
            "target_met": bool(hardening_summary["target_met"]),
            "target_success_rate": float(hardening_summary["target_success_rate"]),
            "top1": top1,
            "top3": hardening_top3 or [],
        }

    return payload


def write_markdown(path: Path, payload: dict[str, Any]) -> None:
    plant = payload["plant"]
    bacterial = payload["bacterial"]
    recommendation = payload["recommendation"]

    lines: list[str] = []
    lines.append("# Plant vs Bacterial Decision Report")
    lines.append("")
    lines.append("Decision checkpoint for the BlackCotton primary path and fallback path.")
    lines.append("")
    lines.append("## Recommendation")
    lines.append("")
    lines.append(f"- Primary path: `{recommendation['primary_path']}`")
    lines.append(f"- Secondary path: `{recommendation['secondary_path']}`")
    lines.append(f"- Deprioritized path: `{recommendation['deprioritized_path']}`")
    lines.append(f"- Reason: {recommendation['reason']}")
    lines.append("")
    lines.append("## Plant Path")
    lines.append("")
    lines.append(
        f"- Reference candidate: `mat/scw {plant['reference_candidate']['mat_activation_dpa']:.0f}/"
        f"{plant['reference_candidate']['scw_activation_dpa']:.0f}`, "
        f"`k {plant['reference_candidate']['k_competition']:.2f}`, "
        f"`eff {plant['reference_candidate']['melanin_efficiency']:.2f}`, "
        f"`ret {plant['reference_candidate']['late_retention_factor']:.2f}`"
    )
    lines.append(
        f"- Reference performance: success `{plant['reference_candidate']['success_rate']:.3f}`, "
        f"robust `{plant['reference_candidate']['robust_score']:.3f}`, "
        f"strength `{plant['reference_candidate']['strength_g_tex']:.2f} g/tex`, "
        f"yield `{plant['reference_candidate']['yield_index']:.3f}`"
    )
    lines.append(
        f"- Chemical gate: `{plant['chemical_gate']['num_refined_candidates_chemical_safe']}` / "
        f"`{plant['chemical_gate']['num_refined_candidates']}` refined survivors remain chemically safe"
    )
    lines.append(
        f"- Chemical ROS: best refined ratio `{plant['chemical_gate']['best_chemical_peak_ros_ratio']:.3f}`, "
        f"reference ratio `{plant['chemical_gate']['reference_ros_ratio']:.3f}`, "
        f"first forced-kill point `melA {plant['chemical_gate']['unsafe_melA_multiplier']:.2f}x`"
    )
    lines.append(
        f"- Worst correlated scenario: `{plant['worst_case_correlated']['worst_scenario']}` with "
        f"best min success `{plant['worst_case_correlated']['best_min_success_rate']:.3f}`"
    )
    lines.append(
        f"- Current search frontier: strength up to `{plant['current_search_frontier']['best_strength_g_tex']:.2f} g/tex`, "
        f"yield up to `{plant['current_search_frontier']['best_yield_index']:.3f}`, "
        f"ROS ratio down to `{plant['current_search_frontier']['best_ros_ratio_seen']:.3f}`"
    )
    if "hardening" in plant:
        hardening = plant["hardening"]
        lines.append(
            f"- Chemical-gated hardening: strict+high-noise top success "
            f"`{hardening['baseline_strict_high_noise_top_success_rate']:.3f}` -> "
            f"`{hardening['hardened_strict_high_noise_top_success_rate']:.3f}` "
            f"(delta `{hardening['delta_strict_high_noise_top_success_rate']:+.3f}`)"
        )
        lines.append(
            f"- Hardening target status: `{hardening['target_met']}` "
            f"against target `{hardening['target_success_rate']:.3f}`"
        )
        lines.append(
            f"- Chemical-safe seed retention during hardening: `{hardening['seed_candidates_after_chemical_gate']}` / "
            f"`{hardening['seed_candidates_before_chemical_gate']}`"
        )
        if hardening["top1"] is not None:
            top1 = hardening["top1"]["base_candidate"]
            lines.append(
                f"- Best hardened candidate: `mat/scw {top1['mat_activation_dpa']:.0f}/{top1['scw_activation_dpa']:.0f}`, "
                f"`k {top1['k_competition']:.2f}`, `eff {top1['melanin_efficiency']:.2f}`, "
                f"`ret {top1['late_retention_factor']:.2f}`"
            )
    lines.append("")
    lines.append("## Bacterial Fallback")
    lines.append("")
    lines.append(f"- Architecture: `{bacterial['architecture']}`")
    lines.append(
        f"- Required binder at target ratio `{bacterial['target_ratio']:.2f}`: "
        f"`{bacterial['required_binding_factor_at_target']:.2f}`"
    )
    lines.append(
        f"- Best focused-window design: ratio `{bacterial['recommended_design']['ratio']:.2f}`, "
        f"binder `{bacterial['recommended_design']['binding_protein_factor']:.2f}`, "
        f"hbond `{bacterial['recommended_design']['hbond_retention']:.3f}`, "
        f"wet strength `{bacterial['recommended_design']['wet_matrix_strength']:.2f}`"
    )
    lines.append(
        f"- Best safe strength seen in focused window: `{bacterial['best_safe_strength_in_focus_window']:.2f}`"
    )
    lines.append("")
    lines.append("## Decision")
    lines.append("")
    lines.append("- Keep the plant chemical path as the main program.")
    lines.append("- Keep the bacterial path as a backup branch focused on binder engineering, not as the lead path.")
    lines.append("- Keep the CSC structural path deprioritized under the current screen.")
    lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append("- Plant chemical safety is tied to the current ODE-based ROS proxy and enforced as a hard gate.")
    lines.append("- Bacterial results are coarse-grained surrogate outputs, not fermentation-scale wet-lab predictions.")

    path.write_text("\n".join(lines) + "\n")


def run(out_prefix: str, hardening_prefix: str | None) -> None:
    payload = build_payload(hardening_prefix)
    summary_path = RESULTS_DIR / f"{out_prefix}_summary.json"
    report_path = RESULTS_DIR / f"{out_prefix}_report.md"

    summary_path.write_text(json.dumps(payload, indent=2) + "\n")
    write_markdown(report_path, payload)

    print("Plant vs bacterial decision report complete")
    print(f"Saved: results/{summary_path.name}")
    print(f"Saved: results/{report_path.name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a plant-vs-bacterial decision report.")
    parser.add_argument(
        "--out-prefix",
        default="plant_vs_bacterial_decision_2026_03_07",
        help="Output prefix under results/.",
    )
    parser.add_argument(
        "--hardening-prefix",
        default=None,
        help="Optional hardening output prefix to include in the report.",
    )
    args = parser.parse_args()
    run(
        out_prefix=str(args.out_prefix),
        hardening_prefix=None if args.hardening_prefix in (None, "", "none") else str(args.hardening_prefix),
    )
