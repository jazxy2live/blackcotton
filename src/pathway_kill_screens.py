#!/usr/bin/env python3
"""
pathway_kill_screens.py — Kill-screen simulations for BlackCotton pathways.

Simulation 1 reuses the existing cotton expression + melanin ODEs and adds a
ROS/H2O2 proxy tied to the repo's toxicity threshold. Simulations 2 and 3 are
coarse-grained surrogate screens meant to answer kill questions quickly; they
are not atomistic MD or full bioreactor process models.
"""

from __future__ import annotations

import argparse
import copy
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from src.config_loader import load_config
from src.expression_model import run_simulation
from src.failure_risk_model import resolved_failure_risks
from src.melanin_pathway import create_enzyme_interpolator, michaelis_menten, run_melanin_simulation
from src.tradeoff_optimizer import build_candidate_params

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"

CHEMICAL_MELA_MULTIPLIERS = [round(x, 3) for x in np.linspace(0.5, 3.0, 21)]
CHEMICAL_SUPPORT_ARMS = [
    {
        "label": "balanced_current_pathway",
        "support_multiplier": 1.0,
        "description": "Current balanced melA + TYRP1 + DCT pathway.",
    },
    {
        "label": "support_starved_pathway",
        "support_multiplier": 0.10,
        "description": "melA-heavy pathway with weak downstream drain capacity.",
    },
]
CHEMICAL_H2O2_YIELD = 4.0
CHEMICAL_H2O2_CLEAR_BASE = 0.10
CHEMICAL_H2O2_CLEAR_ROS_COEF = 0.22
CHEMICAL_OPACITY_TARGET = 0.90

STRUCTURAL_MIN_STRENGTH_G_TEX = 28.0
STRUCTURAL_TOTAL_ABSORPTION_TARGET = 0.99
STRUCTURAL_NEAR_TOTAL_TARGETS = [0.90, 0.95, 0.97, 0.99]

BACTERIAL_TARGET_MELANIN_RATIO = 0.24
BACTERIAL_MIN_HBOND_RETENTION = 0.80
BACTERIAL_BASE_WET_STRENGTH = 22.0
BACTERIAL_BINDER_GRID = [round(x, 2) for x in np.linspace(0.0, 0.8, 9)]
BACTERIAL_RATIO_GRID = [round(x, 3) for x in np.linspace(0.08, 0.30, 12)]
BACTERIAL_ARCHITECTURES = [
    {
        "label": "co_culture",
        "co_localization": 0.55,
        "description": "Separate cellulose and melanin producers sharing one broth.",
    },
    {
        "label": "single_strain",
        "co_localization": 0.72,
        "description": "Single engineered Komagataeibacter strain with tighter co-localization.",
    },
]


def _clip(value: float, lo: float, hi: float) -> float:
    return float(np.clip(value, lo, hi))


def _logistic(value: float) -> float:
    return 1.0 / (1.0 + math.exp(-value))


def load_reference_candidate() -> tuple[dict[str, float] | None, str]:
    top3_path = RESULTS_DIR / "final_lab_top3.json"
    if not top3_path.exists():
        return None, "active config only (results/final_lab_top3.json missing)"

    rows = json.loads(top3_path.read_text())
    if not rows:
        return None, "active config only (results/final_lab_top3.json empty)"

    rank = 1
    worst_case_path = RESULTS_DIR / "adversarial_robustness_correlated_2026_03_07_candidate_summary.json"
    if worst_case_path.exists():
        candidate_rows = json.loads(worst_case_path.read_text())
        if candidate_rows:
            rank = max(1, int(candidate_rows[0].get("candidate_input_rank", 1)))

    idx = min(rank - 1, len(rows) - 1)
    chosen = rows[idx].get("base_candidate", rows[idx])
    return chosen, f"results/final_lab_top3.json candidate {idx + 1}"


def apply_reference_candidate(base_params: dict[str, Any], candidate: dict[str, float] | None) -> dict[str, Any]:
    if candidate is None:
        return copy.deepcopy(base_params)
    return build_candidate_params(base_params, candidate)


def chemical_arm_params(
    base_params: dict[str, Any],
    candidate: dict[str, float] | None,
    melA_multiplier: float,
    support_multiplier: float,
) -> dict[str, Any]:
    params = apply_reference_candidate(base_params, candidate)
    params["promoters"]["pGhMat1"]["strength_relative"] *= float(melA_multiplier)

    traffic = params.setdefault("melanin_pathway", {}).setdefault("traffic_control", {})
    traffic["melA_entry_scale"] = float(traffic.get("melA_entry_scale", 1.0)) * float(melA_multiplier)

    for key in ("pGhSCW_late", "pGhSCW_TYRP1", "pGhSCW_DCT"):
        if key in params["promoters"]:
            params["promoters"][key]["strength_relative"] *= float(support_multiplier)

    if "DCT" in params["melanin_pathway"]:
        params["melanin_pathway"]["DCT"]["Vmax"] *= float(support_multiplier)
    if "TYRP1" in params["melanin_pathway"]:
        params["melanin_pathway"]["TYRP1"]["Vmax"] *= float(support_multiplier)
    return params


def simulate_chemical_case(params: dict[str, Any]) -> dict[str, Any]:
    expr = run_odes_silent(run_simulation, params)
    mel = run_odes_silent(run_melanin_simulation, params, expr)

    enzyme_func = create_enzyme_interpolator(expr)
    risk = resolved_failure_risks(params)
    comp = compartmentalization_profile(params)
    threshold = float(toxicity_threshold(params))
    mp = params["melanin_pathway"]
    melA_params = mp["tyrosinase_melA"]

    t_hours = np.array(mel["t_hours"], dtype=float)
    dt = np.diff(t_hours, prepend=t_hours[0])
    h2o2_proxy = np.zeros_like(t_hours)
    ros_proxy = np.zeros_like(t_hours)
    melA_flux_v2 = np.zeros_like(t_hours)

    expression_scale = max(
        0.05,
        (1.0 - float(risk["silencing_probability"])) * math.exp(-0.5 * (float(risk["event_expression_cv"]) ** 2)),
    )
    melA_entry_scale = float(mp.get("traffic_control", {}).get("melA_entry_scale", 1.0))
    h2o2_clear = CHEMICAL_H2O2_CLEAR_BASE + CHEMICAL_H2O2_CLEAR_ROS_COEF * float(risk["ros_buffer_capacity"])

    for i, (t, dopa, dq, dc, iq) in enumerate(
        zip(
            t_hours,
            mel["L_DOPA"],
            mel["dopaquinone"],
            mel["dopachrome"],
            mel["indole_quinone"],
        )
    ):
        melA_level, _, _ = enzyme_func(float(t))
        melA_active = (
            float(melA_level)
            * float(risk["copper_loading_fraction"])
            * float(risk["tyrosinase_activation_fraction"])
            * expression_scale
            * melA_entry_scale
        )
        v2 = melA_active * michaelis_menten(
            float(dopa),
            float(melA_params["Vmax"]) * 0.8,
            float(melA_params["Km_tyrosine"]) * 0.5,
        )
        melA_flux_v2[i] = float(v2)
        if i > 0:
            h2o2_proxy[i] = max(
                0.0,
                h2o2_proxy[i - 1] + float(dt[i]) * (CHEMICAL_H2O2_YIELD * v2 - h2o2_clear * h2o2_proxy[i - 1]),
            )

        quinone_burden = float(comp["cytosolic_quinone_scale"]) * (
            float(dq) + 0.55 * float(dc) + 0.20 * float(iq)
        )
        ros_proxy[i] = h2o2_proxy[i] + quinone_burden

    final_melanin = max(float(mel["melanin"][-1]), 1e-9)
    opacity_fraction = np.array(mel["melanin"], dtype=float) / final_melanin
    opacity_90_idx = int(np.argmax(opacity_fraction >= CHEMICAL_OPACITY_TARGET))
    kill_indices = np.where(ros_proxy > threshold)[0]
    first_kill_idx = int(kill_indices[0]) if len(kill_indices) else None
    kill_before_opacity = bool(first_kill_idx is not None and first_kill_idx < opacity_90_idx)

    out = {
        "final_melanin_mM_equiv": float(final_melanin),
        "opacity_target_fraction": float(CHEMICAL_OPACITY_TARGET),
        "opacity_90_dpa": float(mel["t_dpa"][opacity_90_idx]),
        "peak_ros_proxy": float(np.max(ros_proxy)),
        "peak_ros_ratio": float(np.max(ros_proxy) / max(threshold, 1e-9)),
        "peak_h2o2_proxy": float(np.max(h2o2_proxy)),
        "peak_dopaquinone": float(np.max(mel["dopaquinone"])),
        "peak_dopachrome": float(np.max(mel["dopachrome"])),
        "peak_melA_flux_v2": float(np.max(melA_flux_v2)),
        "ros_threshold": float(threshold),
        "kill_before_opacity_90": bool(kill_before_opacity),
        "kill_dpa": None if first_kill_idx is None else float(mel["t_dpa"][first_kill_idx]),
        "kill_melanin_mM_equiv": None if first_kill_idx is None else float(mel["melanin"][first_kill_idx]),
        "kill_opacity_fraction": None if first_kill_idx is None else float(opacity_fraction[first_kill_idx]),
        "use_vacuolar_transit_peptides": bool(comp["use_vacuolar_transit_peptides"]),
        "cytosolic_quinone_scale": float(comp["cytosolic_quinone_scale"]),
    }
    return out


def summarize_chemical_arm(rows: list[dict[str, Any]]) -> dict[str, Any]:
    reference_row = next((row for row in rows if abs(float(row["melA_multiplier"]) - 1.0) < 1e-9), rows[0])
    safe_rows = [row for row in rows if not bool(row["kill_before_opacity_90"])]
    dead_rows = [row for row in rows if bool(row["kill_before_opacity_90"])]
    first_dead = dead_rows[0] if dead_rows else None
    safe_ceiling = safe_rows[-1] if safe_rows else None

    return {
        "reference_pass": not bool(reference_row["kill_before_opacity_90"]),
        "reference_peak_ros_ratio": float(reference_row["peak_ros_ratio"]),
        "reference_final_melanin_mM_equiv": float(reference_row["final_melanin_mM_equiv"]),
        "max_safe_melA_multiplier": None if safe_ceiling is None else float(safe_ceiling["melA_multiplier"]),
        "first_kill_melA_multiplier": None if first_dead is None else float(first_dead["melA_multiplier"]),
        "first_kill_dpa": None if first_dead is None else first_dead["kill_dpa"],
        "first_kill_melanin_mM_equiv": None if first_dead is None else first_dead["kill_melanin_mM_equiv"],
        "first_kill_opacity_fraction": None if first_dead is None else first_dead["kill_opacity_fraction"],
        "rows": rows,
    }


def run_chemical_screen(base_params: dict[str, Any], candidate: dict[str, float] | None) -> dict[str, Any]:
    arm_summaries = []
    for arm in CHEMICAL_SUPPORT_ARMS:
        rows = []
        for melA_multiplier in CHEMICAL_MELA_MULTIPLIERS:
            params = chemical_arm_params(
                base_params=base_params,
                candidate=candidate,
                melA_multiplier=float(melA_multiplier),
                support_multiplier=float(arm["support_multiplier"]),
            )
            metrics = simulate_chemical_case(params)
            rows.append(
                {
                    "arm_label": arm["label"],
                    "description": arm["description"],
                    "melA_multiplier": float(melA_multiplier),
                    "support_multiplier": float(arm["support_multiplier"]),
                    **metrics,
                }
            )
        arm_summaries.append(
            {
                "label": arm["label"],
                "description": arm["description"],
                "support_multiplier": float(arm["support_multiplier"]),
                **summarize_chemical_arm(rows),
            }
        )

    balanced = arm_summaries[0]
    return {
        "kill_metric": "ROS proxy exceeds antioxidant capacity before 90% of simulated final pigment load.",
        "current_pathway_dead": not bool(balanced["reference_pass"]),
        "arms": arm_summaries,
    }


def simulate_structural_case(base_strength_g_tex: float, porosity: float, jaggedness: float, branching: float) -> dict[str, float]:
    alignment = max(0.60, 1.0 - 0.08 * jaggedness - 0.28 * porosity - 0.10 * branching)
    crystallinity = max(0.55, 1.0 - 0.16 * porosity - 0.12 * branching - 0.04 * jaggedness)
    hydrogen_bond_retention = max(0.35, 1.0 - 0.40 * porosity - 0.12 * jaggedness - 0.15 * branching)
    tensile_strength = float(base_strength_g_tex * alignment * math.sqrt(crystallinity) * hydrogen_bond_retention)
    absorption = float(
        1.0 - math.exp(-(2.7 * porosity + 1.0 * jaggedness + 0.75 * branching + 0.9 * (1.0 - alignment)))
    )
    return {
        "porosity_fraction": float(porosity),
        "surface_jaggedness": float(jaggedness),
        "branching_factor": float(branching),
        "alignment_retention": float(alignment),
        "crystallinity_retention": float(crystallinity),
        "hydrogen_bond_retention": float(hydrogen_bond_retention),
        "absorption_fraction": float(absorption),
        "tensile_strength_g_tex": float(tensile_strength),
        "passes_strength_floor": bool(tensile_strength >= STRUCTURAL_MIN_STRENGTH_G_TEX),
    }


def run_structural_screen(base_params: dict[str, Any]) -> dict[str, Any]:
    base_strength = float(base_params["fiber_development"]["target_strength_g_tex"])
    rows = []
    for porosity in np.linspace(0.05, 0.70, 14):
        for jaggedness in np.linspace(0.0, 1.0, 11):
            for branching in np.linspace(0.0, 0.8, 9):
                rows.append(simulate_structural_case(base_strength, float(porosity), float(jaggedness), float(branching)))

    rows.sort(key=lambda row: (row["absorption_fraction"], row["tensile_strength_g_tex"]), reverse=True)
    safe_rows = [row for row in rows if bool(row["passes_strength_floor"])]
    feasible_total = [
        row
        for row in rows
        if row["absorption_fraction"] >= STRUCTURAL_TOTAL_ABSORPTION_TARGET and row["tensile_strength_g_tex"] >= STRUCTURAL_MIN_STRENGTH_G_TEX
    ]
    best_any = rows[0]
    best_safe = safe_rows[0] if safe_rows else None

    strength_at_targets = {}
    for target in STRUCTURAL_NEAR_TOTAL_TARGETS:
        eligible = [row["tensile_strength_g_tex"] for row in rows if row["absorption_fraction"] >= float(target)]
        strength_at_targets[f"{target:.2f}"] = None if not eligible else float(max(eligible))

    return {
        "kill_metric": "No morphology that reaches total absorption keeps tensile strength above 28 g/tex.",
        "structural_pathway_dead": len(feasible_total) == 0,
        "best_absorber": best_any,
        "best_safe_absorber": best_safe,
        "best_strength_by_absorption_target": strength_at_targets,
        "rows": rows[:20],
    }


def simulate_bacterial_case(architecture: dict[str, Any], melanin_ratio: float, binder_factor: float) -> dict[str, float]:
    co_localization = float(architecture["co_localization"])
    cellulose_scale = max(0.65, 1.0 - 0.12 * melanin_ratio)
    binding_affinity = _clip(0.20 + 0.55 * binder_factor + 0.40 * co_localization - 0.70 * melanin_ratio, 0.0, 1.0)
    bound_fraction = float(_logistic(7.0 * (binding_affinity - 0.48)))
    hbond_retention = _clip(
        1.0 - 1.25 * melanin_ratio * (1.0 - bound_fraction) - 0.22 * melanin_ratio * bound_fraction + 0.04 * binder_factor,
        0.0,
        1.05,
    )
    wet_strength = float(BACTERIAL_BASE_WET_STRENGTH * cellulose_scale * hbond_retention)
    return {
        "architecture": str(architecture["label"]),
        "melanin_to_cellulose_ratio": float(melanin_ratio),
        "binding_protein_factor": float(binder_factor),
        "co_localization": float(co_localization),
        "cellulose_scale": float(cellulose_scale),
        "binding_affinity": float(binding_affinity),
        "bound_fraction": float(bound_fraction),
        "hbond_retention": float(hbond_retention),
        "wet_matrix_strength": float(wet_strength),
        "passes_hbond_gate": bool(hbond_retention >= BACTERIAL_MIN_HBOND_RETENTION),
    }


def summarize_bacterial_architecture(rows: list[dict[str, Any]]) -> dict[str, Any]:
    target_rows = [row for row in rows if abs(float(row["melanin_to_cellulose_ratio"]) - BACTERIAL_TARGET_MELANIN_RATIO) < 1e-9]
    target_rows.sort(key=lambda row: float(row["binding_protein_factor"]))
    no_linker = next((row for row in target_rows if abs(float(row["binding_protein_factor"])) < 1e-9), target_rows[0])
    rescued = next((row for row in target_rows if bool(row["passes_hbond_gate"])), None)

    no_linker_pass_ratios = [
        float(row["melanin_to_cellulose_ratio"])
        for row in rows
        if abs(float(row["binding_protein_factor"])) < 1e-9 and bool(row["passes_hbond_gate"])
    ]

    return {
        "no_linker_target_row": no_linker,
        "required_binding_protein_factor": None if rescued is None else float(rescued["binding_protein_factor"]),
        "target_ratio_rescued": rescued is not None,
        "max_no_linker_safe_melanin_ratio": None if not no_linker_pass_ratios else float(max(no_linker_pass_ratios)),
    }


def run_bacterial_screen() -> dict[str, Any]:
    summaries = []
    for architecture in BACTERIAL_ARCHITECTURES:
        rows = []
        for ratio in BACTERIAL_RATIO_GRID:
            for binder in BACTERIAL_BINDER_GRID:
                rows.append(simulate_bacterial_case(architecture, float(ratio), float(binder)))
        summaries.append(
            {
                "label": str(architecture["label"]),
                "description": str(architecture["description"]),
                **summarize_bacterial_architecture(rows),
                "rows": rows,
            }
        )

    return {
        "kill_metric": "If hydrogen-bond retention falls below 0.80 at the target melanin load, the architecture needs binding-protein rescue.",
        "target_melanin_to_cellulose_ratio": float(BACTERIAL_TARGET_MELANIN_RATIO),
        "architectures": summaries,
    }


def build_report(summary: dict[str, Any]) -> str:
    chemical = summary["chemical"]
    structural = summary["structural"]
    bacterial = summary["bacterial"]

    lines: list[str] = []
    lines.append("# Pathway Kill Screens")
    lines.append("")
    lines.append("These screens answer three kill questions using the current BlackCotton model stack.")
    lines.append("")
    lines.append("## Reference")
    lines.append("")
    lines.append(f"- Config: `{summary['config_path']}`")
    lines.append(f"- Reference candidate source: `{summary['reference_candidate_source']}`")
    if summary["reference_candidate"] is not None:
        ref = summary["reference_candidate"]
        lines.append(
            f"- Reference candidate: mat/scw `{ref['mat_activation_dpa']:.0f}/{ref['scw_activation_dpa']:.0f}`, "
            f"k `{ref['k_competition']:.2f}`, eff `{ref['melanin_efficiency']:.2f}`, ret `{ref['late_retention_factor']:.2f}`"
        )
    lines.append("")
    lines.append("## Simulation 1: Cotton Toxicity Threshold")
    lines.append("")
    lines.append(f"- Kill metric: {chemical['kill_metric']}")
    lines.append(
        f"- Current pathway status: `{'DEAD' if chemical['current_pathway_dead'] else 'PASS'}`"
    )
    for arm in chemical["arms"]:
        lines.append(
            f"- `{arm['label']}`: reference ROS ratio `{arm['reference_peak_ros_ratio']:.3f}`, "
            f"safe through melA `{arm['max_safe_melA_multiplier']}`, first kill `{arm['first_kill_melA_multiplier']}`, "
            f"kill melanin `{arm['first_kill_melanin_mM_equiv']}`"
        )
    lines.append("")
    lines.append("## Simulation 2: CSC Metamaterial Mutation")
    lines.append("")
    lines.append(f"- Kill metric: {structural['kill_metric']}")
    lines.append(
        f"- Structural pathway status: `{'DEAD' if structural['structural_pathway_dead'] else 'PASS'}`"
    )
    best_any = structural["best_absorber"]
    lines.append(
        f"- Best absorber found: absorption `{best_any['absorption_fraction']:.3f}`, strength `{best_any['tensile_strength_g_tex']:.2f} g/tex`"
    )
    best_safe = structural["best_safe_absorber"]
    if best_safe is not None:
        lines.append(
            f"- Best safe absorber: absorption `{best_safe['absorption_fraction']:.3f}`, strength `{best_safe['tensile_strength_g_tex']:.2f} g/tex`"
        )
    for target, strength in structural["best_strength_by_absorption_target"].items():
        lines.append(f"- Best strength at absorption >= `{target}`: `{strength}`")
    lines.append("")
    lines.append("## Simulation 3: Bacterial Co-Culture Alternative")
    lines.append("")
    lines.append(f"- Kill metric: {bacterial['kill_metric']}")
    lines.append(f"- Target melanin/cellulose ratio: `{bacterial['target_melanin_to_cellulose_ratio']:.2f}`")
    for architecture in bacterial["architectures"]:
        no_linker = architecture["no_linker_target_row"]
        lines.append(
            f"- `{architecture['label']}` no-linker: affinity `{no_linker['binding_affinity']:.3f}`, "
            f"hbond `{no_linker['hbond_retention']:.3f}`, wet strength `{no_linker['wet_matrix_strength']:.2f}`"
        )
        lines.append(
            f"- `{architecture['label']}` rescue requirement: binding protein factor `{architecture['required_binding_protein_factor']}`"
        )
    lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append("- Simulation 1 is grounded in the existing ODE stack.")
    lines.append("- Simulations 2 and 3 are coarse-grained surrogate screens, not atomistic MD or full fermentation process models.")
    return "\n".join(lines) + "\n"


def run(config_path: str | Path | None, out_prefix: str) -> dict[str, Any]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    base_params = load_config(config_path)
    reference_candidate, reference_source = load_reference_candidate()

    chemical = run_chemical_screen(base_params, reference_candidate)
    structural = run_structural_screen(base_params)
    bacterial = run_bacterial_screen()

    summary = {
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "config_path": str(Path(config_path).resolve()) if config_path is not None else "active production config",
        "reference_candidate_source": reference_source,
        "reference_candidate": reference_candidate,
        "chemical": chemical,
        "structural": structural,
        "bacterial": bacterial,
    }

    summary_path = RESULTS_DIR / f"{out_prefix}_summary.json"
    report_path = RESULTS_DIR / f"{out_prefix}_report.md"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    report_path.write_text(build_report(summary))
    return {
        "summary_path": summary_path,
        "report_path": report_path,
        "summary": summary,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run BlackCotton pathway kill screens.")
    parser.add_argument(
        "--config",
        default=None,
        help="Optional config path. Defaults to active production config.",
    )
    parser.add_argument(
        "--out-prefix",
        default="pathway_kill_screens_2026_03_07",
        help="Output prefix under results/.",
    )
    args = parser.parse_args()

    output = run(config_path=args.config, out_prefix=str(args.out_prefix))
    summary = output["summary"]
    print("\nPathway kill screens complete")
    print("=" * 32)
    print(f"Chemical current pathway dead: {summary['chemical']['current_pathway_dead']}")
    print(f"Structural pathway dead: {summary['structural']['structural_pathway_dead']}")
    print(f"Saved: {output['summary_path'].relative_to(BASE_DIR)}")
    print(f"Saved: {output['report_path'].relative_to(BASE_DIR)}")
