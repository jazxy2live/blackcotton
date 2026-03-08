#!/usr/bin/env python3
"""
pathway_safety_models.py — Shared kill-screen safety models.

These are coarse-grained safety models used to screen pathway concepts.
The cotton chemical model reuses the existing ODE outputs; the bacterial
model is a surrogate for broth-phase cellulose/melanin compatibility.
"""

from __future__ import annotations

import io
import math
from contextlib import redirect_stdout
from typing import Any

import numpy as np

from src.expression_model import run_simulation
from src.failure_risk_model import resolved_failure_risks
from src.melanin_pathway import create_enzyme_interpolator, michaelis_menten, run_melanin_simulation

CHEMICAL_H2O2_YIELD = 4.0
CHEMICAL_H2O2_CLEAR_BASE = 0.10
CHEMICAL_H2O2_CLEAR_ROS_COEF = 0.22
CHEMICAL_OPACITY_TARGET = 0.90

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


def _run_silent(func, *args, **kwargs):
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        return func(*args, **kwargs)


def compartmentalization_profile(params: dict[str, Any]) -> dict[str, float]:
    mp = params.get("melanin_pathway", {})
    comp = mp.get("compartmentalization", {})
    construct_cfg = params.get("construct", {})
    has_transit_targeting = bool(construct_cfg.get("use_vacuolar_transit_peptides", False))

    if has_transit_targeting:
        sequestration = float(np.clip(comp.get("vacuolar_sequestration_fraction", 0.0), 0.0, 0.95))
        leak = float(np.clip(comp.get("cytosolic_quinone_leak_fraction", 1.0), 0.05, 1.0))
    else:
        sequestration = 0.0
        leak = 1.0
    cytosolic_scale = float(np.clip((1.0 - sequestration) * leak, 0.02, 1.0))
    return {
        "use_vacuolar_transit_peptides": bool(has_transit_targeting),
        "vacuolar_sequestration_fraction": sequestration,
        "cytosolic_quinone_leak_fraction": leak,
        "cytosolic_quinone_scale": cytosolic_scale,
    }


def chemical_ros_threshold(params: dict[str, Any]) -> float:
    risk = resolved_failure_risks(params)
    return float(max(1e-6, risk["quinone_stress_threshold"] * (0.70 + 0.60 * risk["ros_buffer_capacity"])))


def simulate_chemical_ros_safety(
    params: dict[str, Any],
    expression_results: dict[str, Any] | None = None,
    melanin_results: dict[str, Any] | None = None,
) -> dict[str, Any]:
    expr = expression_results if expression_results is not None else _run_silent(run_simulation, params)
    mel = melanin_results if melanin_results is not None else _run_silent(run_melanin_simulation, params, expr)

    enzyme_func = create_enzyme_interpolator(expr)
    risk = resolved_failure_risks(params)
    comp = compartmentalization_profile(params)
    threshold = float(chemical_ros_threshold(params))
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
    h2o2_clear *= float(params.get("failure_risks", {}).get("catalase_coexpression_multiplier", 1.0))

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
        # Quasi-steady state approximation for highly reactive H2O2 sink avoids numerical stiffness
        h2o2_proxy[i] = (CHEMICAL_H2O2_YIELD * v2) / max(1e-6, h2o2_clear)

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

    return {
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


def sweep_bacterial_designs() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for architecture in BACTERIAL_ARCHITECTURES:
        for ratio in BACTERIAL_RATIO_GRID:
            for binder in BACTERIAL_BINDER_GRID:
                rows.append(simulate_bacterial_case(architecture, float(ratio), float(binder)))
    return rows
