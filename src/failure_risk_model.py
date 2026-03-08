#!/usr/bin/env python3
"""
failure_risk_model.py — Shared failure-risk resolution for BlackCotton.
"""

from typing import Any

import numpy as np

DEFAULT_FAILURE_RISKS = {
    "copper_loading_fraction": 1.0,
    "tyrosinase_activation_fraction": 1.0,
    "quinone_stress_threshold": 0.08,
    "ros_buffer_capacity": 1.0,
    "silencing_probability": 0.0,
    "event_expression_cv": 0.0,
}

DEFAULT_STABILITY_DESIGN = {
    "insulator_effectiveness": 0.0,
    "intron_silencing_relief": 0.0,
    "matrix_attachment_relief": 0.0,
    "repeat_penalty": 0.0,
}


def bounded(value: float, lo: float, hi: float) -> float:
    return float(np.clip(value, lo, hi))


def construct_stability_design(params: dict[str, Any]) -> dict[str, float]:
    construct_cfg = params.get("construct", {})
    cfg = dict(DEFAULT_STABILITY_DESIGN)
    design_cfg = construct_cfg.get("stability_design", {})
    if isinstance(design_cfg, dict):
        cfg.update(design_cfg)
    return {
        "insulator_effectiveness": bounded(float(cfg["insulator_effectiveness"]), 0.0, 0.8),
        "intron_silencing_relief": bounded(float(cfg["intron_silencing_relief"]), 0.0, 0.8),
        "matrix_attachment_relief": bounded(float(cfg["matrix_attachment_relief"]), 0.0, 0.8),
        "repeat_penalty": bounded(float(cfg["repeat_penalty"]), 0.0, 0.8),
    }


def stability_design_metadata(params: dict[str, Any]) -> dict[str, Any]:
    design = construct_stability_design(params)
    silencing_scale = float(np.clip(
        1.0
        - 0.55 * design["insulator_effectiveness"]
        - 0.60 * design["intron_silencing_relief"]
        - 0.45 * design["matrix_attachment_relief"]
        + 0.90 * design["repeat_penalty"],
        0.10,
        1.60,
    ))
    event_cv_scale = float(np.clip(
        1.0
        - 0.28 * design["insulator_effectiveness"]
        - 0.18 * design["intron_silencing_relief"]
        - 0.40 * design["matrix_attachment_relief"]
        + 0.70 * design["repeat_penalty"],
        0.20,
        1.60,
    ))
    return {
        "design": design,
        "silencing_scale": silencing_scale,
        "event_expression_cv_scale": event_cv_scale,
    }


def resolved_failure_risks(params: dict[str, Any]) -> dict[str, float]:
    risk = dict(DEFAULT_FAILURE_RISKS)
    raw_cfg = params.get("failure_risks", {})
    if isinstance(raw_cfg, dict):
        risk.update(raw_cfg)

    meta = stability_design_metadata(params)
    out = {
        "copper_loading_fraction": bounded(float(risk["copper_loading_fraction"]), 0.0, 1.0),
        "tyrosinase_activation_fraction": bounded(float(risk["tyrosinase_activation_fraction"]), 0.0, 1.0),
        "quinone_stress_threshold": float(max(float(risk["quinone_stress_threshold"]), 1e-6)),
        "ros_buffer_capacity": bounded(float(risk["ros_buffer_capacity"]), 0.0, 1.5),
        "silencing_probability": bounded(
            float(risk["silencing_probability"]) * float(meta["silencing_scale"]),
            0.0,
            0.80,
        ),
        "event_expression_cv": bounded(
            float(risk["event_expression_cv"]) * float(meta["event_expression_cv_scale"]),
            0.0,
            0.80,
        ),
    }
    return out
