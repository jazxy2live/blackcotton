#!/usr/bin/env python3
"""
transcriptome_calibrator.py — Fit promoter parameters to expression reference data
==================================================================================

Calibrates pGhSCW_late and pGhMat1 promoter parameters against a reference
time-series table (DPA vs normalized expression markers), then writes a new
calibrated parameter file and fit diagnostics.

Usage:
    python -m src.transcriptome_calibrator
"""

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import yaml
from scipy.optimize import differential_evolution

from src.expression_model import load_params, promoter_activity, run_simulation

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data" / "reference"
CONFIG_DIR = BASE_DIR / "config"
RESULTS_DIR = BASE_DIR / "results"


def clamp(x: float, lo: float, hi: float) -> float:
    return float(np.clip(x, lo, hi))


def read_reference_table(path: Path) -> dict[str, np.ndarray]:
    if not path.exists():
        raise FileNotFoundError(f"Reference expression table not found: {path}")

    dpa = []
    scw = []
    mat = []
    weight = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        required = {"dpa", "scw_marker_norm", "maturation_marker_norm"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(
                f"Reference table must include columns {sorted(required)}; "
                f"found {reader.fieldnames}"
            )
        for row in reader:
            dpa.append(float(row["dpa"]))
            scw.append(clamp(float(row["scw_marker_norm"]), 0.0, 1.2))
            mat.append(clamp(float(row["maturation_marker_norm"]), 0.0, 1.2))
            w = row.get("weight")
            weight.append(float(w) if w not in (None, "") else 1.0)

    order = np.argsort(np.array(dpa, dtype=float))
    return {
        "dpa": np.array(dpa, dtype=float)[order],
        "scw": np.array(scw, dtype=float)[order],
        "mat": np.array(mat, dtype=float)[order],
        "weight": np.array(weight, dtype=float)[order],
    }


def unpack_promoter_vector(vector: np.ndarray) -> tuple[float, float, float, float, float]:
    activation = float(vector[0])
    delta_peak = max(float(vector[1]), 1.0)
    peak = activation + delta_peak
    hill = max(float(vector[2]), 1.0)
    leakage = clamp(float(vector[3]), 0.0, 0.30)
    strength = clamp(float(vector[4]), 0.4, 2.0)
    return activation, peak, hill, leakage, strength


def predict_promoter_curve(dpa: np.ndarray, vector: np.ndarray) -> np.ndarray:
    activation, peak, hill, leakage, strength = unpack_promoter_vector(vector)
    return np.array(
        [
            promoter_activity(
                t_dpa=float(t),
                activation_dpa=activation,
                peak_dpa=peak,
                hill_n=hill,
                leakage=leakage,
                strength=strength,
            )
            for t in dpa
        ],
        dtype=float,
    )


def weighted_rmse(y_true: np.ndarray, y_pred: np.ndarray, weight: np.ndarray) -> float:
    w = np.array(weight, dtype=float)
    if np.allclose(w.sum(), 0.0):
        w = np.ones_like(w)
    w = w / w.sum()
    err2 = (y_true - y_pred) ** 2
    return float(np.sqrt(np.sum(w * err2)))


def fit_single_promoter(
    dpa: np.ndarray,
    observed: np.ndarray,
    weight: np.ndarray,
    bounds: list[tuple[float, float]],
    seed: int,
    maxiter: int = 140,
) -> dict[str, Any]:
    def objective(vec: np.ndarray) -> float:
        pred = predict_promoter_curve(dpa, vec)
        rmse = weighted_rmse(observed, pred, weight)
        # Soft prior against overly flat curves.
        _, _, hill, leakage, _ = unpack_promoter_vector(vec)
        slope_penalty = 0.01 * max(0.0, 2.0 - hill)
        leak_penalty = 0.05 * max(0.0, leakage - 0.12)
        return rmse + slope_penalty + leak_penalty

    result = differential_evolution(
        objective,
        bounds=bounds,
        seed=seed,
        maxiter=maxiter,
        popsize=12,
        polish=True,
        tol=1e-5,
    )
    vec = np.array(result.x, dtype=float)
    pred = predict_promoter_curve(dpa, vec)
    rmse = weighted_rmse(observed, pred, weight)
    activation, peak, hill, leakage, strength = unpack_promoter_vector(vec)
    return {
        "activation_dpa": float(activation),
        "peak_dpa": float(peak),
        "hill_coefficient": float(hill),
        "leakage_fraction": float(leakage),
        "strength_relative": float(strength),
        "rmse": float(rmse),
        "predicted_curve": pred.tolist(),
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
        "optimizer_fun": float(result.fun),
    }


def promoter_vector_from_params(promoter_params: dict[str, Any]) -> np.ndarray:
    activation = float(promoter_params["activation_dpa"])
    peak = float(promoter_params["peak_dpa"])
    delta = max(1.0, peak - activation)
    return np.array(
        [
            activation,
            delta,
            float(promoter_params["hill_coefficient"]),
            float(promoter_params["leakage_fraction"]),
            float(promoter_params["strength_relative"]),
        ],
        dtype=float,
    )


def apply_calibration_to_params(params: dict[str, Any], scw_fit: dict[str, Any], mat_fit: dict[str, Any]) -> dict[str, Any]:
    out = json.loads(json.dumps(params))
    out["promoters"]["pGhSCW_late"]["activation_dpa"] = float(scw_fit["activation_dpa"])
    out["promoters"]["pGhSCW_late"]["peak_dpa"] = float(scw_fit["peak_dpa"])
    out["promoters"]["pGhSCW_late"]["hill_coefficient"] = float(scw_fit["hill_coefficient"])
    out["promoters"]["pGhSCW_late"]["leakage_fraction"] = float(scw_fit["leakage_fraction"])
    out["promoters"]["pGhSCW_late"]["strength_relative"] = float(scw_fit["strength_relative"])

    out["promoters"]["pGhMat1"]["activation_dpa"] = float(mat_fit["activation_dpa"])
    out["promoters"]["pGhMat1"]["peak_dpa"] = float(mat_fit["peak_dpa"])
    out["promoters"]["pGhMat1"]["hill_coefficient"] = float(mat_fit["hill_coefficient"])
    out["promoters"]["pGhMat1"]["leakage_fraction"] = float(mat_fit["leakage_fraction"])
    out["promoters"]["pGhMat1"]["strength_relative"] = float(mat_fit["strength_relative"])
    return out


def expression_snapshot(results: dict[str, np.ndarray]) -> dict[str, float]:
    t = np.array(results["t_dpa"], dtype=float)
    mel = np.array(results["protein_melA"], dtype=float)
    scw = np.array(results["protein_TYRP1"], dtype=float)
    cel = np.array(results["cellulose"], dtype=float)
    mel_half_idx = int(np.argmax(mel >= 0.5 * np.max(mel))) if np.max(mel) > 0 else 0
    cel_90_idx = int(np.argmax(cel >= 0.9 * np.max(cel))) if np.max(cel) > 0 else 0
    return {
        "melA_peak_dpa": float(t[int(np.argmax(mel))]) if mel.size else 0.0,
        "TYRP1_peak_dpa": float(t[int(np.argmax(scw))]) if scw.size else 0.0,
        "cellulose_90pct_dpa": float(t[cel_90_idx]) if cel.size else 0.0,
        "melA_halfmax_dpa": float(t[mel_half_idx]) if mel.size else 0.0,
        "temporal_gap_days": float(t[mel_half_idx] - t[cel_90_idx]) if mel.size and cel.size else 0.0,
    }


def enforce_min_temporal_gap(
    params: dict[str, Any],
    current_gap: float,
    min_gap: float,
) -> tuple[dict[str, Any], dict[str, Any]]:
    if current_gap >= min_gap:
        return params, {"applied": False, "shift_days": 0.0}

    delta = float(min_gap - current_gap)
    out = json.loads(json.dumps(params))
    total_days = float(out["fiber_development"]["total_development_days"])
    mat = out["promoters"]["pGhMat1"]
    mat["activation_dpa"] = min(total_days - 1.0, float(mat["activation_dpa"]) + delta)
    mat["peak_dpa"] = min(total_days, float(mat["peak_dpa"]) + delta)
    return out, {"applied": True, "shift_days": delta}


def save_fit_curves_csv(
    path: Path,
    dpa: np.ndarray,
    observed_scw: np.ndarray,
    observed_mat: np.ndarray,
    baseline_scw: np.ndarray,
    baseline_mat: np.ndarray,
    calibrated_scw: np.ndarray,
    calibrated_mat: np.ndarray,
) -> None:
    rows = []
    for i in range(len(dpa)):
        rows.append(
            {
                "dpa": float(dpa[i]),
                "observed_scw_marker_norm": float(observed_scw[i]),
                "baseline_scw_fit": float(baseline_scw[i]),
                "calibrated_scw_fit": float(calibrated_scw[i]),
                "observed_maturation_marker_norm": float(observed_mat[i]),
                "baseline_maturation_fit": float(baseline_mat[i]),
                "calibrated_maturation_fit": float(calibrated_mat[i]),
            }
        )
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def run_calibration(
    reference_path: Path,
    output_config_path: Path,
    out_prefix: str,
    seed: int,
    min_temporal_gap_days: float,
) -> dict[str, Any]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output_config_path.parent.mkdir(parents=True, exist_ok=True)

    params = load_params()
    table = read_reference_table(reference_path)
    dpa = table["dpa"]
    obs_scw = table["scw"]
    obs_mat = table["mat"]
    weight = table["weight"]

    base_scw_vec = promoter_vector_from_params(params["promoters"]["pGhSCW_late"])
    base_mat_vec = promoter_vector_from_params(params["promoters"]["pGhMat1"])

    baseline_scw = predict_promoter_curve(dpa, base_scw_vec)
    baseline_mat = predict_promoter_curve(dpa, base_mat_vec)
    baseline_scw_rmse = weighted_rmse(obs_scw, baseline_scw, weight)
    baseline_mat_rmse = weighted_rmse(obs_mat, baseline_mat, weight)

    scw_bounds = [
        (20.0, 40.0),   # activation
        (1.0, 14.0),    # peak delta
        (1.0, 12.0),    # hill
        (0.0, 0.20),    # leakage
        (0.6, 1.8),     # strength
    ]
    mat_bounds = [
        (25.0, 48.0),   # activation
        (1.0, 14.0),    # peak delta
        (1.0, 12.0),    # hill
        (0.0, 0.20),    # leakage
        (0.6, 1.8),     # strength
    ]

    scw_fit = fit_single_promoter(dpa, obs_scw, weight, bounds=scw_bounds, seed=seed + 1)
    mat_fit = fit_single_promoter(dpa, obs_mat, weight, bounds=mat_bounds, seed=seed + 2)
    cal_scw = np.array(scw_fit["predicted_curve"], dtype=float)
    cal_mat = np.array(mat_fit["predicted_curve"], dtype=float)

    raw_calibrated_params = apply_calibration_to_params(params, scw_fit=scw_fit, mat_fit=mat_fit)

    baseline_expr = run_simulation(params)
    raw_calibrated_expr = run_simulation(raw_calibrated_params)
    baseline_snapshot = expression_snapshot(baseline_expr)
    raw_calibrated_snapshot = expression_snapshot(raw_calibrated_expr)

    final_calibrated_params, safety_adjustment = enforce_min_temporal_gap(
        params=raw_calibrated_params,
        current_gap=float(raw_calibrated_snapshot["temporal_gap_days"]),
        min_gap=float(min_temporal_gap_days),
    )
    final_calibrated_expr = run_simulation(final_calibrated_params)
    final_calibrated_snapshot = expression_snapshot(final_calibrated_expr)
    total_shift = float(safety_adjustment["shift_days"])
    # Coarse 1-hour sampling can under-shoot target gap by one bin; refine iteratively.
    for _ in range(12):
        gap = float(final_calibrated_snapshot["temporal_gap_days"])
        if gap + 1e-9 >= float(min_temporal_gap_days):
            break
        extra = float(min_temporal_gap_days - gap + (1.0 / 24.0))
        total_shift += extra
        mat = final_calibrated_params["promoters"]["pGhMat1"]
        total_days = float(final_calibrated_params["fiber_development"]["total_development_days"])
        mat["activation_dpa"] = min(total_days - 1.0, float(mat["activation_dpa"]) + extra)
        mat["peak_dpa"] = min(total_days, float(mat["peak_dpa"]) + extra)
        final_calibrated_expr = run_simulation(final_calibrated_params)
        final_calibrated_snapshot = expression_snapshot(final_calibrated_expr)
    safety_adjustment = {
        "applied": bool(total_shift > 0.0),
        "shift_days": float(total_shift),
        "achieved_temporal_gap_days": float(final_calibrated_snapshot["temporal_gap_days"]),
    }

    with open(output_config_path, "w") as f:
        yaml.safe_dump(final_calibrated_params, f, sort_keys=False)

    curves_path = RESULTS_DIR / f"{out_prefix}_curves.csv"
    save_fit_curves_csv(
        curves_path,
        dpa=dpa,
        observed_scw=obs_scw,
        observed_mat=obs_mat,
        baseline_scw=baseline_scw,
        baseline_mat=baseline_mat,
        calibrated_scw=cal_scw,
        calibrated_mat=cal_mat,
    )

    summary = {
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "reference_table": str(reference_path.relative_to(BASE_DIR)),
        "output_config": str(output_config_path.relative_to(BASE_DIR)),
        "baseline_rmse": {
            "pGhSCW_late": float(baseline_scw_rmse),
            "pGhMat1": float(baseline_mat_rmse),
            "mean": float((baseline_scw_rmse + baseline_mat_rmse) / 2.0),
        },
        "calibrated_rmse": {
            "pGhSCW_late": float(scw_fit["rmse"]),
            "pGhMat1": float(mat_fit["rmse"]),
            "mean": float((scw_fit["rmse"] + mat_fit["rmse"]) / 2.0),
        },
        "promoter_baseline": {
            "pGhSCW_late": {
                k: float(params["promoters"]["pGhSCW_late"][k])
                for k in ("activation_dpa", "peak_dpa", "hill_coefficient", "leakage_fraction", "strength_relative")
            },
            "pGhMat1": {
                k: float(params["promoters"]["pGhMat1"][k])
                for k in ("activation_dpa", "peak_dpa", "hill_coefficient", "leakage_fraction", "strength_relative")
            },
        },
        "promoter_calibrated": {
            "pGhSCW_late": {k: float(scw_fit[k]) for k in ("activation_dpa", "peak_dpa", "hill_coefficient", "leakage_fraction", "strength_relative")},
            "pGhMat1": {k: float(mat_fit[k]) for k in ("activation_dpa", "peak_dpa", "hill_coefficient", "leakage_fraction", "strength_relative")},
        },
        "expression_snapshot_baseline": baseline_snapshot,
        "expression_snapshot_raw_calibrated": raw_calibrated_snapshot,
        "expression_snapshot_calibrated": final_calibrated_snapshot,
        "safety_adjustment": {
            "min_temporal_gap_days": float(min_temporal_gap_days),
            **safety_adjustment,
        },
        "improvements": {
            "scw_rmse_delta": float(scw_fit["rmse"] - baseline_scw_rmse),
            "mat_rmse_delta": float(mat_fit["rmse"] - baseline_mat_rmse),
            "mean_rmse_delta": float(((scw_fit["rmse"] + mat_fit["rmse"]) / 2.0) - ((baseline_scw_rmse + baseline_mat_rmse) / 2.0)),
        },
        "files": {
            "curves_csv": str(curves_path.relative_to(BASE_DIR)),
        },
    }

    summary_path = RESULTS_DIR / f"{out_prefix}_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    expr_cmp_path = RESULTS_DIR / f"{out_prefix}_expression_comparison.json"
    with open(expr_cmp_path, "w") as f:
        json.dump(
            {
                "baseline": baseline_snapshot,
                "raw_calibrated": raw_calibrated_snapshot,
                "calibrated": final_calibrated_snapshot,
                "safety_adjustment": summary["safety_adjustment"],
            },
            f,
            indent=2,
        )

    report_lines = []
    report_lines.append("# Transcriptome Calibration Report")
    report_lines.append("")
    report_lines.append(f"- Reference table: `{summary['reference_table']}`")
    report_lines.append(f"- Output calibrated config: `{summary['output_config']}`")
    report_lines.append("")
    report_lines.append("## Fit Quality")
    report_lines.append("")
    report_lines.append(f"- Baseline mean RMSE: `{summary['baseline_rmse']['mean']:.4f}`")
    report_lines.append(f"- Calibrated mean RMSE: `{summary['calibrated_rmse']['mean']:.4f}`")
    report_lines.append(f"- Delta (cal - base): `{summary['improvements']['mean_rmse_delta']:.4f}`")
    report_lines.append("")
    report_lines.append("## Promoter Updates")
    report_lines.append("")
    report_lines.append("| Promoter | Activation DPA | Peak DPA | Hill | Leakage | Strength |")
    report_lines.append("|---|---:|---:|---:|---:|---:|")
    for name in ("pGhSCW_late", "pGhMat1"):
        row = summary["promoter_calibrated"][name]
        report_lines.append(
            f"| {name} | {row['activation_dpa']:.2f} | {row['peak_dpa']:.2f} | "
            f"{row['hill_coefficient']:.2f} | {row['leakage_fraction']:.4f} | {row['strength_relative']:.3f} |"
        )
    report_lines.append("")
    report_lines.append("## Expression Timing Snapshot")
    report_lines.append("")
    report_lines.append(
        f"- Baseline temporal gap: `{baseline_snapshot['temporal_gap_days']:.3f}` days"
    )
    report_lines.append(
        f"- Raw calibrated temporal gap: `{raw_calibrated_snapshot['temporal_gap_days']:.3f}` days"
    )
    report_lines.append(
        f"- Final calibrated temporal gap: `{final_calibrated_snapshot['temporal_gap_days']:.3f}` days "
        f"(min required `{min_temporal_gap_days:.3f}`)"
    )
    report_lines.append(
        f"- Baseline melA peak DPA: `{baseline_snapshot['melA_peak_dpa']:.2f}`, "
        f"calibrated: `{final_calibrated_snapshot['melA_peak_dpa']:.2f}`"
    )
    if safety_adjustment["applied"]:
        report_lines.append(
            f"- Safety adjustment applied: shifted `pGhMat1` by `{safety_adjustment['shift_days']:.3f}` days"
        )
    else:
        report_lines.append("- Safety adjustment applied: none")
    report_lines.append("")

    report_path = RESULTS_DIR / f"{out_prefix}_report.md"
    report_path.write_text("\n".join(report_lines) + "\n")

    print("\n🧠 BlackCotton Transcriptome Calibrator")
    print("=" * 44)
    print(f"Reference table: {summary['reference_table']}")
    print(f"Baseline mean RMSE:  {summary['baseline_rmse']['mean']:.4f}")
    print(f"Calibrated mean RMSE:{summary['calibrated_rmse']['mean']:.4f}")
    print(f"RMSE delta:          {summary['improvements']['mean_rmse_delta']:.4f}")
    print(f"Output config:       {summary['output_config']}")
    print("Saved:")
    print(f"  - {summary_path.relative_to(BASE_DIR)}")
    print(f"  - {expr_cmp_path.relative_to(BASE_DIR)}")
    print(f"  - {curves_path.relative_to(BASE_DIR)}")
    print(f"  - {report_path.relative_to(BASE_DIR)}")
    print("=" * 44)

    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calibrate promoter parameters from transcriptome reference table.")
    parser.add_argument(
        "--reference",
        default=str(DATA_DIR / "cotton_expression_reference.csv"),
        help="CSV with columns: dpa, scw_marker_norm, maturation_marker_norm[, weight].",
    )
    parser.add_argument(
        "--output-config",
        default=str(CONFIG_DIR / "parameters_calibrated.yaml"),
        help="Path to write calibrated YAML parameters.",
    )
    parser.add_argument("--out-prefix", default="transcriptome_calibration", help="Output prefix under results/.")
    parser.add_argument("--seed", type=int, default=20260302, help="RNG seed for calibration search.")
    parser.add_argument(
        "--min-temporal-gap-days",
        type=float,
        default=0.5,
        help="Safety floor; if raw calibration violates this, pGhMat1 is shifted later.",
    )
    args = parser.parse_args()

    run_calibration(
        reference_path=Path(args.reference),
        output_config_path=Path(args.output_config),
        out_prefix=str(args.out_prefix),
        seed=int(args.seed),
        min_temporal_gap_days=float(args.min_temporal_gap_days),
    )
