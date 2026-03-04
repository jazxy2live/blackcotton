#!/usr/bin/env python3
"""
latest_graphs_builder.py — Build readable graphs for latest calibration impact
==============================================================================

Generates a compact graph pack from the newest baseline-vs-calibrated analysis.

Usage:
    python -m src.latest_graphs_builder
"""

import argparse
import csv
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent.parent / "results" / ".mplconfig"))

import matplotlib.pyplot as plt
import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Missing required JSON: {path}")
    with open(path) as f:
        return json.load(f)


def load_curves(path: Path) -> list[dict[str, float]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing required CSV: {path}")
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({k: float(v) for k, v in row.items()})
    return rows


def save_figure(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=170)
    plt.close()


def graph_optimization_footprint(summary: dict[str, Any], out_dir: Path) -> str:
    b = summary["baseline"]
    c = summary["calibrated"]
    labels = ["coarse", "seed", "refined", "pareto", "top"]
    base_vals = [b["n_coarse"], b["n_seed"], b["n_refined"], b["n_pareto"], b["n_top"]]
    cal_vals = [c["n_coarse"], c["n_seed"], c["n_refined"], c["n_pareto"], c["n_top"]]

    x = np.arange(len(labels))
    w = 0.38
    plt.figure(figsize=(10, 5))
    plt.bar(x - w / 2, base_vals, width=w, label="Baseline", color="#C45E3E")
    plt.bar(x + w / 2, cal_vals, width=w, label="Calibrated", color="#2D7F84")
    plt.xticks(x, labels)
    plt.ylabel("Count")
    plt.title("Optimization Footprint (Baseline vs Calibrated)")
    plt.legend(frameon=False)
    plt.grid(axis="y", alpha=0.25)
    out = out_dir / "01_optimization_footprint.png"
    save_figure(out)
    return out.name


def graph_expression_timing(summary: dict[str, Any], out_dir: Path) -> str:
    b = summary["baseline"]["expression_snapshot"]
    c = summary["calibrated"]["expression_snapshot"]
    labels = ["Cellulose 90% DPA", "melA half-max DPA", "Temporal gap (days)"]
    base_vals = [b["cellulose_90pct_dpa"], b["melA_halfmax_dpa"], b["temporal_gap_days"]]
    cal_vals = [c["cellulose_90pct_dpa"], c["melA_halfmax_dpa"], c["temporal_gap_days"]]

    x = np.arange(len(labels))
    w = 0.38
    plt.figure(figsize=(10, 5))
    plt.bar(x - w / 2, base_vals, width=w, label="Baseline", color="#C45E3E")
    plt.bar(x + w / 2, cal_vals, width=w, label="Calibrated", color="#2D7F84")
    plt.xticks(x, labels, rotation=10)
    plt.ylabel("Value")
    plt.title("Expression Timing Snapshot")
    plt.legend(frameon=False)
    plt.grid(axis="y", alpha=0.25)
    out = out_dir / "02_expression_timing.png"
    save_figure(out)
    return out.name


def graph_robustness_comparison(summary: dict[str, Any], out_dir: Path) -> str:
    b_def = summary["baseline"]["robust_default_best"]
    c_def = summary["calibrated"]["robust_default_best"]
    b_wc = summary["baseline"]["robust_worst_case_best"]
    c_wc = summary["calibrated"]["robust_worst_case_best"]

    labels = [
        "Default success",
        "Default robust",
        "Worst success",
        "Worst robust",
    ]
    base_vals = [b_def["success_rate"], b_def["robust_score"], b_wc["success_rate"], b_wc["robust_score"]]
    cal_vals = [c_def["success_rate"], c_def["robust_score"], c_wc["success_rate"], c_wc["robust_score"]]

    x = np.arange(len(labels))
    w = 0.38
    plt.figure(figsize=(10, 5))
    plt.bar(x - w / 2, base_vals, width=w, label="Baseline", color="#C45E3E")
    plt.bar(x + w / 2, cal_vals, width=w, label="Calibrated", color="#2D7F84")
    plt.xticks(x, labels)
    plt.ylabel("Score / Probability")
    plt.ylim(0.0, 1.0)
    plt.title("Robustness Comparison (Higher is Better)")
    plt.legend(frameon=False)
    plt.grid(axis="y", alpha=0.25)
    out = out_dir / "03_robustness_comparison.png"
    save_figure(out)
    return out.name


def graph_deterministic_metrics(summary: dict[str, Any], out_dir: Path) -> str:
    b = summary["baseline"]["deterministic_best"]
    c = summary["calibrated"]["deterministic_best"]
    labels = ["Color L* (lower better)", "Strength", "Yield", "Gap"]
    base_vals = [b["color_L"], b["strength_g_tex"], b["yield_index"], b["temporal_gap_days"]]
    cal_vals = [c["color_L"], c["strength_g_tex"], c["yield_index"], c["temporal_gap_days"]]

    x = np.arange(len(labels))
    w = 0.38
    plt.figure(figsize=(11, 5))
    plt.bar(x - w / 2, base_vals, width=w, label="Baseline", color="#C45E3E")
    plt.bar(x + w / 2, cal_vals, width=w, label="Calibrated", color="#2D7F84")
    plt.xticks(x, labels, rotation=8)
    plt.ylabel("Metric value")
    plt.title("Deterministic Best-Candidate Metrics")
    plt.legend(frameon=False)
    plt.grid(axis="y", alpha=0.25)
    out = out_dir / "04_deterministic_metrics.png"
    save_figure(out)
    return out.name


def graph_calibration_fit_curves(curves_rows: list[dict[str, float]], out_dir: Path) -> str:
    dpa = np.array([r["dpa"] for r in curves_rows], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)

    axes[0].plot(dpa, [r["observed_scw_marker_norm"] for r in curves_rows], "ko", ms=4, label="Reference")
    axes[0].plot(dpa, [r["baseline_scw_fit"] for r in curves_rows], color="#C45E3E", lw=2, label="Baseline fit")
    axes[0].plot(dpa, [r["calibrated_scw_fit"] for r in curves_rows], color="#2D7F84", lw=2, label="Calibrated fit")
    axes[0].set_title("SCW Marker Fit")
    axes[0].set_xlabel("DPA")
    axes[0].set_ylabel("Normalized expression")
    axes[0].grid(alpha=0.2)
    axes[0].legend(frameon=False)

    axes[1].plot(dpa, [r["observed_maturation_marker_norm"] for r in curves_rows], "ko", ms=4, label="Reference")
    axes[1].plot(dpa, [r["baseline_maturation_fit"] for r in curves_rows], color="#C45E3E", lw=2, label="Baseline fit")
    axes[1].plot(dpa, [r["calibrated_maturation_fit"] for r in curves_rows], color="#2D7F84", lw=2, label="Calibrated fit")
    axes[1].set_title("Maturation Marker Fit")
    axes[1].set_xlabel("DPA")
    axes[1].grid(alpha=0.2)
    axes[1].legend(frameon=False)

    out = out_dir / "05_calibration_fit_curves.png"
    save_figure(out)
    return out.name


def write_index(
    out_dir: Path,
    summary: dict[str, Any],
    generated_files: list[str],
) -> None:
    decision = summary.get("decision", {})
    b_wc = summary["baseline"]["robust_worst_case_best"]
    c_wc = summary["calibrated"]["robust_worst_case_best"]

    lines = []
    lines.append("# Latest Graphs Index")
    lines.append("")
    lines.append(f"- Generated: `{datetime.now(timezone.utc).replace(microsecond=0).isoformat()}`")
    lines.append(f"- Winner from latest impact analysis: `{decision.get('winner', 'n/a')}`")
    lines.append(f"- Reason: `{decision.get('reason', 'n/a')}`")
    lines.append("")
    lines.append("## Fast Read")
    lines.append("")
    lines.append(
        f"- Worst-case success improved from `{b_wc['success_rate']:.3f}` "
        f"to `{c_wc['success_rate']:.3f}`."
    )
    lines.append(
        f"- Worst-case robust score improved from `{b_wc['robust_score']:.3f}` "
        f"to `{c_wc['robust_score']:.3f}`."
    )
    lines.append("")
    lines.append("## Graph Files")
    lines.append("")
    for fname in generated_files:
        lines.append(f"- `{fname}`")
    lines.append("")
    (out_dir / "INDEX.md").write_text("\n".join(lines) + "\n")


def run(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    impact_summary = load_json(RESULTS_DIR / "calibration_impact_summary.json")
    curves_rows = load_curves(RESULTS_DIR / "transcriptome_calibration_curves.csv")

    generated = []
    generated.append(graph_optimization_footprint(impact_summary, out_dir))
    generated.append(graph_expression_timing(impact_summary, out_dir))
    generated.append(graph_robustness_comparison(impact_summary, out_dir))
    generated.append(graph_deterministic_metrics(impact_summary, out_dir))
    generated.append(graph_calibration_fit_curves(curves_rows, out_dir))

    manifest = {
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_summary": "results/calibration_impact_summary.json",
        "source_curves": "results/transcriptome_calibration_curves.csv",
        "output_dir": str(out_dir.relative_to(BASE_DIR)),
        "files": generated,
    }
    with open(out_dir / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)

    write_index(out_dir, impact_summary, generated)

    print("\n📈 BlackCotton Latest Graph Pack")
    print("=" * 40)
    print(f"Output folder: {out_dir.relative_to(BASE_DIR)}")
    print("Files:")
    for f in generated:
        print(f"  - {f}")
    print("  - INDEX.md")
    print("  - manifest.json")
    print("=" * 40)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build latest easy-read graph folder.")
    parser.add_argument(
        "--out-dir",
        default=str(RESULTS_DIR / "latest_graphs"),
        help="Folder to place generated graph pack.",
    )
    args = parser.parse_args()
    run(Path(args.out_dir))
