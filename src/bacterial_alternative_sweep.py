#!/usr/bin/env python3
"""
bacterial_alternative_sweep.py — Focused sweep for the bacterial fallback path.

Turns binding-protein factor into an explicit design variable and ranks
broth-phase cellulose/melanin designs by safe melanin load, binder burden,
and retained matrix strength.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from src.pathway_safety_models import (
    BACTERIAL_ARCHITECTURES,
    BACTERIAL_TARGET_MELANIN_RATIO,
    simulate_bacterial_case,
)

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"


def _architecture_meta() -> dict[str, dict[str, Any]]:
    return {str(row["label"]): row for row in BACTERIAL_ARCHITECTURES}


def _frange(start: float, stop: float, step: float) -> list[float]:
    if step <= 0:
        raise ValueError("step must be > 0")
    out: list[float] = []
    value = float(start)
    guard = 0
    while value <= float(stop) + 1e-9:
        out.append(round(value, 6))
        value += float(step)
        guard += 1
        if guard > 10000:
            raise RuntimeError("range generation exceeded guard limit")
    return out


def _target_rows(rows: list[dict[str, Any]], target_ratio: float) -> list[dict[str, Any]]:
    return [row for row in rows if abs(float(row["melanin_to_cellulose_ratio"]) - float(target_ratio)) < 1e-9]


def _safe_sort_key(row: dict[str, Any]) -> tuple[float, float, float]:
    return (
        float(row["melanin_to_cellulose_ratio"]),
        -float(row["binding_protein_factor"]),
        float(row["wet_matrix_strength"]),
    )


def top_safe_designs(rows: list[dict[str, Any]], n: int = 5) -> list[dict[str, Any]]:
    safe_rows = [row for row in rows if bool(row["passes_hbond_gate"])]
    safe_rows.sort(
        key=lambda row: (
            float(row["melanin_to_cellulose_ratio"]),
            -float(row["binding_protein_factor"]),
            float(row["wet_matrix_strength"]),
            float(row["binding_affinity"]),
        ),
        reverse=True,
    )

    top: list[dict[str, Any]] = []
    seen: set[tuple[float, float]] = set()
    for row in safe_rows:
        key = (
            round(float(row["melanin_to_cellulose_ratio"]), 4),
            round(float(row["binding_protein_factor"]), 4),
        )
        if key in seen:
            continue
        seen.add(key)
        top.append(row)
        if len(top) >= int(n):
            break
    return top


def summarize_architecture(rows: list[dict[str, Any]], target_ratio: float) -> dict[str, Any]:
    meta = _architecture_meta()[str(rows[0]["architecture"])]
    safe_rows = [row for row in rows if bool(row["passes_hbond_gate"])]
    safe_rows.sort(
        key=lambda row: (
            float(row["melanin_to_cellulose_ratio"]),
            -float(row["binding_protein_factor"]),
            float(row["wet_matrix_strength"]),
        ),
        reverse=True,
    )

    target_rows = sorted(_target_rows(rows, target_ratio), key=lambda row: float(row["binding_protein_factor"]))
    target_safe = [row for row in target_rows if bool(row["passes_hbond_gate"])]
    no_linker_target = next((row for row in target_rows if abs(float(row["binding_protein_factor"])) < 1e-9), None)
    recommended = safe_rows[0] if safe_rows else None

    return {
        "label": str(meta["label"]),
        "description": str(meta["description"]),
        "target_ratio": float(target_ratio),
        "required_binding_factor_at_target": None if not target_safe else float(target_safe[0]["binding_protein_factor"]),
        "no_linker_target_row": no_linker_target,
        "recommended_design": recommended,
        "top_safe_designs": top_safe_designs(rows, n=5),
    }


def overall_recommendation(summaries: list[dict[str, Any]]) -> dict[str, Any] | None:
    candidates = [row for row in summaries if row["recommended_design"] is not None]
    if not candidates:
        return None
    candidates.sort(
        key=lambda row: (
            float(row["recommended_design"]["melanin_to_cellulose_ratio"]),
            -float(row["required_binding_factor_at_target"] or 999.0),
            float(row["recommended_design"]["wet_matrix_strength"]),
        ),
        reverse=True,
    )
    return candidates[0]


def write_report(path: Path, summary: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("# Bacterial Alternative Sweep")
    lines.append("")
    lines.append("Focused sweep over broth-phase melanin/cellulose compatibility with binding-protein factor as a design variable.")
    lines.append("")
    lines.append("## Headline")
    lines.append("")
    lines.append(f"- Target melanin/cellulose ratio: `{summary['target_ratio']:.2f}`")
    lines.append(
        f"- Architectures scanned: `{', '.join(summary['architectures_filter'])}`"
    )
    lines.append(
        f"- Ratio window: `{summary['ratio_window']['min']:.2f}` to `{summary['ratio_window']['max']:.2f}` "
        f"(step `{summary['ratio_window']['step']:.2f}`)"
    )
    lines.append(
        f"- Binder window: `{summary['binder_window']['min']:.2f}` to `{summary['binder_window']['max']:.2f}` "
        f"(step `{summary['binder_window']['step']:.2f}`)"
    )
    if summary["overall_recommendation"] is not None:
        rec = summary["overall_recommendation"]
        design = rec["recommended_design"]
        lines.append(
            f"- Recommended architecture: `{rec['label']}` at ratio `{design['melanin_to_cellulose_ratio']:.2f}`, "
            f"binder `{design['binding_protein_factor']:.2f}`, hbond `{design['hbond_retention']:.3f}`, "
            f"wet strength `{design['wet_matrix_strength']:.2f}`"
        )
    lines.append("")
    lines.append("## Architecture Summary")
    lines.append("")
    for architecture in summary["architectures"]:
        no_linker = architecture["no_linker_target_row"]
        lines.append(f"### {architecture['label']}")
        lines.append("")
        lines.append(
            f"- Required binding factor at target ratio: `{architecture['required_binding_factor_at_target']}`"
        )
        if no_linker is not None:
            lines.append(
                f"- No-linker target row: affinity `{no_linker['binding_affinity']:.3f}`, "
                f"hbond `{no_linker['hbond_retention']:.3f}`, wet strength `{no_linker['wet_matrix_strength']:.2f}`"
            )
        rec = architecture["recommended_design"]
        if rec is not None:
            lines.append(
                f"- Recommended design: ratio `{rec['melanin_to_cellulose_ratio']:.2f}`, "
                f"binder `{rec['binding_protein_factor']:.2f}`, hbond `{rec['hbond_retention']:.3f}`, "
                f"wet strength `{rec['wet_matrix_strength']:.2f}`"
            )
        lines.append("")
        lines.append("| Rank | Ratio | Binder | Affinity | Hbond | Wet Strength |")
        lines.append("|---:|---:|---:|---:|---:|---:|")
        for idx, row in enumerate(architecture["top_safe_designs"], 1):
            lines.append(
                f"| {idx} | {row['melanin_to_cellulose_ratio']:.2f} | {row['binding_protein_factor']:.2f} | "
                f"{row['binding_affinity']:.3f} | {row['hbond_retention']:.3f} | {row['wet_matrix_strength']:.2f} |"
            )
        lines.append("")

    path.write_text("\n".join(lines) + "\n")


def run(
    target_ratio: float,
    out_prefix: str,
    architectures_filter: list[str] | None = None,
    ratio_min: float = 0.08,
    ratio_max: float = 0.30,
    ratio_step: float = 0.02,
    binder_min: float = 0.0,
    binder_max: float = 0.8,
    binder_step: float = 0.1,
) -> dict[str, Any]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    meta = _architecture_meta()
    if architectures_filter:
        selected_architectures = [meta[name] for name in architectures_filter if name in meta]
    else:
        selected_architectures = [meta[name] for name in sorted(meta.keys())]

    ratio_values = _frange(ratio_min, ratio_max, ratio_step)
    binder_values = _frange(binder_min, binder_max, binder_step)

    all_rows: list[dict[str, Any]] = []
    for architecture in selected_architectures:
        for ratio in ratio_values:
            for binder in binder_values:
                all_rows.append(simulate_bacterial_case(architecture, float(ratio), float(binder)))

    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in all_rows:
        grouped.setdefault(str(row["architecture"]), []).append(row)

    architectures = [summarize_architecture(rows, target_ratio=target_ratio) for _, rows in sorted(grouped.items())]
    summary = {
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "target_ratio": float(target_ratio),
        "architectures_filter": [str(row["label"]) for row in selected_architectures],
        "ratio_window": {
            "min": float(ratio_min),
            "max": float(ratio_max),
            "step": float(ratio_step),
        },
        "binder_window": {
            "min": float(binder_min),
            "max": float(binder_max),
            "step": float(binder_step),
        },
        "architectures": architectures,
        "overall_recommendation": overall_recommendation(architectures),
    }

    summary_path = RESULTS_DIR / f"{out_prefix}_summary.json"
    report_path = RESULTS_DIR / f"{out_prefix}_report.md"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    write_report(report_path, summary)
    return {
        "summary_path": summary_path,
        "report_path": report_path,
        "summary": summary,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run focused bacterial alternative sweep.")
    parser.add_argument(
        "--target-ratio",
        type=float,
        default=BACTERIAL_TARGET_MELANIN_RATIO,
        help="Target melanin/cellulose ratio to test rescue requirements at.",
    )
    parser.add_argument(
        "--out-prefix",
        default="bacterial_alternative_sweep_2026_03_07",
        help="Output prefix under results/.",
    )
    parser.add_argument(
        "--architectures",
        default="",
        help="Comma-separated architecture labels to include (e.g. single_strain).",
    )
    parser.add_argument("--ratio-min", type=float, default=0.08, help="Minimum melanin/cellulose ratio.")
    parser.add_argument("--ratio-max", type=float, default=0.30, help="Maximum melanin/cellulose ratio.")
    parser.add_argument("--ratio-step", type=float, default=0.02, help="Step for melanin/cellulose ratio.")
    parser.add_argument("--binder-min", type=float, default=0.0, help="Minimum binding-protein factor.")
    parser.add_argument("--binder-max", type=float, default=0.8, help="Maximum binding-protein factor.")
    parser.add_argument("--binder-step", type=float, default=0.1, help="Step for binding-protein factor.")
    args = parser.parse_args()

    architectures_filter = [item.strip() for item in str(args.architectures).split(",") if item.strip()]
    output = run(
        target_ratio=float(args.target_ratio),
        out_prefix=str(args.out_prefix),
        architectures_filter=architectures_filter or None,
        ratio_min=float(args.ratio_min),
        ratio_max=float(args.ratio_max),
        ratio_step=float(args.ratio_step),
        binder_min=float(args.binder_min),
        binder_max=float(args.binder_max),
        binder_step=float(args.binder_step),
    )
    summary = output["summary"]
    print("\nBacterial alternative sweep complete")
    print("=" * 37)
    if summary["overall_recommendation"] is not None:
        rec = summary["overall_recommendation"]
        design = rec["recommended_design"]
        print(
            f"Recommended: {rec['label']} ratio={design['melanin_to_cellulose_ratio']:.2f} "
            f"binder={design['binding_protein_factor']:.2f} hbond={design['hbond_retention']:.3f}"
        )
    print(f"Saved: {output['summary_path'].relative_to(BASE_DIR)}")
    print(f"Saved: {output['report_path'].relative_to(BASE_DIR)}")
