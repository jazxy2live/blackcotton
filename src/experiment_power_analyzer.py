#!/usr/bin/env python3
"""
experiment_power_analyzer.py — Replicate and timepoint planning from trial data
===============================================================================

Uses robust Monte Carlo trial records to estimate how many replicates are needed
to pass one-sided confidence gates for each endpoint.

Usage:
    python -m src.experiment_power_analyzer
"""

import argparse
import csv
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from statistics import NormalDist
from typing import Any

import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"

DEFAULT_ALPHA = 0.05
DEFAULT_TARGET_POWER = 0.80
DEFAULT_MAX_REPLICATES = 24
DEFAULT_SIMULATIONS = 4000


@dataclass(frozen=True)
class EndpointSpec:
    key: str
    threshold_key: str
    direction: str  # "higher" means >= threshold, "lower" means <= threshold
    unit: str


ENDPOINT_SPECS = [
    EndpointSpec("color_L", "max_color_L", "lower", "L*"),
    EndpointSpec("strength_g_tex", "min_strength_g_tex", "higher", "g/tex"),
    EndpointSpec("yield_index", "min_yield_index", "higher", "index"),
    EndpointSpec("temporal_gap_days", "min_temporal_gap_days", "higher", "days"),
]

# Additional assay noise to avoid overconfident n from simulation-only variance.
MEASUREMENT_NOISE_SD = {
    "color_L": 0.30,
    "strength_g_tex": 0.10,
    "yield_index": 0.006,
    "temporal_gap_days": 0.15,
}


def load_json(path: Path, default: Any) -> Any:
    if not path.exists():
        return default
    with open(path) as f:
        return json.load(f)


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    rows = []
    if not path.exists():
        return rows
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def one_sided_confident_pass(
    sample: np.ndarray,
    threshold: float,
    direction: str,
    z_crit: float,
) -> bool:
    n = sample.size
    if n <= 1:
        return False
    mean = float(np.mean(sample))
    sd = float(np.std(sample, ddof=1))
    se = sd / float(np.sqrt(n))
    if direction == "higher":
        lower_bound = mean - z_crit * se
        return lower_bound >= threshold
    upper_bound = mean + z_crit * se
    return upper_bound <= threshold


def estimate_power_for_n(
    values: np.ndarray,
    threshold: float,
    direction: str,
    n_replicates: int,
    n_sims: int,
    alpha: float,
    measurement_sd: float,
    rng: np.random.Generator,
) -> float:
    if values.size < 2 or n_replicates < 2:
        return 0.0
    z_crit = NormalDist().inv_cdf(1.0 - alpha)
    passes = 0
    for _ in range(n_sims):
        sample = rng.choice(values, size=n_replicates, replace=True)
        if measurement_sd > 0:
            sample = sample + rng.normal(0.0, measurement_sd, size=n_replicates)
        if one_sided_confident_pass(sample, threshold, direction, z_crit):
            passes += 1
    return float(passes) / float(n_sims)


def required_n_for_target_power(
    values: np.ndarray,
    threshold: float,
    direction: str,
    target_power: float,
    max_replicates: int,
    n_sims: int,
    alpha: float,
    measurement_sd: float,
    seed: int,
) -> tuple[int | None, float]:
    rng = np.random.default_rng(seed)
    best_power = 0.0
    for n in range(2, max_replicates + 1):
        power = estimate_power_for_n(
            values=values,
            threshold=threshold,
            direction=direction,
            n_replicates=n,
            n_sims=n_sims,
            alpha=alpha,
            measurement_sd=measurement_sd,
            rng=rng,
        )
        best_power = max(best_power, power)
        if power >= target_power:
            return n, power
    return None, best_power


def summarize_controls(fiber_quality_rows: list[dict[str, Any]]) -> dict[str, Any]:
    lookup = {str(row.get("name", "")).strip(): row for row in fiber_quality_rows}
    return {
        "white_cotton": lookup.get("White Cotton (Coker-312)", {}),
        "dyed_black_cotton": lookup.get("Dyed Black Cotton", {}),
    }


def recommend_timepoints(final_top3: list[dict[str, Any]]) -> list[int]:
    points: set[int] = set()
    for row in final_top3:
        base = row.get("base_candidate", {})
        scw = int(round(float(base.get("scw_activation_dpa", 0.0))))
        mat = int(round(float(base.get("mat_activation_dpa", 0.0))))
        for p in (scw - 2, scw, scw + 2, mat - 2, mat, mat + 2, mat + 4):
            if p >= 0:
                points.add(p)
    return sorted(points)


def parse_current_replicates(path: Path) -> int:
    if not path.exists():
        return 8
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                return max(2, int(row.get("replicates", "8")))
            except ValueError:
                continue
    return 8


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    summary: dict[str, Any],
    detailed_rows: list[dict[str, Any]],
    out_prefix: str,
) -> None:
    lines: list[str] = []
    lines.append("# Experiment Power Analysis Report")
    lines.append("")
    lines.append("Replicate planning from robust trial distributions and one-sided confidence gates.")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Generated: `{summary['generated_at']}`")
    lines.append(f"- Trial records used: `{summary['n_trial_records']}`")
    lines.append(f"- Candidates analyzed: `{summary['n_candidates']}`")
    lines.append(f"- Current planned replicates: `{summary['current_replicates']}`")
    lines.append(f"- Target power: `{summary['target_power']:.2f}`")
    lines.append(f"- Confidence alpha (one-sided): `{summary['alpha']:.2f}`")
    lines.append(f"- Statistical minimum replicates/candidate: `{summary['recommended_replicates_by_candidate']}`")
    lines.append(f"- Statistical minimum overall: `{summary['recommended_replicates_overall']}`")
    lines.append(f"- Practical minimum floor: `{summary['practical_min_replicates']}`")
    lines.append(f"- Practical recommendation/candidate: `{summary['recommended_replicates_by_candidate_practical']}`")
    lines.append(f"- Practical recommendation overall: `{summary['recommended_replicates_overall_practical']}`")
    lines.append(f"- Recommended DPA timepoints: `{summary['recommended_timepoints_dpa']}`")
    lines.append("")
    lines.append("## Endpoint-Level Requirements")
    lines.append("")
    lines.append("| Candidate | Endpoint | Threshold | Mean | SD | Power @ Current n | Required n | Achieved Power |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|")
    for row in detailed_rows:
        req_n = row["required_n_target_power"] if row["required_n_target_power"] is not None else "n/a"
        lines.append(
            f"| {row['candidate_arm']} | {row['endpoint']} | {row['threshold']:.3f} | {row['mean']:.3f} | "
            f"{row['sd']:.3f} | {row['power_at_current_n']:.3f} | {req_n} | {row['power_at_required_n_or_max']:.3f} |"
        )
    lines.append("")

    (RESULTS_DIR / f"{out_prefix}_report.md").write_text("\n".join(lines) + "\n")


def run(
    trials_path: Path,
    final_top_path: Path,
    robust_summary_path: Path,
    fiber_quality_path: Path,
    target_power: float,
    alpha: float,
    max_replicates: int,
    n_sims: int,
    seed: int,
    out_prefix: str,
    practical_min_replicates: int,
) -> None:
    trials = load_jsonl(trials_path)
    if not trials:
        raise FileNotFoundError(f"Missing or empty trial records: {trials_path}")

    final_top3 = load_json(final_top_path, [])
    if not final_top3:
        raise FileNotFoundError(f"Missing final top candidates: {final_top_path}")

    robust_summary = load_json(robust_summary_path, {})
    thresholds = robust_summary.get("thresholds", {})
    if not thresholds:
        raise ValueError(f"Thresholds not found in {robust_summary_path}")

    controls = summarize_controls(load_json(fiber_quality_path, []))
    current_replicates = parse_current_replicates(RESULTS_DIR / "lab_validation_matrix.csv")

    by_rank: dict[int, list[dict[str, Any]]] = {}
    for row in trials:
        by_rank.setdefault(int(row["candidate_input_rank"]), []).append(row)

    detailed_rows: list[dict[str, Any]] = []
    recommended_by_candidate: dict[str, int] = {}

    for idx, top in enumerate(final_top3, 1):
        input_rank = int(top.get("input_rank", idx))
        arm_name = f"CAND_{idx}"
        rows = by_rank.get(input_rank, [])
        if len(rows) < 3:
            continue

        reqs = []
        for spec in ENDPOINT_SPECS:
            values = np.array([float(r[spec.key]) for r in rows], dtype=float)
            threshold = float(thresholds[spec.threshold_key])
            measurement_sd = float(MEASUREMENT_NOISE_SD.get(spec.key, 0.0))

            power_current = estimate_power_for_n(
                values=values,
                threshold=threshold,
                direction=spec.direction,
                n_replicates=current_replicates,
                n_sims=n_sims,
                alpha=alpha,
                measurement_sd=measurement_sd,
                rng=np.random.default_rng(seed + input_rank * 100 + len(spec.key)),
            )
            required_n, achieved = required_n_for_target_power(
                values=values,
                threshold=threshold,
                direction=spec.direction,
                target_power=target_power,
                max_replicates=max_replicates,
                n_sims=n_sims,
                alpha=alpha,
                measurement_sd=measurement_sd,
                seed=seed + input_rank * 1000 + len(spec.key),
            )
            reqs.append(required_n if required_n is not None else max_replicates)

            detailed_rows.append(
                {
                    "candidate_arm": arm_name,
                    "candidate_input_rank": input_rank,
                    "endpoint": spec.key,
                    "unit": spec.unit,
                    "direction": spec.direction,
                    "threshold": threshold,
                    "mean": float(np.mean(values)),
                    "sd": float(np.std(values, ddof=1)),
                    "n_trial_records": int(values.size),
                    "power_at_current_n": float(power_current),
                    "required_n_target_power": required_n,
                    "power_at_required_n_or_max": float(achieved),
                }
            )

        recommended_by_candidate[arm_name] = int(max(reqs)) if reqs else current_replicates

    recommended_overall = max(recommended_by_candidate.values()) if recommended_by_candidate else current_replicates
    recommended_by_candidate_practical = {
        k: max(v, practical_min_replicates)
        for k, v in recommended_by_candidate.items()
    }
    recommended_overall_practical = (
        max(recommended_by_candidate_practical.values())
        if recommended_by_candidate_practical
        else max(current_replicates, practical_min_replicates)
    )
    timepoints = recommend_timepoints(final_top3)

    summary = {
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "trials_path": str(trials_path.relative_to(BASE_DIR)),
        "final_top_path": str(final_top_path.relative_to(BASE_DIR)),
        "n_trial_records": int(len(trials)),
        "n_candidates": int(len(recommended_by_candidate)),
        "current_replicates": int(current_replicates),
        "target_power": float(target_power),
        "alpha": float(alpha),
        "max_replicates_checked": int(max_replicates),
        "simulations_per_n": int(n_sims),
        "recommended_replicates_by_candidate": recommended_by_candidate,
        "recommended_replicates_overall": int(recommended_overall),
        "practical_min_replicates": int(practical_min_replicates),
        "recommended_replicates_by_candidate_practical": recommended_by_candidate_practical,
        "recommended_replicates_overall_practical": int(recommended_overall_practical),
        "recommended_timepoints_dpa": timepoints,
        "control_references": controls,
    }

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_DIR / f"{out_prefix}_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    with open(RESULTS_DIR / f"{out_prefix}_details.json", "w") as f:
        json.dump(detailed_rows, f, indent=2)
    write_csv(RESULTS_DIR / f"{out_prefix}_details.csv", detailed_rows)
    write_report(summary=summary, detailed_rows=detailed_rows, out_prefix=out_prefix)

    print("\n📈 BlackCotton Experiment Power Analyzer")
    print("=" * 47)
    print(f"Trial records: {len(trials)}")
    print(f"Candidates analyzed: {len(recommended_by_candidate)}")
    print(f"Current planned replicates: {current_replicates}")
    print(f"Target power: {target_power:.2f} (alpha={alpha:.2f}, one-sided)")
    print(f"Statistical min replicates by candidate: {recommended_by_candidate}")
    print(f"Statistical min overall replicates: {recommended_overall}")
    print(f"Practical min floor: {practical_min_replicates}")
    print(f"Practical replicates by candidate: {recommended_by_candidate_practical}")
    print(f"Practical overall replicates: {recommended_overall_practical}")
    print(f"Recommended DPA timepoints: {timepoints}")
    print("Saved:")
    print(f"  - results/{out_prefix}_summary.json")
    print(f"  - results/{out_prefix}_details.json")
    print(f"  - results/{out_prefix}_details.csv")
    print(f"  - results/{out_prefix}_report.md")
    print("=" * 47)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Power analysis for lab experiment design.")
    parser.add_argument(
        "--trials",
        default=str(RESULTS_DIR / "adaptive_ode_robust_trial_records.jsonl"),
        help="Path to robust trial records JSONL.",
    )
    parser.add_argument(
        "--final-top",
        default=str(RESULTS_DIR / "final_lab_top3.json"),
        help="Path to final Top 3 JSON.",
    )
    parser.add_argument(
        "--robust-summary",
        default=str(RESULTS_DIR / "adaptive_ode_robust_summary.json"),
        help="Path to robust summary (for thresholds).",
    )
    parser.add_argument(
        "--fiber-quality",
        default=str(RESULTS_DIR / "fiber_quality_comparison.json"),
        help="Path to reference cotton quality comparison JSON.",
    )
    parser.add_argument("--target-power", type=float, default=DEFAULT_TARGET_POWER, help="Target statistical power.")
    parser.add_argument("--alpha", type=float, default=DEFAULT_ALPHA, help="One-sided alpha.")
    parser.add_argument("--max-replicates", type=int, default=DEFAULT_MAX_REPLICATES, help="Max replicates to scan.")
    parser.add_argument("--sims", type=int, default=DEFAULT_SIMULATIONS, help="Simulations per replicate count.")
    parser.add_argument("--seed", type=int, default=31415, help="RNG seed.")
    parser.add_argument(
        "--practical-min-replicates",
        type=int,
        default=6,
        help="Practical biological floor even if statistical minimum is lower.",
    )
    parser.add_argument("--out-prefix", default="power_analysis", help="Output prefix under results/.")
    args = parser.parse_args()

    run(
        trials_path=Path(args.trials),
        final_top_path=Path(args.final_top),
        robust_summary_path=Path(args.robust_summary),
        fiber_quality_path=Path(args.fiber_quality),
        target_power=float(args.target_power),
        alpha=float(args.alpha),
        max_replicates=max(2, int(args.max_replicates)),
        n_sims=max(500, int(args.sims)),
        seed=int(args.seed),
        out_prefix=str(args.out_prefix),
        practical_min_replicates=max(2, int(args.practical_min_replicates)),
    )
