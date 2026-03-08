#!/usr/bin/env python3
"""
worst_case_hardening_sprint.py — Worst-case robust hardening sprint
====================================================================

Purpose:
  - Improve strict-high-noise robustness instead of average-case robustness.
  - Produce a clean before/after package without mutating frozen baseline files.

Usage:
  python -m src.worst_case_hardening_sprint
"""

import argparse
import json
import os
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent.parent / "results" / ".mplconfig"))

import matplotlib.pyplot as plt
import numpy as np
import yaml

from src.adversarial_robustness_suite import (
    aggregate_worst_case,
    normalize_top_candidate_rows,
    scenario_configs,
)
from src.config_loader import resolve_config_path
from src.pathway_safety_models import simulate_chemical_ros_safety
from src.robustness_analyzer import run_robustness_analysis
from src.tradeoff_optimizer import build_candidate_params, evaluate_candidate_proxy, load_params

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"

DEFAULT_GRID = {
    "mat_activation_dpa": [37.0, 38.0, 39.0, 40.0, 41.0, 42.0],
    "scw_activation_dpa": [31.0, 32.0, 33.0, 34.0, 35.0],
    "mat_strength": [1.0, 1.1, 1.2, 1.3, 1.4],
    "scw_strength": [1.0, 1.1, 1.2, 1.3, 1.4],
    "k_competition": [0.03, 0.05, 0.07, 0.10, 0.12],
    "melanin_efficiency": [1.4, 1.5, 1.6, 1.7],
    "late_retention_factor": [0.65, 0.75, 0.85],
}


def find_named_scenario(
    name: str,
    correlated_profile: str | None = None,
    include_correlated: bool = False,
) -> dict[str, Any]:
    for scenario in scenario_configs(
        correlated_profile=correlated_profile,
        include_correlated=include_correlated,
    ):
        if str(scenario["name"]) == str(name):
            return scenario
    raise ValueError(f"Scenario not found: {name}")


def scenario_display_name(name: str) -> str:
    return str(name).replace("_", " ")


def _candidate_key(row: dict[str, Any]) -> tuple[float, float, float, float, float, float, float]:
    return (
        round(float(row["mat_activation_dpa"]), 4),
        round(float(row["scw_activation_dpa"]), 4),
        round(float(row["mat_strength"]), 4),
        round(float(row["scw_strength"]), 4),
        round(float(row["k_competition"]), 4),
        round(float(row["melanin_efficiency"]), 4),
        round(float(row["late_retention_factor"]), 4),
    )


def strict_margin_score(row: dict[str, Any], thresholds: dict[str, float]) -> tuple[float, bool]:
    color_margin = float(thresholds["max_color_L"]) - float(row["color_L"])
    strength_margin = float(row["strength_g_tex"]) - float(thresholds["min_strength_g_tex"])
    yield_margin = float(row["yield_index"]) - float(thresholds["min_yield_index"])
    gap_margin = float(row["temporal_gap_days"]) - float(thresholds["min_temporal_gap_days"])
    tox_margin = float(thresholds.get("max_toxicity_pre_cellulose", 1e9)) - float(
        row.get("toxicity_pre_cellulose", 0.0)
    )

    strict_pass = (
        color_margin >= 0.0
        and strength_margin >= 0.0
        and yield_margin >= 0.0
        and gap_margin >= 0.0
        and tox_margin >= 0.0
        and bool(row.get("toxicity_gate_pass", True))
    )

    # Yield is narrow-range, so it gets scaled to avoid being drowned by color.
    score = (
        color_margin
        + 1.3 * strength_margin
        + 20.0 * yield_margin
        + 0.2 * gap_margin
        + 1.8 * min(tox_margin, 1.0)
    )
    return float(score), bool(strict_pass)


def build_deterministic_pool(
    params: dict[str, Any],
    strict_thresholds: dict[str, float],
    min_temporal_gap_days: float,
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    grid = DEFAULT_GRID
    for values in product(
        grid["mat_activation_dpa"],
        grid["scw_activation_dpa"],
        grid["mat_strength"],
        grid["scw_strength"],
        grid["k_competition"],
        grid["melanin_efficiency"],
        grid["late_retention_factor"],
    ):
        row = evaluate_candidate_proxy(
            params=params,
            mat_activation_dpa=float(values[0]),
            scw_activation_dpa=float(values[1]),
            mat_strength=float(values[2]),
            scw_strength=float(values[3]),
            k_competition=float(values[4]),
            melanin_efficiency=float(values[5]),
            late_retention_factor=float(values[6]),
        )
        if float(row["temporal_gap_days"]) < float(min_temporal_gap_days):
            continue
        if not bool(row.get("toxicity_gate_pass", True)):
            continue
        margin_score, strict_pass = strict_margin_score(row, strict_thresholds)
        row["_strict_margin_score"] = float(margin_score)
        row["_strict_pass"] = bool(strict_pass)
        out.append(row)
    return out


def select_seed_candidates(pool: list[dict[str, Any]], limit: int) -> list[dict[str, Any]]:
    sorted_rows = sorted(
        pool,
        key=lambda r: (
            1 if bool(r.get("_strict_pass", False)) else 0,
            float(r.get("_strict_margin_score", -999.0)),
            float(r.get("composite_score", 0.0)),
        ),
        reverse=True,
    )
    deduped: list[dict[str, Any]] = []
    seen = set()
    for row in sorted_rows:
        key = _candidate_key(row)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(row)
        if len(deduped) >= int(limit):
            break
    return deduped


def filter_chemical_safe_seed_candidates(
    params: dict[str, Any],
    candidates: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for row in candidates:
        candidate_params = build_candidate_params(params, row)
        chemical = simulate_chemical_ros_safety(candidate_params)
        if bool(chemical["kill_before_opacity_90"]):
            continue
        item = dict(row)
        item["_chemical_peak_ros_ratio"] = float(chemical["peak_ros_ratio"])
        item["_chemical_peak_ros_proxy"] = float(chemical["peak_ros_proxy"])
        item["_chemical_opacity_90_dpa"] = float(chemical["opacity_90_dpa"])
        out.append(item)
    return out


def strip_internal_fields(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for row in rows:
        item = dict(row)
        for key in list(item.keys()):
            if key.startswith("_"):
                item.pop(key, None)
        out.append(item)
    return out


def pick_top_unique_robust_rows(
    robust_rows: list[dict[str, Any]],
    n: int = 3,
    min_success_rate: float = 0.20,
    max_fragility_index: float = 0.40,
) -> list[dict[str, Any]]:
    pool = [
        r
        for r in robust_rows
        if (
            float(r["success_rate"]) >= float(min_success_rate)
            and float(r.get("fragility_index", 0.0)) <= float(max_fragility_index)
        )
    ]
    if not pool:
        pool = robust_rows
    selected: list[dict[str, Any]] = []
    seen = set()
    for row in pool:
        key = _candidate_key(row["base_candidate"])
        if key in seen:
            continue
        seen.add(key)
        selected.append(row)
        if len(selected) >= int(n):
            break
    return selected


def run_top3_adversarial_pack(
    params: dict[str, Any],
    top3_candidates: list[dict[str, Any]],
    n_trials: int,
    seed: int,
    correlated_profile: str | None = None,
    include_correlated: bool = False,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    scenario_summary: list[dict[str, Any]] = []
    all_rows: list[dict[str, Any]] = []
    for i, scenario in enumerate(
        scenario_configs(
            correlated_profile=correlated_profile,
            include_correlated=include_correlated,
        )
    ):
        robust_rows = run_robustness_analysis(
            params=params,
            candidates=top3_candidates,
            n_trials=int(n_trials),
            seed=int(seed) + i * 101,
            thresholds=dict(scenario["thresholds"]),
            noise=dict(scenario["noise"]),
            correlated_profile=scenario.get("correlated_profile"),
            collect_trial_records=False,
        )
        top = robust_rows[0] if robust_rows else {}
        scenario_summary.append(
            {
                "scenario": scenario["name"],
                "description": scenario["description"],
                "top_candidate_input_rank": int(top.get("input_rank", 0)) if top else 0,
                "top_success_rate": float(top.get("success_rate", 0.0)) if top else 0.0,
                "top_robust_score": float(top.get("robust_score", 0.0)) if top else 0.0,
            }
        )
        for row in robust_rows:
            all_rows.append(
                {
                    "scenario": scenario["name"],
                    "candidate_input_rank": int(row["input_rank"]),
                    "success_rate": float(row["success_rate"]),
                    "robust_score": float(row["robust_score"]),
                }
            )
    candidate_summary = aggregate_worst_case(all_rows)
    return scenario_summary, candidate_summary, all_rows


def _scenario_map(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {str(r["scenario"]): r for r in rows}


def build_scenario_delta(
    baseline: list[dict[str, Any]],
    hardened: list[dict[str, Any]],
    scenario_names: list[str],
) -> list[dict[str, Any]]:
    b_map = _scenario_map(baseline)
    h_map = _scenario_map(hardened)
    out = []
    for scenario in scenario_names:
        b = b_map.get(str(scenario), {})
        h = h_map.get(str(scenario), {})
        b_s = float(b.get("top_success_rate", 0.0))
        h_s = float(h.get("top_success_rate", 0.0))
        b_r = float(b.get("top_robust_score", 0.0))
        h_r = float(h.get("top_robust_score", 0.0))
        out.append(
            {
                "scenario": str(scenario),
                "baseline_top_success_rate": b_s,
                "hardened_top_success_rate": h_s,
                "delta_top_success_rate": h_s - b_s,
                "baseline_top_robust_score": b_r,
                "hardened_top_robust_score": h_r,
                "delta_top_robust_score": h_r - b_r,
            }
        )
    return out


def save_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)


def save_fig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=170)
    plt.close()


def build_graph_pack(
    out_dir: Path,
    baseline_objective_rows: list[dict[str, Any]],
    hardened_objective_rows: list[dict[str, Any]],
    scenario_delta: list[dict[str, Any]],
    objective_scenario_name: str,
) -> list[str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    files: list[str] = []
    objective_label = scenario_display_name(objective_scenario_name)

    b_top = baseline_objective_rows[0] if baseline_objective_rows else {}
    h_top = hardened_objective_rows[0] if hardened_objective_rows else {}

    labels = ["Top success", "Top robust"]
    b_vals = [float(b_top.get("success_rate", 0.0)), float(b_top.get("robust_score", 0.0))]
    h_vals = [float(h_top.get("success_rate", 0.0)), float(h_top.get("robust_score", 0.0))]
    x = np.arange(len(labels))
    w = 0.38
    plt.figure(figsize=(8, 4.5))
    plt.bar(x - w / 2, b_vals, width=w, color="#C45E3E", label="Baseline Top3")
    plt.bar(x + w / 2, h_vals, width=w, color="#2D7F84", label="Hardened Top3")
    plt.xticks(x, labels)
    plt.ylim(0.0, 1.0)
    plt.ylabel("Score")
    plt.title(f"{objective_label} best candidate")
    plt.grid(axis="y", alpha=0.25)
    plt.legend(frameon=False)
    p1 = out_dir / "01_objective_top_gain.png"
    save_fig(p1)
    files.append(p1.name)

    scenarios = [str(r["scenario"]) for r in scenario_delta]
    b_success = [float(r["baseline_top_success_rate"]) for r in scenario_delta]
    h_success = [float(r["hardened_top_success_rate"]) for r in scenario_delta]
    x2 = np.arange(len(scenarios))
    plt.figure(figsize=(9, 4.8))
    plt.plot(x2, b_success, marker="o", lw=2.2, color="#C45E3E", label="Baseline Top3")
    plt.plot(x2, h_success, marker="o", lw=2.2, color="#2D7F84", label="Hardened Top3")
    plt.xticks(x2, scenarios)
    plt.ylim(0.0, 1.0)
    plt.ylabel("Top success rate")
    plt.title("Scenario-wise top success")
    plt.grid(alpha=0.25)
    plt.legend(frameon=False)
    p2 = out_dir / "02_scenario_top_success.png"
    save_fig(p2)
    files.append(p2.name)

    b3 = baseline_objective_rows[:3]
    h3 = hardened_objective_rows[:3]
    labels3 = ["Rank-1", "Rank-2", "Rank-3"]
    b3_s = [float(r["success_rate"]) for r in b3] + [0.0] * (3 - len(b3))
    h3_s = [float(r["success_rate"]) for r in h3] + [0.0] * (3 - len(h3))
    x3 = np.arange(len(labels3))
    plt.figure(figsize=(8.5, 4.8))
    plt.bar(x3 - w / 2, b3_s, width=w, color="#C45E3E", label="Baseline Top3")
    plt.bar(x3 + w / 2, h3_s, width=w, color="#2D7F84", label="Hardened Top3")
    plt.xticks(x3, labels3)
    plt.ylim(0.0, 1.0)
    plt.ylabel("Success rate")
    plt.title(f"{objective_label} top-3 success distribution")
    plt.grid(axis="y", alpha=0.25)
    plt.legend(frameon=False)
    p3 = out_dir / "03_top3_success_distribution.png"
    save_fig(p3)
    files.append(p3.name)

    lines = []
    lines.append("# Worst-Case Hardening Graph Pack")
    lines.append("")
    lines.append(f"- Generated: `{datetime.now(timezone.utc).replace(microsecond=0).isoformat()}`")
    lines.append(f"- Objective scenario: `{objective_scenario_name}`")
    lines.append(
        f"- Objective top success: `{float(b_top.get('success_rate', 0.0)):.3f}` -> "
        f"`{float(h_top.get('success_rate', 0.0)):.3f}`"
    )
    lines.append("")
    lines.append("## Files")
    lines.append("")
    for name in files:
        lines.append(f"- `{name}`")
    lines.append("")
    (out_dir / "INDEX.md").write_text("\n".join(lines) + "\n")

    save_json(
        out_dir / "manifest.json",
        {
            "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            "output_dir": str(out_dir.relative_to(BASE_DIR)),
            "files": files,
        },
    )
    return files


def write_report(
    path: Path,
    summary: dict[str, Any],
    scenario_delta: list[dict[str, Any]],
    hardened_top3: list[dict[str, Any]],
) -> None:
    objective_label = scenario_display_name(summary["objective_scenario"])
    lines: list[str] = []
    lines.append("# Worst-Case Hardening Sprint Report")
    lines.append("")
    lines.append(
        f"Objective: maximize `{summary['objective_scenario']}` success instead of average-case robustness."
    )
    lines.append("")
    lines.append("## Headline")
    lines.append("")
    lines.append(
        f"- {objective_label} top success: `{summary['baseline_objective_top_success_rate']:.3f}` -> "
        f"`{summary['hardened_objective_top_success_rate']:.3f}` "
        f"(delta `{summary['delta_objective_top_success_rate']:+.3f}`)"
    )
    lines.append(
        f"- {objective_label} top robust score: `{summary['baseline_objective_top_robust_score']:.3f}` -> "
        f"`{summary['hardened_objective_top_robust_score']:.3f}` "
        f"(delta `{summary['delta_objective_top_robust_score']:+.3f}`)"
    )
    lines.append(f"- Target success (`>= {summary['target_success_rate']:.3f}`): `{summary['target_met']}`")
    lines.append(f"- Correlated profile: `{summary['correlated_profile']}`")
    lines.append(
        f"- Chemical-safe seeds: `{summary['seed_candidates_after_chemical_gate']}` / "
        f"`{summary['seed_candidates_before_chemical_gate']}`"
    )
    lines.append("")
    lines.append("## Scenario Delta")
    lines.append("")
    lines.append("| Scenario | Baseline Top Success | Hardened Top Success | Delta |")
    lines.append("|---|---:|---:|---:|")
    for row in scenario_delta:
        lines.append(
            f"| {row['scenario']} | {row['baseline_top_success_rate']:.3f} | "
            f"{row['hardened_top_success_rate']:.3f} | {row['delta_top_success_rate']:+.3f} |"
        )
    lines.append("")
    lines.append(f"## Hardened Top 3 ({objective_label} ranked)")
    lines.append("")
    lines.append("| Rank | Success | Robust | Fragility | p50 L* | p50 Strength | p50 Yield | p50 Gap | Params |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for i, row in enumerate(hardened_top3, 1):
        b = row["base_candidate"]
        lines.append(
            f"| {i} | {row['success_rate']:.3f} | {row['robust_score']:.3f} | {row.get('fragility_index', 0.0):.3f} | {row['p50_color_L']:.2f} | "
            f"{row['p50_strength_g_tex']:.2f} | {row['p50_yield_index']:.3f} | {row['p50_temporal_gap_days']:.2f} | "
            f"{b['mat_activation_dpa']:.0f}/{b['scw_activation_dpa']:.0f}, "
            f"k={b['k_competition']:.2f}, eff={b['melanin_efficiency']:.2f}, ret={b['late_retention_factor']:.2f} |"
        )
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


def load_params_from_path(config_path: Path) -> dict[str, Any]:
    default_path = resolve_config_path()
    if config_path.resolve() == default_path.resolve():
        return load_params()
    with open(config_path) as f:
        return yaml.safe_load(f)


def run(
    config_path: Path,
    n_trials: int,
    seed: int,
    seed_candidates_limit: int,
    out_prefix: str,
    min_temporal_gap_days: float,
    target_success_rate: float,
    include_graphs: bool,
    objective_scenario_name: str,
    correlated_profile: str | None,
) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    params = load_params_from_path(config_path)
    include_correlated = bool(correlated_profile and str(correlated_profile) != "none") or str(
        objective_scenario_name
    ).endswith("_correlated")
    objective_scenario = find_named_scenario(
        objective_scenario_name,
        correlated_profile=correlated_profile,
        include_correlated=include_correlated,
    )
    objective_thr = dict(objective_scenario["thresholds"])
    objective_noise = dict(objective_scenario["noise"])
    objective_correlated = objective_scenario.get("correlated_profile")
    scenario_names = [
        str(s["name"])
        for s in scenario_configs(
            correlated_profile=correlated_profile,
            include_correlated=include_correlated,
        )
    ]

    current_path = RESULTS_DIR / "final_lab_top3.json"
    if not current_path.exists():
        raise FileNotFoundError(
            "Missing results/final_lab_top3.json. Run `python -m src.adaptive_ode_robust_pipeline` first."
        )
    baseline_top3 = normalize_top_candidate_rows(json.loads(current_path.read_text()))

    baseline_objective_rows = run_robustness_analysis(
        params=params,
        candidates=baseline_top3,
        n_trials=int(n_trials),
        seed=int(seed),
        thresholds=objective_thr,
        noise=objective_noise,
        correlated_profile=objective_correlated,
        collect_trial_records=False,
    )

    deterministic_pool = build_deterministic_pool(
        params=params,
        strict_thresholds=objective_thr,
        min_temporal_gap_days=float(min_temporal_gap_days),
    )
    seed_candidates_pre_chemical = select_seed_candidates(deterministic_pool, limit=seed_candidates_limit)
    seed_candidates = filter_chemical_safe_seed_candidates(params, seed_candidates_pre_chemical)
    if not seed_candidates:
        raise RuntimeError("Chemical kill gate removed every hardening seed candidate.")

    hardened_objective_rows = run_robustness_analysis(
        params=params,
        candidates=seed_candidates,
        n_trials=int(n_trials),
        seed=int(seed) + 77,
        thresholds=objective_thr,
        noise=objective_noise,
        correlated_profile=objective_correlated,
        collect_trial_records=False,
    )
    hardened_top3 = pick_top_unique_robust_rows(
        hardened_objective_rows,
        n=3,
        min_success_rate=0.20,
        max_fragility_index=0.40,
    )
    hardened_top3_candidates = [dict(r["base_candidate"]) for r in hardened_top3]

    baseline_scenario_summary, baseline_candidate_summary, baseline_all_rows = run_top3_adversarial_pack(
        params=params,
        top3_candidates=baseline_top3,
        n_trials=int(n_trials),
        seed=int(seed) + 900,
        correlated_profile=correlated_profile,
        include_correlated=include_correlated,
    )
    hardened_scenario_summary, hardened_candidate_summary, hardened_all_rows = run_top3_adversarial_pack(
        params=params,
        top3_candidates=hardened_top3_candidates,
        n_trials=int(n_trials),
        seed=int(seed) + 1200,
        correlated_profile=correlated_profile,
        include_correlated=include_correlated,
    )
    scenario_delta = build_scenario_delta(
        baseline_scenario_summary,
        hardened_scenario_summary,
        scenario_names=scenario_names,
    )

    b_best = baseline_objective_rows[0] if baseline_objective_rows else {}
    h_best = hardened_objective_rows[0] if hardened_objective_rows else {}
    b_success = float(b_best.get("success_rate", 0.0))
    h_success = float(h_best.get("success_rate", 0.0))
    b_robust = float(b_best.get("robust_score", 0.0))
    h_robust = float(h_best.get("robust_score", 0.0))

    summary = {
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "config_path": str(config_path.relative_to(BASE_DIR)),
        "objective_scenario": str(objective_scenario_name),
        "objective_description": str(objective_scenario.get("description", "")),
        "correlated_profile": str(objective_correlated or "none"),
        "n_trials": int(n_trials),
        "seed": int(seed),
        "grid_size_total": int(len(DEFAULT_GRID["mat_activation_dpa"]) * len(DEFAULT_GRID["scw_activation_dpa"]) * len(DEFAULT_GRID["mat_strength"]) * len(DEFAULT_GRID["scw_strength"]) * len(DEFAULT_GRID["k_competition"]) * len(DEFAULT_GRID["melanin_efficiency"]) * len(DEFAULT_GRID["late_retention_factor"])),
        "deterministic_pool_size": int(len(deterministic_pool)),
        "seed_candidates_before_chemical_gate": int(len(seed_candidates_pre_chemical)),
        "seed_candidates_after_chemical_gate": int(len(seed_candidates)),
        "baseline_objective_top_success_rate": b_success,
        "hardened_objective_top_success_rate": h_success,
        "delta_objective_top_success_rate": h_success - b_success,
        "baseline_objective_top_robust_score": b_robust,
        "hardened_objective_top_robust_score": h_robust,
        "delta_objective_top_robust_score": h_robust - b_robust,
        "target_success_rate": float(target_success_rate),
        "target_met": bool(h_success >= float(target_success_rate)),
    }

    save_json(RESULTS_DIR / f"{out_prefix}_summary.json", summary)
    save_json(RESULTS_DIR / f"{out_prefix}_seed_candidates.json", strip_internal_fields(seed_candidates))
    save_json(RESULTS_DIR / f"{out_prefix}_objective_candidates.json", hardened_objective_rows)
    save_json(RESULTS_DIR / f"{out_prefix}_top3.json", hardened_top3)
    save_json(RESULTS_DIR / f"{out_prefix}_baseline_objective_rows.json", baseline_objective_rows)
    save_json(RESULTS_DIR / f"{out_prefix}_scenario_baseline.json", baseline_scenario_summary)
    save_json(RESULTS_DIR / f"{out_prefix}_scenario_hardened.json", hardened_scenario_summary)
    save_json(RESULTS_DIR / f"{out_prefix}_scenario_delta.json", scenario_delta)
    save_json(RESULTS_DIR / f"{out_prefix}_candidate_baseline_summary.json", baseline_candidate_summary)
    save_json(RESULTS_DIR / f"{out_prefix}_candidate_hardened_summary.json", hardened_candidate_summary)
    save_json(RESULTS_DIR / f"{out_prefix}_rows_baseline.json", baseline_all_rows)
    save_json(RESULTS_DIR / f"{out_prefix}_rows_hardened.json", hardened_all_rows)

    write_report(
        RESULTS_DIR / f"{out_prefix}_report.md",
        summary=summary,
        scenario_delta=scenario_delta,
        hardened_top3=hardened_top3,
    )

    graph_files: list[str] = []
    if include_graphs:
        graph_dir = RESULTS_DIR / f"{out_prefix}_graphs"
        graph_files = build_graph_pack(
            out_dir=graph_dir,
            baseline_objective_rows=baseline_objective_rows,
            hardened_objective_rows=hardened_objective_rows,
            scenario_delta=scenario_delta,
            objective_scenario_name=str(objective_scenario_name),
        )

    print("\n🚀 BlackCotton Worst-Case Hardening Sprint")
    print("=" * 46)
    print(f"Config: {config_path.relative_to(BASE_DIR)}")
    print(f"Objective scenario:            {objective_scenario_name}")
    print(f"Correlated profile:           {objective_correlated or 'none'}")
    print(f"Objective top success:        {b_success:.3f} -> {h_success:.3f} ({h_success - b_success:+.3f})")
    print(f"Objective top robust:         {b_robust:.3f} -> {h_robust:.3f} ({h_robust - b_robust:+.3f})")
    print(
        f"Chemical-safe seeds:           {summary['seed_candidates_after_chemical_gate']} / "
        f"{summary['seed_candidates_before_chemical_gate']}"
    )
    print(f"Target met (>= {target_success_rate:.3f}): {summary['target_met']}")
    print("Saved:")
    print(f"  - results/{out_prefix}_summary.json")
    print(f"  - results/{out_prefix}_top3.json")
    print(f"  - results/{out_prefix}_report.md")
    if graph_files:
        print(f"  - results/{out_prefix}_graphs ({len(graph_files)} PNGs + index)")
    print("=" * 46)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run strict-high-noise hardening sprint.")
    parser.add_argument(
        "--config-path",
        default=str(resolve_config_path()),
        help="Path to YAML parameter config.",
    )
    parser.add_argument("--n-trials", type=int, default=220, help="Monte Carlo trials per candidate.")
    parser.add_argument("--seed", type=int, default=20260302, help="RNG seed.")
    parser.add_argument("--seed-candidates", type=int, default=260, help="Number of deterministic seeds to stress-test.")
    parser.add_argument("--out-prefix", default="worst_case_hardening", help="Output prefix under results/.")
    parser.add_argument("--min-temporal-gap-days", type=float, default=0.5, help="Deterministic minimum gap filter.")
    parser.add_argument("--target-success-rate", type=float, default=0.85, help="Sprint target for strict-high-noise top success.")
    parser.add_argument(
        "--objective-scenario",
        default="strict_high_noise",
        help="Scenario name from adversarial_robustness_suite to optimize against.",
    )
    parser.add_argument(
        "--correlated-profile",
        default="none",
        help="Optional correlated profile name to enable correlated scenarios.",
    )
    parser.add_argument(
        "--no-graphs",
        action="store_true",
        help="Disable graph generation.",
    )
    args = parser.parse_args()

    run(
        config_path=Path(args.config_path),
        n_trials=max(int(args.n_trials), 60),
        seed=int(args.seed),
        seed_candidates_limit=max(40, int(args.seed_candidates)),
        out_prefix=str(args.out_prefix),
        min_temporal_gap_days=float(args.min_temporal_gap_days),
        target_success_rate=float(args.target_success_rate),
        include_graphs=not bool(args.no_graphs),
        objective_scenario_name=str(args.objective_scenario),
        correlated_profile=None if str(args.correlated_profile) == "none" else str(args.correlated_profile),
    )
