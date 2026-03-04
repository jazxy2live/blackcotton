#!/usr/bin/env python3
"""
proof_packager.py — Rebuild showable proof artifacts from latest results
=======================================================================

Creates a consistent reporting package from the most recent optimization and
robustness outputs so external stakeholders see the current project status.

Usage:
    python -m src.proof_packager
"""

import argparse
import hashlib
import json
import tarfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"
CONFIG_DIR = BASE_DIR / "config"


def first_present(d: dict[str, Any], keys: list[str], default: Any = "n/a") -> Any:
    for key in keys:
        if key in d:
            return d[key]
    return default


def load_json(path: Path, default: Any) -> Any:
    if not path.exists():
        return default
    with open(path) as f:
        return json.load(f)


def sha256_file(path: Path) -> str | None:
    if not path.exists() or not path.is_file():
        return None
    hasher = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def parse_test_report(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"status": "missing", "ran_tests": 0}

    text = path.read_text()
    ran_tests = 0
    status = "unknown"
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("Ran ") and " tests" in line:
            tokens = line.split()
            if len(tokens) >= 2 and tokens[1].isdigit():
                ran_tests = int(tokens[1])
        if line == "OK":
            status = "pass"
        if line.startswith("FAILED"):
            status = "fail"
    return {"status": status, "ran_tests": ran_tests}


def _normalize_final_lab_row(row: dict[str, Any]) -> dict[str, Any]:
    base = row.get("base_candidate", row)
    return {
        "source": "final_lab_top3",
        "model_stage": "adaptive_ode_robust",
        "success_rate": float(row.get("success_rate", 0.0)),
        "robust_score": float(row.get("robust_score", 0.0)),
        "color_L": float(base.get("color_L", row.get("p50_color_L", row.get("color_L", 0.0)))),
        "strength_g_tex": float(base.get("strength_g_tex", row.get("p50_strength_g_tex", row.get("strength_g_tex", 0.0)))),
        "yield_index": float(base.get("yield_index", row.get("p50_yield_index", row.get("yield_index", 0.0)))),
        "temporal_gap_days": float(base.get("temporal_gap_days", row.get("p50_temporal_gap_days", row.get("temporal_gap_days", 0.0)))),
        "mat_activation_dpa": float(base.get("mat_activation_dpa", 0.0)),
        "scw_activation_dpa": float(base.get("scw_activation_dpa", 0.0)),
        "mat_strength": float(base.get("mat_strength", 0.0)),
        "scw_strength": float(base.get("scw_strength", 0.0)),
        "k_competition": float(base.get("k_competition", 0.0)),
        "melanin_efficiency": float(base.get("melanin_efficiency", 0.0)),
        "late_retention_factor": float(base.get("late_retention_factor", 0.0)),
        "p50_color_L": float(row.get("p50_color_L", 0.0)),
        "p50_strength_g_tex": float(row.get("p50_strength_g_tex", 0.0)),
        "p50_yield_index": float(row.get("p50_yield_index", 0.0)),
        "p50_temporal_gap_days": float(row.get("p50_temporal_gap_days", 0.0)),
    }


def _normalize_ode_row(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "source": "top_candidates",
        "model_stage": "ode_refined",
        "success_rate": None,
        "robust_score": float(row.get("composite_score", 0.0)),
        "color_L": float(row.get("color_L", 0.0)),
        "strength_g_tex": float(row.get("strength_g_tex", 0.0)),
        "yield_index": float(row.get("yield_index", 0.0)),
        "temporal_gap_days": float(row.get("temporal_gap_days", 0.0)),
        "mat_activation_dpa": float(row.get("mat_activation_dpa", 0.0)),
        "scw_activation_dpa": float(row.get("scw_activation_dpa", 0.0)),
        "mat_strength": float(row.get("mat_strength", 0.0)),
        "scw_strength": float(row.get("scw_strength", 0.0)),
        "k_competition": float(row.get("k_competition", 0.0)),
        "melanin_efficiency": float(row.get("melanin_efficiency", 0.0)),
        "late_retention_factor": float(row.get("late_retention_factor", 0.0)),
        "p50_color_L": None,
        "p50_strength_g_tex": None,
        "p50_yield_index": None,
        "p50_temporal_gap_days": None,
    }


def select_lead_candidate(
    final_lab_top3: list[dict[str, Any]],
    top_candidates: list[dict[str, Any]],
) -> dict[str, Any]:
    if final_lab_top3:
        return _normalize_final_lab_row(final_lab_top3[0])
    if top_candidates:
        return _normalize_ode_row(top_candidates[0])
    return {
        "source": "none",
        "model_stage": "none",
        "success_rate": None,
        "robust_score": 0.0,
        "color_L": 0.0,
        "strength_g_tex": 0.0,
        "yield_index": 0.0,
        "temporal_gap_days": 0.0,
        "mat_activation_dpa": 0.0,
        "scw_activation_dpa": 0.0,
        "mat_strength": 0.0,
        "scw_strength": 0.0,
        "k_competition": 0.0,
        "melanin_efficiency": 0.0,
        "late_retention_factor": 0.0,
        "p50_color_L": None,
        "p50_strength_g_tex": None,
        "p50_yield_index": None,
        "p50_temporal_gap_days": None,
    }


def build_artifact_inventory(paths: list[Path]) -> list[dict[str, Any]]:
    inventory = []
    for path in paths:
        inventory.append(
            {
                "path": str(path.relative_to(BASE_DIR)),
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else 0,
                "sha256": sha256_file(path),
            }
        )
    return inventory


def write_canonical_report(
    generated_at: str,
    construct_summary: dict[str, Any],
    expression_summary: dict[str, Any],
    melanin_summary: dict[str, Any],
    optimization_summary: dict[str, Any],
    adaptive_ode_summary: dict[str, Any],
    lead: dict[str, Any],
    final_lab_top3: list[dict[str, Any]],
) -> None:
    lines: list[str] = []
    lines.append("# Canonical Pipeline Report")
    lines.append("")
    lines.append(f"- Generated at: {generated_at}")
    lines.append(f"- Workspace: `{BASE_DIR}`")
    lines.append("- Interpreter baseline: `./venv/bin/python` for all pipeline modules")
    lines.append("- Lead candidate source priority: `final_lab_top3 > top_candidates`")
    lines.append("")
    lines.append("## Key Metrics")
    lines.append("")
    lines.append(
        f"- Construct: `{construct_summary.get('name', 'n/a')}`, length `{construct_summary.get('total_length_bp', 'n/a')}` bp, "
        f"cassettes `{len(construct_summary.get('cassettes', []))}`"
    )
    lines.append(
        f"- Expression: cellulose 90% at `{first_present(expression_summary, ['cellulose_90pct_dpa'])}` DPA, "
        f"melA peak at `{first_present(expression_summary, ['melA_peak_dpa', 'peak_melA_dpa'])}` DPA"
    )
    lines.append(
        f"- Melanin pathway: final melanin `{melanin_summary.get('final_melanin_mM', 'n/a')}` mM eq, "
        f"tyrosine consumed `{melanin_summary.get('tyrosine_consumed_pct', 'n/a')}`%"
    )
    lines.append(
        f"- Optimization: coarse `{optimization_summary.get('num_candidates_coarse', 'n/a')}`, "
        f"refined `{optimization_summary.get('num_candidates_refined', 'n/a')}`, "
        f"pareto `{optimization_summary.get('num_pareto_refined', 'n/a')}`"
    )
    lines.append(
        f"- Final robust stage: ODE top `{adaptive_ode_summary.get('num_ode_top_for_robustness', 'n/a')}`, "
        f"trials/candidate `{adaptive_ode_summary.get('n_trials_per_candidate', 'n/a')}`, "
        f"best success `{adaptive_ode_summary.get('best_success_rate', 0.0):.3f}`"
    )
    lines.append("")
    lines.append("## Lead Candidate Snapshot")
    lines.append("")
    lines.append(
        f"- Stage `{lead['model_stage']}` from `{lead['source']}`: "
        f"L* `{lead['color_L']:.3f}`, strength `{lead['strength_g_tex']:.3f}` g/tex, "
        f"yield `{lead['yield_index']:.3f}`, gap `{lead['temporal_gap_days']:.3f}` days"
    )
    if lead["success_rate"] is not None:
        lines.append(
            f"- Robustness: success `{lead['success_rate']:.3f}`, robust score `{lead['robust_score']:.3f}`, "
            f"p50 L* `{lead['p50_color_L']:.3f}`, p50 strength `{lead['p50_strength_g_tex']:.3f}`, "
            f"p50 yield `{lead['p50_yield_index']:.3f}`"
        )
    else:
        lines.append(f"- Composite score `{lead['robust_score']:.3f}` (deterministic stage)")
    lines.append(
        f"- Parameters: mat/scw `{lead['mat_activation_dpa']:.0f}/{lead['scw_activation_dpa']:.0f}`, "
        f"k `{lead['k_competition']:.2f}`, eff `{lead['melanin_efficiency']:.2f}`, "
        f"ret `{lead['late_retention_factor']:.2f}`"
    )
    lines.append("")
    lines.append("## Final Lab Top 3")
    lines.append("")
    if final_lab_top3:
        lines.append("| Rank | Success | Robust | p50 L* | p50 Str | p50 Yield | Params |")
        lines.append("|---:|---:|---:|---:|---:|---:|---|")
        for i, row in enumerate(final_lab_top3[:3], 1):
            base = row.get("base_candidate", {})
            lines.append(
                f"| {i} | {float(row.get('success_rate', 0.0)):.3f} | {float(row.get('robust_score', 0.0)):.3f} | "
                f"{float(row.get('p50_color_L', 0.0)):.2f} | {float(row.get('p50_strength_g_tex', 0.0)):.2f} | "
                f"{float(row.get('p50_yield_index', 0.0)):.3f} | "
                f"{float(base.get('mat_activation_dpa', 0.0)):.0f}/{float(base.get('scw_activation_dpa', 0.0)):.0f}, "
                f"k={float(base.get('k_competition', 0.0)):.2f}, eff={float(base.get('melanin_efficiency', 0.0)):.2f}, "
                f"ret={float(base.get('late_retention_factor', 0.0)):.2f} |"
            )
    else:
        lines.append("- `final_lab_top3.json` not available.")
    lines.append("")

    (RESULTS_DIR / "canonical_report.md").write_text("\n".join(lines) + "\n")


def write_proof_package(
    lead: dict[str, Any],
    test_info: dict[str, Any],
    adaptive_ode_summary: dict[str, Any],
    has_lab_package: bool = False,
    has_power_package: bool = False,
    has_adversarial_package: bool = False,
    has_calibration_package: bool = False,
    has_calibration_impact_package: bool = False,
) -> None:
    lines: list[str] = []
    lines.append("# BlackCotton Proof Package")
    lines.append("")
    lines.append("## What This Proves")
    lines.append("")
    lines.append("This run produces timing-safe dark-fiber candidates with robust uncertainty stress-testing.")
    lines.append("- Hard safety constraint used: `temporal_gap_days >= 0 (hard)`")
    lines.append("- Final shortlist is ODE-refined and Monte Carlo-validated")
    lines.append(f"- Automated tests: `{test_info['status']}` (`{test_info['ran_tests']}` tests)")
    lines.append("")
    lines.append("## Current Lead Candidate")
    lines.append("")
    lines.append(f"- Stage/source: `{lead['model_stage']}` / `{lead['source']}`")
    lines.append(f"- Color (L*): `{lead['color_L']:.3f}`")
    lines.append(f"- Strength: `{lead['strength_g_tex']:.3f} g/tex`")
    lines.append(f"- Yield index: `{lead['yield_index']:.3f}`")
    lines.append(f"- Temporal gap: `{lead['temporal_gap_days']:.3f} days`")
    if lead["success_rate"] is not None:
        lines.append(f"- Robust success rate: `{lead['success_rate']:.3f}`")
        lines.append(f"- Robust score: `{lead['robust_score']:.3f}`")
    lines.append("")
    lines.append("## Robust Stage Summary")
    lines.append("")
    lines.append(f"- ODE candidates evaluated under uncertainty: `{adaptive_ode_summary.get('num_robust_candidates', 'n/a')}`")
    lines.append(f"- Trials per candidate: `{adaptive_ode_summary.get('n_trials_per_candidate', 'n/a')}`")
    lines.append(f"- Best robust success rate: `{adaptive_ode_summary.get('best_success_rate', 0.0):.3f}`")
    lines.append("")
    if has_lab_package:
        lines.append("## Lab Validation Package")
        lines.append("")
        lines.append("- Pilot matrix and data-capture templates were generated from final Top 3.")
        lines.append("- Phase-1 and stretch pass/fail gates are frozen in JSON for objective go/no-go.")
        lines.append("")
    if has_power_package:
        lines.append("## Power Planning Package")
        lines.append("")
        lines.append("- Replicate counts were estimated per endpoint/candidate from robust trial distributions.")
        lines.append("- DPA timepoints were generated from Top 3 activation windows.")
        lines.append("")
    if has_adversarial_package:
        lines.append("## Adversarial Robustness Package")
        lines.append("")
        lines.append("- Top candidates were red-teamed under strict gates and high-noise scenarios.")
        lines.append("- Worst-case candidate ranking and scenario leaderboard were generated.")
        lines.append("")
    if has_calibration_package:
        lines.append("## Transcriptome Calibration Package")
        lines.append("")
        lines.append("- Promoter timing/shape parameters were fit to reference expression trajectories.")
        lines.append("- Calibrated parameter file and before/after fit diagnostics were generated.")
        lines.append("")
    if has_calibration_impact_package:
        lines.append("## Calibration Impact Package")
        lines.append("")
        lines.append("- Baseline vs calibrated configs were compared on optimization and robustness outcomes.")
        lines.append("- A decision summary reports which config is better under worst-case stress.")
        lines.append("")
    lines.append("## Evidence Files")
    lines.append("")
    evidence_lines = [
        "- `results/canonical_report.md`",
        "- `results/proof_manifest.json`",
        "- `results/one_page_proof_brief.md`",
        "- `results/final_lab_top3.json`",
        "- `results/adaptive_ode_robust_summary.json`",
        "- `results/adaptive_ode_robust_report.md`",
        "- `results/test_report.txt`",
        "- `results/full_pipeline_run.log`",
    ]
    if has_lab_package:
        evidence_lines.extend(
            [
                "- `results/lab_validation_plan.md`",
                "- `results/lab_validation_matrix.csv`",
                "- `results/lab_measurements_template.csv`",
                "- `results/lab_pass_fail_rules.json`",
            ]
        )
    if has_power_package:
        evidence_lines.extend(
            [
                "- `results/power_analysis_summary.json`",
                "- `results/power_analysis_details.json`",
                "- `results/power_analysis_details.csv`",
                "- `results/power_analysis_report.md`",
            ]
        )
    if has_adversarial_package:
        evidence_lines.extend(
            [
                "- `results/adversarial_robustness_summary.json`",
                "- `results/adversarial_robustness_scenarios.json`",
                "- `results/adversarial_robustness_rows.json`",
                "- `results/adversarial_robustness_candidate_summary.json`",
                "- `results/adversarial_robustness_report.md`",
            ]
        )
    if has_calibration_package:
        evidence_lines.extend(
            [
                "- `results/transcriptome_calibration_summary.json`",
                "- `results/transcriptome_calibration_expression_comparison.json`",
                "- `results/transcriptome_calibration_curves.csv`",
                "- `results/transcriptome_calibration_report.md`",
                "- `config/parameters_calibrated.yaml`",
            ]
        )
    if has_calibration_impact_package:
        evidence_lines.extend(
            [
                "- `results/calibration_impact_summary.json`",
                "- `results/calibration_impact_baseline.json`",
                "- `results/calibration_impact_calibrated.json`",
                "- `results/calibration_impact_report.md`",
            ]
        )
    lines.extend(evidence_lines)
    lines.append("")

    (RESULTS_DIR / "proof_package.md").write_text("\n".join(lines) + "\n")


def write_one_page_brief(
    generated_at: str,
    lead: dict[str, Any],
    optimization_summary: dict[str, Any],
    adaptive_ode_summary: dict[str, Any],
    test_info: dict[str, Any],
) -> None:
    lines: list[str] = []
    lines.append("# BlackCotton One-Page Proof Brief")
    lines.append("")
    lines.append(f"Date: {generated_at[:10]}")
    lines.append("Project: In-silico design for non-fading black cotton via delayed melanin biosynthesis")
    lines.append("")
    lines.append("## 1) Core Problem")
    lines.append("Dark cotton typically loses fiber quality because pigment and cellulose formation overlap.")
    lines.append("")
    lines.append("## 2) Hypothesis")
    lines.append("Delay pigment pathway activation until late fiber development to keep strength while achieving deep black color.")
    lines.append("")
    lines.append("## 3) What Exists Now")
    lines.append("- End-to-end simulation pipeline from construct design to robust candidate ranking")
    lines.append("- Hard non-negative temporal gap constraint in optimizer")
    lines.append("- ODE refinement + uncertainty stress-test + frozen final lab Top 3")
    lines.append("")
    lines.append("## 4) Current Evidence (Reproducible)")
    lines.append(f"- Coarse candidates scanned: `{optimization_summary.get('num_candidates_coarse', 'n/a')}`")
    lines.append(f"- ODE refined candidates: `{optimization_summary.get('num_candidates_refined', 'n/a')}`")
    lines.append(f"- ODE top candidates stress-tested: `{adaptive_ode_summary.get('num_ode_top_for_robustness', 'n/a')}`")
    lines.append(f"- Trials per candidate: `{adaptive_ode_summary.get('n_trials_per_candidate', 'n/a')}`")
    lines.append(f"- Best robust success rate: `{adaptive_ode_summary.get('best_success_rate', 0.0):.3f}`")
    lines.append(
        f"- Lead candidate: L* `{lead['color_L']:.3f}`, strength `{lead['strength_g_tex']:.3f} g/tex`, "
        f"yield `{lead['yield_index']:.3f}`, gap `{lead['temporal_gap_days']:.3f} d`"
    )
    lines.append(f"- Tests: `{test_info['ran_tests']}` run, status `{test_info['status']}`")
    lines.append("")
    lines.append("## 5) Why This Matters")
    lines.append("This converts trial-and-error into ranked, uncertainty-aware experiments so lab cycles are faster and cheaper.")
    lines.append("")
    lines.append("## 6) Next Step")
    lines.append("Validate the final Top 3 in lab for expression timing, darkness durability, and fiber quality.")
    lines.append("")
    lines.append("## 7) Proof Files")
    lines.extend(
        [
            "- `results/canonical_report.md`",
            "- `results/proof_manifest.json`",
            "- `results/proof_package.md`",
            "- `results/final_lab_top3.json`",
            "- `results/adaptive_ode_robust_report.md`",
        ]
    )
    lines.append("")

    (RESULTS_DIR / "one_page_proof_brief.md").write_text("\n".join(lines) + "\n")


def create_bundle(bundle_name: str, files: list[Path]) -> Path:
    bundle_path = RESULTS_DIR / bundle_name
    with tarfile.open(bundle_path, "w:gz") as tar:
        for path in files:
            if path.exists():
                tar.add(path, arcname=str(path.relative_to(BASE_DIR)))
    return bundle_path


def write_manifest(
    generated_at: str,
    optimization_summary: dict[str, Any],
    adaptive_ode_summary: dict[str, Any],
    lead: dict[str, Any],
    test_info: dict[str, Any],
    inventory: list[dict[str, Any]],
    bundle_path: Path,
) -> None:
    manifest = {
        "generated_at": generated_at,
        "workspace": str(BASE_DIR),
        "optimizer_summary": optimization_summary,
        "adaptive_ode_robust_summary": adaptive_ode_summary,
        "lead_candidate": lead,
        "tests": test_info,
        "artifacts": inventory,
        "proof_bundle": str(bundle_path.relative_to(BASE_DIR)),
    }
    with open(RESULTS_DIR / "proof_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)


def run() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    generated_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()

    construct_summary = load_json(RESULTS_DIR / "construct_summary.json", {})
    expression_summary = load_json(RESULTS_DIR / "expression_summary.json", {})
    melanin_summary = load_json(RESULTS_DIR / "melanin_summary.json", {})
    optimization_summary = load_json(RESULTS_DIR / "optimization_summary.json", {})
    adaptive_ode_summary = load_json(RESULTS_DIR / "adaptive_ode_robust_summary.json", {})
    final_lab_top3 = load_json(RESULTS_DIR / "final_lab_top3.json", [])
    top_candidates = load_json(RESULTS_DIR / "top_candidates.json", [])
    test_info = parse_test_report(RESULTS_DIR / "test_report.txt")
    lead = select_lead_candidate(final_lab_top3=final_lab_top3, top_candidates=top_candidates)
    lab_package_files = [
        RESULTS_DIR / "lab_validation_plan.md",
        RESULTS_DIR / "lab_validation_matrix.csv",
        RESULTS_DIR / "lab_measurements_template.csv",
        RESULTS_DIR / "lab_pass_fail_rules.json",
    ]
    has_lab_package = all(path.exists() for path in lab_package_files)
    power_package_files = [
        RESULTS_DIR / "power_analysis_summary.json",
        RESULTS_DIR / "power_analysis_details.json",
        RESULTS_DIR / "power_analysis_details.csv",
        RESULTS_DIR / "power_analysis_report.md",
    ]
    has_power_package = all(path.exists() for path in power_package_files)
    adversarial_package_files = [
        RESULTS_DIR / "adversarial_robustness_summary.json",
        RESULTS_DIR / "adversarial_robustness_scenarios.json",
        RESULTS_DIR / "adversarial_robustness_rows.json",
        RESULTS_DIR / "adversarial_robustness_candidate_summary.json",
        RESULTS_DIR / "adversarial_robustness_report.md",
    ]
    has_adversarial_package = all(path.exists() for path in adversarial_package_files)
    calibration_package_files = [
        RESULTS_DIR / "transcriptome_calibration_summary.json",
        RESULTS_DIR / "transcriptome_calibration_expression_comparison.json",
        RESULTS_DIR / "transcriptome_calibration_curves.csv",
        RESULTS_DIR / "transcriptome_calibration_report.md",
        CONFIG_DIR / "parameters_calibrated.yaml",
    ]
    has_calibration_package = all(path.exists() for path in calibration_package_files)
    calibration_impact_package_files = [
        RESULTS_DIR / "calibration_impact_summary.json",
        RESULTS_DIR / "calibration_impact_baseline.json",
        RESULTS_DIR / "calibration_impact_calibrated.json",
        RESULTS_DIR / "calibration_impact_report.md",
    ]
    has_calibration_impact_package = all(path.exists() for path in calibration_impact_package_files)

    write_canonical_report(
        generated_at=generated_at,
        construct_summary=construct_summary,
        expression_summary=expression_summary,
        melanin_summary=melanin_summary,
        optimization_summary=optimization_summary,
        adaptive_ode_summary=adaptive_ode_summary,
        lead=lead,
        final_lab_top3=final_lab_top3,
    )
    write_proof_package(
        lead=lead,
        test_info=test_info,
        adaptive_ode_summary=adaptive_ode_summary,
        has_lab_package=has_lab_package,
        has_power_package=has_power_package,
        has_adversarial_package=has_adversarial_package,
        has_calibration_package=has_calibration_package,
        has_calibration_impact_package=has_calibration_impact_package,
    )
    write_one_page_brief(
        generated_at=generated_at,
        lead=lead,
        optimization_summary=optimization_summary,
        adaptive_ode_summary=adaptive_ode_summary,
        test_info=test_info,
    )

    bundle_day = generated_at[:10]
    bundle_name = f"proof_bundle_{bundle_day}.tgz"
    evidence_files = [
        RESULTS_DIR / "canonical_report.md",
        RESULTS_DIR / "proof_package.md",
        RESULTS_DIR / "one_page_proof_brief.md",
        RESULTS_DIR / "final_lab_top3.json",
        RESULTS_DIR / "adaptive_ode_robust_summary.json",
        RESULTS_DIR / "adaptive_ode_robust_report.md",
        RESULTS_DIR / "test_report.txt",
        RESULTS_DIR / "full_pipeline_run.log",
    ]
    if has_lab_package:
        evidence_files.extend(lab_package_files)
    if has_power_package:
        evidence_files.extend(power_package_files)
    if has_adversarial_package:
        evidence_files.extend(adversarial_package_files)
    if has_calibration_package:
        evidence_files.extend(calibration_package_files)
    if has_calibration_impact_package:
        evidence_files.extend(calibration_impact_package_files)
    bundle_path = RESULTS_DIR / bundle_name
    inventory = build_artifact_inventory(evidence_files)

    write_manifest(
        generated_at=generated_at,
        optimization_summary=optimization_summary,
        adaptive_ode_summary=adaptive_ode_summary,
        lead=lead,
        test_info=test_info,
        inventory=inventory,
        bundle_path=bundle_path,
    )
    bundle_files = evidence_files + [RESULTS_DIR / "proof_manifest.json"]
    create_bundle(bundle_name=bundle_name, files=bundle_files)

    print("\n📦 BlackCotton Proof Packager")
    print("=" * 42)
    print(f"Generated: {generated_at}")
    print(f"Lead stage: {lead['model_stage']} ({lead['source']})")
    print(f"Lead metrics: L*={lead['color_L']:.3f}, strength={lead['strength_g_tex']:.3f}, "
          f"yield={lead['yield_index']:.3f}, gap={lead['temporal_gap_days']:.3f}")
    if lead["success_rate"] is not None:
        print(f"Robust: success={lead['success_rate']:.3f}, score={lead['robust_score']:.3f}")
    print(f"Tests: {test_info['status']} ({test_info['ran_tests']} tests)")
    print("Saved:")
    print("  - results/canonical_report.md")
    print("  - results/proof_package.md")
    print("  - results/one_page_proof_brief.md")
    print("  - results/proof_manifest.json")
    print(f"  - results/{bundle_name}")
    print("=" * 42)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build proof artifacts from latest pipeline outputs.")
    parser.parse_args()
    run()
