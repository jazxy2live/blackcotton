# BlackCotton Proof Package

## What This Proves

This run produces timing-safe dark-fiber candidates with robust uncertainty stress-testing.
- Hard safety constraint used: `temporal_gap_days >= 0 (hard)`
- Final shortlist is ODE-refined and Monte Carlo-validated
- Automated tests: `pass` (`28` tests)

## Current Lead Candidate

- Stage/source: `adaptive_ode_robust` / `final_lab_top3`
- Color (L*): `21.587`
- Strength: `29.350 g/tex`
- Yield index: `0.934`
- Temporal gap: `5.500 days`
- Robust success rate: `0.267`
- Robust score: `0.455`

## Robust Stage Summary

- ODE candidates evaluated under uncertainty: `3`
- Trials per candidate: `120`
- Best robust success rate: `0.833`

## Lab Validation Package

- Pilot matrix and data-capture templates were generated from final Top 3.
- Phase-1 and stretch pass/fail gates are frozen in JSON for objective go/no-go.

## Power Planning Package

- Replicate counts were estimated per endpoint/candidate from robust trial distributions.
- DPA timepoints were generated from Top 3 activation windows.

## Adversarial Robustness Package

- Top candidates were red-teamed under strict gates and high-noise scenarios.
- Worst-case candidate ranking and scenario leaderboard were generated.

## Transcriptome Calibration Package

- Promoter timing/shape parameters were fit to reference expression trajectories.
- Calibrated parameter file and before/after fit diagnostics were generated.

## Calibration Impact Package

- Baseline vs calibrated configs were compared on optimization and robustness outcomes.
- A decision summary reports which config is better under worst-case stress.

## Evidence Files

- `results/canonical_report.md`
- `results/proof_manifest.json`
- `results/one_page_proof_brief.md`
- `results/final_lab_top3.json`
- `results/adaptive_ode_robust_summary.json`
- `results/adaptive_ode_robust_report.md`
- `results/test_report.txt`
- `results/full_pipeline_run.log`
- `results/lab_validation_plan.md`
- `results/lab_validation_matrix.csv`
- `results/lab_measurements_template.csv`
- `results/lab_pass_fail_rules.json`
- `results/power_analysis_summary.json`
- `results/power_analysis_details.json`
- `results/power_analysis_details.csv`
- `results/power_analysis_report.md`
- `results/adversarial_robustness_summary.json`
- `results/adversarial_robustness_scenarios.json`
- `results/adversarial_robustness_rows.json`
- `results/adversarial_robustness_candidate_summary.json`
- `results/adversarial_robustness_report.md`
- `results/transcriptome_calibration_summary.json`
- `results/transcriptome_calibration_expression_comparison.json`
- `results/transcriptome_calibration_curves.csv`
- `results/transcriptome_calibration_report.md`
- `config/parameters_calibrated.yaml`
- `results/calibration_impact_summary.json`
- `results/calibration_impact_baseline.json`
- `results/calibration_impact_calibrated.json`
- `results/calibration_impact_report.md`

