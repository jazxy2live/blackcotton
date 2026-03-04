# Sensitivity Analysis Report

Global sensitivity derived from robust trial records.

## Summary

- Trial records analyzed: `16000`
- Candidate ranks analyzed: `160`
- Overall success rate: `0.000`
- Main failure mode: `darkness_failure`

## Top Feature Impacts

| Rank | Feature | Impact | Corr Success | Corr Dark Fail | Corr Gap Fail | Recommended Direction |
|---:|---|---:|---:|---:|---:|---|
| 1 | sampled_k_competition | 0.686 | 0.000 | -0.002 | 0.016 | neutral |
| 2 | sampled_mat_activation_dpa | 0.634 | 0.000 | 0.016 | -0.305 | increase |
| 3 | sampled_melanin_efficiency | 0.397 | 0.000 | -0.089 | -0.086 | increase |
| 4 | sampled_scw_activation_dpa | 0.207 | 0.000 | 0.009 | -0.012 | neutral |
| 5 | sampled_late_retention_factor | 0.189 | 0.000 | -0.017 | 0.019 | neutral |
| 6 | sampled_copper_loading_fraction | 0.117 | 0.000 | -0.033 | 0.007 | neutral |
| 7 | sampled_tyrosinase_activation_fraction | 0.108 | 0.000 | -0.027 | -0.004 | neutral |
| 8 | sampled_hill_melA | 0.107 | 0.000 | -0.006 | 0.025 | neutral |
| 9 | sampled_scw_strength | 0.086 | 0.000 | 0.008 | -0.005 | neutral |
| 10 | sampled_silencing_probability | 0.080 | 0.000 | 0.026 | 0.002 | neutral |

## Actionable Tuning

- `sampled_k_competition`: `neutral` (impact 0.686, dark-fail corr -0.002, gap-fail corr 0.016)
- `sampled_mat_activation_dpa`: `increase` (impact 0.634, dark-fail corr 0.016, gap-fail corr -0.305)
- `sampled_melanin_efficiency`: `increase` (impact 0.397, dark-fail corr -0.089, gap-fail corr -0.086)
- `sampled_scw_activation_dpa`: `neutral` (impact 0.207, dark-fail corr 0.009, gap-fail corr -0.012)
- `sampled_late_retention_factor`: `neutral` (impact 0.189, dark-fail corr -0.017, gap-fail corr 0.019)

