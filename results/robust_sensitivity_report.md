# Sensitivity Analysis Report

Global sensitivity derived from robust trial records.

## Summary

- Trial records analyzed: `16000`
- Candidate ranks analyzed: `160`
- Overall success rate: `0.042`
- Main failure mode: `darkness_failure`

## Top Feature Impacts

| Rank | Feature | Impact | Corr Success | Corr Dark Fail | Corr Gap Fail | Recommended Direction |
|---:|---|---:|---:|---:|---:|---|
| 1 | sampled_melanin_efficiency | 1.292 | 0.318 | -0.481 | 0.033 | increase |
| 2 | sampled_k_competition | 0.978 | -0.086 | -0.091 | 0.051 | increase |
| 3 | sampled_mat_activation_dpa | 0.743 | -0.014 | 0.120 | -0.404 | increase |
| 4 | sampled_scw_activation_dpa | 0.478 | 0.041 | 0.023 | -0.075 | increase |
| 5 | sampled_late_retention_factor | 0.444 | 0.097 | -0.140 | -0.018 | increase |
| 6 | sampled_hill_melA | 0.246 | 0.064 | -0.079 | -0.016 | increase |
| 7 | sampled_mat_strength | 0.207 | 0.017 | -0.039 | 0.069 | decrease |
| 8 | sampled_hill_scw | 0.175 | 0.024 | -0.051 | 0.052 | decrease |
| 9 | sampled_scw_strength | 0.097 | 0.007 | -0.034 | 0.024 | neutral |
| 10 | sampled_leak_melA | 0.073 | -0.020 | 0.006 | 0.014 | neutral |

## Actionable Tuning

- `sampled_melanin_efficiency`: `increase` (impact 1.292, dark-fail corr -0.481, gap-fail corr 0.033)
- `sampled_k_competition`: `increase` (impact 0.978, dark-fail corr -0.091, gap-fail corr 0.051)
- `sampled_mat_activation_dpa`: `increase` (impact 0.743, dark-fail corr 0.120, gap-fail corr -0.404)
- `sampled_scw_activation_dpa`: `increase` (impact 0.478, dark-fail corr 0.023, gap-fail corr -0.075)
- `sampled_late_retention_factor`: `increase` (impact 0.444, dark-fail corr -0.140, gap-fail corr -0.018)

