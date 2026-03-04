# Sensitivity Analysis Report

Global sensitivity derived from robust trial records.

## Summary

- Trial records analyzed: `360`
- Candidate ranks analyzed: `3`
- Overall success rate: `0.992`
- Main failure mode: `darkness_failure`

## Top Feature Impacts

| Rank | Feature | Impact | Corr Success | Corr Dark Fail | Corr Gap Fail | Recommended Direction |
|---:|---|---:|---:|---:|---:|---|
| 1 | sampled_mat_activation_dpa | 0.813 | 0.033 | 0.033 | -0.104 | increase |
| 2 | sampled_scw_activation_dpa | 0.673 | 0.074 | -0.011 | -0.112 | increase |
| 3 | sampled_k_competition | 0.609 | -0.047 | 0.014 | 0.060 | decrease |
| 4 | sampled_hill_melA | 0.480 | 0.108 | -0.092 | -0.057 | increase |
| 5 | sampled_melanin_efficiency | 0.473 | 0.105 | -0.135 | 0.010 | increase |
| 6 | sampled_scw_strength | 0.438 | 0.106 | -0.098 | -0.046 | increase |
| 7 | sampled_mat_strength | 0.410 | -0.022 | 0.052 | -0.035 | decrease |
| 8 | sampled_late_retention_factor | 0.232 | 0.064 | -0.089 | 0.016 | increase |
| 9 | sampled_leak_melA | 0.228 | -0.017 | 0.075 | -0.076 | increase |
| 10 | sampled_hill_scw | 0.221 | 0.066 | -0.093 | 0.017 | increase |

## Actionable Tuning

- `sampled_mat_activation_dpa`: `increase` (impact 0.813, dark-fail corr 0.033, gap-fail corr -0.104)
- `sampled_scw_activation_dpa`: `increase` (impact 0.673, dark-fail corr -0.011, gap-fail corr -0.112)
- `sampled_k_competition`: `decrease` (impact 0.609, dark-fail corr 0.014, gap-fail corr 0.060)
- `sampled_hill_melA`: `increase` (impact 0.480, dark-fail corr -0.092, gap-fail corr -0.057)
- `sampled_melanin_efficiency`: `increase` (impact 0.473, dark-fail corr -0.135, gap-fail corr 0.010)

