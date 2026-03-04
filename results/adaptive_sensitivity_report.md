# Sensitivity Analysis Report

Global sensitivity derived from robust trial records.

## Summary

- Trial records analyzed: `31200`
- Candidate ranks analyzed: `260`
- Overall success rate: `0.963`
- Main failure mode: `temporal_overlap_failure`

## Top Feature Impacts

| Rank | Feature | Impact | Corr Success | Corr Dark Fail | Corr Gap Fail | Recommended Direction |
|---:|---|---:|---:|---:|---:|---|
| 1 | sampled_mat_activation_dpa | 0.825 | 0.116 | 0.062 | -0.197 | increase |
| 2 | sampled_k_competition | 0.645 | -0.009 | -0.043 | 0.004 | neutral |
| 3 | sampled_scw_activation_dpa | 0.510 | 0.103 | 0.013 | -0.135 | increase |
| 4 | sampled_melanin_efficiency | 0.502 | 0.112 | -0.184 | -0.009 | increase |
| 5 | sampled_late_retention_factor | 0.428 | 0.126 | -0.168 | -0.029 | increase |
| 6 | sampled_hill_melA | 0.205 | 0.045 | -0.029 | -0.027 | increase |
| 7 | sampled_scw_strength | 0.141 | -0.026 | 0.014 | 0.022 | neutral |
| 8 | sampled_hill_scw | 0.122 | 0.005 | -0.030 | 0.027 | neutral |
| 9 | sampled_leak_scw | 0.047 | -0.014 | 0.007 | 0.008 | neutral |
| 10 | sampled_leak_melA | 0.044 | -0.009 | 0.008 | 0.004 | neutral |

## Actionable Tuning

- `sampled_mat_activation_dpa`: `increase` (impact 0.825, dark-fail corr 0.062, gap-fail corr -0.197)
- `sampled_k_competition`: `neutral` (impact 0.645, dark-fail corr -0.043, gap-fail corr 0.004)
- `sampled_scw_activation_dpa`: `increase` (impact 0.510, dark-fail corr 0.013, gap-fail corr -0.135)
- `sampled_melanin_efficiency`: `increase` (impact 0.502, dark-fail corr -0.184, gap-fail corr -0.009)
- `sampled_late_retention_factor`: `increase` (impact 0.428, dark-fail corr -0.168, gap-fail corr -0.029)

