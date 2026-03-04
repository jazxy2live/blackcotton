# Calibration Impact Report

End-to-end comparison of baseline vs calibrated parameters.

## Decision

- Winner: `calibrated`
- Reason: `better worst-case success/robust score`

## Optimization Footprint

| Metric | Baseline | Calibrated | Delta (cal-base) |
|---|---:|---:|---:|
| n_coarse | 14080 | 21952 | 7872.000 |
| n_seed | 96 | 96 | 0.000 |
| n_refined | 95 | 96 | 1.000 |
| n_pareto | 3 | 3 | 0.000 |
| n_top | 12 | 12 | 0.000 |

## Robustness Comparison

| Metric | Baseline | Calibrated | Delta (cal-base) |
|---|---:|---:|---:|
| default.success_rate | 0.142 | 0.192 | +0.050 |
| default.robust_score | 0.346 | 0.382 | +0.036 |
| default.p50_color_L | 27.361 | 27.743 | +0.382 |
| default.p50_strength_g_tex | 28.824 | 28.865 | +0.041 |
| default.p50_yield_index | 0.890 | 0.894 | +0.004 |
| default.p50_temporal_gap_days | 0.250 | 1.000 | +0.750 |
| worst.success_rate | 0.033 | 0.075 | +0.042 |
| worst.robust_score | 0.267 | 0.298 | +0.032 |
| worst.p50_color_L | 30.775 | 28.877 | -1.898 |
| worst.p50_strength_g_tex | 29.004 | 28.904 | -0.100 |
| worst.p50_yield_index | 0.906 | 0.897 | -0.009 |
| worst.p50_temporal_gap_days | 2.125 | 1.375 | -0.750 |

## Expression Timing Snapshot

- `cellulose_90pct_dpa` baseline `38.375` vs calibrated `38.375` (delta `+0.000`)
- `melA_halfmax_dpa` baseline `41.250` vs calibrated `38.875` (delta `-2.375`)
- `temporal_gap_days` baseline `2.875` vs calibrated `0.500` (delta `-2.375`)

