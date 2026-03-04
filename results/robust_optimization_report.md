# Robust Optimization Report

Uncertainty-aware ranking of black cotton candidates.

## Run Summary

- Coarse safe candidates: `12848`
- Robust seeds evaluated: `160`
- Trials per seed: `100`
- Best robust score: `0.244`
- Best success rate: `0.000`

## Top 12 (Robust)

| Rank | Success | Robust Score | p50 L* | p50 Strength | p50 Yield | Risk Gap | Risk Darkness | Params (mat/scw, k, eff, ret) |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 0.000 | 0.244 | 36.52 | 29.13 | 0.918 | 0.000 | 1.000 | 38/34, 0.10, 1.30, 0.45 |
| 2 | 0.000 | 0.244 | 37.18 | 29.16 | 0.920 | 0.000 | 0.990 | 40/34, 0.10, 1.30, 0.45 |
| 3 | 0.000 | 0.243 | 36.78 | 29.12 | 0.917 | 0.000 | 1.000 | 40/34, 0.10, 1.30, 0.45 |
| 4 | 0.000 | 0.243 | 36.90 | 29.15 | 0.920 | 0.000 | 1.000 | 40/34, 0.10, 1.30, 0.45 |
| 5 | 0.000 | 0.242 | 35.84 | 28.95 | 0.906 | 0.000 | 1.000 | 38/34, 0.20, 1.30, 0.45 |
| 6 | 0.000 | 0.242 | 36.56 | 28.98 | 0.907 | 0.000 | 1.000 | 40/34, 0.20, 1.30, 0.45 |
| 7 | 0.000 | 0.242 | 36.66 | 29.10 | 0.915 | 0.000 | 1.000 | 38/34, 0.10, 1.30, 0.45 |
| 8 | 0.000 | 0.242 | 36.28 | 29.11 | 0.916 | 0.000 | 1.000 | 38/34, 0.10, 1.30, 0.45 |
| 9 | 0.000 | 0.241 | 36.37 | 28.84 | 0.898 | 0.000 | 1.000 | 40/34, 0.30, 1.30, 0.45 |
| 10 | 0.000 | 0.241 | 37.18 | 29.08 | 0.914 | 0.000 | 1.000 | 36/34, 0.10, 1.30, 0.45 |
| 11 | 0.000 | 0.241 | 35.93 | 29.03 | 0.910 | 0.000 | 1.000 | 40/32, 0.10, 1.30, 0.45 |
| 12 | 0.000 | 0.241 | 35.91 | 29.07 | 0.913 | 0.000 | 1.000 | 40/32, 0.10, 1.30, 0.45 |

## Interpretation

- Robust winners are less dark in median conditions but survive uncertainty better.
- Main failure mode is still darkness variability; strength/yield are comparatively stable.
- Next optimization cycles should focus on reducing darkness-risk without creating timing overlap.

