# Worst-Case Hardening Sprint Report

Objective: maximize strict+high-noise success instead of average-case robustness.

## Headline

- Strict+high-noise top success: `0.318` -> `0.727` (delta `+0.409`)
- Strict+high-noise top robust score: `0.487` -> `0.778` (delta `+0.291`)
- Target success (`>= 0.120`): `True`

## Scenario Delta

| Scenario | Baseline Top Success | Hardened Top Success | Delta |
|---|---:|---:|---:|
| baseline | 0.850 | 1.000 | +0.150 |
| strict_gate | 0.423 | 0.905 | +0.482 |
| high_noise | 0.686 | 0.882 | +0.195 |
| strict_high_noise | 0.341 | 0.705 | +0.364 |

## Hardened Top 3 (strict-high-noise ranked)

| Rank | Success | Robust | Fragility | p50 L* | p50 Strength | p50 Yield | p50 Gap | Params |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 0.727 | 0.778 | 0.085 | 19.72 | 29.08 | 0.911 | 2.19 | 38/35, k=0.07, eff=1.70, ret=0.85 |
| 2 | 0.718 | 0.772 | 0.089 | 19.95 | 29.19 | 0.919 | 2.38 | 38/35, k=0.03, eff=1.70, ret=0.85 |
| 3 | 0.709 | 0.767 | 0.088 | 19.60 | 29.16 | 0.916 | 2.12 | 38/35, k=0.03, eff=1.70, ret=0.85 |

