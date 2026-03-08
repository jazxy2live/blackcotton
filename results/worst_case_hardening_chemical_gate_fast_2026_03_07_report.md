# Worst-Case Hardening Sprint Report

Objective: maximize strict+high-noise success instead of average-case robustness.

## Headline

- Strict+high-noise top success: `0.433` -> `0.850` (delta `+0.417`)
- Strict+high-noise top robust score: `0.569` -> `0.865` (delta `+0.296`)
- Target success (`>= 0.850`): `True`
- Chemical-safe seeds: `180` / `180`

## Scenario Delta

| Scenario | Baseline Top Success | Hardened Top Success | Delta |
|---|---:|---:|---:|
| baseline | 0.958 | 1.000 | +0.042 |
| strict_gate | 0.667 | 0.975 | +0.308 |
| high_noise | 0.800 | 0.942 | +0.142 |
| strict_high_noise | 0.433 | 0.833 | +0.400 |

## Hardened Top 3 (strict-high-noise ranked)

| Rank | Success | Robust | Fragility | p50 L* | p50 Strength | p50 Yield | p50 Gap | Params |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 0.850 | 0.865 | 0.047 | 19.42 | 29.22 | 0.921 | 2.88 | 39/35, k=0.03, eff=1.70, ret=0.85 |
| 2 | 0.833 | 0.854 | 0.049 | 19.20 | 29.18 | 0.918 | 3.19 | 39/35, k=0.05, eff=1.70, ret=0.85 |
| 3 | 0.825 | 0.849 | 0.053 | 18.31 | 29.14 | 0.914 | 2.00 | 38/35, k=0.03, eff=1.70, ret=0.85 |

