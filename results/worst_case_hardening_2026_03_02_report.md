# Worst-Case Hardening Sprint Report

Objective: maximize strict+high-noise success instead of average-case robustness.

## Headline

- Strict+high-noise top success: `0.712` -> `0.942` (delta `+0.231`)
- Strict+high-noise top robust score: `0.764` -> `0.930` (delta `+0.166`)
- Target success (`>= 0.120`): `True`

## Scenario Delta

| Scenario | Baseline Top Success | Hardened Top Success | Delta |
|---|---:|---:|---:|
| baseline | 0.996 | 1.000 | +0.004 |
| strict_gate | 0.892 | 1.000 | +0.108 |
| high_noise | 0.885 | 0.988 | +0.104 |
| strict_high_noise | 0.712 | 0.954 | +0.242 |

## Hardened Top 3 (strict-high-noise ranked)

| Rank | Success | Robust | p50 L* | p50 Strength | p50 Yield | p50 Gap | Params |
|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 0.942 | 0.930 | 16.94 | 29.19 | 0.916 | 3.00 | 39/35, k=0.03, eff=1.70, ret=0.85 |
| 2 | 0.935 | 0.925 | 16.68 | 29.21 | 0.919 | 3.12 | 39/35, k=0.03, eff=1.70, ret=0.85 |
| 3 | 0.935 | 0.925 | 16.29 | 29.10 | 0.910 | 2.38 | 38/35, k=0.05, eff=1.70, ret=0.85 |

