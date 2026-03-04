# Worst-Case Hardening Sprint Report

Objective: maximize strict+high-noise success instead of average-case robustness.

## Headline

- Strict+high-noise top success: `0.017` -> `0.267` (delta `+0.250`)
- Strict+high-noise top robust score: `0.271` -> `0.455` (delta `+0.185`)
- Target success (`>= 0.120`): `True`

## Scenario Delta

| Scenario | Baseline Top Success | Hardened Top Success | Delta |
|---|---:|---:|---:|
| baseline | 0.113 | 0.779 | +0.667 |
| strict_gate | 0.000 | 0.342 | +0.342 |
| high_noise | 0.150 | 0.600 | +0.450 |
| strict_high_noise | 0.025 | 0.217 | +0.192 |

## Hardened Top 3 (strict-high-noise ranked)

| Rank | Success | Robust | Fragility | p50 L* | p50 Strength | p50 Yield | p50 Gap | Params |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 0.267 | 0.455 | 0.227 | 21.91 | 29.36 | 0.936 | 5.94 | 42/35, k=0.03, eff=1.70, ret=0.85 |
| 2 | 0.267 | 0.454 | 0.224 | 21.82 | 29.25 | 0.928 | 5.62 | 42/34, k=0.12, eff=1.70, ret=0.85 |
| 3 | 0.263 | 0.452 | 0.231 | 22.44 | 29.31 | 0.932 | 5.62 | 42/35, k=0.07, eff=1.60, ret=0.85 |

