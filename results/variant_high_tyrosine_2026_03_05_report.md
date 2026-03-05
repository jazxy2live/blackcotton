# Worst-Case Hardening Sprint Report

Objective: maximize strict+high-noise success instead of average-case robustness.

## Headline

- Strict+high-noise top success: `0.392` -> `0.733` (delta `+0.342`)
- Strict+high-noise top robust score: `0.540` -> `0.783` (delta `+0.243`)
- Target success (`>= 0.850`): `False`

## Scenario Delta

| Scenario | Baseline Top Success | Hardened Top Success | Delta |
|---|---:|---:|---:|
| baseline | 0.875 | 1.000 | +0.125 |
| strict_gate | 0.425 | 0.892 | +0.467 |
| high_noise | 0.608 | 0.892 | +0.283 |
| strict_high_noise | 0.317 | 0.650 | +0.333 |

## Hardened Top 3 (strict-high-noise ranked)

| Rank | Success | Robust | Fragility | p50 L* | p50 Strength | p50 Yield | p50 Gap | Params |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 0.733 | 0.783 | 0.080 | 19.81 | 29.16 | 0.916 | 2.12 | 38/35, k=0.03, eff=1.70, ret=0.85 |
| 2 | 0.717 | 0.772 | 0.087 | 19.53 | 29.16 | 0.917 | 2.31 | 38/35, k=0.03, eff=1.70, ret=0.85 |
| 3 | 0.717 | 0.771 | 0.085 | 19.63 | 29.16 | 0.917 | 1.88 | 38/35, k=0.03, eff=1.70, ret=0.85 |

