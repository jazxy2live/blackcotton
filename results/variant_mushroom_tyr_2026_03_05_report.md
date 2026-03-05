# Worst-Case Hardening Sprint Report

Objective: maximize strict+high-noise success instead of average-case robustness.

## Headline

- Strict+high-noise top success: `0.433` -> `0.825` (delta `+0.392`)
- Strict+high-noise top robust score: `0.570` -> `0.847` (delta `+0.278`)
- Target success (`>= 0.850`): `False`

## Scenario Delta

| Scenario | Baseline Top Success | Hardened Top Success | Delta |
|---|---:|---:|---:|
| baseline | 0.933 | 1.000 | +0.067 |
| strict_gate | 0.525 | 0.933 | +0.408 |
| high_noise | 0.658 | 0.933 | +0.275 |
| strict_high_noise | 0.333 | 0.700 | +0.367 |

## Hardened Top 3 (strict-high-noise ranked)

| Rank | Success | Robust | Fragility | p50 L* | p50 Strength | p50 Yield | p50 Gap | Params |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 0.825 | 0.847 | 0.050 | 19.43 | 29.12 | 0.913 | 2.38 | 38/35, k=0.05, eff=1.70, ret=0.85 |
| 2 | 0.792 | 0.824 | 0.060 | 19.28 | 29.15 | 0.916 | 2.12 | 38/35, k=0.03, eff=1.70, ret=0.85 |
| 3 | 0.767 | 0.807 | 0.071 | 19.05 | 29.13 | 0.915 | 2.44 | 39/34, k=0.05, eff=1.70, ret=0.85 |

