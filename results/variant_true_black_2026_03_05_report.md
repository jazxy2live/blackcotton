# Worst-Case Hardening Sprint Report

Objective: maximize strict+high-noise success instead of average-case robustness.

## Headline

- Strict+high-noise top success: `0.525` -> `0.825` (delta `+0.300`)
- Strict+high-noise top robust score: `0.635` -> `0.848` (delta `+0.213`)
- Target success (`>= 0.850`): `False`

## Scenario Delta

| Scenario | Baseline Top Success | Hardened Top Success | Delta |
|---|---:|---:|---:|
| baseline | 0.967 | 1.000 | +0.033 |
| strict_gate | 0.708 | 0.983 | +0.275 |
| high_noise | 0.692 | 0.967 | +0.275 |
| strict_high_noise | 0.417 | 0.758 | +0.342 |

## Hardened Top 3 (strict-high-noise ranked)

| Rank | Success | Robust | Fragility | p50 L* | p50 Strength | p50 Yield | p50 Gap | Params |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 0.825 | 0.848 | 0.053 | 18.81 | 29.17 | 0.917 | 2.38 | 38/35, k=0.03, eff=1.70, ret=0.85 |
| 2 | 0.817 | 0.841 | 0.054 | 19.35 | 29.04 | 0.908 | 2.12 | 38/35, k=0.07, eff=1.70, ret=0.85 |
| 3 | 0.800 | 0.830 | 0.061 | 19.56 | 29.15 | 0.915 | 2.69 | 39/35, k=0.05, eff=1.70, ret=0.85 |

