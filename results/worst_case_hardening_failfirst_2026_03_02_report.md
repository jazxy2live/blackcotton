# Worst-Case Hardening Sprint Report

Objective: maximize strict+high-noise success instead of average-case robustness.

## Headline

- Strict+high-noise top success: `0.000` -> `0.000` (delta `+0.000`)
- Strict+high-noise top robust score: `0.221` -> `0.000` (delta `-0.221`)
- Target success (`>= 0.120`): `False`

## Scenario Delta

| Scenario | Baseline Top Success | Hardened Top Success | Delta |
|---|---:|---:|---:|
| baseline | 0.000 | 0.000 | +0.000 |
| strict_gate | 0.000 | 0.000 | +0.000 |
| high_noise | 0.000 | 0.000 | +0.000 |
| strict_high_noise | 0.000 | 0.000 | +0.000 |

## Hardened Top 3 (strict-high-noise ranked)

| Rank | Success | Robust | Fragility | p50 L* | p50 Strength | p50 Yield | p50 Gap | Params |
|---:|---:|---:|---:|---:|---:|---:|---:|---|

