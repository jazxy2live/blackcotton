# Worst-Case Hardening Sprint Report

Objective: maximize strict+high-noise success instead of average-case robustness.

## Headline

- Strict+high-noise top success: `0.445` -> `0.827` (delta `+0.382`)
- Strict+high-noise top robust score: `0.578` -> `0.849` (delta `+0.271`)
- Target success (`>= 0.850`): `False`
- Chemical-safe seeds: `260` / `260`

## Scenario Delta

| Scenario | Baseline Top Success | Hardened Top Success | Delta |
|---|---:|---:|---:|
| baseline | 0.964 | 1.000 | +0.036 |
| strict_gate | 0.700 | 0.973 | +0.273 |
| high_noise | 0.750 | 0.941 | +0.191 |
| strict_high_noise | 0.455 | 0.805 | +0.350 |

## Hardened Top 3 (strict-high-noise ranked)

| Rank | Success | Robust | Fragility | p50 L* | p50 Strength | p50 Yield | p50 Gap | Params |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 0.827 | 0.849 | 0.050 | 19.65 | 29.19 | 0.919 | 2.75 | 39/34, k=0.03, eff=1.70, ret=0.85 |
| 2 | 0.814 | 0.840 | 0.057 | 19.27 | 29.22 | 0.921 | 3.12 | 39/35, k=0.03, eff=1.70, ret=0.85 |
| 3 | 0.809 | 0.838 | 0.056 | 18.53 | 29.08 | 0.910 | 2.06 | 38/35, k=0.05, eff=1.70, ret=0.85 |

