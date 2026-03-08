# Correlated Construct Audit

Objective: compare frozen production against an anti-silencing construct under strict high-noise stress with correlated failure bundles.

## Headline

- Baseline correlated top success: `0.708`
- Anti-silencing correlated top success: `0.800`
- Delta: `+0.092`
- Correlated profile: `construct_bundle_v1`

## Comparison

| Arm | Independent Top Success | Correlated Top Success | Correlated Top Robust | Fragility |
|---|---:|---:|---:|---:|
| Production | 0.842 | 0.708 | 0.766 | 0.086 |
| Anti-Silencing | 0.883 | 0.800 | 0.830 | 0.061 |

## Anti-Silencing Top 3

| Rank | Success | Robust | Fragility | p50 L* | p50 Strength | p50 Yield | Params |
|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 0.800 | 0.830 | 0.061 | 18.88 | 29.10 | 0.910 | 38/35, k=0.05, eff=1.70, ret=0.85 |
| 2 | 0.767 | 0.808 | 0.076 | 18.65 | 29.04 | 0.907 | 38/34, k=0.05, eff=1.70, ret=0.85 |
| 3 | 0.742 | 0.790 | 0.074 | 18.93 | 29.13 | 0.913 | 38/35, k=0.03, eff=1.70, ret=0.85 |

## Construct Stability Assumptions

- Anti-silencing silencing scale: `0.662`
- Anti-silencing event-CV scale: `0.830`

