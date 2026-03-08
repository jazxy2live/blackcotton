# Correlated Construct Audit

Objective: compare Anti-Silencing against Split-Support under strict high-noise stress with correlated failure bundles.

## Headline

- Baseline correlated top success: `0.764`
- Variant correlated top success: `0.771`
- Delta: `+0.007`
- Correlated profile: `construct_bundle_v1`

## Comparison

| Arm | Independent Top Success | Correlated Top Success | Correlated Top Robust | Fragility |
|---|---:|---:|---:|---:|
| Anti-Silencing | 0.871 | 0.764 | 0.804 | 0.072 |
| Split-Support | 0.857 | 0.771 | 0.810 | 0.067 |

## Split-Support Top 3

| Rank | Success | Robust | Fragility | p50 L* | p50 Strength | p50 Yield | Params |
|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 0.771 | 0.810 | 0.067 | 19.30 | 29.10 | 0.913 | 39/35, k=0.10, eff=1.70, ret=0.85 |
| 2 | 0.750 | 0.794 | 0.075 | 20.22 | 29.23 | 0.923 | 40/35, k=0.05, eff=1.70, ret=0.85 |
| 3 | 0.736 | 0.786 | 0.079 | 19.49 | 29.24 | 0.923 | 39/35, k=0.03, eff=1.70, ret=0.85 |

## Construct Stability Assumptions

- Variant silencing scale: `0.662`
- Variant event-CV scale: `0.830`

