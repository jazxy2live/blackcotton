# Adaptive ODE Robust Final Report

Finalization pipeline: adaptive robust shortlist -> ODE refinement -> robustness -> lab Top 3.

## Summary

- Input adaptive candidates: `12`
- ODE refined candidates: `12`
- ODE pareto candidates: `1`
- ODE top candidates for robustness: `6`
- Robust trials per candidate: `120`
- Best robust success rate: `0.883`
- Best robust score: `0.880`

## Final Lab Top 3

| Rank | Success | Robust | p50 L* | p50 Str | p50 Yield | Risk Gap | Risk Dark | Params |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 0.883 | 0.880 | 22.16 | 28.72 | 0.886 | 0.000 | 0.117 | 38/34, k=0.20, eff=1.60, ret=0.75 |
| 2 | 0.858 | 0.864 | 22.36 | 28.95 | 0.902 | 0.008 | 0.133 | 38/34, k=0.10, eff=1.60, ret=0.75 |
| 3 | 0.858 | 0.864 | 22.28 | 28.90 | 0.899 | 0.000 | 0.142 | 38/34, k=0.15, eff=1.60, ret=0.75 |

