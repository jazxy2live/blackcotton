# Canonical Pipeline Report

- Generated at: 2026-03-02T14:11:08+00:00
- Workspace: `/Users/chaitanyajongra/Desktop/blackcotton`
- Interpreter baseline: `./venv/bin/python` for all pipeline modules
- Lead candidate source priority: `final_lab_top3 > top_candidates`

## Key Metrics

- Construct: `pBC-MelaninCotton-v2`, length `8194` bp, cassettes `4`
- Expression: cellulose 90% at `38.375` DPA, melA peak at `50.0` DPA
- Melanin pathway: final melanin `0.49999080355024245` mM eq, tyrosine consumed `100.0`%
- Optimization: coarse `12848`, refined `112`, pareto `2`
- Final robust stage: ODE top `3`, trials/candidate `120`, best success `0.833`

## Lead Candidate Snapshot

- Stage `adaptive_ode_robust` from `final_lab_top3`: L* `21.587`, strength `29.350` g/tex, yield `0.934`, gap `5.500` days
- Robustness: success `0.267`, robust score `0.455`, p50 L* `21.911`, p50 strength `29.361`, p50 yield `0.936`
- Parameters: mat/scw `42/35`, k `0.03`, eff `1.70`, ret `0.85`

## Final Lab Top 3

| Rank | Success | Robust | p50 L* | p50 Str | p50 Yield | Params |
|---:|---:|---:|---:|---:|---:|---|
| 1 | 0.267 | 0.455 | 21.91 | 29.36 | 0.936 | 42/35, k=0.03, eff=1.70, ret=0.85 |
| 2 | 0.267 | 0.454 | 21.82 | 29.25 | 0.928 | 42/34, k=0.12, eff=1.70, ret=0.85 |
| 3 | 0.263 | 0.452 | 22.44 | 29.31 | 0.932 | 42/35, k=0.07, eff=1.60, ret=0.85 |

