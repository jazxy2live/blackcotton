# Lab Validation Plan — Final Top 3

- Run ID: `BC-LAB-2026-03-02`
- Generated (UTC): `2026-03-02T14:10:35+00:00`
- Candidate replicates per arm: `8`

## Objective

Validate whether the simulated Top 3 preserve quality while reaching dark fiber under timing-safe conditions.

## Phase 1 Pass/Fail Gate

- `color_L <= 25.00`
- `strength_g_tex >= 28.00`
- `yield_index >= 0.850`
- `temporal_gap_days >= 0.00`
- `toxicity_pre_cellulose <= 0.090`

## Candidate Arms

| Arm | Pred Success | Pred Robust | Pred p50 L* | Pred p50 Strength | Pred p50 Yield | Timing (mat/scw) | Key Params |
|---|---:|---:|---:|---:|---:|---|---|
| CAND_1 | 0.267 | 0.455 | 21.91 | 29.36 | 0.936 | 42/35 DPA | k=0.03, eff=1.70, ret=0.85 |
| CAND_2 | 0.267 | 0.454 | 21.82 | 29.25 | 0.928 | 42/34 DPA | k=0.12, eff=1.70, ret=0.85 |
| CAND_3 | 0.263 | 0.452 | 22.44 | 29.31 | 0.932 | 42/35 DPA | k=0.07, eff=1.60, ret=0.85 |

## Required Measurements

- Expression timing: melA, TYRP1, DCT, cellulose markers across DPA stages
- Pigment endpoint: melanin content and color L*
- Quality endpoint: strength (g/tex), yield index
- Durability endpoint: color L* shift after wash cycles

## Generated Files

- `results/lab_validation_matrix.csv`
- `results/lab_measurements_template.csv`
- `results/lab_pass_fail_rules.json`

