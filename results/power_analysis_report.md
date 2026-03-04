# Experiment Power Analysis Report

Replicate planning from robust trial distributions and one-sided confidence gates.

## Summary

- Generated: `2026-03-02T14:10:58+00:00`
- Trial records used: `360`
- Candidates analyzed: `3`
- Current planned replicates: `8`
- Target power: `0.80`
- Confidence alpha (one-sided): `0.05`
- Statistical minimum replicates/candidate: `{'CAND_1': 2, 'CAND_2': 3, 'CAND_3': 3}`
- Statistical minimum overall: `3`
- Practical minimum floor: `6`
- Practical recommendation/candidate: `{'CAND_1': 6, 'CAND_2': 6, 'CAND_3': 6}`
- Practical recommendation overall: `6`
- Recommended DPA timepoints: `[32, 33, 34, 35, 36, 37, 40, 42, 44, 46]`

## Endpoint-Level Requirements

| Candidate | Endpoint | Threshold | Mean | SD | Power @ Current n | Required n | Achieved Power |
|---|---|---:|---:|---:|---:|---:|---:|
| CAND_1 | color_L | 25.000 | 21.531 | 1.978 | 0.991 | 2 | 0.859 |
| CAND_1 | strength_g_tex | 28.000 | 29.343 | 0.047 | 1.000 | 2 | 1.000 |
| CAND_1 | yield_index | 0.850 | 0.934 | 0.005 | 1.000 | 2 | 1.000 |
| CAND_1 | temporal_gap_days | 0.000 | 5.569 | 1.029 | 1.000 | 2 | 1.000 |
| CAND_2 | color_L | 25.000 | 21.712 | 1.876 | 0.999 | 3 | 0.857 |
| CAND_2 | strength_g_tex | 28.000 | 29.252 | 0.072 | 1.000 | 2 | 1.000 |
| CAND_2 | yield_index | 0.850 | 0.928 | 0.006 | 1.000 | 2 | 1.000 |
| CAND_2 | temporal_gap_days | 0.000 | 5.693 | 0.872 | 1.000 | 2 | 1.000 |
| CAND_3 | color_L | 25.000 | 22.188 | 1.802 | 0.997 | 3 | 0.827 |
| CAND_3 | strength_g_tex | 28.000 | 29.317 | 0.067 | 1.000 | 2 | 1.000 |
| CAND_3 | yield_index | 0.850 | 0.933 | 0.006 | 1.000 | 2 | 1.000 |
| CAND_3 | temporal_gap_days | 0.000 | 5.467 | 0.875 | 1.000 | 2 | 1.000 |

