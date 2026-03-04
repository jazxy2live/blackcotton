# Transcriptome Calibration Report

- Reference table: `data/reference/cotton_expression_reference.csv`
- Output calibrated config: `config/parameters_calibrated.yaml`

## Fit Quality

- Baseline mean RMSE: `0.2245`
- Calibrated mean RMSE: `0.0117`
- Delta (cal - base): `-0.2128`

## Promoter Updates

| Promoter | Activation DPA | Peak DPA | Hill | Leakage | Strength |
|---|---:|---:|---:|---:|---:|
| pGhSCW_late | 23.72 | 26.44 | 6.52 | 0.0339 | 1.016 |
| pGhMat1 | 31.13 | 43.85 | 9.75 | 0.0214 | 1.072 |

## Expression Timing Snapshot

- Baseline temporal gap: `2.875` days
- Raw calibrated temporal gap: `0.000` days
- Final calibrated temporal gap: `0.500` days (min required `0.500`)
- Baseline melA peak DPA: `50.00`, calibrated: `50.00`
- Safety adjustment applied: shifted `pGhMat1` by `0.625` days

