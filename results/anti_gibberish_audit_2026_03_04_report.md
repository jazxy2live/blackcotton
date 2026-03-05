# Anti-Gibberish Audit

Focused audit: ablation + unseen holdouts + confidence scoring.

## Run Setup

- Generated: `2026-03-04T09:08:41+00:00`
- Strict trials/candidate: `120`
- Holdout trials/candidate: `140`
- Seed candidates per variant: `120`
- Min temporal gap filter: `0.50` days
- RNG seed: `20260304`
- Strict scenario: `strict_high_noise`
- Holdouts: `2`

## Ablation Results

| Variant | Traffic | Storage | Transit | Strict Top Success | Strict Top Robust | Holdout Worst Top Success |
|---|---:|---:|---:|---:|---:|---:|
| All Off (control) | 0 | 0 | 0 | 0.300 | 0.477 | 0.093 |
| Transit Only | 0 | 0 | 1 | 0.300 | 0.478 | 0.086 |
| Traffic Only | 1 | 0 | 0 | 0.292 | 0.473 | 0.114 |
| Storage Only (Transit+Compartment) | 0 | 1 | 1 | 0.767 | 0.805 | 0.464 |
| All On | 1 | 1 | 1 | 0.767 | 0.806 | 0.507 |

## Holdout Leaderboard (Top Success by Variant)

### holdout_ultra_noise_gate

| Variant | Top Success | Top Robust |
|---|---:|---:|
| All On | 0.507 | 0.624 |
| Storage Only (Transit+Compartment) | 0.486 | 0.609 |
| Traffic Only | 0.129 | 0.359 |
| All Off (control) | 0.107 | 0.342 |
| Transit Only | 0.086 | 0.329 |

### holdout_dark_yield_clamp

| Variant | Top Success | Top Robust |
|---|---:|---:|
| All On | 0.529 | 0.640 |
| Storage Only (Transit+Compartment) | 0.464 | 0.595 |
| Transit Only | 0.129 | 0.357 |
| Traffic Only | 0.114 | 0.347 |
| All Off (control) | 0.093 | 0.332 |

## Confidence

- Confidence rating: `HIGH`
- Combined strict gain vs control: `+0.467`
- Combined holdout worst-top gain vs control: `+0.414`
- Marginal traffic gain: `-0.008`
- Marginal storage gain (vs transit-only): `+0.467`
- Marginal transit-only gain: `+0.000`

### Strong Evidence

- Combined upgrades improve strict-high-noise top success by +0.467 versus control.
- Combined upgrades improve worst holdout top success by +0.414 versus control.
- Storage targeting (with transit) contributes positive gain (+0.467).

### Caveats

- Transit-only change has near-zero direct robustness impact unless paired with compartment/toxicity assumptions.

