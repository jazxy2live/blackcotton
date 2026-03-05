# Anti-Gibberish Audit

Focused audit: ablation + unseen holdouts + confidence scoring.

## Run Setup

- Generated: `2026-03-04T09:14:45+00:00`
- Strict trials/candidate: `320`
- Holdout trials/candidate: `420`
- Seed candidates per variant: `140`
- Min temporal gap filter: `0.50` days
- RNG seed: `20260304`
- Strict scenario: `strict_high_noise`
- Holdouts: `2`

## Ablation Results

| Variant | Traffic | Storage | Transit | Strict Top Success | Strict Top Robust | Holdout Worst Top Success |
|---|---:|---:|---:|---:|---:|---:|
| All Off (control) | 0 | 0 | 0 | 0.266 | 0.454 | 0.100 |
| Transit Only | 0 | 0 | 1 | 0.247 | 0.442 | 0.114 |
| Traffic Only | 1 | 0 | 0 | 0.275 | 0.460 | 0.117 |
| Storage Only (Transit+Compartment) | 0 | 1 | 1 | 0.703 | 0.761 | 0.493 |
| All On | 1 | 1 | 1 | 0.709 | 0.765 | 0.483 |

## Holdout Leaderboard (Top Success by Variant)

### holdout_ultra_noise_gate

| Variant | Top Success | Top Robust |
|---|---:|---:|
| Storage Only (Transit+Compartment) | 0.498 | 0.617 |
| All On | 0.483 | 0.607 |
| Traffic Only | 0.157 | 0.377 |
| All Off (control) | 0.129 | 0.357 |
| Transit Only | 0.121 | 0.351 |

### holdout_dark_yield_clamp

| Variant | Top Success | Top Robust |
|---|---:|---:|
| All On | 0.498 | 0.618 |
| Storage Only (Transit+Compartment) | 0.493 | 0.614 |
| Traffic Only | 0.117 | 0.350 |
| Transit Only | 0.114 | 0.348 |
| All Off (control) | 0.100 | 0.337 |

## Confidence

- Confidence rating: `HIGH`
- Combined strict gain vs control: `+0.444`
- Combined holdout worst-top gain vs control: `+0.383`
- Marginal traffic gain: `+0.009`
- Marginal storage gain (vs transit-only): `+0.456`
- Marginal transit-only gain: `-0.019`

### Strong Evidence

- Combined upgrades improve strict-high-noise top success by +0.444 versus control.
- Combined upgrades improve worst holdout top success by +0.383 versus control.
- Traffic tuning alone contributes positive gain (+0.009).
- Storage targeting (with transit) contributes positive gain (+0.456).

### Caveats

- Transit-only change has near-zero direct robustness impact unless paired with compartment/toxicity assumptions.

