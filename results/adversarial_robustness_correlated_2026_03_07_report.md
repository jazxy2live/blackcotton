# Adversarial Robustness Report

Red-team robustness check across strict thresholds and high-noise scenarios.

## Run Summary

- Generated: `2026-03-07T09:45:25+00:00`
- Candidates tested: `3`
- Scenarios: `8`
- Trials per scenario/candidate: `240`
- Correlated profile: `construct_bundle_v1`
- Worst-case best candidate input rank: `2`

## Scenario Leaderboard

| Scenario | Top Candidate Rank | Top Success | Top Robust |
|---|---:|---:|---:|
| baseline | 2 | 0.971 | 0.945 |
| strict_gate | 2 | 0.662 | 0.732 |
| high_noise | 1 | 0.729 | 0.775 |
| strict_high_noise | 3 | 0.433 | 0.570 |
| baseline_correlated | 2 | 0.838 | 0.851 |
| strict_gate_correlated | 2 | 0.487 | 0.608 |
| high_noise_correlated | 3 | 0.608 | 0.689 |
| strict_high_noise_correlated | 2 | 0.362 | 0.518 |

## Worst-Case Candidate Ranking

| Candidate Rank | Min Success | Mean Success | Min Robust | Worst Scenario |
|---:|---:|---:|---:|---|
| 2 | 0.362 | 0.626 | 0.518 | strict_high_noise_correlated |
| 3 | 0.287 | 0.610 | 0.465 | strict_high_noise_correlated |
| 1 | 0.233 | 0.517 | 0.426 | strict_high_noise_correlated |

