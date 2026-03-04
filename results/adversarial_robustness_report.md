# Adversarial Robustness Report

Red-team robustness check across strict thresholds and high-noise scenarios.

## Run Summary

- Generated: `2026-03-02T14:11:10+00:00`
- Candidates tested: `3`
- Scenarios: `4`
- Trials per scenario/candidate: `240`
- Worst-case best candidate input rank: `3`

## Scenario Leaderboard

| Scenario | Top Candidate Rank | Top Success | Top Robust |
|---|---:|---:|---:|
| baseline | 3 | 0.796 | 0.824 |
| strict_gate | 3 | 0.292 | 0.474 |
| high_noise | 1 | 0.600 | 0.687 |
| strict_high_noise | 3 | 0.242 | 0.438 |

## Worst-Case Candidate Ranking

| Candidate Rank | Min Success | Mean Success | Min Robust | Worst Scenario |
|---:|---:|---:|---:|---|
| 3 | 0.242 | 0.480 | 0.438 | strict_high_noise |
| 2 | 0.229 | 0.454 | 0.429 | strict_high_noise |
| 1 | 0.171 | 0.410 | 0.388 | strict_high_noise |

