# BlackCotton One-Page Proof Brief

Date: 2026-03-02
Project: In-silico design for non-fading black cotton via delayed melanin biosynthesis

## 1) Core Problem
Dark cotton typically loses fiber quality because pigment and cellulose formation overlap.

## 2) Hypothesis
Delay pigment pathway activation until late fiber development to keep strength while achieving deep black color.

## 3) What Exists Now
- End-to-end simulation pipeline from construct design to robust candidate ranking
- Hard non-negative temporal gap constraint in optimizer
- ODE refinement + uncertainty stress-test + frozen final lab Top 3

## 4) Current Evidence (Reproducible)
- Coarse candidates scanned: `12848`
- ODE refined candidates: `112`
- ODE top candidates stress-tested: `3`
- Trials per candidate: `120`
- Best robust success rate: `0.833`
- Lead candidate: L* `21.587`, strength `29.350 g/tex`, yield `0.934`, gap `5.500 d`
- Tests: `28` run, status `pass`

## 5) Why This Matters
This converts trial-and-error into ranked, uncertainty-aware experiments so lab cycles are faster and cheaper.

## 6) Next Step
Validate the final Top 3 in lab for expression timing, darkness durability, and fiber quality.

## 7) Proof Files
- `results/canonical_report.md`
- `results/proof_manifest.json`
- `results/proof_package.md`
- `results/final_lab_top3.json`
- `results/adaptive_ode_robust_report.md`

