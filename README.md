# 🧬 BlackCotton — In Silico Genetic Engineering Pipeline

**Genetically engineering *Gossypium hirsutum* to grow natively black cotton fibers.**

## The Strategy

Instead of surface-dyeing white cotton (which fades), we engineer the plant to flood its fiber lumen with melanin **after** structural cellulose has fully formed. This is achieved via stage-specific genetic promoters that act as biological "if/then" switches.

```
0 DPA ────── 21 DPA ────── 35 DPA ────── 50 DPA
  │              │              │              │
  │  ELONGATION  │   CELLULOSE  │   MELANIN    │
  │  (grow long) │  (grow strong)│  (go black)  │
  └──────────────┴──────────────┴──────────────┘
       NO pigment genes active     Switch ON!
```

## Pipeline Components

| Module | Purpose |
|--------|---------|
| `construct_designer.py` | Assemble the T-DNA genetic construct |
| `codon_optimizer.py` | Optimize gene codons for *G. hirsutum* |
| `expression_model.py` | ODE-based gene expression simulation |
| `melanin_pathway.py` | Melanin biosynthesis kinetics |
| `fiber_model.py` | Predict fiber quality metrics |
| `tradeoff_optimizer.py` | Multi-objective search for blackness-strength-yield balance |
| `robustness_analyzer.py` | Uncertainty stress-test and robust candidate ranking |
| `robust_tradeoff_optimizer.py` | Robust search for high-success candidates under uncertainty |
| `sensitivity_analyzer.py` | Quantify which parameters drive failures/success under uncertainty |
| `adaptive_robust_optimizer.py` | Sensitivity-guided retuned robust optimization |
| `adaptive_ode_robust_pipeline.py` | ODE-refine adaptive shortlist and freeze final lab Top 3 |
| `lab_validation_planner.py` | Generate first-cycle lab matrix, measurement template, and pass/fail gates |
| `experiment_power_analyzer.py` | Estimate replicate counts + DPA timepoints needed to hit target confidence power |
| `adversarial_robustness_suite.py` | Red-team Top 3 against strict gates and high-noise worst-case scenarios |
| `transcriptome_calibrator.py` | Fit promoter timing/shape parameters from reference expression time-series |
| `calibration_impact_analyzer.py` | Compare baseline vs calibrated configs on full optimization + robustness metrics |
| `proof_packager.py` | Rebuild canonical proof docs + manifest + shareable bundle |
| `visualization.py` | Generate analysis plots |

## Quick Start

```bash
# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Run the full pipeline
python -m src.construct_designer    # Design the genetic construct
python -m src.expression_model      # Simulate gene expression
python -m src.melanin_pathway       # Model melanin biosynthesis
python -m src.fiber_model           # Predict fiber quality
python -m src.tradeoff_optimizer    # Safe-timing constrained sweep + ODE-refined candidates
python -m src.robustness_analyzer   # Monte Carlo robustness ranking under uncertainty
python -m src.robust_tradeoff_optimizer  # Robust optimization from coarse sweep to shortlist
python -m src.sensitivity_analyzer  # Feature sensitivity from robust trial records
python -m src.adaptive_robust_optimizer  # Retuned robust search based on sensitivity
python -m src.adaptive_ode_robust_pipeline  # ODE refine + robust rerank + final lab Top 3
python -m src.lab_validation_planner  # Build practical lab validation package from Top 3
python -m src.experiment_power_analyzer  # Compute recommended replicates/timepoints from trial distributions
python -m src.adversarial_robustness_suite  # Stress-test Top 3 under strict/worst-case scenarios
python -m src.transcriptome_calibrator  # Calibrate promoter parameters from reference expression data
python -m src.calibration_impact_analyzer  # Check if calibration improves end-to-end robustness
python -m src.proof_packager        # Refresh proof docs and evidence bundle
python -m src.visualization         # Generate all plots

# Run automated checks
python -m unittest discover -s tests -v
```

## Key Outputs

- `results/construct_map.txt` — Annotated genetic construct
- `results/expression_kinetics.png` — Gene expression over fiber development
- `results/melanin_accumulation.png` — Melanin yield prediction
- `results/fiber_quality_comparison.png` — Quality vs. conventional cotton
- `results/construct_summary.json` — Machine-readable construct data
- `results/optimization_refined.json` — ODE-refined candidate metrics
- `results/top_candidates.json` — Best blackness-strength-yield parameter regimes
- `results/robustness_top_candidates.json` — Uncertainty-aware ranking (success probability + risk)
- `results/robust_top_candidates.json` — Robust-optimized shortlist for lab prioritization
- `results/robust_optimization_report.md` — Human-readable robust optimization summary
- `results/sensitivity_report.md` — Ranked feature impacts + tuning directions
- `results/adaptive_robust_report.md` — Retuned robust optimization report
- `results/final_lab_top3.json` — Frozen ODE-refined, robustness-reranked lab shortlist
- `results/lab_validation_plan.md` — Phase-1 lab protocol brief + measurement checklist
- `results/lab_validation_matrix.csv` — Candidate/control arm matrix for the pilot
- `results/lab_measurements_template.csv` — Structured template for experimental data capture
- `results/lab_pass_fail_rules.json` — Explicit pass/fail gates for go/no-go decisions
- `results/power_analysis_summary.json` — Recommended replicates + DPA timepoints for target power
- `results/power_analysis_details.csv` — Endpoint-wise required replicate counts per candidate
- `results/power_analysis_report.md` — Human-readable experimental power planning report
- `results/adversarial_robustness_summary.json` — Run-level summary for red-team robustness checks
- `results/adversarial_robustness_candidate_summary.json` — Worst-case candidate ranking across scenarios
- `results/adversarial_robustness_report.md` — Scenario-wise stress-test report
- `results/transcriptome_calibration_summary.json` — Baseline vs calibrated promoter fit quality
- `results/transcriptome_calibration_report.md` — Human-readable calibration summary
- `config/parameters_calibrated.yaml` — Calibrated parameter file (baseline config untouched)
- `results/calibration_impact_summary.json` — End-to-end baseline vs calibrated comparison
- `results/calibration_impact_report.md` — Human-readable decision report (which config wins)
- `results/canonical_report.md` — Current snapshot report (auto-refreshed by proof packager)
- `results/proof_manifest.json` — Artifact checksums + metadata for verification
- `results/proof_bundle_YYYY-MM-DD.tgz` — Shareable evidence package

## The Genetic Construct

```
LB ─── [pGhMat1 → melA] ─── [pGhSCW-late → TYRP1] ─── [pGhSCW-late → DCT] ─── [p35S → nptII] ─── RB
        Maturation promoter    Late SCW promoter         Late SCW promoter        Selection marker
        fires at ~35 DPA       fires at ~30 DPA          fires at ~30 DPA         (Kanamycin resistance)
```

## License

Research use only. Patent pending on temporal decoupling strategy.
