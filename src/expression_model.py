#!/usr/bin/env python3
"""
expression_model.py — ODE-based Gene Expression Simulation
============================================================

Models the temporal expression of transgenes in cotton fibers using
ordinary differential equations (ODEs). Simulates:

  1. Promoter activation kinetics (stage-specific switches)
  2. mRNA transcription and degradation
  3. Protein translation and degradation
  4. All three melanin genes: melA, TYRP1, DCT
  5. Cellulose deposition for comparison

The key question this answers:
  "When do the melanin enzymes actually accumulate, and is it truly
   AFTER cellulose deposition is complete?"
"""

import numpy as np
from scipy.integrate import solve_ivp
import yaml
import json
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
CONFIG_DIR = BASE_DIR / "config"
RESULTS_DIR = BASE_DIR / "results"


def load_params() -> dict:
    """Load biological parameters from YAML config."""
    with open(CONFIG_DIR / "parameters.yaml") as f:
        return yaml.safe_load(f)


def promoter_activity(t_dpa: float, activation_dpa: float, peak_dpa: float,
                       hill_n: float, leakage: float, strength: float) -> float:
    """
    Model promoter activity as a function of developmental time.
    
    Uses a Hill function centered at activation_dpa:
        activity = leakage + (strength - leakage) * (t^n / (K^n + t^n))
    
    where t is DPA and K is the half-max point.
    
    Parameters:
        t_dpa: Current days post anthesis
        activation_dpa: DPA when promoter starts firing
        peak_dpa: DPA of peak expression
        hill_n: Steepness of the activation curve
        leakage: Basal expression level (fraction of max)
        strength: Maximum relative expression level
    
    Returns:
        Promoter activity (0 to strength)
    """
    if t_dpa < 0:
        return leakage * strength
    
    # Half-max point between activation and peak
    K = (activation_dpa + peak_dpa) / 2.0
    
    if K <= 0:
        return strength  # Constitutive
    
    activity = leakage * strength + (strength - leakage * strength) * \
               (t_dpa ** hill_n / (K ** hill_n + t_dpa ** hill_n))
    
    return activity


def cellulose_accumulation(t_dpa: float, params: dict) -> float:
    """
    Model cellulose accumulation in cotton fiber.
    
    Uses a sigmoidal (Hill) function based on published CesA expression data:
    - Cellulose deposition begins ~16 DPA
    - Peaks at ~28 DPA
    - Plateaus by ~40 DPA
    """
    cell_params = params['cellulose']
    K = cell_params['half_max_dpa']
    n = cell_params['hill_coefficient']
    C_max = cell_params['total_cellulose_per_fiber']
    
    if t_dpa <= 0:
        return 0.0
    
    return C_max * (t_dpa ** n / (K ** n + t_dpa ** n))


def expression_odes(t, y, params):
    """
    System of ODEs for gene expression dynamics.
    
    State variables (y):
        y[0] = mRNA_melA       (melA transcript)
        y[1] = protein_melA    (Tyrosinase protein)
        y[2] = mRNA_TYRP1
        y[3] = protein_TYRP1
        y[4] = mRNA_DCT
        y[5] = protein_DCT
        y[6] = mRNA_nptII      (selection marker)
        y[7] = protein_nptII
    
    Parameters from config:
        - Promoter kinetics
        - mRNA half-life (~2-6 hours for plant mRNAs)
        - Protein half-life (~12-48 hours)
        - Transcription rate
        - Translation rate
    """
    # Unpack state
    mRNA_melA, prot_melA = y[0], y[1]
    mRNA_TYRP1, prot_TYRP1 = y[2], y[3]
    mRNA_DCT, prot_DCT = y[4], y[5]
    mRNA_nptII, prot_nptII = y[6], y[7]
    
    # Kinetic constants (reasonable estimates for plant systems)
    # mRNA half-life: ~4 hours → degradation rate = ln(2)/4 = 0.173 hr^-1
    # Protein half-life: ~24 hours → degradation rate = ln(2)/24 = 0.029 hr^-1
    k_mrna_deg = 0.173     # hr^-1
    k_prot_deg = 0.029     # hr^-1
    k_transcription = 5.0  # mRNA molecules/hr at full promoter activity
    k_translation = 2.0    # protein molecules/mRNA/hr
    
    # Time in DPA (convert hours to days for promoter lookup)
    t_dpa = t / 24.0  # simulation runs in hours, promoters defined in days
    
    prom = params['promoters']
    
    # ── Promoter activities ───────────────────────────────────────────
    
    # melA: maturation promoter (pGhMat1)
    p_melA = promoter_activity(
        t_dpa,
        prom['pGhMat1']['activation_dpa'],
        prom['pGhMat1']['peak_dpa'],
        prom['pGhMat1']['hill_coefficient'],
        prom['pGhMat1']['leakage_fraction'],
        prom['pGhMat1']['strength_relative']
    )
    
    # TYRP1 & DCT: late SCW promoter
    p_SCW = promoter_activity(
        t_dpa,
        prom['pGhSCW_late']['activation_dpa'],
        prom['pGhSCW_late']['peak_dpa'],
        prom['pGhSCW_late']['hill_coefficient'],
        prom['pGhSCW_late']['leakage_fraction'],
        prom['pGhSCW_late']['strength_relative']
    )
    
    # nptII: constitutive 35S promoter
    p_35S = promoter_activity(
        t_dpa,
        prom['p35S']['activation_dpa'],
        prom['p35S']['peak_dpa'],
        prom['p35S']['hill_coefficient'],
        prom['p35S']['leakage_fraction'],
        prom['p35S']['strength_relative']
    )
    
    # ── ODEs ──────────────────────────────────────────────────────────
    
    # melA
    d_mRNA_melA = k_transcription * p_melA - k_mrna_deg * mRNA_melA
    d_prot_melA = k_translation * mRNA_melA - k_prot_deg * prot_melA
    
    # TYRP1
    d_mRNA_TYRP1 = k_transcription * p_SCW - k_mrna_deg * mRNA_TYRP1
    d_prot_TYRP1 = k_translation * mRNA_TYRP1 - k_prot_deg * prot_TYRP1
    
    # DCT
    d_mRNA_DCT = k_transcription * p_SCW - k_mrna_deg * mRNA_DCT
    d_prot_DCT = k_translation * mRNA_DCT - k_prot_deg * prot_DCT
    
    # nptII
    d_mRNA_nptII = k_transcription * p_35S - k_mrna_deg * mRNA_nptII
    d_prot_nptII = k_translation * mRNA_nptII - k_prot_deg * prot_nptII
    
    return [d_mRNA_melA, d_prot_melA,
            d_mRNA_TYRP1, d_prot_TYRP1,
            d_mRNA_DCT, d_prot_DCT,
            d_mRNA_nptII, d_prot_nptII]


def run_simulation(params: dict) -> dict:
    """
    Run the full gene expression simulation over the cotton fiber lifecycle.
    
    Returns:
        Dictionary with time arrays and all state variables.
    """
    total_hours = params['fiber_development']['total_development_days'] * 24
    t_span = (0, total_hours)
    t_eval = np.linspace(0, total_hours, total_hours + 1)
    
    # Initial conditions (all zero — no expression before anthesis)
    y0 = [0.0] * 8
    
    print("  Running ODE solver...")
    sol = solve_ivp(
        expression_odes,
        t_span,
        y0,
        args=(params,),
        t_eval=t_eval,
        method='RK45',
        rtol=1e-8,
        atol=1e-10
    )
    
    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")
    
    t_dpa = sol.t / 24.0
    
    # Compute cellulose accumulation for comparison
    cellulose = np.array([cellulose_accumulation(t, params) for t in t_dpa])
    
    # Compute promoter activities
    prom = params['promoters']
    prom_melA = np.array([promoter_activity(
        t, prom['pGhMat1']['activation_dpa'], prom['pGhMat1']['peak_dpa'],
        prom['pGhMat1']['hill_coefficient'], prom['pGhMat1']['leakage_fraction'],
        prom['pGhMat1']['strength_relative']
    ) for t in t_dpa])
    
    prom_SCW = np.array([promoter_activity(
        t, prom['pGhSCW_late']['activation_dpa'], prom['pGhSCW_late']['peak_dpa'],
        prom['pGhSCW_late']['hill_coefficient'], prom['pGhSCW_late']['leakage_fraction'],
        prom['pGhSCW_late']['strength_relative']
    ) for t in t_dpa])
    
    results = {
        't_dpa': t_dpa,
        't_hours': sol.t,
        'mRNA_melA': sol.y[0],
        'protein_melA': sol.y[1],
        'mRNA_TYRP1': sol.y[2],
        'protein_TYRP1': sol.y[3],
        'mRNA_DCT': sol.y[4],
        'protein_DCT': sol.y[5],
        'mRNA_nptII': sol.y[6],
        'protein_nptII': sol.y[7],
        'cellulose': cellulose,
        'promoter_melA': prom_melA,
        'promoter_SCW_late': prom_SCW,
    }
    
    return results


def print_summary(results: dict, params: dict):
    """Print key simulation results."""
    t = results['t_dpa']
    
    print("\n" + "=" * 60)
    print("  GENE EXPRESSION SIMULATION RESULTS")
    print("=" * 60)
    
    # Find when each protein reaches 50% of max
    for name, key in [('melA (Tyrosinase)', 'protein_melA'),
                       ('TYRP1', 'protein_TYRP1'),
                       ('DCT', 'protein_DCT')]:
        protein = results[key]
        max_val = protein.max()
        if max_val > 0:
            half_max_idx = np.argmax(protein >= max_val * 0.5)
            half_max_dpa = t[half_max_idx]
            peak_dpa = t[np.argmax(protein)]
            print(f"\n  {name}:")
            print(f"    50% max at:  {half_max_dpa:.1f} DPA")
            print(f"    Peak at:     {peak_dpa:.1f} DPA")
            print(f"    Max level:   {max_val:.1f} (relative units)")
    
    # Cellulose completion
    cellulose = results['cellulose']
    cel_90_idx = np.argmax(cellulose >= 0.9 * cellulose.max())
    print(f"\n  Cellulose Deposition:")
    print(f"    90% complete: {t[cel_90_idx]:.1f} DPA")
    print(f"    Max:          {cellulose.max():.2f} mg/fiber")
    
    # Critical check: Is melanin production truly AFTER cellulose?
    melA_half = t[np.argmax(results['protein_melA'] >= results['protein_melA'].max() * 0.5)]
    print(f"\n  ✓ TEMPORAL SEPARATION CHECK:")
    print(f"    Cellulose 90% done:    {t[cel_90_idx]:.1f} DPA")
    print(f"    melA 50% expression:   {melA_half:.1f} DPA")
    gap = melA_half - t[cel_90_idx]
    if gap > 0:
        print(f"    Gap:                   {gap:.1f} days ✓ SAFE — melanin is after cellulose")
    else:
        print(f"    Gap:                   {gap:.1f} days ⚠ OVERLAP — consider adjusting promoter")
    
    print("\n" + "=" * 60)


def save_results(results: dict):
    """Save simulation results as NumPy arrays and JSON summary."""
    import os
    os.makedirs(RESULTS_DIR, exist_ok=True)
    
    # Save as compressed NumPy archive
    np_path = RESULTS_DIR / "expression_data.npz"
    np.savez_compressed(str(np_path), **{k: v for k, v in results.items()})
    print(f"  📊 Saved: {np_path}")
    
    # Save summary JSON
    summary = {
        'simulation_days': float(results['t_dpa'][-1]),
        'peak_melA_dpa': float(results['t_dpa'][np.argmax(results['protein_melA'])]),
        'peak_TYRP1_dpa': float(results['t_dpa'][np.argmax(results['protein_TYRP1'])]),
        'peak_DCT_dpa': float(results['t_dpa'][np.argmax(results['protein_DCT'])]),
        'cellulose_90pct_dpa': float(results['t_dpa'][
            np.argmax(results['cellulose'] >= 0.9 * results['cellulose'].max())
        ]),
        'max_protein_melA': float(results['protein_melA'].max()),
        'max_protein_TYRP1': float(results['protein_TYRP1'].max()),
        'max_protein_DCT': float(results['protein_DCT'].max()),
    }
    
    json_path = RESULTS_DIR / "expression_summary.json"
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"  📊 Saved: {json_path}")


# ── Entry Point ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("\n🧬 BlackCotton Gene Expression Simulator")
    print("=" * 50)
    
    params = load_params()
    print(f"  Simulation: 0 → {params['fiber_development']['total_development_days']} DPA")
    print(f"  Resolution: {params['simulation']['time_resolution_hours']} hour timesteps")
    
    results = run_simulation(params)
    print_summary(results, params)
    
    print("\nSaving results...")
    save_results(results)
    
    print("\n✅ Expression simulation complete!")
