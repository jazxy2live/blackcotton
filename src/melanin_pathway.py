#!/usr/bin/env python3
"""
melanin_pathway.py — Melanin Biosynthesis Kinetics Model
=========================================================

Models the enzymatic melanin biosynthesis pathway in cotton fiber lumen:

    L-Tyrosine  ──[melA]──▸  L-DOPA  ──[melA]──▸  Dopaquinone
                                                        │
                                      ┌─────────────────┘
                                      ▼
                                 Leucodopachrome
                                      │
                              ──[DCT]─┤
                                      ▼
                                   DHICA
                                      │
                             ──[TYRP1]─┤
                                      ▼
                              Indole-quinone
                                      │
                             (polymerization)
                                      ▼
                                   MELANIN

Uses Michaelis-Menten kinetics for each enzymatic step and
first-order kinetics for spontaneous reactions.
"""

import numpy as np
from scipy.integrate import solve_ivp
import json
from pathlib import Path

from src.config_loader import load_config
from src.failure_risk_model import resolved_failure_risks

BASE_DIR = Path(__file__).resolve().parent.parent
CONFIG_DIR = BASE_DIR / "config"
RESULTS_DIR = BASE_DIR / "results"


def load_params() -> dict:
    return load_config()


def michaelis_menten(substrate: float, Vmax: float, Km: float) -> float:
    """Classical Michaelis-Menten rate equation."""
    if substrate <= 0:
        return 0.0
    return Vmax * substrate / (Km + substrate)


def melanin_odes(t, y, enzyme_levels, params):
    """
    System of ODEs for the melanin biosynthesis pathway.
    
    State variables:
        y[0] = [L-Tyrosine]        (mM)
        y[1] = [L-DOPA]            (mM)
        y[2] = [Dopaquinone]       (mM)
        y[3] = [Leucodopachrome]   (mM) — cyclization product
        y[4] = [Dopachrome]        (mM)
        y[5] = [DHICA]             (mM) — DCT product
        y[6] = [Indole-quinone]    (mM) — TYRP1 product
        y[7] = [Melanin]           (mM equivalent, polymer)
    
    enzyme_levels: callable(t) → (melA_level, TYRP1_level, DCT_level)
        Returns relative enzyme concentrations at time t (from expression model)
    """
    Tyr, DOPA, DQ, LDC, DC, DHICA, IQ, Melanin = y
    
    mp = params['melanin_pathway']
    melA_params = mp['tyrosinase_melA']
    tyrp1_params = mp['TYRP1']
    dct_params = mp['DCT']
    
    # Get enzyme levels at current time
    melA_e, tyrp1_e, dct_e = enzyme_levels(t)

    # Fail-first activity scaling: cofactor loading + activation state +
    # transgene stability affect catalytically active enzyme availability.
    risk = resolved_failure_risks(params)
    copper_loading = float(risk["copper_loading_fraction"])
    tyrosinase_activation = float(risk["tyrosinase_activation_fraction"])
    silencing_prob = float(risk["silencing_probability"])
    event_cv = float(risk["event_expression_cv"])
    ros_capacity = float(risk["ros_buffer_capacity"])

    tc = mp.get("traffic_control", {})
    melA_entry_scale = float(np.clip(tc.get("melA_entry_scale", 1.0), 0.05, 2.0))
    dopaquinone_drain_scale = float(np.clip(tc.get("dopaquinone_drain_scale", 1.0), 0.20, 3.0))
    dopachrome_drain_scale = float(np.clip(tc.get("dopachrome_drain_scale", 1.0), 0.20, 3.0))
    dct_flux_scale = float(np.clip(tc.get("dct_flux_scale", 1.0), 0.20, 3.0))
    tyrp1_flux_scale = float(np.clip(tc.get("tyrp1_flux_scale", 1.0), 0.20, 3.0))
    eumelanin_push_scale = float(np.clip(tc.get("eumelanin_push_scale", 1.0), 0.20, 3.0))

    expression_scale = max(0.05, (1.0 - silencing_prob) * np.exp(-0.5 * (event_cv ** 2)))
    
    # ── THE DETONATOR: melC1 Chaperone Switch ──
    # Instead of a static activation fraction, melA is locked in an inactive state.
    # It requires the melC1 chaperone to fold and load copper. We model melC1 
    # as being blasted late by the maturation promoter (e.g. peaking around day 45).
    t_dpa = t / 24.0
    melC1_level = 0.0
    if t_dpa > 38.0:
        # Logistic release of the detonator
        melC1_level = 1.0 / (1.0 + np.exp(-2.0 * (t_dpa - 40.0)))
    
    # True active melA is now conditionally dependent on the detonator
    true_activation_fraction = tyrosinase_activation * melC1_level
    melA_active = melA_e * copper_loading * true_activation_fraction * expression_scale * melA_entry_scale
    
    tyrp1_active = tyrp1_e * expression_scale
    dct_active = dct_e * expression_scale
    detox_scale = float(np.clip(0.80 + 0.30 * ros_capacity, 0.60, 1.25))
    
    # ── Enzymatic rates ───────────────────────────────────────────────
    
    # Step 1: L-Tyrosine → L-DOPA (melA hydroxylase activity)
    v1 = melA_active * michaelis_menten(Tyr, melA_params['Vmax'], melA_params['Km_tyrosine'])
    
    # Step 2: L-DOPA → Dopaquinone (melA oxidase activity)
    v2 = melA_active * michaelis_menten(DOPA, melA_params['Vmax'] * 0.8, melA_params['Km_tyrosine'] * 0.5)
    
    # Step 3: Dopaquinone → Leucodopachrome (spontaneous cyclization)
    v3 = mp['dopaquinone_cyclization_rate'] * dopaquinone_drain_scale * DQ
    
    # Step 4: Leucodopachrome → Dopachrome (spontaneous oxidation)
    v4 = 0.3 * dopachrome_drain_scale * LDC  # Spontaneous oxidation + downstream pull

    # Step 5: Dopachrome → DHICA (DCT-catalyzed)
    v5 = dct_active * dct_flux_scale * michaelis_menten(DC, dct_params['Vmax'], dct_params['Km_dopachrome'])

    # Step 6: DHICA → Indole-5,6-quinone (TYRP1-catalyzed oxidation)
    v6 = tyrp1_active * tyrp1_flux_scale * michaelis_menten(DHICA, tyrp1_params['Vmax'], tyrp1_params['Km_DHICA'])

    # Step 7: Indole-quinone → Melanin (polymerization)
    v7 = mp['indole_polymerization_rate'] * detox_scale * eumelanin_push_scale * IQ
    
    # ── THE FUEL PUMP: Mutant TyrA (Feedback-Insensitive) ──
    # Instead of a static 0.5mM pool, we use a feedback-insensitive prephenate 
    # dehydrogenase (TyrA) to constantly convert shikimate pathway carbon into 
    # fresh L-Tyrosine. This pump is driven by the late promoter (proxied by tyrp1_e).
    pump_rate = float(mp.get('tyra_mutant_pump_rate', 0.0))
    tyra_pump = pump_rate * tyrp1_e * expression_scale

    # ── ODEs ──────────────────────────────────────────────────────────
    
    dTyr = tyra_pump - v1               # Synthesized/Imported by TyrA pump, consumed by melA
    dDOPA = v1 - v2                     # Produced by step 1, consumed by step 2
    dDQ = v2 - v3                       # Produced by step 2, consumed by cyclization
    dLDC = v3 - v4                      # Produced by cyclization, consumed by oxidation
    dDC = v4 - v5                       # Produced by oxidation, consumed by DCT
    dDHICA = v5 - v6                    # Produced by DCT, consumed by TYRP1
    dIQ = v6 - v7                       # Produced by TYRP1, consumed by polymerization
    dMelanin = v7                       # Accumulates (polymer, irreversible)
    
    return [dTyr, dDOPA, dDQ, dLDC, dDC, dDHICA, dIQ, dMelanin]


def create_enzyme_interpolator(expression_results: dict = None):
    """
    Create a function that returns enzyme levels at any time point.
    
    If expression_results are provided (from expression_model.py),
    use those. Otherwise, use a simplified model based on promoter timing.
    """
    if expression_results is not None:
        t_hours = expression_results['t_hours']
        
        # Normalize to 0-1 range for use as enzyme scaling factors
        melA_max = expression_results['protein_melA'].max()
        tyrp1_max = expression_results['protein_TYRP1'].max()
        dct_max = expression_results['protein_DCT'].max()
        
        melA_norm = expression_results['protein_melA'] / melA_max if melA_max > 0 else expression_results['protein_melA']
        tyrp1_norm = expression_results['protein_TYRP1'] / tyrp1_max if tyrp1_max > 0 else expression_results['protein_TYRP1']
        dct_norm = expression_results['protein_DCT'] / dct_max if dct_max > 0 else expression_results['protein_DCT']
        
        # SciPy interp1d is preferred, but NumPy fallback keeps the pipeline
        # usable in environments with partial SciPy builds.
        try:
            from scipy.interpolate import interp1d

            f_melA = interp1d(t_hours, melA_norm, fill_value=(0, melA_norm[-1]), bounds_error=False)
            f_tyrp1 = interp1d(t_hours, tyrp1_norm, fill_value=(0, tyrp1_norm[-1]), bounds_error=False)
            f_dct = interp1d(t_hours, dct_norm, fill_value=(0, dct_norm[-1]), bounds_error=False)

            def enzyme_levels(t):
                return float(f_melA(t)), float(f_tyrp1(t)), float(f_dct(t))
        except Exception:
            x_min = float(t_hours[0]) if len(t_hours) else 0.0
            x_max = float(t_hours[-1]) if len(t_hours) else 0.0

            def _interp(arr, t):
                if t <= x_min:
                    return 0.0
                if t >= x_max:
                    return float(arr[-1]) if len(arr) else 0.0
                return float(np.interp(t, t_hours, arr))

            def enzyme_levels(t):
                return _interp(melA_norm, t), _interp(tyrp1_norm, t), _interp(dct_norm, t)

        return enzyme_levels
    else:
        # Simplified model: step function activation
        def enzyme_levels(t):
            t_dpa = t / 24.0
            melA = 1.0 if t_dpa >= 35.0 else 0.01
            tyrp1 = 1.0 if t_dpa >= 30.0 else 0.02
            dct = 1.0 if t_dpa >= 30.0 else 0.02
            return melA, tyrp1, dct
        
        return enzyme_levels


def run_melanin_simulation(params: dict, expression_results: dict = None) -> dict:
    """
    Run the melanin biosynthesis simulation.
    
    Parameters:
        params: Configuration parameters
        expression_results: Optional results from expression_model.py
    
    Returns:
        Dictionary with time arrays and metabolite concentrations.
    """
    mp = params['melanin_pathway']
    total_hours = params['fiber_development']['total_development_days'] * 24
    
    # Initial conditions
    y0 = [
        mp['initial_tyrosine_mM'],  # L-Tyrosine
        0.0,                         # L-DOPA
        0.0,                         # Dopaquinone
        0.0,                         # Leucodopachrome
        0.0,                         # Dopachrome
        0.0,                         # DHICA
        0.0,                         # Indole-quinone
        0.0,                         # Melanin
    ]
    
    t_span = (0, total_hours)
    t_eval = np.linspace(0, total_hours, 2000)
    
    enzyme_func = create_enzyme_interpolator(expression_results)
    
    print("  Running melanin pathway simulation...")
    sol = solve_ivp(
        melanin_odes,
        t_span,
        y0,
        args=(enzyme_func, params),
        t_eval=t_eval,
        method='RK45',
        rtol=1e-8,
        atol=1e-10
    )
    
    if not sol.success:
        raise RuntimeError(f"Melanin ODE solver failed: {sol.message}")
    
    results = {
        't_hours': sol.t,
        't_dpa': sol.t / 24.0,
        'tyrosine': sol.y[0],
        'L_DOPA': sol.y[1],
        'dopaquinone': sol.y[2],
        'leucodopachrome': sol.y[3],
        'dopachrome': sol.y[4],
        'DHICA': sol.y[5],
        'indole_quinone': sol.y[6],
        'melanin': sol.y[7],
    }
    
    return results


def print_melanin_summary(results: dict, params: dict):
    """Print pathway simulation results."""
    t = results['t_dpa']
    melanin = results['melanin']
    
    print("\n" + "=" * 60)
    print("  MELANIN BIOSYNTHESIS SIMULATION RESULTS")
    print("=" * 60)
    
    # Final melanin concentration
    final_melanin = melanin[-1]
    target = params['melanin_pathway']['melanin_target_mg_per_fiber']
    
    print(f"\n  Final melanin:      {final_melanin:.4f} mM equivalent")
    print(f"  Target:             {target} mg/fiber")
    
    # When does melanin accumulation begin?
    threshold = 0.01 * final_melanin  # 1% of max
    onset_idx = np.argmax(melanin >= threshold) if final_melanin > 0 else -1
    if onset_idx > 0:
        print(f"  Onset (1% max):     {t[onset_idx]:.1f} DPA")
    
    # When does it reach 50%?
    half_idx = np.argmax(melanin >= 0.5 * final_melanin) if final_melanin > 0 else -1
    if half_idx > 0:
        print(f"  50% accumulation:   {t[half_idx]:.1f} DPA")
    
    # Substrate depletion
    tyr_remaining = results['tyrosine'][-1] / results['tyrosine'][0] * 100
    print(f"\n  Tyrosine remaining: {tyr_remaining:.1f}%")
    print(f"  Tyrosine consumed:  {100 - tyr_remaining:.1f}%")
    
    # Intermediate accumulation (check for toxic buildup)
    print(f"\n  Intermediate Levels (final, should be low):")
    for name, key in [('L-DOPA', 'L_DOPA'), ('Dopaquinone', 'dopaquinone'),
                       ('Dopachrome', 'dopachrome'), ('DHICA', 'DHICA')]:
        val = results[key][-1]
        status = "✓ LOW" if val < 0.01 else "⚠ ELEVATED"
        print(f"    {name:<18}: {val:.6f} mM  {status}")
    
    print("\n" + "=" * 60)


def save_melanin_results(results: dict):
    """Save melanin simulation results."""
    import os
    os.makedirs(RESULTS_DIR, exist_ok=True)
    
    np_path = RESULTS_DIR / "melanin_data.npz"
    np.savez_compressed(str(np_path), **results)
    print(f"  📊 Saved: {np_path}")
    
    summary = {
        'final_melanin_mM': float(results['melanin'][-1]),
        'tyrosine_consumed_pct': float(
            (1 - results['tyrosine'][-1] / results['tyrosine'][0]) * 100
        ),
        'simulation_days': float(results['t_dpa'][-1]),
    }
    
    json_path = RESULTS_DIR / "melanin_summary.json"
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"  📊 Saved: {json_path}")


# ── Entry Point ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("\n🧬 BlackCotton Melanin Pathway Simulator")
    print("=" * 50)
    
    params = load_params()
    
    # Try to load expression model results
    expr_data_path = RESULTS_DIR / "expression_data.npz"
    expression_results = None
    if expr_data_path.exists():
        print("  Loading expression data from previous simulation...")
        data = np.load(str(expr_data_path))
        expression_results = {k: data[k] for k in data.files}
        print("  ✓ Loaded expression profiles")
    else:
        print("  ⚠ No expression data found — using simplified enzyme model")
        print("    (Run expression_model.py first for full simulation)")
    
    results = run_melanin_simulation(params, expression_results)
    print_melanin_summary(results, params)
    
    print("\nSaving results...")
    save_melanin_results(results)
    
    print("\n✅ Melanin pathway simulation complete!")
