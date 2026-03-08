#!/usr/bin/env python3
import numpy as np
from src.config_loader import load_config
from src.expression_model import run_simulation
from src.melanin_pathway import run_melanin_simulation
import json
import copy
import sys

def calculate_silencing_probability(params, expression_results, activation_dpa=38.0):
    """
    Rule 4 Engine: The Antivirus Problem
    Calculates the RNA Immunogenicity Score and the resulting probability of gene silencing (RNAi).
    """
    t_dpa = expression_results["t_dpa"]
    mRNA_melA = expression_results["mRNA_melA"]
    mRNA_TYRP1 = expression_results["mRNA_TYRP1"]
    mRNA_DCT = expression_results["mRNA_DCT"]
    
    # 1. Codon Optimization Stealth Factor
    codon_config = params.get("codon_optimization", {})
    if codon_config.get("target_organism") == "Gossypium hirsutum" and 0.40 <= codon_config.get("target_gc_content", 0) <= 0.45:
        codon_stealth_factor = 0.85  # Highly optimized, looks like native cotton mRNA
    else:
        codon_stealth_factor = 0.20  # Looks like viral/bacterial garbage, highly immunogenic
        
    probabilities = []
    max_scores = []
    
    # RNAi machinery shuts down late in fiber development
    # At DPA 10, RNAi is 1.0 (maximum). At DPA 45, RNAi is ~0.1.
    rnai_activity = 1.0 - 1.0 / (1.0 + np.exp(-0.5 * (t_dpa - 32.0)))
    
    # Total foreign mRNA load over time
    total_foreign_mrna = mRNA_melA + mRNA_TYRP1 + mRNA_DCT
    
    # Immunogenicity Score: Volume of foreign mRNA * Plant immune state * Codon visibility
    immunogenicity_scores = (total_foreign_mrna / 10.0) * rnai_activity * (1.0 - codon_stealth_factor)
    
    peak_immunogenicity = np.max(immunogenicity_scores)
    
    # Sigmoidal response to immunogenicity peak
    # If the score goes above an arbitrary threshold (e.g. 50), silencing triggers hard.
    silencing_trigger_threshold = 30.0
    prob_silencing = 1.0 / (1.0 + np.exp(-0.5 * (peak_immunogenicity - silencing_trigger_threshold)))
    
    # Bound the probability between base risk and 1.0
    base_risk = params.get("failure_risks", {}).get("silencing_probability", 0.0)
    prob_silencing = max(base_risk, prob_silencing)
    
    return {
        "peak_immunogenicity_score": float(peak_immunogenicity),
        "silencing_probability": float(prob_silencing),
        "codon_stealth_factor": codon_stealth_factor,
        "peak_total_foreign_mrna": float(np.max(total_foreign_mrna))
    }

def print_audit_report():
    print("🛡️  Rule 4 Audit: Antiviral Gene Silencing (RNAi)")
    print("==================================================")
    
    params = load_config()
    expr = run_simulation(params)
    
    # Baseline with current optimized construct (Late Promoter @ 38 DPA)
    res = calculate_silencing_probability(params, expr)
    
    print("\n[Current Construct: melC1 Late Detonator (Day 38)]")
    print(f"  Codon Stealth Rating : {res['codon_stealth_factor']*100:.1f}% (Gossypium hirsutum optimal)")
    print(f"  Peak Foreign mRNA    : {res['peak_total_foreign_mrna']:.1f} transcripts/hr")
    print(f"  Immunogenicity Score : {res['peak_immunogenicity_score']:.2f} (Threshold to trigger: 30.0)")
    print(f"  Silencing Risk       : {res['silencing_probability']*100:.1f}%\n")
    
    if res['silencing_probability'] < 0.10: # Under 10%
        print("  ✓ PASSED: The plant's immune system will ignore the code.")
    else:
        print("  ❌ FAILED: The code is too loud. The plant will delete the instructions.")
        
    print("\n--------------------------------------------------")
    print("🧪 Silencing Sensitivity Sweep (What if we ran the code earlier?)")
    print("--------------------------------------------------")
    
    # Sweep through hypothetical activation timings
    test_timings = [10.0, 15.0, 20.0, 25.0, 30.0, 38.0]
    
    for timing in test_timings:
        # Simulate an earlier promoter
        test_params = copy.deepcopy(params)
        test_params["promoters"]["pGhSCW_late"]["activation_dpa"] = timing
        test_params["promoters"]["pGhSCW_late"]["peak_dpa"] = timing + 6.0
        
        # We need to run the heavy ODE solver silently to get the new mRNA peak
        test_expr = run_simulation(test_params)
        test_res = calculate_silencing_probability(test_params, test_expr, activation_dpa=timing)
        
        status = "💀 SILENCED" if test_res['silencing_probability'] > 0.50 else "✅ SAFE"
        print(f"  Trigger @ DPA {timing:04.1f} | Score: {test_res['peak_immunogenicity_score']:05.1f} | Risk: {test_res['silencing_probability']*100:05.1f}% | {status}")


if __name__ == "__main__":
    print_audit_report()
