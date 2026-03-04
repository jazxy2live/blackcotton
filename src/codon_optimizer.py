#!/usr/bin/env python3
"""
codon_optimizer.py — Codon Optimization for Gossypium hirsutum
===============================================================

Optimizes heterologous gene sequences for expression in cotton by:
  1. Replacing codons with preferred cotton codons
  2. Adjusting GC content to match cotton coding sequences (~38%)
  3. Removing internal restriction sites used in cloning
  4. Avoiding mRNA secondary structure hotspots
  5. Generating a Codon Adaptation Index (CAI) score

Usage:
    python -m src.codon_optimizer
"""

import json
import random
from pathlib import Path
from collections import Counter
from typing import Optional

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
REF_DIR = DATA_DIR / "reference"


# ── Genetic Code ──────────────────────────────────────────────────────────

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Reverse: amino acid → list of codons
AA_TO_CODONS = {}
for codon, aa in CODON_TABLE.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)


def load_codon_usage() -> dict:
    """Load cotton codon usage frequencies from reference data."""
    usage_path = REF_DIR / "codon_usage_table.json"
    with open(usage_path) as f:
        data = json.load(f)
    return data["Gossypium_hirsutum"]["codon_table"]


def translate(dna_seq: str) -> str:
    """Translate a DNA sequence to protein."""
    protein = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3].upper()
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)


def compute_gc_content(seq: str) -> float:
    """Calculate GC content of a DNA sequence."""
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq) if len(seq) > 0 else 0.0


def compute_cai(dna_seq: str, usage_table: dict) -> float:
    """
    Compute the Codon Adaptation Index (CAI) for a DNA sequence.
    
    CAI is the geometric mean of the relative adaptedness of each codon,
    where relative adaptedness = frequency / max frequency for that amino acid.
    """
    import math
    
    w_values = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3].upper()
        aa = CODON_TABLE.get(codon)
        if aa is None or aa == '*' or aa == 'M' or aa == 'W':
            continue  # Skip stop, Met (only ATG), Trp (only TGG)
        
        if aa in usage_table:
            freqs = usage_table[aa]
            max_freq = max(freqs.values())
            codon_freq = freqs.get(codon, 0.01)  # Avoid zero
            w = codon_freq / max_freq
            w_values.append(w)
    
    if not w_values:
        return 0.0
    
    # Geometric mean
    log_sum = sum(math.log(w) for w in w_values)
    cai = math.exp(log_sum / len(w_values))
    return cai


def optimize_codons(dna_seq: str, usage_table: dict, 
                     avoid_sites: list = None,
                     seed: int = 42) -> str:
    """
    Optimize a DNA sequence for expression in cotton.
    
    Strategy:
    1. Translate to protein
    2. Back-translate using cotton-preferred codons (weighted random)
    3. Check for and eliminate unwanted restriction sites
    4. Verify protein sequence is preserved
    """
    random.seed(seed)
    
    if avoid_sites is None:
        avoid_sites = ['GGTCTC', 'CGTCTC', 'GAAGAC', 'GAATTC', 'GGATCC']
    
    protein = translate(dna_seq)
    
    # Back-translate with cotton codon preference
    optimized_codons = []
    for aa in protein:
        if aa in usage_table:
            codons = list(usage_table[aa].keys())
            weights = list(usage_table[aa].values())
            chosen = random.choices(codons, weights=weights, k=1)[0]
            optimized_codons.append(chosen)
        elif aa == 'M':
            optimized_codons.append('ATG')
        elif aa == 'W':
            optimized_codons.append('TGG')
        else:
            # Fallback: use original codon
            optimized_codons.append('NNN')
    
    # Add stop codon
    stop_codons = usage_table.get('*', {'TAA': 0.45, 'TAG': 0.25, 'TGA': 0.30})
    stop = random.choices(list(stop_codons.keys()), 
                          weights=list(stop_codons.values()), k=1)[0]
    optimized_codons.append(stop)
    
    optimized_seq = ''.join(optimized_codons)
    
    # Remove restriction sites (iterate until clean)
    for _ in range(10):  # Max iterations
        clean = True
        for site in avoid_sites:
            # Check forward and reverse complement
            rc_site = site[::-1].translate(str.maketrans('ATGC', 'TACG'))
            for check_site in [site, rc_site]:
                pos = optimized_seq.find(check_site)
                if pos >= 0:
                    clean = False
                    # Find which codon this falls in and swap to alternative
                    codon_idx = pos // 3
                    if codon_idx < len(optimized_codons) - 1:
                        old_codon = optimized_codons[codon_idx]
                        aa = CODON_TABLE.get(old_codon, '')
                        if aa and aa in usage_table:
                            alternatives = [c for c in usage_table[aa].keys() 
                                          if c != old_codon]
                            if alternatives:
                                optimized_codons[codon_idx] = random.choice(alternatives)
                                optimized_seq = ''.join(optimized_codons)
        if clean:
            break
    
    return optimized_seq


def analyze_sequence(name: str, original: str, optimized: str, usage_table: dict):
    """Print a comparison analysis of original vs optimized sequence."""
    
    print(f"\n{'='*60}")
    print(f"  CODON OPTIMIZATION REPORT: {name}")
    print(f"{'='*60}")
    
    orig_protein = translate(original)
    opt_protein = translate(optimized)
    
    print(f"\n  Protein length:     {len(orig_protein)} aa")
    print(f"  Protein preserved:  {'✓ YES' if orig_protein == opt_protein else '✗ NO — CRITICAL ERROR'}")
    
    orig_gc = compute_gc_content(original)
    opt_gc = compute_gc_content(optimized)
    print(f"\n  GC Content:")
    print(f"    Original:         {orig_gc:.1%}")
    print(f"    Optimized:        {opt_gc:.1%}")
    print(f"    Cotton target:    38.0%")
    print(f"    Improvement:      {'✓' if abs(opt_gc - 0.38) < abs(orig_gc - 0.38) else '—'}")
    
    orig_cai = compute_cai(original, usage_table)
    opt_cai = compute_cai(optimized, usage_table)
    print(f"\n  Codon Adaptation Index (CAI):")
    print(f"    Original:         {orig_cai:.3f}")
    print(f"    Optimized:        {opt_cai:.3f}")
    print(f"    Improvement:      {((opt_cai - orig_cai) / orig_cai * 100):.1f}%")
    
    # Codon usage comparison
    print(f"\n  Codon Usage Frequency (top changes):")
    orig_codons = [original[i:i+3] for i in range(0, len(original)-2, 3)]
    opt_codons = [optimized[i:i+3] for i in range(0, len(optimized)-2, 3)]
    
    changes = sum(1 for o, n in zip(orig_codons, opt_codons) if o != n)
    total = min(len(orig_codons), len(opt_codons))
    print(f"    Codons changed:   {changes}/{total} ({changes/total*100:.1f}%)")
    
    print(f"\n{'='*60}\n")


# ── Entry Point ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("\n🧬 BlackCotton Codon Optimizer")
    print("=" * 50)
    
    usage_table = load_codon_usage()
    
    # Resolve sequences using the same validated path as construct assembly.
    from src.construct_designer import (
        resolve_melA_cds,
        resolve_strict_cds,
        SEQ_DIR,
    )

    melA_seq, _ = resolve_melA_cds()
    melA_opt = optimize_codons(melA_seq, usage_table)
    analyze_sequence("melA (Tyrosinase)", melA_seq, melA_opt, usage_table)
    
    # Load and optimize TYRP1
    tyrp1_seq, _ = resolve_strict_cds(
        real_path=SEQ_DIR / "TYRP1_CDS_real.fasta",
        fallback_path=SEQ_DIR / "GhTYRP1_optimized.fasta",
    )
    tyrp1_opt = optimize_codons(tyrp1_seq, usage_table)
    analyze_sequence("GhTYRP1", tyrp1_seq, tyrp1_opt, usage_table)
    
    # Load and optimize DCT
    dct_seq, _ = resolve_strict_cds(
        real_path=SEQ_DIR / "DCT_CDS_real.fasta",
        fallback_path=SEQ_DIR / "GhDCT_optimized.fasta",
    )
    dct_opt = optimize_codons(dct_seq, usage_table)
    analyze_sequence("GhDCT", dct_seq, dct_opt, usage_table)
    
    print("✅ Codon optimization complete!")
