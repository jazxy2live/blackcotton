#!/usr/bin/env python3
"""
construct_designer.py — T-DNA Genetic Construct Assembly
=========================================================

Assembles the complete T-DNA construct for melanin production in cotton fibers.
Reads gene/promoter/terminator sequences and outputs:
  - Annotated construct map
  - Complete construct sequence (FASTA)
  - Machine-readable construct summary (JSON)

The construct architecture:
  LB — [pGhMat1 → melA → tNOS] — [pGhSCW-late → TYRP1 → tNOS] — 
       [pGhSCW-late → DCT → tNOS] — [p35S → nptII → tNOS] — RB
"""

import json
import os
import sys
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional

from src.config_loader import load_config

# ── Configuration ──────────────────────────────────────────────────────────

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
SEQ_DIR = DATA_DIR / "sequences"
RESULTS_DIR = BASE_DIR / "results"

# T-DNA border sequences (25 bp conserved repeats)
LEFT_BORDER = "TGGCAGGATATATTGTGGTGTAAAC"   # LB
RIGHT_BORDER = "TGGCAGGATATATTGTGGTGTAAAC"  # RB

STOP_CODONS = {"TAA", "TAG", "TGA"}
STRICT_START_CODONS = {"ATG"}
BACTERIAL_START_CODONS = {"ATG", "GTG", "TTG"}


@dataclass
class GeneticElement:
    """A single genetic element (promoter, gene, or terminator)."""
    name: str
    element_type: str  # 'promoter', 'gene', 'terminator', 'border'
    sequence: str
    description: str = ""
    length: int = 0
    
    def __post_init__(self):
        self.sequence = self.sequence.upper().replace('\n', '').replace(' ', '')
        self.length = len(self.sequence)


@dataclass
class Cassette:
    """A transcription unit: promoter → gene → terminator."""
    name: str
    promoter: GeneticElement
    gene: GeneticElement
    terminator: GeneticElement
    description: str = ""
    
    @property
    def sequence(self) -> str:
        return self.promoter.sequence + self.gene.sequence + self.terminator.sequence
    
    @property
    def length(self) -> int:
        return len(self.sequence)
    
    def get_feature_map(self, offset: int = 0) -> list:
        """Return annotation features with positions relative to construct."""
        features = []
        pos = offset
        
        features.append({
            'name': self.promoter.name,
            'type': 'promoter',
            'start': pos,
            'end': pos + self.promoter.length,
            'strand': '+'
        })
        pos += self.promoter.length
        
        features.append({
            'name': self.gene.name,
            'type': 'CDS',
            'start': pos,
            'end': pos + self.gene.length,
            'strand': '+'
        })
        pos += self.gene.length
        
        features.append({
            'name': self.terminator.name,
            'type': 'terminator',
            'start': pos,
            'end': pos + self.terminator.length,
            'strand': '+'
        })
        
        return features


@dataclass
class TDNAConstruct:
    """Complete T-DNA construct with all cassettes."""
    name: str
    cassettes: list = field(default_factory=list)
    left_border: str = LEFT_BORDER
    right_border: str = RIGHT_BORDER
    
    @property
    def sequence(self) -> str:
        """Full construct sequence from LB to RB."""
        parts = [self.left_border]
        for cassette in self.cassettes:
            parts.append(cassette.sequence)
        parts.append(self.right_border)
        return ''.join(parts)
    
    @property
    def total_length(self) -> int:
        return len(self.sequence)
    
    def get_all_features(self) -> list:
        """Get complete feature annotation map."""
        features = []
        pos = 0
        
        # Left border
        features.append({
            'name': 'LB',
            'type': 'T-DNA_border',
            'start': pos,
            'end': pos + len(self.left_border),
            'strand': '+'
        })
        pos += len(self.left_border)
        
        # Cassettes
        for cassette in self.cassettes:
            cassette_features = cassette.get_feature_map(offset=pos)
            features.extend(cassette_features)
            pos += cassette.length
        
        # Right border
        features.append({
            'name': 'RB',
            'type': 'T-DNA_border',
            'start': pos,
            'end': pos + len(self.right_border),
            'strand': '+'
        })
        
        return features
    
    def get_construct_map(self) -> str:
        """Generate a visual text map of the construct."""
        lines = []
        lines.append("=" * 80)
        lines.append(f"  CONSTRUCT: {self.name}")
        lines.append(f"  Total Length: {self.total_length:,} bp")
        lines.append("=" * 80)
        lines.append("")
        
        # Visual diagram
        diagram_parts = ["LB"]
        for cassette in self.cassettes:
            diagram_parts.append(
                f"[{cassette.promoter.name} → {cassette.gene.name} → {cassette.terminator.name}]"
            )
        diagram_parts.append("RB")
        lines.append("  " + " ── ".join(diagram_parts))
        lines.append("")
        
        # Detailed feature table
        lines.append("  FEATURE TABLE")
        lines.append("  " + "-" * 76)
        lines.append(f"  {'Feature':<25} {'Type':<15} {'Start':>8} {'End':>8} {'Length':>8}")
        lines.append("  " + "-" * 76)
        
        for feat in self.get_all_features():
            length = feat['end'] - feat['start']
            lines.append(
                f"  {feat['name']:<25} {feat['type']:<15} {feat['start']:>8} {feat['end']:>8} {length:>8}"
            )
        
        lines.append("  " + "-" * 76)
        lines.append("")
        
        # Cassette summaries
        lines.append("  CASSETTE DETAILS")
        lines.append("  " + "-" * 76)
        for i, cassette in enumerate(self.cassettes, 1):
            lines.append(f"  Cassette {i}: {cassette.name}")
            lines.append(f"    Promoter:   {cassette.promoter.name} ({cassette.promoter.length} bp)")
            lines.append(f"    Gene:       {cassette.gene.name} ({cassette.gene.length} bp)")
            lines.append(f"    Terminator: {cassette.terminator.name} ({cassette.terminator.length} bp)")
            lines.append(f"    Total:      {cassette.length} bp")
            if cassette.description:
                lines.append(f"    Function:   {cassette.description}")
            lines.append("")
        
        # GC content analysis
        seq = self.sequence
        gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
        at = (seq.count('A') + seq.count('T')) / len(seq) * 100
        lines.append(f"  SEQUENCE COMPOSITION")
        lines.append(f"  " + "-" * 76)
        lines.append(f"  GC Content: {gc:.1f}%")
        lines.append(f"  AT Content: {at:.1f}%")
        lines.append(f"  A: {seq.count('A'):,}  T: {seq.count('T'):,}  "
                     f"G: {seq.count('G'):,}  C: {seq.count('C'):,}")
        lines.append("")
        
        # Restriction site scan
        lines.append("  RESTRICTION SITE SCAN (sites to avoid)")
        lines.append("  " + "-" * 76)
        restriction_sites = {
            'EcoRI': 'GAATTC',
            'BamHI': 'GGATCC',
            'HindIII': 'AAGCTT',
            'BsaI': 'GGTCTC',
            'BsmBI': 'CGTCTC',
            'NotI': 'GCGGCCGC',
            'XbaI': 'TCTAGA',
            'SacI': 'GAGCTC'
        }
        for name, site in restriction_sites.items():
            count = seq.count(site) + seq.count(
                site[::-1].translate(str.maketrans('ATGC', 'TACG'))  # reverse complement
            )
            status = "✓ CLEAR" if count == 0 else f"⚠ FOUND ({count}x)"
            lines.append(f"  {name:<10} ({site}): {status}")
        
        lines.append("")
        lines.append("=" * 80)
        
        return '\n'.join(lines)
    
    def to_json(self) -> dict:
        """Export construct data as JSON-serializable dict."""
        return {
            'name': self.name,
            'total_length_bp': self.total_length,
            'gc_content': (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence),
            'num_cassettes': len(self.cassettes),
            'cassettes': [
                {
                    'name': c.name,
                    'length_bp': c.length,
                    'promoter': {'name': c.promoter.name, 'length': c.promoter.length},
                    'gene': {
                        'name': c.gene.name,
                        'length': c.gene.length,
                        'transit_peptide_fused': bool("_vTP" in c.gene.name),
                        'description': c.gene.description,
                    },
                    'terminator': {'name': c.terminator.name, 'length': c.terminator.length},
                    'description': c.description
                }
                for c in self.cassettes
            ],
            'features': self.get_all_features()
        }


# ── Sequence I/O ──────────────────────────────────────────────────────────

def read_fasta(filepath: str) -> str:
    """Read a FASTA file and return the sequence (no header)."""
    sequence_lines = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            sequence_lines.append(line)
    return ''.join(sequence_lines).upper().replace(' ', '')


def load_params() -> dict:
    """Load active project parameters for construct-level options."""
    return load_config()


def read_genbank_sequence(filepath: str) -> str:
    """Read a GenBank file and extract the ORIGIN sequence."""
    in_origin = False
    sequence_parts = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('ORIGIN'):
                in_origin = True
                continue
            if line.startswith('//'):
                break
            if in_origin:
                # Remove line numbers and spaces
                cleaned = ''.join(c for c in line if c.isalpha())
                sequence_parts.append(cleaned)
    return ''.join(sequence_parts).upper()


def write_fasta(filepath: str, header: str, sequence: str, line_width: int = 60):
    """Write a sequence to FASTA format."""
    with open(filepath, 'w') as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i:i+line_width] + '\n')


def clean_dna(sequence: str) -> str:
    """Normalize to uppercase DNA alphabet while preserving frame length."""
    normalized = []
    for base in sequence.upper():
        if base in {'A', 'T', 'G', 'C', 'N'}:
            normalized.append(base)
        elif base.isalpha():
            normalized.append('N')
    return ''.join(normalized)


def is_valid_cds(sequence: str, start_codons: set[str]) -> bool:
    """Basic CDS validation: frame, start codon, terminal stop, no internal stops."""
    seq = clean_dna(sequence)
    if len(seq) < 6 or len(seq) % 3 != 0:
        return False
    if seq[:3] not in start_codons:
        return False
    if seq[-3:] not in STOP_CODONS:
        return False
    for i in range(3, len(seq) - 3, 3):
        if seq[i:i+3] in STOP_CODONS:
            return False
    return True


def prepare_transit_peptide(sequence: str, peptide_name: str) -> str:
    """
    Validate and normalize a transit peptide coding sequence for N-terminal fusion.

    Rules:
      - DNA only, frame-aligned length.
      - Must start with ATG.
      - No internal stop codons.
      - Terminal stop codon is removed automatically (if present).
    """
    seq = clean_dna(sequence)
    if len(seq) < 9 or len(seq) % 3 != 0:
        raise ValueError(f"Transit peptide {peptide_name} must be >=9 bp and divisible by 3.")
    if seq[:3] != "ATG":
        raise ValueError(f"Transit peptide {peptide_name} must start with ATG.")

    for i in range(0, len(seq) - 3, 3):
        codon = seq[i:i+3]
        if codon in STOP_CODONS:
            raise ValueError(f"Transit peptide {peptide_name} has an internal stop codon at {i+1}.")

    if seq[-3:] in STOP_CODONS:
        seq = seq[:-3]

    if len(seq) < 9 or len(seq) % 3 != 0:
        raise ValueError(f"Transit peptide {peptide_name} became invalid after stop-codon trimming.")
    return seq


def fuse_transit_peptide(gene: GeneticElement, transit_dna: str, transit_name: str) -> GeneticElement:
    """
    Fuse N-terminal transit peptide CDS to a gene CDS in frame.

    The peptide start codon is kept; the native gene start codon is removed.
    The gene terminal stop codon is retained.
    """
    cds = clean_dna(gene.sequence)
    if len(cds) < 6 or len(cds) % 3 != 0:
        raise ValueError(f"Gene {gene.name} CDS must be frame-aligned before transit fusion.")
    if cds[-3:] not in STOP_CODONS:
        raise ValueError(f"Gene {gene.name} CDS must have a terminal stop codon before transit fusion.")

    fused = transit_dna + cds[3:]
    if len(fused) % 3 != 0:
        raise ValueError(f"Fused sequence for {gene.name} is out-of-frame.")
    if fused[:3] != "ATG":
        raise ValueError(f"Fused sequence for {gene.name} must start with ATG.")
    if fused[-3:] not in STOP_CODONS:
        raise ValueError(f"Fused sequence for {gene.name} must end with a stop codon.")
    for i in range(3, len(fused) - 3, 3):
        if fused[i:i+3] in STOP_CODONS:
            raise ValueError(f"Fused sequence for {gene.name} has an internal stop codon at {i+1}.")

    return GeneticElement(
        name=f"{gene.name}_vTP",
        element_type=gene.element_type,
        sequence=fused,
        description=f"{gene.description} | N-terminal vacuolar transit peptide fusion ({transit_name})",
    )


def maybe_apply_transit_peptides(
    params: dict,
    melA: GeneticElement,
    tyrp1: GeneticElement,
    dct: GeneticElement,
) -> tuple[GeneticElement, GeneticElement, GeneticElement]:
    """Apply configured transit peptide fusions to melA/TYRP1/DCT when enabled."""
    construct_cfg = params.get("construct", {})
    if not bool(construct_cfg.get("use_vacuolar_transit_peptides", False)):
        return melA, tyrp1, dct

    transit_map = construct_cfg.get("transit_peptides", {})
    if not isinstance(transit_map, dict):
        raise ValueError("construct.transit_peptides must be a mapping in the active config")

    gene_plan = [
        ("melA", melA),
        ("TYRP1", tyrp1),
        ("DCT", dct),
    ]
    fused_genes = []
    print("\nApplying vacuolar transit peptide fusions...")
    for gene_name, gene in gene_plan:
        fasta_name = transit_map.get(gene_name)
        if not fasta_name:
            raise ValueError(f"Missing construct.transit_peptides.{gene_name} in the active config")
        fasta_path = SEQ_DIR / str(fasta_name)
        if not fasta_path.exists():
            raise FileNotFoundError(f"Transit peptide FASTA missing for {gene_name}: {fasta_path}")

        transit_raw = read_fasta(str(fasta_path))
        transit_clean = prepare_transit_peptide(transit_raw, gene_name)
        fused = fuse_transit_peptide(gene, transit_clean, fasta_path.name)
        fused_genes.append(fused)
        print(
            f"  ✓ {gene_name}: transit {fasta_path.name} ({len(transit_clean)} bp) "
            f"+ mature CDS ({len(gene.sequence) - 3} bp) -> {len(fused.sequence)} bp"
        )

    return tuple(fused_genes)


def find_best_orf(
    dna_sequence: str,
    start_codons: set[str],
    target_len: Optional[int] = None,
    min_len: int = 600,
    max_len: int = 2000,
) -> Optional[str]:
    """
    Find a plausible ORF in the forward strand using first in-frame stop codon.
    Preference order:
      1) Closest to target length (if provided)
      2) Longer ORF
      3) Earlier start position
    """
    seq = clean_dna(dna_sequence)
    candidates = []
    for start in range(0, len(seq) - 2):
        codon = seq[start:start+3]
        if codon not in start_codons:
            continue
        for end in range(start + 3, len(seq) - 2, 3):
            stop = seq[end:end+3]
            if stop in STOP_CODONS:
                cds = seq[start:end+3]
                if min_len <= len(cds) <= max_len:
                    candidates.append((start, cds))
                break

    if not candidates:
        return None

    if target_len is None:
        target_len = max(len(cds) for _, cds in candidates)

    ranked = sorted(
        candidates,
        key=lambda item: (abs(len(item[1]) - target_len), -len(item[1]), item[0])
    )
    return ranked[0][1]


def reverse_complement(sequence: str) -> str:
    """Return reverse complement for DNA (N preserved)."""
    return sequence.translate(str.maketrans("ATGCN", "TACGN"))[::-1]


def resolve_melA_cds() -> tuple[str, str]:
    """
    Resolve melA CDS with validation/fallback.

    Priority:
      1) melC2_CDS_real.fasta if valid CDS
      2) Annotated range 331..1278 from GenBank if valid
      3) Best ORF scan from GenBank with bacterial start codons
    """
    mel_real_path = SEQ_DIR / "melC2_CDS_real.fasta"
    mel_gb_path = SEQ_DIR / "melA_tyrosinase.gb"

    if mel_real_path.exists():
        seq_real = clean_dna(read_fasta(str(mel_real_path)))
        if is_valid_cds(seq_real, BACTERIAL_START_CODONS):
            return seq_real, f"{mel_real_path.name} (validated)"
        print("  ⚠ melC2_CDS_real.fasta is not a valid CDS; falling back to GenBank resolution")

    if mel_gb_path.exists():
        gb_seq = clean_dna(read_genbank_sequence(str(mel_gb_path)))

        # Historical annotation in this project: CDS 331..1278 (1-based inclusive)
        if len(gb_seq) >= 1278:
            annotated = gb_seq[330:1278]
            if is_valid_cds(annotated, BACTERIAL_START_CODONS):
                return annotated, f"{mel_gb_path.name}:CDS(331..1278)"

        orf = find_best_orf(
            gb_seq,
            start_codons=BACTERIAL_START_CODONS,
            target_len=948,
            min_len=750,
            max_len=1200,
        )
        if orf is None and "N" in gb_seq:
            # Fallback for ambiguous loci where a single unknown base disrupts frame.
            orf = find_best_orf(
                gb_seq.replace("N", ""),
                start_codons=BACTERIAL_START_CODONS,
                target_len=948,
                min_len=750,
                max_len=1200,
            )

        if orf and is_valid_cds(orf, BACTERIAL_START_CODONS):
            return orf, f"{mel_gb_path.name}:auto_orf_forward"

        # If forward strand fails, try reverse complement.
        rc_seq = reverse_complement(gb_seq)
        orf_rc = find_best_orf(
            rc_seq,
            start_codons=BACTERIAL_START_CODONS,
            target_len=948,
            min_len=750,
            max_len=1200,
        )
        if orf_rc is None and "N" in rc_seq:
            orf_rc = find_best_orf(
                rc_seq.replace("N", ""),
                start_codons=BACTERIAL_START_CODONS,
                target_len=948,
                min_len=750,
                max_len=1200,
            )
        if orf_rc and is_valid_cds(orf_rc, BACTERIAL_START_CODONS):
            return orf_rc, f"{mel_gb_path.name}:auto_orf_reverse"

    raise ValueError(
        "Unable to resolve a valid melA CDS. "
        "Provide a valid data/sequences/melC2_CDS_real.fasta or fix melA_tyrosinase.gb."
    )


def resolve_strict_cds(real_path: Path, fallback_path: Path) -> tuple[str, str]:
    """
    Resolve CDS requiring canonical ATG...STOP structure.
    Used for eukaryotic TYRP1/DCT sequences.
    """
    if real_path.exists():
        seq_real = clean_dna(read_fasta(str(real_path)))
        if is_valid_cds(seq_real, STRICT_START_CODONS):
            return seq_real, f"{real_path.name} (validated)"
        print(f"  ⚠ {real_path.name} is not a valid CDS; falling back to {fallback_path.name}")

    seq_fallback = clean_dna(read_fasta(str(fallback_path)))
    if is_valid_cds(seq_fallback, STRICT_START_CODONS):
        return seq_fallback, f"{fallback_path.name} (fallback)"

    raise ValueError(
        f"Unable to resolve valid CDS from {real_path.name} or {fallback_path.name}."
    )


# ── CaMV 35S Promoter (for selection marker) ─────────────────────────────

# Minimal CaMV 35S promoter (-343 to +8)
P35S_SEQUENCE = (
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    "ATCGATCGATCGATCGATCAAGAGGATCTCAGAAGAATTTGTCCACTTATGTTTCAGAATC"
    "TAAACTCAGTCAATTTTCACTGTCAACAGAAATAAACAATGGCATTCTGATATCCTTTTTC"
    "AAATTCATATTTGTGCTAGGGAAACTTAGTTTTGAGAGAAGATTTCGTATTTTGTCAGAAT"
    "TCGTTGATCATGCAAAAGTCCCAATTTTGTAATCAAAGAAAAGGATCTTCAAGGTCTCTTT"
    "TGGAGATATAGGTAAGTTCTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGA"
    "GCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCAT"
    "GGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGC"
)

# nptII gene (Kanamycin resistance) - synthetic, cotton-codon-compatible
NPTII_SEQUENCE = (
    "ATGATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTCG"
    "GCTATGACTGGGCACAACAGACAATCGGCTGCTCTGATGCCGCCGTGTTCCGGCTGTCAGC"
    "GCAGGGGCGCCCGGTTCTTTTTGTCAAGACCGACCTGTCCGGTGCCCTGAATGAACTGCAG"
    "GACGAGGCAGCGCGGCTATCGTGGCTGGCCACGACGGGCGTTCCTTGCGCAGCTGTGCTCG"
    "ACGTTGTCACTGAAGCGGGAAGGGACTGGCTGCTATTGGGCGAAGTGCCGGGGCAGGATCT"
    "CCTGTCATCTCACCTTGCTCCTGCCGAGAAAGTATCCATCATGGCTGATGCAATGCGGCGG"
    "CTGCATACGCTTGATCCGGCTACCTGCCCATTCGACCACCAAGCGAAACATCGCATCGAGC"
    "GAGCACGTACTCGGATGGAAGCCGGTCTTGTCGATCAGGATGATCTGGACGAAGAGCATCA"
    "GGGGCTCGCGCCAGCCGAACTGTTCGCCAGGCTCAAGGCGCGCATGCCCGACGGCGATGAT"
    "CTCGTCGTGACCCATGGCGATGCCTGCTTGCCGAATATCATGGTGGAAAATGGCCGCTTTT"
    "CTGGATTCATCGACTGTGGCCGGCTGGGTGTGGCGGACCGCTATCAGGACATAGCGTTGGC"
    "TACCCGTGATATTGCTGAAGAGCTTGGCGGCGAATGGGCTGACCGCTTCCTCGTGCTTTAC"
    "GGTATCGCCGCTCCCGATTCGCAGCGCATCGCCTTCTATCGCCTTCTTGACGAGTTCTTCT"
    "GA"
)


# ── Main Assembly ─────────────────────────────────────────────────────────

def build_construct() -> TDNAConstruct:
    """Assemble the complete T-DNA construct for black cotton."""
    params = load_params()

    print("\n🧬 BlackCotton Construct Designer")
    print("=" * 50)
    print("Loading genetic elements...\n")
    
    # ── Load sequences ────────────────────────────────────────────────
    
    # Promoters
    pGhMat1 = GeneticElement(
        name="pGhMat1",
        element_type="promoter",
        sequence=read_fasta(str(SEQ_DIR / "pGhMat1_promoter.fasta")),
        description="Maturation-specific promoter, active at ~35 DPA"
    )
    print(f"  ✓ Loaded {pGhMat1.name}: {pGhMat1.length} bp")
    
    pGhSCW = GeneticElement(
        name="pGhSCW-late",
        element_type="promoter",
        sequence=read_fasta(str(SEQ_DIR / "pGhSCW_late_promoter.fasta")),
        description="Late SCW promoter, active at ~30 DPA"
    )
    print(f"  ✓ Loaded {pGhSCW.name}: {pGhSCW.length} bp")
    
    p35S = GeneticElement(
        name="p35S",
        element_type="promoter",
        sequence=P35S_SEQUENCE.upper().replace('\n', ''),
        description="CaMV 35S constitutive promoter"
    )
    print(f"  ✓ Loaded {p35S.name}: {p35S.length} bp")
    
    # Genes
    # melA resolver validates CDS integrity and falls back to GenBank ORF extraction.
    melA_seq, melA_source = resolve_melA_cds()
    melA = GeneticElement(
        name="melA",
        element_type="gene",
        sequence=melA_seq,
        description="Tyrosinase CDS from S. antibioticus (NCBI M11582). Tyr → DOPA → Dopaquinone"
    )
    print(f"  ✓ Loaded {melA.name}: {melA.length} bp (source: {melA_source})")

    # Use real TYRP1/DCT CDS when valid; otherwise fallback to legacy optimized versions.
    tyrp1_seq, tyrp1_source = resolve_strict_cds(
        real_path=SEQ_DIR / "TYRP1_CDS_real.fasta",
        fallback_path=SEQ_DIR / "GhTYRP1_optimized.fasta",
    )
    tyrp1 = GeneticElement(
        name="TYRP1",
        element_type="gene",
        sequence=tyrp1_seq,
        description="Tyrosinase-related protein 1 CDS (NCBI NM_000550, H. sapiens)"
    )
    print(f"  ✓ Loaded {tyrp1.name}: {tyrp1.length} bp (source: {tyrp1_source})")

    dct_seq, dct_source = resolve_strict_cds(
        real_path=SEQ_DIR / "DCT_CDS_real.fasta",
        fallback_path=SEQ_DIR / "GhDCT_optimized.fasta",
    )
    dct = GeneticElement(
        name="DCT",
        element_type="gene",
        sequence=dct_seq,
        description="Dopachrome tautomerase CDS (NCBI NM_001922, H. sapiens)"
    )
    print(f"  ✓ Loaded {dct.name}: {dct.length} bp (source: {dct_source})")

    melA, tyrp1, dct = maybe_apply_transit_peptides(params, melA, tyrp1, dct)
    if any(g.name.endswith("_vTP") for g in (melA, tyrp1, dct)):
        print("  ✓ Transit peptide compartmentalization enabled for melA/TYRP1/DCT")
    
    nptII = GeneticElement(
        name="nptII",
        element_type="gene",
        sequence=NPTII_SEQUENCE.upper().replace('\n', ''),
        description="Neomycin phosphotransferase II (Kanamycin resistance)"
    )
    print(f"  ✓ Loaded {nptII.name}: {nptII.length} bp")
    
    # Terminator
    tNOS = GeneticElement(
        name="tNOS",
        element_type="terminator",
        sequence=read_fasta(str(SEQ_DIR / "tNOS_terminator.fasta")),
        description="Nopaline synthase terminator"
    )
    print(f"  ✓ Loaded {tNOS.name}: {tNOS.length} bp")
    
    # ── Assemble Cassettes ────────────────────────────────────────────
    print("\nAssembling cassettes...")
    
    cassette1 = Cassette(
        name="Melanin Primary (melA)",
        promoter=pGhMat1,
        gene=melA,
        terminator=tNOS,
        description="Primary melanin enzyme. Converts L-Tyrosine → L-DOPA → Dopaquinone. "
                    "Driven by maturation promoter to fire AFTER cellulose deposition."
    )
    print(f"  ✓ Cassette 1: {cassette1.name} ({cassette1.length} bp)")
    
    cassette2 = Cassette(
        name="TYRP1 Support",
        promoter=pGhSCW,
        gene=tyrp1,
        terminator=tNOS,
        description="Stabilizes melanin intermediates. Prevents toxic accumulation "
                    "of dopaquinone by channeling it toward eumelanin."
    )
    print(f"  ✓ Cassette 2: {cassette2.name} ({cassette2.length} bp)")
    
    cassette3 = Cassette(
        name="DCT Support",
        promoter=pGhSCW,
        gene=dct,
        terminator=tNOS,
        description="Converts dopachrome to DHICA, improving melanin polymer quality "
                    "and ensuring deep black color (vs. brown)."
    )
    print(f"  ✓ Cassette 3: {cassette3.name} ({cassette3.length} bp)")
    
    cassette4 = Cassette(
        name="Selection Marker (nptII)",
        promoter=p35S,
        gene=nptII,
        terminator=tNOS,
        description="Constitutive kanamycin resistance for transgenic selection."
    )
    print(f"  ✓ Cassette 4: {cassette4.name} ({cassette4.length} bp)")
    
    # ── Build T-DNA ───────────────────────────────────────────────────
    print("\nBuilding T-DNA construct...")
    
    construct = TDNAConstruct(
        name="pBC-MelaninCotton-v2",
        cassettes=[cassette1, cassette2, cassette3, cassette4]
    )
    
    print(f"  ✓ Construct assembled: {construct.total_length:,} bp\n")
    
    return construct


def save_results(construct: TDNAConstruct):
    """Save all construct outputs to the results directory."""
    
    os.makedirs(RESULTS_DIR, exist_ok=True)
    
    # ── Construct map (text) ──────────────────────────────────────────
    map_path = RESULTS_DIR / "construct_map.txt"
    with open(map_path, 'w') as f:
        f.write(construct.get_construct_map())
    print(f"  📄 Saved: {map_path}")
    
    # ── Full sequence (FASTA) ─────────────────────────────────────────
    fasta_path = RESULTS_DIR / "pBC_MelaninCotton_v1.fasta"
    write_fasta(
        str(fasta_path),
        f"{construct.name} | {construct.total_length} bp | BlackCotton T-DNA",
        construct.sequence
    )
    print(f"  📄 Saved: {fasta_path}")
    
    # ── Machine-readable summary (JSON) ───────────────────────────────
    json_path = RESULTS_DIR / "construct_summary.json"
    with open(json_path, 'w') as f:
        json.dump(construct.to_json(), f, indent=2)
    print(f"  📄 Saved: {json_path}")


# ── Entry Point ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    construct = build_construct()
    
    # Print the construct map to console
    print(construct.get_construct_map())
    
    # Save all outputs
    print("\nSaving results...")
    save_results(construct)
    print("\n✅ Construct design complete!")
