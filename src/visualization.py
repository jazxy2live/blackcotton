#!/usr/bin/env python3
"""
visualization.py — Generate all analysis plots for BlackCotton
================================================================
Produces publication-quality figures for construct design, expression
kinetics, melanin accumulation, fiber quality, and wash fastness.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import json
from pathlib import Path

from src.config_loader import load_config

BASE_DIR = Path(__file__).resolve().parent.parent
CONFIG_DIR = BASE_DIR / "config"
RESULTS_DIR = BASE_DIR / "results"

# Style
plt.rcParams.update({
    'figure.facecolor': '#0a0a0a',
    'axes.facecolor': '#111111',
    'axes.edgecolor': '#333333',
    'text.color': '#e0e0e0',
    'axes.labelcolor': '#e0e0e0',
    'xtick.color': '#999999',
    'ytick.color': '#999999',
    'grid.color': '#222222',
    'font.family': 'sans-serif',
    'font.size': 11,
})

COLORS = {
    'melanin': '#1a1a2e',
    'cellulose': '#16a085',
    'melA': '#e74c3c',
    'TYRP1': '#e67e22',
    'DCT': '#9b59b6',
    'nptII': '#3498db',
    'white': '#f5f5dc',
    'brown': '#8B4513',
    'dyed': '#1a1a1a',
    'engineered': '#0d0d0d',
    'accent': '#00d4ff',
    'warning': '#ff6b6b',
}


def load_params():
    return load_config()


def plot_expression_kinetics():
    """Plot gene expression over fiber development timeline."""
    data_path = RESULTS_DIR / "expression_data.npz"
    if not data_path.exists():
        print("  ⚠ No expression data — skipping expression plot")
        return

    data = np.load(str(data_path))
    t = data['t_dpa']
    
    fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=True)
    fig.suptitle('Gene Expression Kinetics in Cotton Fiber Development',
                 fontsize=16, fontweight='bold', color='#00d4ff', y=0.98)

    # Panel 1: Promoter activity
    ax = axes[0]
    ax.plot(t, data['promoter_melA'], color=COLORS['melA'], lw=2.5, label='pGhMat1 (melA)')
    ax.plot(t, data['promoter_SCW_late'], color=COLORS['TYRP1'], lw=2.5, label='pGhSCW-late (TYRP1/DCT)')
    ax.set_ylabel('Promoter Activity\n(relative)')
    ax.legend(framealpha=0.3, loc='upper left')
    ax.set_title('Promoter Activation', fontsize=12, color='#cccccc')
    ax.axvspan(35, 50, alpha=0.1, color=COLORS['melA'], label='Melanin window')
    ax.axvline(35, color=COLORS['warning'], ls='--', alpha=0.5, lw=1)
    ax.annotate('Switch ON', xy=(35, 0.5), fontsize=9, color=COLORS['warning'])
    ax.grid(True, alpha=0.3)

    # Panel 2: mRNA levels
    ax = axes[1]
    ax.plot(t, data['mRNA_melA'], color=COLORS['melA'], lw=2, label='melA mRNA')
    ax.plot(t, data['mRNA_TYRP1'], color=COLORS['TYRP1'], lw=2, label='TYRP1 mRNA')
    ax.plot(t, data['mRNA_DCT'], color=COLORS['DCT'], lw=2, label='DCT mRNA')
    ax.plot(t, data['mRNA_nptII'], color=COLORS['nptII'], lw=1.5, ls='--', label='nptII mRNA', alpha=0.6)
    ax.set_ylabel('mRNA Level\n(molecules)')
    ax.legend(framealpha=0.3, loc='upper left')
    ax.set_title('Transcript Accumulation', fontsize=12, color='#cccccc')
    ax.grid(True, alpha=0.3)

    # Panel 3: Protein levels + cellulose
    ax = axes[2]
    ax2 = ax.twinx()
    
    ax.plot(t, data['protein_melA'], color=COLORS['melA'], lw=2.5, label='melA protein')
    ax.plot(t, data['protein_TYRP1'], color=COLORS['TYRP1'], lw=2.5, label='TYRP1 protein')
    ax.plot(t, data['protein_DCT'], color=COLORS['DCT'], lw=2.5, label='DCT protein')
    ax.set_ylabel('Protein Level\n(molecules)', color='#e0e0e0')
    ax.set_xlabel('Days Post Anthesis (DPA)')
    
    ax2.fill_between(t, 0, data['cellulose'], alpha=0.2, color=COLORS['cellulose'])
    ax2.plot(t, data['cellulose'], color=COLORS['cellulose'], lw=2, ls='--', label='Cellulose')
    ax2.set_ylabel('Cellulose\n(mg/fiber)', color=COLORS['cellulose'])
    ax2.tick_params(axis='y', labelcolor=COLORS['cellulose'])
    
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1+lines2, labels1+labels2, framealpha=0.3, loc='center left')
    ax.set_title('Protein Accumulation vs Cellulose Deposition', fontsize=12, color='#cccccc')
    ax.grid(True, alpha=0.3)

    # Stage annotations
    for ax_i in axes:
        ax_i.axvspan(0, 3, alpha=0.05, color='cyan')
        ax_i.axvspan(3, 21, alpha=0.05, color='green')
        ax_i.axvspan(16, 40, alpha=0.05, color='yellow')
        ax_i.axvspan(35, 50, alpha=0.08, color='red')

    plt.tight_layout()
    out = RESULTS_DIR / "expression_kinetics.png"
    plt.savefig(str(out), dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  📊 Saved: {out}")


def plot_melanin_accumulation():
    """Plot melanin pathway metabolite kinetics."""
    data_path = RESULTS_DIR / "melanin_data.npz"
    if not data_path.exists():
        print("  ⚠ No melanin data — skipping melanin plot")
        return

    data = np.load(str(data_path))
    t = data['t_dpa']

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 9), sharex=True)
    fig.suptitle('Melanin Biosynthesis Pathway Kinetics',
                 fontsize=16, fontweight='bold', color='#00d4ff', y=0.98)

    # Panel 1: Substrate and product
    ax1.plot(t, data['tyrosine'], color='#f39c12', lw=2.5, label='L-Tyrosine (substrate)')
    ax1.plot(t, data['melanin'], color='#1a1a2e', lw=3, label='Melanin (product)')
    ax1.fill_between(t, 0, data['melanin'], alpha=0.3, color='#1a1a2e')
    ax1.set_ylabel('Concentration (mM)')
    ax1.legend(framealpha=0.3)
    ax1.set_title('Substrate Consumption → Melanin Production', fontsize=12, color='#cccccc')
    ax1.axvline(35, color=COLORS['warning'], ls='--', alpha=0.5)
    ax1.grid(True, alpha=0.3)

    # Panel 2: Intermediates
    ax2.plot(t, data['L_DOPA'], lw=1.5, label='L-DOPA', color='#2ecc71')
    ax2.plot(t, data['dopaquinone'], lw=1.5, label='Dopaquinone', color='#e74c3c')
    ax2.plot(t, data['dopachrome'], lw=1.5, label='Dopachrome', color='#9b59b6')
    ax2.plot(t, data['DHICA'], lw=1.5, label='DHICA', color='#3498db')
    ax2.plot(t, data['indole_quinone'], lw=1.5, label='Indole-quinone', color='#e67e22')
    ax2.set_ylabel('Concentration (mM)')
    ax2.set_xlabel('Days Post Anthesis (DPA)')
    ax2.legend(framealpha=0.3, ncol=3)
    ax2.set_title('Pathway Intermediates (should be low → no toxic buildup)', fontsize=12, color='#cccccc')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    out = RESULTS_DIR / "melanin_accumulation.png"
    plt.savefig(str(out), dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  📊 Saved: {out}")


def plot_fiber_comparison():
    """Bar chart comparing fiber quality across cotton types."""
    json_path = RESULTS_DIR / "fiber_quality_comparison.json"
    if not json_path.exists():
        print("  ⚠ No fiber data — skipping fiber comparison plot")
        return

    with open(json_path) as f:
        fibers = json.load(f)

    names = [fb['name'][:25] for fb in fibers]
    metrics = ['uhml_mm', 'strength_g_tex', 'micronaire', 'color_L']
    labels = ['Length (mm)', 'Strength (g/tex)', 'Micronaire', 'L* (lightness)']
    bar_colors = ['#f5f5dc', '#8B4513', '#1a1a1a', '#0d0d0d', '#222222', '#333333']

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    fig.suptitle('Fiber Quality Comparison — BlackCotton vs Conventional',
                 fontsize=16, fontweight='bold', color='#00d4ff', y=0.98)

    for idx, (ax, metric, label) in enumerate(zip(axes.flat, metrics, labels)):
        vals = [fb[metric] for fb in fibers]
        bars = ax.barh(range(len(names)), vals,
                      color=bar_colors[:len(names)], edgecolor='#444', height=0.6)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names, fontsize=9)
        ax.set_xlabel(label)
        ax.set_title(label, fontsize=12, color='#cccccc')
        ax.grid(True, axis='x', alpha=0.3)
        for bar, val in zip(bars, vals):
            ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height()/2,
                   f'{val:.1f}', va='center', fontsize=9, color='#aaa')

    plt.tight_layout()
    out = RESULTS_DIR / "fiber_quality_comparison.png"
    plt.savefig(str(out), dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  📊 Saved: {out}")


def plot_wash_fastness():
    """Plot color fade over 50 washes."""
    json_path = RESULTS_DIR / "wash_fastness_data.json"
    if not json_path.exists():
        print("  ⚠ No wash data — skipping wash fastness plot")
        return

    with open(json_path) as f:
        wash_data = json.load(f)

    fig, ax = plt.subplots(figsize=(12, 7))
    fig.suptitle('Wash Fastness — Color Stability Over 50 Washes',
                 fontsize=16, fontweight='bold', color='#00d4ff', y=0.98)

    colors_map = {'Dyed': '#e74c3c', 'Engineered': '#00d4ff', 'White': '#f5f5dc', 'Brown': '#8B4513'}
    for name, washes in wash_data.items():
        ws = [w['wash'] for w in washes]
        ls = [w['L_star'] for w in washes]
        c = '#888888'
        lw = 1.5
        for key, col in colors_map.items():
            if key in name:
                c = col
                lw = 2.5 if key in ['Dyed', 'Engineered'] else 1.5
                break
        short_name = name[:35]
        ax.plot(ws, ls, color=c, lw=lw, label=short_name)

    ax.set_xlabel('Number of Washes')
    ax.set_ylabel('L* (Lightness — lower = darker)')
    ax.legend(framealpha=0.3, fontsize=9, loc='center right')
    ax.grid(True, alpha=0.3)
    ax.annotate('Dyed black fades →', xy=(30, 30), fontsize=10, color=COLORS['warning'])
    ax.annotate('Engineered black stays ●', xy=(30, 10), fontsize=10, color=COLORS['accent'])

    plt.tight_layout()
    out = RESULTS_DIR / "wash_fastness.png"
    plt.savefig(str(out), dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  📊 Saved: {out}")


def plot_construct_map():
    """Visual map of the genetic construct."""
    json_path = RESULTS_DIR / "construct_summary.json"
    if not json_path.exists():
        print("  ⚠ No construct data — skipping construct map")
        return

    with open(json_path) as f:
        construct = json.load(f)

    fig, ax = plt.subplots(figsize=(16, 4))
    fig.suptitle(f"T-DNA Construct: {construct['name']} ({construct['total_length_bp']:,} bp)",
                 fontsize=14, fontweight='bold', color='#00d4ff')

    type_colors = {
        'T-DNA_border': '#666666', 'promoter': '#2ecc71',
        'CDS': '#e74c3c', 'terminator': '#3498db'
    }
    total = construct['total_length_bp']
    y_center = 0.5

    for feat in construct['features']:
        x_start = feat['start'] / total
        width = (feat['end'] - feat['start']) / total
        color = type_colors.get(feat['type'], '#888888')
        rect = mpatches.FancyBboxPatch((x_start, y_center-0.15), width, 0.3,
            boxstyle="round,pad=0.01", facecolor=color, edgecolor='white',
            linewidth=0.5, alpha=0.85)
        ax.add_patch(rect)
        if width > 0.03:
            ax.text(x_start + width/2, y_center, feat['name'],
                   ha='center', va='center', fontsize=7, fontweight='bold', color='white')

    # Legend
    legend_elements = [mpatches.Patch(facecolor=c, label=t) for t, c in type_colors.items()]
    ax.legend(handles=legend_elements, loc='upper right', framealpha=0.3, fontsize=9)

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Relative Position in T-DNA')
    ax.set_yticks([])
    plt.tight_layout()
    out = RESULTS_DIR / "construct_visual_map.png"
    plt.savefig(str(out), dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  📊 Saved: {out}")


if __name__ == "__main__":
    print("\n🧬 BlackCotton Visualization Engine")
    print("="*50)
    import os; os.makedirs(RESULTS_DIR, exist_ok=True)
    plot_construct_map()
    plot_expression_kinetics()
    plot_melanin_accumulation()
    plot_fiber_comparison()
    plot_wash_fastness()
    print("\n✅ All visualizations generated!")
