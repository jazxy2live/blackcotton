#!/usr/bin/env python3
"""
dashboard_generator.py — Interactive Scientific Dashboard for BlackCotton
=========================================================================

Generates a self-contained HTML dashboard with:
  1. 3D protein structure viewers (3Dmol.js) for melA, TYRP1, DCT
  2. Gene expression dynamics over time (Chart.js)
  3. Melanin pathway metabolite profiles
  4. Fiber quality comparison
  5. T-DNA construct map

Usage:
    python -m src.dashboard_generator
"""

import json
import webbrowser
from pathlib import Path

import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"


def load_expression_data():
    """Load expression simulation NPZ and return JSON-serializable dict."""
    path = RESULTS_DIR / "expression_data.npz"
    if not path.exists():
        return None
    data = np.load(str(path))
    # Downsample to ~200 points for the browser
    step = max(1, len(data["t_dpa"]) // 200)
    return {
        "t_dpa": data["t_dpa"][::step].tolist(),
        "protein_melA": data["protein_melA"][::step].tolist(),
        "protein_TYRP1": data["protein_TYRP1"][::step].tolist(),
        "protein_DCT": data["protein_DCT"][::step].tolist(),
        "cellulose": data["cellulose"][::step].tolist(),
        "promoter_melA": data["promoter_melA"][::step].tolist(),
        "promoter_SCW_late": data["promoter_SCW_late"][::step].tolist(),
    }


def load_melanin_data():
    """Load melanin pathway NPZ and return JSON-serializable dict."""
    path = RESULTS_DIR / "melanin_data.npz"
    if not path.exists():
        return None
    data = np.load(str(path))
    step = max(1, len(data["t_dpa"]) // 200)
    return {
        "t_dpa": data["t_dpa"][::step].tolist(),
        "tyrosine": data["tyrosine"][::step].tolist(),
        "L_DOPA": data["L_DOPA"][::step].tolist(),
        "dopaquinone": data["dopaquinone"][::step].tolist(),
        "dopachrome": data["dopachrome"][::step].tolist(),
        "DHICA": data["DHICA"][::step].tolist(),
        "indole_quinone": data["indole_quinone"][::step].tolist(),
        "melanin": data["melanin"][::step].tolist(),
    }


def load_fiber_data():
    path = RESULTS_DIR / "fiber_quality_comparison.json"
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


def load_construct_data():
    path = RESULTS_DIR / "construct_summary.json"
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


def load_wash_data():
    path = RESULTS_DIR / "wash_fastness_data.json"
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


def generate_html(expr_data, melanin_data, fiber_data, construct_data, wash_data):
    """Generate complete self-contained HTML dashboard."""
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>BlackCotton — Scientific Dashboard</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/4.4.1/chart.umd.min.js"></script>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
:root {{
  --bg: #ffffff;
  --card: #fdfdfd;
  --card-border: #e0e0e0;
  --text: #111111;
  --text-dim: #555555;
  --accent: #000000;
  --accent2: #222222;
  --accent3: #444444;
  --accent4: #666666;
  --danger: #d32f2f;
  --success: #388e3c;
}}
body {{
  font-family: 'Inter', -apple-system, sans-serif;
  background: var(--bg);
  color: var(--text);
  line-height: 1.6;
  overflow-x: hidden;
}}
.hero {{
  text-align: center;
  padding: 60px 20px 40px;
  background: var(--bg);
  border-bottom: 1px solid var(--card-border);
}}
.hero h1 {{
  font-size: 2.8rem;
  font-weight: 700;
  letter-spacing: -1px;
  color: var(--text);
  margin-bottom: 8px;
}}
.hero p {{
  font-size: 1.1rem;
  color: var(--text-dim);
  max-width: 700px;
  margin: 0 auto;
}}
.hero .tag {{
  display: inline-block;
  padding: 4px 14px;
  border-radius: 4px;
  font-size: 0.78rem;
  font-weight: 500;
  margin: 16px 6px 0;
  background: #f0f0f0;
  color: var(--text);
  border: 1px solid #d0d0d0;
}}
nav {{
  display: flex;
  justify-content: center;
  gap: 6px;
  padding: 16px 20px;
  background: var(--bg);
  border-bottom: 1px solid var(--card-border);
  position: sticky;
  top: 0;
  z-index: 100;
}}
nav button {{
  background: #f9f9f9;
  color: var(--text-dim);
  border: 1px solid var(--card-border);
  padding: 8px 18px;
  border-radius: 4px;
  font-family: 'Inter', sans-serif;
  font-size: 0.85rem;
  font-weight: 500;
  cursor: pointer;
}}
nav button.active {{
  background: #111111;
  color: #ffffff;
  border-color: #111111;
}}
.container {{
  max-width: 1400px;
  margin: 0 auto;
  padding: 24px;
}}
.section {{
  display: none;
}}
.section.active {{ display: block; }}
.grid-2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-bottom: 20px; }}
.grid-3 {{ display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 20px; margin-bottom: 20px; }}
@media(max-width: 900px) {{ .grid-2, .grid-3 {{ grid-template-columns: 1fr; }} }}
.card {{
  background: var(--card);
  border: 1px solid var(--card-border);
  border-radius: 8px;
  padding: 24px;
}}
.card h3 {{
  font-size: 1rem;
  font-weight: 600;
  margin-bottom: 16px;
  display: flex;
  align-items: center;
  gap: 8px;
}}
.card h3 .icon {{ font-size: 1.2rem; }}
.card.full {{ grid-column: 1 / -1; }}
.protein-viewer {{
  width: 100%;
  height: 380px;
  border-radius: 8px;
  overflow: hidden;
  position: relative;
  background: #f7f7f7;
  border: 1px solid var(--card-border);
}}
.protein-label {{
  position: absolute;
  bottom: 12px;
  left: 12px;
  background: #ffffff;
  padding: 6px 14px;
  border-radius: 4px;
  font-size: 0.78rem;
  color: var(--text-dim);
  z-index: 10;
  border: 1px solid var(--card-border);
}}
.protein-label strong {{ color: var(--text); }}
.section-title {{
  font-size: 1.5rem;
  font-weight: 700;
  margin-bottom: 6px;
  display: flex;
  align-items: center;
  gap: 10px;
}}
.section-sub {{
  color: var(--text-dim);
  font-size: 0.9rem;
  margin-bottom: 24px;
}}
.chart-container {{ position: relative; height: 320px; }}
.construct-map {{
  display: flex;
  align-items: center;
  gap: 0;
  padding: 24px 0;
  overflow-x: auto;
}}
.construct-element {{
  display: flex;
  flex-direction: column;
  align-items: center;
  min-width: 0;
  flex-shrink: 0;
  position: relative;
}}
.construct-element .bar {{
  height: 44px;
  border-radius: 2px;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 0.7rem;
  font-weight: 600;
  color: #111;
  white-space: nowrap;
  padding: 0 8px;
  letter-spacing: 0.3px;
  cursor: default;
  border: 1px solid #111;
}}
.construct-element .size {{
  font-size: 0.65rem;
  color: var(--text-dim);
  margin-top: 4px;
  font-family: 'JetBrains Mono', monospace;
}}
.bar.border {{ background: #f0f0f0; width: 30px; }}
.bar.promoter {{ background: #e8e8e8; }}
.bar.cds {{ background: #d0d0d0; }}
.bar.terminator {{ background: #f9f9f9; }}
.arrow {{
  width: 12px;
  height: 2px;
  background: #111111;
  flex-shrink: 0;
}}
.stat-grid {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
  gap: 12px;
  margin-bottom: 20px;
}}
.stat {{
  background: var(--card);
  border: 1px solid var(--card-border);
  border-radius: 6px;
  padding: 16px 18px;
  text-align: center;
}}
.stat .value {{
  font-size: 1.8rem;
  font-weight: 700;
  font-family: 'JetBrains Mono', monospace;
  color: var(--text);
}}
.stat .label {{
  font-size: 0.78rem;
  color: var(--text-dim);
  margin-top: 2px;
}}
.fiber-table {{
  width: 100%;
  border-collapse: collapse;
  font-size: 0.85rem;
}}
.fiber-table th {{
  text-align: left;
  padding: 10px 14px;
  border-bottom: 2px solid var(--text);
  color: var(--text);
  font-weight: 600;
  font-size: 0.78rem;
  text-transform: uppercase;
  letter-spacing: 0.5px;
}}
.fiber-table td {{
  padding: 10px 14px;
  border-bottom: 1px solid var(--card-border);
  font-family: 'JetBrains Mono', monospace;
  font-size: 0.82rem;
}}
.fiber-table .name {{ font-family: 'Inter', sans-serif; font-weight: 500; }}
.grade {{
  display: inline-block;
  padding: 2px 10px;
  border-radius: 4px;
  font-size: 0.72rem;
  font-weight: 600;
  border: 1px solid var(--card-border);
  background: #f0f0f0;
  color: #111;
}}
.color-swatch {{
  display: inline-block;
  width: 28px;
  height: 28px;
  border-radius: 4px;
  vertical-align: middle;
  border: 1px solid #ccc;
}}
</style>
</head>
<body>

<div class="hero">
  <h1>🧬 BlackCotton</h1>
  <p>Computational Design for Naturally Black Cotton Fibers via Engineered Melanin Biosynthesis</p>
  <div>
    <span class="tag purple">ODE Simulation</span>
    <span class="tag teal">Melanin Pathway</span>
    <span class="tag pink">Gossypium hirsutum</span>
  </div>
</div>

<nav id="nav">
  <button class="active" onclick="showSection('proteins')">🔬 Protein Structures</button>
  <button onclick="showSection('expression')">📈 Expression Dynamics</button>
  <button onclick="showSection('pathway')">⚗️ Melanin Pathway</button>
  <button onclick="showSection('fiber')">🧵 Fiber Quality</button>
  <button onclick="showSection('construct')">🧬 T-DNA Construct</button>
</nav>

<div class="container">

<!-- ═══════════ PROTEIN STRUCTURES ═══════════ -->
<div id="sec-proteins" class="section active">
  <h2 class="section-title"><span>🔬</span> Enzyme Structures</h2>
  <p class="section-sub">Interactive 3D structures of the three melanin pathway enzymes used in the BlackCotton construct. Drag to rotate, scroll to zoom.</p>

  <div class="grid-3">
    <div class="card">
      <h3><span class="icon">🟣</span> melA — Tyrosinase</h3>
      <div class="protein-viewer" id="viewer-melA">
        <div class="protein-label"><strong>PDB: 1WX2</strong> · Streptomyces tyrosinase</div>
      </div>
    </div>
    <div class="card">
      <h3><span class="icon">🟢</span> TYRP1</h3>
      <div class="protein-viewer" id="viewer-tyrp1">
        <div class="protein-label"><strong>PDB: 5M8Q</strong> · Human TYRP1 intramelanosomal domain</div>
      </div>
    </div>
    <div class="card">
      <h3><span class="icon">🔵</span> DCT (Dopachrome Tautomerase)</h3>
      <div class="protein-viewer" id="viewer-dct">
        <div class="protein-label"><strong>PDB: 3NQ1</strong> · Dopachrome tautomerase</div>
      </div>
    </div>
  </div>

  <div class="card full">
    <h3><span class="icon">📋</span> Enzyme Overview</h3>
    <table class="fiber-table">
      <thead>
        <tr><th>Enzyme</th><th>Gene</th><th>Source</th><th>Function</th><th>MW (kDa)</th><th>PDB</th></tr>
      </thead>
      <tbody>
        <tr>
          <td class="name">Tyrosinase</td><td><code>melA (melC2)</code></td>
          <td>S. antibioticus</td><td>L-Tyr → L-DOPA → Dopaquinone</td>
          <td>30.9</td><td>1WX2</td>
        </tr>
        <tr>
          <td class="name">TYRP1</td><td><code>TYRP1</code></td>
          <td>H. sapiens</td><td>DHICA → Indole-5,6-quinone</td>
          <td>~60</td><td>5M8L</td>
        </tr>
        <tr>
          <td class="name">DCT</td><td><code>DCT</code></td>
          <td>H. sapiens</td><td>Dopachrome → DHICA</td>
          <td>~58</td><td>3NQ1</td>
        </tr>
      </tbody>
    </table>
  </div>
</div>

<!-- ═══════════ EXPRESSION DYNAMICS ═══════════ -->
<div id="sec-expression" class="section">
  <h2 class="section-title"><span>📈</span> Gene Expression Dynamics</h2>
  <p class="section-sub">ODE simulation of transgene expression over 50 days of cotton fiber development (Days Post Anthesis). The key question: do melanin enzymes arrive AFTER cellulose is deposited?</p>

  <div class="stat-grid">
    <div class="stat"><div class="value">38.4</div><div class="label">Cellulose 90% DPA</div></div>
    <div class="stat"><div class="value">50.0</div><div class="label">Peak melA DPA</div></div>
    <div class="stat"><div class="value">50.0</div><div class="label">Peak TYRP1 DPA</div></div>
    <div class="stat"><div class="value">11.6</div><div class="label">Temporal Gap (days)</div></div>
  </div>

  <div class="grid-2">
    <div class="card">
      <h3><span class="icon">🧪</span> Protein Accumulation vs Cellulose</h3>
      <div class="chart-container"><canvas id="chart-expression"></canvas></div>
    </div>
    <div class="card">
      <h3><span class="icon">⚡</span> Promoter Activity Profiles</h3>
      <div class="chart-container"><canvas id="chart-promoters"></canvas></div>
    </div>
  </div>
</div>

<!-- ═══════════ MELANIN PATHWAY ═══════════ -->
<div id="sec-pathway" class="section">
  <h2 class="section-title"><span>⚗️</span> Melanin Biosynthesis Pathway</h2>
  <p class="section-sub">Michaelis-Menten kinetics simulation showing substrate consumption, intermediate flux, and final melanin polymer accumulation.</p>

  <div class="stat-grid">
    <div class="stat"><div class="value">0.50</div><div class="label">Final Melanin (mM eq.)</div></div>
    <div class="stat"><div class="value">100%</div><div class="label">Tyrosine Consumed</div></div>
    <div class="stat"><div class="value">50</div><div class="label">Simulation Days</div></div>
  </div>

  <div class="grid-2">
    <div class="card">
      <h3><span class="icon">📊</span> Substrate → Product Dynamics</h3>
      <div class="chart-container"><canvas id="chart-melanin"></canvas></div>
    </div>
    <div class="card">
      <h3><span class="icon">⚠️</span> Toxic Intermediate Monitoring</h3>
      <div class="chart-container"><canvas id="chart-intermediates"></canvas></div>
    </div>
  </div>

  <div class="card full">
    <h3><span class="icon">🔄</span> Pathway Architecture</h3>
    <div style="text-align:center; padding: 20px; font-family: 'JetBrains Mono', monospace; font-size: 0.9rem; line-height: 2.2; color: var(--text-dim);">
      <span style="color:#fdcb6e; font-weight:600;">L-Tyrosine</span>
      <span style="color:#666;">──[</span><span style="color:#a29bfe;">melA</span><span style="color:#666;">]──▸</span>
      <span style="color:#fdcb6e;">L-DOPA</span>
      <span style="color:#666;">──[</span><span style="color:#a29bfe;">melA</span><span style="color:#666;">]──▸</span>
      <span style="color:#e17055;">Dopaquinone</span>
      <span style="color:#666;">──(spontaneous)──▸</span>
      <span style="color:#e17055;">Leucodopachrome</span>
      <span style="color:#666;">──▸</span>
      <span style="color:#e17055;">Dopachrome</span>
      <span style="color:#666;">──[</span><span style="color:#55efc4;">DCT</span><span style="color:#666;">]──▸</span>
      <span style="color:#81ecec;">DHICA</span>
      <span style="color:#666;">──[</span><span style="color:#55efc4;">TYRP1</span><span style="color:#666;">]──▸</span>
      <span style="color:#81ecec;">Indole-quinone</span>
      <span style="color:#666;">──(polymerization)──▸</span>
      <span style="color:#dfe6e9; font-weight:700; font-size: 1.1rem;">⬛ EUMELANIN</span>
    </div>
  </div>
</div>

<!-- ═══════════ FIBER QUALITY ═══════════ -->
<div id="sec-fiber" class="section">
  <h2 class="section-title"><span>🧵</span> Fiber Quality Comparison</h2>
  <p class="section-sub">USDA HVI metrics comparing engineered black cotton against white, natural brown, and chemically dyed black cotton.</p>

  <div class="grid-2">
    <div class="card">
      <h3><span class="icon">🎨</span> Color Lightness (L*) — Lower = Darker</h3>
      <div class="chart-container"><canvas id="chart-color"></canvas></div>
    </div>
    <div class="card">
      <h3><span class="icon">💪</span> Fiber Strength (g/tex)</h3>
      <div class="chart-container"><canvas id="chart-strength"></canvas></div>
    </div>
  </div>

  <div class="card full">
    <h3><span class="icon">📋</span> Full Comparison Table</h3>
    <div style="overflow-x:auto">
      <table class="fiber-table" id="fiber-table"></table>
    </div>
  </div>

  <div class="card full">
    <h3><span class="icon">🧼</span> Wash Fastness — L* After 50 Washes</h3>
    <div class="chart-container"><canvas id="chart-wash"></canvas></div>
  </div>
</div>

<!-- ═══════════ T-DNA CONSTRUCT ═══════════ -->
<div id="sec-construct" class="section">
  <h2 class="section-title"><span>🧬</span> T-DNA Construct Map</h2>
  <p class="section-sub">pBC-MelaninCotton-v2 — Complete genetic construct for Agrobacterium-mediated cotton transformation.</p>

  <div class="stat-grid">
    <div class="stat"><div class="value">8,437</div><div class="label">Total Length (bp)</div></div>
    <div class="stat"><div class="value">4</div><div class="label">Expression Cassettes</div></div>
    <div class="stat"><div class="value">46.6%</div><div class="label">GC Content</div></div>
    <div class="stat"><div class="value">pCAMBIA2300</div><div class="label">Backbone</div></div>
  </div>

  <div class="card full">
    <h3><span class="icon">🗺️</span> Linear Construct Map</h3>
    <div id="construct-map" style="overflow-x:auto; padding: 16px 0;"></div>
  </div>

  <div class="card full">
    <h3><span class="icon">📐</span> Cassette Details</h3>
    <table class="fiber-table" id="cassette-table"></table>
  </div>
</div>

</div>

<div style="text-align:center; padding: 40px; color: var(--text-dim); font-size: 0.8rem;">
  BlackCotton Scientific Dashboard · Generated {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M')} · Computational Design Platform
</div>

<script>
// ── Data ──
const exprData = {json.dumps(expr_data) if expr_data else 'null'};
const melData = {json.dumps(melanin_data) if melanin_data else 'null'};
const fiberData = {json.dumps(fiber_data) if fiber_data else '[]'};
const constructData = {json.dumps(construct_data) if construct_data else 'null'};
const washData = {json.dumps(wash_data) if wash_data else 'null'};

// ── Navigation ──
function showSection(name) {{
  document.querySelectorAll('.section').forEach(s => s.classList.remove('active'));
  document.querySelectorAll('nav button').forEach(b => b.classList.remove('active'));
  document.getElementById('sec-' + name).classList.add('active');
  event.target.classList.add('active');

  if (name === 'proteins' && !window._proteinsLoaded) {{
    loadProteins();
    window._proteinsLoaded = true;
  }}
}}

// ── Chart.js defaults ──
Chart.defaults.color = '#555555';
Chart.defaults.borderColor = 'rgba(0,0,0,0.1)';
Chart.defaults.font.family = "'Inter', sans-serif";
Chart.defaults.font.size = 11;
Chart.defaults.plugins.legend.labels.usePointStyle = true;
Chart.defaults.plugins.legend.labels.padding = 14;
Chart.defaults.elements.point.radius = 0;
Chart.defaults.elements.point.hoverRadius = 4;
Chart.defaults.elements.line.tension = 0.4;
Chart.defaults.elements.line.borderWidth = 2.2;

// ── Expression Charts ──
if (exprData) {{
  const maxMelA = Math.max(...exprData.protein_melA);
  const maxTYRP1 = Math.max(...exprData.protein_TYRP1);
  const maxDCT = Math.max(...exprData.protein_DCT);
  const maxCel = Math.max(...exprData.cellulose);

  new Chart(document.getElementById('chart-expression'), {{
    type: 'line',
    data: {{
      labels: exprData.t_dpa.map(v => v.toFixed(1)),
      datasets: [
        {{ label: 'Cellulose (norm)', data: exprData.cellulose.map(v => v/maxCel), borderColor: '#888', backgroundColor: 'rgba(0,0,0,0.05)', fill: true }},
        {{ label: 'melA protein (norm)', data: exprData.protein_melA.map(v => v/maxMelA), borderColor: '#222' }},
        {{ label: 'TYRP1 protein (norm)', data: exprData.protein_TYRP1.map(v => v/maxTYRP1), borderColor: '#666' }},
        {{ label: 'DCT protein (norm)', data: exprData.protein_DCT.map(v => v/maxDCT), borderColor: '#999' }},
      ]
    }},
    options: {{
      responsive: true,
      maintainAspectRatio: false,
      scales: {{
        x: {{ title: {{ display: true, text: 'Days Post Anthesis (DPA)' }}, ticks: {{ maxTicksLimit: 12 }} }},
        y: {{ title: {{ display: true, text: 'Normalized Level' }}, min: 0, max: 1.05 }}
      }},
      plugins: {{ tooltip: {{ mode: 'index', intersect: false }} }}
    }}
  }});

  new Chart(document.getElementById('chart-promoters'), {{
    type: 'line',
    data: {{
      labels: exprData.t_dpa.map(v => v.toFixed(1)),
      datasets: [
        {{ label: 'pGhMat1 (melA)', data: exprData.promoter_melA, borderColor: '#333', borderDash: [6,3] }},
        {{ label: 'pGhSCW-late (TYRP1/DCT)', data: exprData.promoter_SCW_late, borderColor: '#888', borderDash: [6,3] }},
      ]
    }},
    options: {{
      responsive: true,
      maintainAspectRatio: false,
      scales: {{
        x: {{ title: {{ display: true, text: 'Days Post Anthesis (DPA)' }}, ticks: {{ maxTicksLimit: 12 }} }},
        y: {{ title: {{ display: true, text: 'Promoter Activity (relative)' }}, min: 0 }}
      }}
    }}
  }});
}}

// ── Melanin Charts ──
if (melData) {{
  new Chart(document.getElementById('chart-melanin'), {{
    type: 'line',
    data: {{
      labels: melData.t_dpa.map(v => v.toFixed(1)),
      datasets: [
        {{ label: 'L-Tyrosine', data: melData.tyrosine, borderColor: '#888' }},
        {{ label: 'L-DOPA', data: melData.L_DOPA, borderColor: '#555' }},
        {{ label: 'Melanin', data: melData.melanin, borderColor: '#111', borderWidth: 3 }},
      ]
    }},
    options: {{
      responsive: true,
      maintainAspectRatio: false,
      scales: {{
        x: {{ title: {{ display: true, text: 'Days Post Anthesis' }}, ticks: {{ maxTicksLimit: 12 }} }},
        y: {{ title: {{ display: true, text: 'Concentration (mM)' }}, min: 0 }}
      }}
    }}
  }});

  new Chart(document.getElementById('chart-intermediates'), {{
    type: 'line',
    data: {{
      labels: melData.t_dpa.map(v => v.toFixed(1)),
      datasets: [
        {{ label: 'Dopaquinone', data: melData.dopaquinone, borderColor: '#333' }},
        {{ label: 'Dopachrome', data: melData.dopachrome, borderColor: '#666' }},
        {{ label: 'DHICA', data: melData.DHICA, borderColor: '#999' }},
        {{ label: 'Indole-quinone', data: melData.indole_quinone, borderColor: '#aaa' }},
        {{ label: 'L-DOPA', data: melData.L_DOPA, borderColor: '#ccc' }},
      ]
    }},
    options: {{
      responsive: true,
      maintainAspectRatio: false,
      scales: {{
        x: {{ title: {{ display: true, text: 'Days Post Anthesis' }}, ticks: {{ maxTicksLimit: 12 }} }},
        y: {{ title: {{ display: true, text: 'Concentration (mM)' }}, min: 0 }}
      }}
    }}
  }});
}}

// ── Fiber Charts ──
if (fiberData.length) {{
  const names = fiberData.map(f => f.name.replace(' (Coker-312)', '').replace('Engineered Black ','Eng '));
  const colors = ['#dfe6e9','#b2bec3','#2d3436','#6c5ce7','#a29bfe','#5f27cd'];

  new Chart(document.getElementById('chart-color'), {{
    type: 'bar',
    data: {{
      labels: names,
      datasets: [{{ label: 'L* (lightness)', data: fiberData.map(f => f.color_L),
        backgroundColor: fiberData.map((f,i) => {{
          const L = f.color_L;
          const gray = Math.round(L * 2.55);
          return `rgb(${{gray}},${{gray}},${{gray}})`;
        }}),
        borderColor: '#111', borderWidth: 1, borderRadius: 2 }}]
    }},
    options: {{
      responsive: true, maintainAspectRatio: false,
      indexAxis: 'y',
      scales: {{ x: {{ min: 0, max: 85, title: {{ display: true, text: 'L* (0=black, 100=white)' }} }} }},
      plugins: {{ legend: {{ display: false }} }}
    }}
  }});

  new Chart(document.getElementById('chart-strength'), {{
    type: 'bar',
    data: {{
      labels: names,
      datasets: [{{ label: 'Strength (g/tex)', data: fiberData.map(f => f.strength_g_tex),
        backgroundColor: '#444',
        borderColor: '#111', borderWidth: 1, borderRadius: 2 }}]
    }},
    options: {{
      responsive: true, maintainAspectRatio: false,
      indexAxis: 'y',
      scales: {{ x: {{ min: 0, max: 35, title: {{ display: true, text: 'Strength (g/tex)' }} }} }},
      plugins: {{ legend: {{ display: false }} }}
    }}
  }});

  // Fiber table
  const tbl = document.getElementById('fiber-table');
  let html = '<thead><tr><th>Color</th><th>Type</th><th>UHML (mm)</th><th>Strength</th><th>Mic</th><th>L*</th><th>Cellulose %</th><th>Grade</th></tr></thead><tbody>';
  fiberData.forEach(f => {{
    const gray = Math.round(f.color_L * 2.55);
    const gc = f.grade === 'PREMIUM' ? 'premium' : f.grade === 'GOOD' ? 'good' : f.grade === 'BASE' ? 'base' : 'below';
    html += `<tr><td><span class="color-swatch" style="background:rgb(${{gray}},${{gray}},${{gray}})"></span></td>
      <td class="name">${{f.name}}</td><td>${{f.uhml_mm.toFixed(1)}}</td><td>${{f.strength_g_tex.toFixed(1)}}</td>
      <td>${{f.micronaire.toFixed(1)}}</td><td>${{f.color_L.toFixed(1)}}</td><td>${{f.cellulose_pct.toFixed(1)}}</td>
      <td><span class="grade ${{gc}}">${{f.grade}}</span></td></tr>`;
  }});
  html += '</tbody>';
  tbl.innerHTML = html;
}}

// ── Wash Fastness ──
if (washData) {{
  const keys = Object.keys(washData);
  const datasets = keys.map((k, i) => ({{
    label: k.replace(' (Coker-312)','').replace('Engineered Black ','Eng '),
    data: washData[k].map(w => w.L_star),
    borderWidth: i >= 3 ? 2.5 : 1.5,
    borderDash: i < 3 ? [4,4] : [],
    borderColor: '#333'
  }}));
  new Chart(document.getElementById('chart-wash'), {{
    type: 'line',
    data: {{ labels: washData[keys[0]].map(w => w.wash), datasets }},
    options: {{
      responsive: true, maintainAspectRatio: false,
      scales: {{
        x: {{ title: {{ display: true, text: 'Number of Washes' }}, ticks: {{ maxTicksLimit: 10 }} }},
        y: {{ title: {{ display: true, text: 'L* (lightness)' }}, min: 0, max: 85 }}
      }}
    }}
  }});
}}

// ── Construct Map ──
if (constructData) {{
  const container = document.getElementById('construct-map');
  const features = constructData.features;
  const total = constructData.total_length_bp;
  const typeColors = {{
    'T-DNA_border': '#dfdfdf', 'promoter': '#f0f0f0',
    'CDS': '#d0d0d0', 'terminator': '#f9f9f9'
  }};
  let html = '<div class="construct-map">';
  features.forEach((f, i) => {{
    const len = f.end - f.start;
    const width = Math.max(30, (len / total) * 900);
    const isGene = f.type === 'CDS';
    const textColor = isGene ? '#0a0a0f' : 'white';
    html += `<div class="construct-element">
      <div class="bar" style="background:${{typeColors[f.type]}};width:${{width}}px;color:${{textColor}}">${{f.name}}</div>
      <div class="size">${{len}} bp</div>
    </div>`;
    if (i < features.length - 1) html += '<div class="arrow"></div>';
  }});
  html += '</div>';
  container.innerHTML = html;

  // Cassette table
  const ct = document.getElementById('cassette-table');
  let chtml = '<thead><tr><th>#</th><th>Cassette</th><th>Promoter</th><th>Gene</th><th>vTP</th><th>Terminator</th><th>Total (bp)</th></tr></thead><tbody>';
  constructData.cassettes.forEach((c, i) => {{
    chtml += `<tr><td>${{i+1}}</td><td class="name">${{c.name}}</td>
      <td>${{c.promoter.name}} (${{c.promoter.length}})</td>
      <td>${{c.gene.name}} (${{c.gene.length}})</td>
      <td>${{c.gene.transit_peptide_fused ? '✅' : '—'}}</td>
      <td>${{c.terminator.name}} (${{c.terminator.length}})</td>
      <td>${{c.length_bp}}</td></tr>`;
  }});
  chtml += '</tbody>';
  ct.innerHTML = chtml;
}}

// ── 3D Protein Viewers ──
function loadProteins() {{
  const configs = [
    {{ id: 'viewer-melA', pdb: '1WX2', color: 'spectrum', style: 'cartoon' }},
    {{ id: 'viewer-tyrp1', pdb: '5M8Q', color: 'spectrum', style: 'cartoon' }},
    {{ id: 'viewer-dct', pdb: '3NQ1', color: 'spectrum', style: 'cartoon' }},
  ];
  configs.forEach(cfg => {{
    const el = document.getElementById(cfg.id);
    const viewer = $3Dmol.createViewer(el, {{
      backgroundColor: '0xf7f7f7',
      antialias: true,
    }});
    $3Dmol.download('pdb:' + cfg.pdb, viewer, {{}}, function() {{
      viewer.setStyle({{}}, {{
        cartoon: {{
          color: 'silver',
          opacity: 0.92,
          thickness: 0.4,
        }}
      }});
      // Highlight metal centers (copper for tyrosinase, zinc for TYRP1)
      viewer.setStyle({{elem: 'CU'}}, {{sphere: {{radius: 0.8, color: '#333'}}}});
      viewer.setStyle({{elem: 'ZN'}}, {{sphere: {{radius: 0.8, color: '#666'}}}});
      viewer.setStyle({{elem: 'FE'}}, {{sphere: {{radius: 0.8, color: '#111'}}}});
      viewer.zoomTo();
      viewer.spin('y', 0.5);
      viewer.render();
    }});
  }});
}}

// Load proteins on first view
window._proteinsLoaded = false;
setTimeout(() => {{
  if (!window._proteinsLoaded) {{
    loadProteins();
    window._proteinsLoaded = true;
  }}
}}, 500);
</script>
</body>
</html>"""


if __name__ == "__main__":
    print("\n🧬 BlackCotton Dashboard Generator")
    print("=" * 50)

    print("  Loading simulation data...")
    expr = load_expression_data()
    melanin = load_melanin_data()
    fiber = load_fiber_data()
    construct = load_construct_data()
    wash = load_wash_data()

    status = []
    if expr:
        status.append("✓ Expression dynamics")
    if melanin:
        status.append("✓ Melanin pathway")
    if fiber:
        status.append("✓ Fiber quality")
    if construct:
        status.append("✓ Construct map")
    if wash:
        status.append("✓ Wash fastness")
    for s in status:
        print(f"  {s}")

    print("\n  Generating HTML dashboard...")
    html = generate_html(expr, melanin, fiber, construct, wash)

    out_path = RESULTS_DIR / "blackcotton_dashboard.html"
    with open(out_path, "w") as f:
        f.write(html)
    print(f"  📊 Saved: {out_path}")

    print("\n  Opening in browser...")
    webbrowser.open(f"file://{out_path.resolve()}")
    print("\n✅ Dashboard ready!")
