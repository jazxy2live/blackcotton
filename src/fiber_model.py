#!/usr/bin/env python3
"""
fiber_model.py — Cotton Fiber Quality Prediction Model
========================================================
Predicts metrics and compares engineered black cotton against conventional.
"""

import numpy as np
import json
from pathlib import Path
from dataclasses import dataclass

from src.config_loader import load_config

BASE_DIR = Path(__file__).resolve().parent.parent
CONFIG_DIR = BASE_DIR / "config"
RESULTS_DIR = BASE_DIR / "results"


def load_params():
    return load_config()


@dataclass
class FiberQuality:
    name: str
    uhml_mm: float
    strength_g_tex: float
    micronaire: float
    uniformity_pct: float
    color_L: float
    color_a: float
    color_b: float
    cellulose_pct: float

    def grade(self):
        if self.uhml_mm >= 32 and self.strength_g_tex >= 30:
            return "PREMIUM"
        elif self.uhml_mm >= 28 and self.strength_g_tex >= 26:
            return "GOOD"
        elif self.uhml_mm >= 25 and self.strength_g_tex >= 23:
            return "BASE"
        return "BELOW-GRADE"

    def to_dict(self):
        return {
            'name': self.name, 'uhml_mm': self.uhml_mm,
            'strength_g_tex': self.strength_g_tex, 'micronaire': self.micronaire,
            'uniformity_pct': self.uniformity_pct, 'color_L': self.color_L,
            'color_a': self.color_a, 'color_b': self.color_b,
            'cellulose_pct': self.cellulose_pct, 'grade': self.grade()
        }


def model_white_cotton(p):
    fd = p['fiber_development']
    return FiberQuality("White Cotton (Coker-312)", fd['target_length_mm'],
        fd['target_strength_g_tex'], fd['target_micronaire'],
        fd['target_uniformity'], 82.0, 1.5, 8.0, 95.0)


def model_natural_brown(p):
    fd = p['fiber_development']
    return FiberQuality("Natural Brown Cotton",
        fd['target_length_mm']*0.75, fd['target_strength_g_tex']*0.70,
        3.2, fd['target_uniformity']*0.90, 45.0, 8.0, 18.0, 82.0)


def model_dyed_black(p):
    w = model_white_cotton(p)
    return FiberQuality("Dyed Black Cotton",
        w.uhml_mm*0.98, w.strength_g_tex*0.95, w.micronaire,
        w.uniformity_pct, 12.0, -0.5, -1.0, 93.0)


def model_engineered_black(
    p,
    melanin=1.0,
    overlap=0.0,
    k_competition=0.25,
    melanin_load=None,
):
    fd = p['fiber_development']
    melanin = float(np.clip(melanin, 0.0, 1.0))
    if melanin_load is None:
        melanin_load = melanin
    melanin_load = float(np.clip(melanin_load, 0.0, 1.0))
    overlap = float(np.clip(overlap, 0.0, 1.0))
    k_competition = max(float(k_competition), 0.0)

    # Competition term captures biochemical/resource conflict between
    # cellulose deposition and pigment synthesis when both run together.
    competition = k_competition * melanin_load * overlap
    quality_factor = max(1.0 - 0.20 * overlap - competition, 0.65)
    cel = max(95.0 - 2.0 * melanin_load - 13.0 * overlap - 18.0 * competition, 70.0)
    L = max(82.0 - 72.0 * melanin, 5.0)

    return FiberQuality(
        f"Engineered Black (mel={melanin:.0%}, ovlp={overlap:.0%})",
        fd['target_length_mm'] * quality_factor,
        fd['target_strength_g_tex'] * quality_factor,
        fd['target_micronaire'] * (1 - 0.15 * overlap),
        fd['target_uniformity'] * (1 - 0.10 * overlap),
        L, 1.5-2.0*melanin, 8.0-9.0*melanin, cel)


def wash_fastness(fiber, washes=50):
    data = []
    for w in range(washes+1):
        if "Dyed" in fiber.name:
            L = fiber.color_L + (82.0-fiber.color_L)*(1-0.988**w)
        elif "Engineered" in fiber.name:
            L = max(fiber.color_L*(1-0.001*w), fiber.color_L*0.95)
        else:
            L = fiber.color_L
        data.append({'wash': w, 'L_star': round(L, 2)})
    return data


def run_comparison(p):
    return [model_white_cotton(p), model_natural_brown(p), model_dyed_black(p),
            model_engineered_black(p, 1.0, 0.0), model_engineered_black(p, 0.7, 0.0),
            model_engineered_black(p, 1.0, 0.2)]


def print_comparison(fibers):
    print("\n" + "="*90)
    print("  FIBER QUALITY COMPARISON — USDA HVI METRICS")
    print("="*90)
    print(f"\n  {'Type':<45} {'UHML':>6} {'Str':>6} {'Mic':>5} {'L*':>6} {'Grade':>8}")
    print("  " + "-"*86)
    for f in fibers:
        print(f"  {f.name:<45} {f.uhml_mm:>6.1f} {f.strength_g_tex:>6.1f} "
              f"{f.micronaire:>5.1f} {f.color_L:>6.1f} {f.grade():>8}")
    print("  " + "-"*86)
    eng = [f for f in fibers if "mel=100%" in f.name and "ovlp=0%" in f.name]
    wht = [f for f in fibers if "White" in f.name]
    if eng and wht:
        e, w = eng[0], wht[0]
        print(f"\n  KEY: Engineered black retains {e.strength_g_tex/w.strength_g_tex*100:.0f}% "
              f"strength at L*={e.color_L:.0f} (deep black) — NO FADING")
    print("="*90)


def save_results(fibers, p):
    import os; os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(RESULTS_DIR/"fiber_quality_comparison.json", 'w') as f:
        json.dump([fb.to_dict() for fb in fibers], f, indent=2)
    wash = {fb.name: wash_fastness(fb) for fb in fibers}
    with open(RESULTS_DIR/"wash_fastness_data.json", 'w') as f:
        json.dump(wash, f, indent=2)
    print("  📊 Saved fiber quality & wash fastness data")


if __name__ == "__main__":
    print("\n🧬 BlackCotton Fiber Quality Model")
    print("="*50)
    p = load_params()
    fibers = run_comparison(p)
    print_comparison(fibers)
    print("\nSaving results...")
    save_results(fibers, p)
    print("\n✅ Fiber quality analysis complete!")
