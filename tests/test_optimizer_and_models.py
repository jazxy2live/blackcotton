import unittest
from unittest.mock import patch

from src.fiber_model import load_params as load_fiber_params
from src.fiber_model import model_engineered_black
from src.construct_designer import resolve_melA_cds
from src.tradeoff_optimizer import (
    effective_pigment_levels,
    load_params as load_optimizer_params,
    pareto_front,
    run_sweep,
    select_top_candidates,
)
from src.robustness_analyzer import DEFAULT_THRESHOLDS, summarize_trials
from src.robust_tradeoff_optimizer import select_robust_top
from src.sensitivity_analyzer import rank_feature_impacts
from src.adaptive_robust_optimizer import targeted_parameter_grid
from src.adaptive_ode_robust_pipeline import normalize_candidate_row, select_final_top_candidates


class FiberModelTests(unittest.TestCase):
    def setUp(self):
        self.params = load_fiber_params()

    def test_melanin_load_decouples_color_from_structural_penalty(self):
        low_load = model_engineered_black(
            self.params,
            melanin=0.8,
            overlap=0.2,
            k_competition=0.2,
            melanin_load=0.4,
        )
        high_load = model_engineered_black(
            self.params,
            melanin=0.8,
            overlap=0.2,
            k_competition=0.2,
            melanin_load=0.9,
        )

        self.assertAlmostEqual(low_load.color_L, high_load.color_L, places=9)
        self.assertGreater(low_load.strength_g_tex, high_load.strength_g_tex)
        self.assertGreater(low_load.cellulose_pct, high_load.cellulose_pct)

    def test_color_tracks_visible_melanin(self):
        lighter = model_engineered_black(
            self.params,
            melanin=0.45,
            overlap=0.1,
            k_competition=0.1,
            melanin_load=0.45,
        )
        darker = model_engineered_black(
            self.params,
            melanin=0.75,
            overlap=0.1,
            k_competition=0.1,
            melanin_load=0.45,
        )
        self.assertLess(darker.color_L, lighter.color_L)


class OptimizerTests(unittest.TestCase):
    def test_effective_pigment_levels_is_bounded_and_retention_sensitive(self):
        load_a, visible_a = effective_pigment_levels(0.70, 0.80, 0.45)
        load_b, visible_b = effective_pigment_levels(0.70, 0.20, 0.45)
        high_load, high_visible = effective_pigment_levels(1.50, 1.00, 1.00)

        self.assertLessEqual(load_a, 0.90)
        self.assertLessEqual(visible_a, 0.98)
        self.assertGreater(visible_a, load_a)
        self.assertGreater(visible_a, visible_b)
        self.assertLessEqual(high_load, 0.90)
        self.assertLessEqual(high_visible, 0.98)

    @patch("src.tradeoff_optimizer.evaluate_candidate_proxy")
    def test_run_sweep_enforces_non_negative_temporal_gap(self, mock_eval):
        def fake_eval(**kwargs):
            # Negative below 36 DPA, non-negative at/after 36 DPA.
            gap = kwargs["mat_activation_dpa"] - 36.0
            return {"temporal_gap_days": gap}

        mock_eval.side_effect = fake_eval
        params = load_optimizer_params()
        candidates = run_sweep(params)

        self.assertTrue(candidates)
        self.assertTrue(all(c["temporal_gap_days"] >= 0.0 for c in candidates))
        # 3 mat activation values (36, 38, 40) survive from 5 total.
        expected = 3 * 5 * 4 * 4 * 4 * 4 * 4
        self.assertEqual(len(candidates), expected)

    def test_pareto_front_deduplicates_and_filters_dominated(self):
        candidates = [
            {"color_L": 20.0, "strength_g_tex": 29.0, "yield_index": 0.90, "temporal_gap_days": 1.0},
            {"color_L": 18.0, "strength_g_tex": 29.0, "yield_index": 0.90, "temporal_gap_days": 1.0},  # dominates first
            {"color_L": 18.0, "strength_g_tex": 29.0, "yield_index": 0.90, "temporal_gap_days": 1.0},  # duplicate
            {"color_L": 19.0, "strength_g_tex": 30.0, "yield_index": 0.88, "temporal_gap_days": 1.2},  # tradeoff
        ]
        front = pareto_front(candidates)

        self.assertEqual(len(front), 2)
        self.assertTrue(any(c["color_L"] == 18.0 and c["strength_g_tex"] == 29.0 for c in front))
        self.assertTrue(any(c["color_L"] == 19.0 and c["strength_g_tex"] == 30.0 for c in front))

    def test_select_top_candidates_removes_objective_duplicates(self):
        candidates = [
            {
                "color_L": 20.0,
                "strength_g_tex": 29.0,
                "yield_index": 0.90,
                "temporal_gap_days": 1.0,
                "composite_score": 0.90,
            },
            {
                "color_L": 20.0,
                "strength_g_tex": 29.0,
                "yield_index": 0.90,
                "temporal_gap_days": 1.0,
                "composite_score": 0.89,
            },
            {
                "color_L": 21.0,
                "strength_g_tex": 28.8,
                "yield_index": 0.89,
                "temporal_gap_days": 1.5,
                "composite_score": 0.88,
            },
        ]
        selected = select_top_candidates(candidates, n=10)

        self.assertEqual(len(selected), 2)


class SequenceResolverTests(unittest.TestCase):
    def test_resolve_mela_returns_valid_cds_structure(self):
        seq, source = resolve_melA_cds()
        self.assertTrue(source)
        self.assertGreaterEqual(len(seq), 750)
        self.assertEqual(len(seq) % 3, 0)
        self.assertIn(seq[:3], {"ATG", "GTG", "TTG"})
        self.assertIn(seq[-3:], {"TAA", "TAG", "TGA"})


class RobustnessTests(unittest.TestCase):
    def test_summarize_trials_tracks_success_and_risks(self):
        trials = [
            {"color_L": 20.0, "strength_g_tex": 29.0, "yield_index": 0.90, "temporal_gap_days": 1.0, "composite_score": 0.85},
            {"color_L": 23.0, "strength_g_tex": 28.4, "yield_index": 0.88, "temporal_gap_days": 0.2, "composite_score": 0.82},
            {"color_L": 27.0, "strength_g_tex": 27.0, "yield_index": 0.80, "temporal_gap_days": -0.3, "composite_score": 0.60},
        ]
        summary = summarize_trials(trials, DEFAULT_THRESHOLDS)

        self.assertEqual(summary["n_trials"], 3)
        self.assertAlmostEqual(summary["success_rate"], 2.0 / 3.0, places=6)
        self.assertAlmostEqual(summary["risk_temporal_overlap"], 1.0 / 3.0, places=6)
        self.assertAlmostEqual(summary["risk_strength_failure"], 1.0 / 3.0, places=6)
        self.assertAlmostEqual(summary["risk_darkness_failure"], 1.0 / 3.0, places=6)
        self.assertAlmostEqual(summary["risk_yield_failure"], 1.0 / 3.0, places=6)
        self.assertGreater(summary["robust_score"], 0.0)

    def test_select_robust_top_prefers_min_success_pool(self):
        rows = [
            {
                "success_rate": 0.03,
                "robust_score": 0.40,
                "p50_color_L": 24.0,
                "p50_strength_g_tex": 28.5,
                "p50_yield_index": 0.88,
                "risk_temporal_overlap": 0.1,
                "risk_darkness_failure": 0.3,
                "base_candidate": {
                    "color_L": 22.0, "strength_g_tex": 28.6, "yield_index": 0.89,
                    "temporal_gap_days": 1.0, "mat_activation_dpa": 38.0,
                    "scw_activation_dpa": 32.0, "k_competition": 0.1,
                    "melanin_efficiency": 1.2, "late_retention_factor": 0.3
                },
            },
            {
                "success_rate": 0.12,
                "robust_score": 0.35,
                "p50_color_L": 23.0,
                "p50_strength_g_tex": 28.9,
                "p50_yield_index": 0.90,
                "risk_temporal_overlap": 0.0,
                "risk_darkness_failure": 0.2,
                "base_candidate": {
                    "color_L": 21.0, "strength_g_tex": 29.0, "yield_index": 0.90,
                    "temporal_gap_days": 1.2, "mat_activation_dpa": 40.0,
                    "scw_activation_dpa": 34.0, "k_competition": 0.1,
                    "melanin_efficiency": 1.3, "late_retention_factor": 0.45
                },
            },
        ]
        top = select_robust_top(rows, n=2, min_success_rate=0.08)
        self.assertEqual(len(top), 1)
        self.assertAlmostEqual(top[0]["success_rate"], 0.12, places=9)


class SensitivityTests(unittest.TestCase):
    def test_rank_feature_impacts_prioritizes_failure_driver(self):
        records = []
        for i in range(40):
            high_comp = 0.8 if i < 20 else 0.1
            color = 31.0 if high_comp > 0.5 else 22.0
            records.append(
                {
                    "candidate_input_rank": 1,
                    "sampled_mat_activation_dpa": 38.0,
                    "sampled_scw_activation_dpa": 32.0,
                    "sampled_mat_strength": 1.0,
                    "sampled_scw_strength": 1.0,
                    "sampled_k_competition": high_comp,
                    "sampled_melanin_efficiency": 1.2,
                    "sampled_late_retention_factor": 0.4,
                    "sampled_hill_melA": 8.0,
                    "sampled_hill_scw": 6.0,
                    "sampled_leak_melA": 0.01,
                    "sampled_leak_scw": 0.02,
                    "color_L": color,
                    "strength_g_tex": 28.8,
                    "yield_index": 0.88,
                    "temporal_gap_days": 1.0,
                    "composite_score": 0.8,
                    "success": color <= 25.0,
                }
            )

        ranked = rank_feature_impacts(records, DEFAULT_THRESHOLDS)
        self.assertTrue(ranked)
        self.assertEqual(ranked[0]["feature"], "sampled_k_competition")
        self.assertGreater(ranked[0]["impact_score"], 0.0)


class AdaptiveOptimizerTests(unittest.TestCase):
    def test_targeted_parameter_grid_has_darkness_focused_ranges(self):
        grid = targeted_parameter_grid()
        self.assertIn(1.60, grid["melanin_efficiency"])
        self.assertIn(0.75, grid["late_retention_factor"])
        self.assertIn(42.0, grid["mat_activation_dpa"])
        self.assertIn(30.0, grid["scw_activation_dpa"])


class AdaptiveODERobustTests(unittest.TestCase):
    def test_normalize_candidate_row_accepts_base_candidate(self):
        row = {
            "base_candidate": {
                "mat_activation_dpa": 38.0,
                "scw_activation_dpa": 34.0,
                "mat_strength": 1.1,
                "scw_strength": 1.0,
                "k_competition": 0.1,
                "melanin_efficiency": 1.6,
                "late_retention_factor": 0.75,
                "color_L": 18.0,
                "strength_g_tex": 29.0,
                "yield_index": 0.90,
                "temporal_gap_days": 1.5,
            },
            "robust_score": 0.9,
        }
        out = normalize_candidate_row(row)
        self.assertEqual(out["k_competition"], 0.1)
        self.assertEqual(out["melanin_efficiency"], 1.6)
        self.assertEqual(out["late_retention_factor"], 0.75)
        self.assertAlmostEqual(out["composite_score"], 0.9)

    def test_select_final_top_candidates_prefers_success_threshold(self):
        rows = [
            {
                "success_rate": 0.10,
                "robust_score": 0.95,
                "base_candidate": {"color_L": 18.0, "strength_g_tex": 29.0, "yield_index": 0.9, "temporal_gap_days": 1.0},
            },
            {
                "success_rate": 0.30,
                "robust_score": 0.80,
                "base_candidate": {"color_L": 19.0, "strength_g_tex": 29.1, "yield_index": 0.9, "temporal_gap_days": 1.0},
            },
        ]
        top = select_final_top_candidates(rows, n=3, min_success_rate=0.20)
        self.assertEqual(len(top), 1)
        self.assertAlmostEqual(top[0]["success_rate"], 0.30)


if __name__ == "__main__":
    unittest.main()
