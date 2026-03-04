import unittest

from src.adversarial_robustness_suite import (
    aggregate_worst_case,
    scale_noise,
    strict_thresholds,
)
from src.robustness_analyzer import DEFAULT_NOISE, DEFAULT_THRESHOLDS


class AdversarialRobustnessSuiteTests(unittest.TestCase):
    def test_scale_noise_multiplies_and_caps_retention_sigma(self):
        noise = dict(DEFAULT_NOISE)
        out = scale_noise(noise, 2.5)
        self.assertGreater(out["timing_sigma_dpa"], noise["timing_sigma_dpa"])
        self.assertLessEqual(out["retention_sigma"], 0.20)

    def test_strict_thresholds_tighten_constraints(self):
        strict = strict_thresholds(DEFAULT_THRESHOLDS)
        self.assertLess(strict["max_color_L"], DEFAULT_THRESHOLDS["max_color_L"])
        self.assertGreater(strict["min_strength_g_tex"], DEFAULT_THRESHOLDS["min_strength_g_tex"])
        self.assertGreater(strict["min_yield_index"], DEFAULT_THRESHOLDS["min_yield_index"])
        self.assertGreater(strict["min_temporal_gap_days"], DEFAULT_THRESHOLDS["min_temporal_gap_days"])

    def test_aggregate_worst_case_orders_by_min_success(self):
        rows = [
            {"scenario": "a", "candidate_input_rank": 1, "success_rate": 0.9, "robust_score": 0.8},
            {"scenario": "b", "candidate_input_rank": 1, "success_rate": 0.7, "robust_score": 0.6},
            {"scenario": "a", "candidate_input_rank": 2, "success_rate": 0.8, "robust_score": 0.7},
            {"scenario": "b", "candidate_input_rank": 2, "success_rate": 0.75, "robust_score": 0.68},
        ]
        out = aggregate_worst_case(rows)
        self.assertEqual(out[0]["candidate_input_rank"], 2)
        self.assertEqual(out[1]["candidate_input_rank"], 1)
        self.assertEqual(out[1]["worst_scenario"], "b")


if __name__ == "__main__":
    unittest.main()
