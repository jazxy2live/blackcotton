import unittest

import numpy as np

from src.experiment_power_analyzer import (
    recommend_timepoints,
    required_n_for_target_power,
)


class ExperimentPowerAnalyzerTests(unittest.TestCase):
    def test_required_n_lower_for_stronger_effect(self):
        strong = np.array([31.0, 30.8, 31.2, 30.9, 31.1, 30.7, 31.0, 31.3], dtype=float)
        weak = np.array([28.1, 28.0, 28.2, 27.9, 28.1, 28.0, 28.2, 28.0], dtype=float)

        n_strong, _ = required_n_for_target_power(
            values=strong,
            threshold=28.0,
            direction="higher",
            target_power=0.80,
            max_replicates=16,
            n_sims=1500,
            alpha=0.05,
            measurement_sd=0.0,
            seed=42,
        )
        n_weak, _ = required_n_for_target_power(
            values=weak,
            threshold=28.0,
            direction="higher",
            target_power=0.80,
            max_replicates=16,
            n_sims=1500,
            alpha=0.05,
            measurement_sd=0.0,
            seed=42,
        )

        self.assertIsNotNone(n_strong)
        self.assertIsNotNone(n_weak)
        self.assertLessEqual(n_strong, n_weak)

    def test_recommend_timepoints_covers_switch_windows(self):
        top3 = [
            {"base_candidate": {"scw_activation_dpa": 34.0, "mat_activation_dpa": 38.0}},
            {"base_candidate": {"scw_activation_dpa": 34.0, "mat_activation_dpa": 40.0}},
            {"base_candidate": {"scw_activation_dpa": 32.0, "mat_activation_dpa": 38.0}},
        ]
        points = recommend_timepoints(top3)
        self.assertEqual(points, sorted(points))
        self.assertIn(32, points)
        self.assertIn(34, points)
        self.assertIn(38, points)
        self.assertIn(40, points)
        self.assertIn(42, points)


if __name__ == "__main__":
    unittest.main()
