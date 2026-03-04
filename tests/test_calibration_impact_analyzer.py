import unittest

from src.calibration_impact_analyzer import delta, pick_winner


class CalibrationImpactAnalyzerTests(unittest.TestCase):
    def test_delta(self):
        self.assertAlmostEqual(delta(3.5, 1.0), 2.5, places=9)
        self.assertAlmostEqual(delta(1.0, 3.5), -2.5, places=9)

    def test_pick_winner_prefers_higher_worst_case(self):
        baseline = {"robust_worst_case_best": {"success_rate": 0.70, "robust_score": 0.76}}
        calibrated = {"robust_worst_case_best": {"success_rate": 0.74, "robust_score": 0.74}}
        out = pick_winner(baseline=baseline, calibrated=calibrated)
        self.assertEqual(out["winner"], "calibrated")

    def test_pick_winner_tie(self):
        baseline = {"robust_worst_case_best": {"success_rate": 0.70, "robust_score": 0.76}}
        calibrated = {"robust_worst_case_best": {"success_rate": 0.70, "robust_score": 0.76}}
        out = pick_winner(baseline=baseline, calibrated=calibrated)
        self.assertEqual(out["winner"], "tie")


if __name__ == "__main__":
    unittest.main()
