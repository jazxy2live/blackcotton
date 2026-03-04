import unittest

from src.lab_validation_planner import build_candidate_rows, evaluate_gate


class LabValidationPlannerTests(unittest.TestCase):
    def test_build_candidate_rows_maps_top3_fields(self):
        final_top3 = [
            {
                "input_rank": 2,
                "success_rate": 0.99,
                "robust_score": 0.95,
                "p50_color_L": 18.2,
                "p50_strength_g_tex": 29.0,
                "p50_yield_index": 0.90,
                "p50_temporal_gap_days": 1.7,
                "base_candidate": {
                    "mat_activation_dpa": 38.0,
                    "scw_activation_dpa": 34.0,
                    "mat_strength": 1.3,
                    "scw_strength": 1.3,
                    "k_competition": 0.05,
                    "melanin_efficiency": 1.6,
                    "late_retention_factor": 0.75,
                },
            }
        ]
        rows = build_candidate_rows(final_top3=final_top3, replicates=8, run_id="BC-LAB-X")
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["arm_id"], "CAND_1")
        self.assertEqual(row["source_input_rank"], 2)
        self.assertAlmostEqual(row["mat_activation_dpa"], 38.0, places=9)
        self.assertAlmostEqual(row["pred_success_rate"], 0.99, places=9)

    def test_evaluate_gate_true_and_false(self):
        thresholds = {
            "max_color_L": 25.0,
            "min_strength_g_tex": 28.0,
            "min_yield_index": 0.85,
            "min_temporal_gap_days": 0.0,
        }
        good = {
            "color_L": 20.0,
            "strength_g_tex": 29.0,
            "yield_index": 0.90,
            "temporal_gap_days": 1.0,
        }
        bad = {
            "color_L": 30.0,
            "strength_g_tex": 27.5,
            "yield_index": 0.82,
            "temporal_gap_days": -0.2,
        }
        self.assertTrue(evaluate_gate(good, thresholds))
        self.assertFalse(evaluate_gate(bad, thresholds))


if __name__ == "__main__":
    unittest.main()
