import tempfile
import unittest
from pathlib import Path

from src.proof_packager import parse_test_report, select_lead_candidate


class ProofPackagerTests(unittest.TestCase):
    def test_select_lead_candidate_prefers_final_lab_top3(self):
        final_lab_top3 = [
            {
                "success_rate": 0.95,
                "robust_score": 0.92,
                "p50_color_L": 18.5,
                "p50_strength_g_tex": 29.0,
                "p50_yield_index": 0.90,
                "p50_temporal_gap_days": 1.8,
                "base_candidate": {
                    "color_L": 13.0,
                    "strength_g_tex": 29.2,
                    "yield_index": 0.91,
                    "temporal_gap_days": 2.5,
                    "mat_activation_dpa": 40.0,
                    "scw_activation_dpa": 34.0,
                    "mat_strength": 1.1,
                    "scw_strength": 1.1,
                    "k_competition": 0.05,
                    "melanin_efficiency": 1.6,
                    "late_retention_factor": 0.75,
                },
            }
        ]
        top_candidates = [
            {
                "color_L": 20.0,
                "strength_g_tex": 29.0,
                "yield_index": 0.89,
                "temporal_gap_days": 1.4,
                "composite_score": 0.88,
            }
        ]

        lead = select_lead_candidate(final_lab_top3=final_lab_top3, top_candidates=top_candidates)
        self.assertEqual(lead["source"], "final_lab_top3")
        self.assertEqual(lead["model_stage"], "adaptive_ode_robust")
        self.assertAlmostEqual(lead["success_rate"], 0.95, places=9)
        self.assertAlmostEqual(lead["color_L"], 13.0, places=9)

    def test_select_lead_candidate_falls_back_to_top_candidates(self):
        lead = select_lead_candidate(
            final_lab_top3=[],
            top_candidates=[
                {
                    "color_L": 19.5,
                    "strength_g_tex": 29.1,
                    "yield_index": 0.90,
                    "temporal_gap_days": 2.0,
                    "mat_activation_dpa": 40.0,
                    "scw_activation_dpa": 32.0,
                    "mat_strength": 1.0,
                    "scw_strength": 1.0,
                    "k_competition": 0.1,
                    "melanin_efficiency": 1.3,
                    "late_retention_factor": 0.45,
                    "composite_score": 0.89,
                }
            ],
        )
        self.assertEqual(lead["source"], "top_candidates")
        self.assertEqual(lead["model_stage"], "ode_refined")
        self.assertIsNone(lead["success_rate"])
        self.assertAlmostEqual(lead["robust_score"], 0.89, places=9)

    def test_parse_test_report_reads_ran_and_status(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "test_report.txt"
            path.write_text(
                "test_a ... ok\n"
                "test_b ... ok\n\n"
                "----------------------------------------------------------------------\n"
                "Ran 13 tests in 0.321s\n\n"
                "OK\n"
            )
            parsed = parse_test_report(path)
            self.assertEqual(parsed["status"], "pass")
            self.assertEqual(parsed["ran_tests"], 13)


if __name__ == "__main__":
    unittest.main()
