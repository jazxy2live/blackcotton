import unittest

import numpy as np

from src.transcriptome_calibrator import (
    fit_single_promoter,
    predict_promoter_curve,
    weighted_rmse,
)


class TranscriptomeCalibratorTests(unittest.TestCase):
    def test_weighted_rmse_zero_when_equal(self):
        y = np.array([0.1, 0.2, 0.3], dtype=float)
        w = np.array([1.0, 1.0, 1.0], dtype=float)
        self.assertAlmostEqual(weighted_rmse(y, y, w), 0.0, places=12)

    def test_fit_single_promoter_improves_over_flat_guess(self):
        dpa = np.linspace(0.0, 50.0, 26)
        true_vec = np.array([30.0, 8.0, 6.0, 0.02, 1.0], dtype=float)
        observed = predict_promoter_curve(dpa, true_vec)
        observed = np.clip(observed + np.linspace(-0.01, 0.01, len(dpa)), 0.0, 1.2)
        weight = np.ones_like(observed)

        flat_guess = np.full_like(observed, 0.5)
        rmse_flat = weighted_rmse(observed, flat_guess, weight)

        fitted = fit_single_promoter(
            dpa=dpa,
            observed=observed,
            weight=weight,
            bounds=[(20.0, 40.0), (1.0, 12.0), (1.0, 12.0), (0.0, 0.2), (0.6, 1.6)],
            seed=123,
            maxiter=45,
        )
        self.assertLess(fitted["rmse"], rmse_flat)
        self.assertGreaterEqual(fitted["activation_dpa"], 20.0)
        self.assertLessEqual(fitted["activation_dpa"], 40.0)


if __name__ == "__main__":
    unittest.main()
