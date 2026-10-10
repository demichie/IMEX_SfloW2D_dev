"""Negative controls for N7-A assertions, separate from numerical model gates."""
import unittest

import numpy as np

import check_cases as checker


class CheckerTests(unittest.TestCase):
    def test_field_components_and_nonfinite(self):
        initial = checker.liquid_state(np.ones((1, 8)))
        checker.field_error("unchanged", initial, initial, initial, checker.ROUND)
        for component in range(5):
            bad = initial.copy()
            bad[0, 3, component] += checker.scales(initial)[component] * 1e-5
            with self.assertRaises(AssertionError):
                checker.field_error("corrupt", bad, initial, initial, checker.ROUND)
        bad[0, 0, 1] = np.nan
        with self.assertRaises(AssertionError):
            checker.field_error("nonfinite", bad, initial, initial, checker.ROUND)

    def test_raw_admissibility_and_budget(self):
        initial = checker.liquid_state(np.ones((1, 8)))
        statistics = np.zeros((2, 15))
        statistics[:, :2] = (0.01, 0.1)
        statistics[:, 2:5] = 1
        statistics[:, 11] = (0.01, 0.02)
        statistics[:, 12:] = (300, 0, 1)
        fields = {"statistics": statistics, "active": np.ones((1, 8)),
                  "known": np.tile(initial, (3, 1, 1, 1)),
                  "solved": np.tile(initial, (3, 1, 1, 1)), "raw_final": initial.copy()}
        checker.validate(fields, initial, 2)
        for column, value in ((2, -0.1), (5, 1e4), (6, 1), (8, 1), (12, 200), (14, -1)):
            bad = {key: value.copy() for key, value in fields.items()}
            bad["statistics"][0, column] = value
            with self.assertRaises(AssertionError):
                checker.validate(bad, initial, 2)
        bad = {key: value.copy() for key, value in fields.items()}
        bad["statistics"][0, 0] = 1
        with self.assertRaises(AssertionError):
            checker.validate(bad, initial, 2)

    def test_ritter_integral(self):
        h, hu = checker.ritter_cell_means(160, 20, 0.2)
        self.assertAlmostEqual(float(h.sum() * 20 / 160), 10, places=12)
        self.assertTrue(np.all(h >= 0))
        self.assertTrue(np.all(hu >= -checker.ROUND))
        # Independent fine midpoint quadrature checks the cell-average formula.
        x = (np.arange(160 * 1000) + 0.5) * 20 / (160 * 1000) - 10
        c = np.sqrt(checker.GRAV * (1 - checker.RHO_A / checker.RHO))
        hf = np.where(x < -c * 0.2, 1, np.where(x < 2*c*0.2, (2*c-x/0.2)**2/(9*c*c), 0))
        uf = np.where((x >= -c*0.2) & (x < 2*c*0.2), 2*(c+x/0.2)/3, 0)
        np.testing.assert_allclose(h, hf.reshape(160, 1000).mean(axis=1), rtol=0, atol=1e-8)
        np.testing.assert_allclose(hu, (hf*uf).reshape(160, 1000).mean(axis=1), rtol=0, atol=4e-8)


if __name__ == "__main__":
    unittest.main()
