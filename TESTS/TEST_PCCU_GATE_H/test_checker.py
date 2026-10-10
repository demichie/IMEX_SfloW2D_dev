"""All-cell comparator negative controls; no solver trajectory is a golden field."""
import unittest

import numpy as np

import compare_fields as c
from build_dt_observer import instrument


class FieldCheckerTests(unittest.TestCase):
    def test_accepted_times_reject_invalid_trajectory(self):
        times = np.array([[0, 0.01], [0.01, 0.02]])
        np.testing.assert_array_equal(c.accepted_times(times, 0.03), [0.01, 0.02])
        for bad in (times[:1], times+np.nan, [[0, 0]], [[0, -1]], [[0, 0.01], [0.02, 0.01]]):
            with self.assertRaises(AssertionError):
                c.accepted_times(bad, 0.03)

    def test_observation_hook_is_unique_and_does_not_replace_advance(self):
        marker = "      CALL simulation%time_integration%advance("
        observed = instrument(marker)
        self.assertEqual(observed.count(marker), 1)
        self.assertIn("WRITE(accepted_time_unit)", observed)
        for bad in ("", marker+marker):
            with self.assertRaises(AssertionError):
                instrument(bad)

    def test_single_cell_corruption(self):
        h = np.ones((8, 10)); bed = np.zeros_like(h)
        self.assertEqual(c.compare(h, h, bed, h, 10)["cells_compared"], 80)
        bad = h.copy(); bad[3, 5] += 1e-5
        with self.assertRaises(AssertionError):
            c.compare(bad, h, bed, h, 10)

    def test_equal_aggregates_do_not_hide_displaced_mass(self):
        h = np.ones((8, 10)); bad = h.copy()
        bad[3, 4] += 0.01; bad[3, 5] -= 0.01
        self.assertEqual(float(bad.sum()), float(h.sum()))
        with self.assertRaises(AssertionError):
            c.compare(bad, h, np.zeros_like(h), h, 10)

    def test_nonfinite_and_shape(self):
        h = np.ones((8, 10)); bad = h.copy(); bad[2, 2] = np.nan
        with self.assertRaises(AssertionError):
            c.compare(bad, h, h, h, 10)
        with self.assertRaises(AssertionError):
            c.compare(h[:, :-1], h, h, h, 10)

    def test_historical_limits_are_separate(self):
        h = np.ones((8, 10)); bed = 5*np.ones_like(h)
        c.historical_comparison(h+1e-4, h, bed)
        for bad in (h+3e-4, np.full_like(h, np.nan)):
            with self.assertRaises(AssertionError):
                c.historical_comparison(bad, h, bed)
        bad = h.copy(); bad[1, 1] += 0.004
        with self.assertRaises(AssertionError):
            c.historical_comparison(bad, h, bed)


if __name__ == "__main__":
    unittest.main()
