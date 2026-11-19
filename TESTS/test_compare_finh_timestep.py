#!/usr/bin/env python3
"""Regression for scientific-notation FINH timestep diagnostic extraction."""

import importlib.util
from pathlib import Path
import unittest


PATH = Path(__file__).with_name("compare_finh_timestep.py")
SPEC = importlib.util.spec_from_file_location("compare_finh_timestep", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class TestFinhDiagnostics(unittest.TestCase):

    def test_reads_full_precision_independently_of_rounded_step_log(self):
        log = (
            "t =     0.0006 dt =  0.00000  mass =  5.537E-002\n"
            "TIMESTEP_DIAG step=100 t=  6.0000000000000000E-04"
            " cfl_dt=  1.2345678901234567E-12"
            " accepted_dt=  1.2345678901234567E-12"
            " max_a_x=  1.9440000000000000E+09"
            " max_a_y=  2.0000000000000000E+00\n"
        )
        values = MODULE.extract_timesteps(log)
        self.assertEqual(len(values), 1)
        self.assertEqual(values[0]["step"], 100)
        self.assertAlmostEqual(values[0]["t"], 0.0006)
        self.assertAlmostEqual(values[0]["cfl_dt"], 1.2345678901234567e-12)
        self.assertEqual(values[0]["cfl_dt"], values[0]["accepted_dt"])
        self.assertAlmostEqual(values[0]["max_a_x"], 1.944e9)

    def test_skips_unrelated_log_lines(self):
        log = (
            "HEARTBEAT: TEST_DYNAMIC_FINH elapsed=60s last_step=t = 0.0006\n"
            "INVALID_TIMESTEP t= 0.0 cfl_dt= 0.0 accepted_dt= 0.0\n"
        )
        self.assertEqual(MODULE.extract_timesteps(log), [])


if __name__ == "__main__":
    unittest.main()
