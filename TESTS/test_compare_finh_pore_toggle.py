#!/usr/bin/env python3
"""Unit checks for the isolated FINH no-pore-pressure input transformation."""

import unittest

from compare_finh_pore_toggle import without_pore_pressure


class TestPoreToggle(unittest.TestCase):

    def test_removes_only_pore_group_and_disables_flag(self):
        original = (
            "&NEWRUN_PARAMETERS\n"
            " PORE_PRESSURE_FLAG = T,\n"
            " RHEOLOGY_FLAG=T,\n"
            " /\n"
            "&RHEOLOGY_PARAMETERS\n"
            " RHEOLOGY%MU_S=0.55D0,\n"
            " /\n"
            "&PORE_PRESSURE_PARAMETERS\n"
            " PORE_PRESSURE%F_INHIBIT_MODE='OFF',\n"
            " PORE_PRESSURE%GAS_LOSS_FLAG=T,\n"
            " /\n"
        )
        changed = without_pore_pressure(original)
        self.assertIn("PORE_PRESSURE_FLAG = F,", changed)
        self.assertIn("RHEOLOGY_FLAG=T", changed)
        self.assertIn("RHEOLOGY%MU_S=0.55D0", changed)
        self.assertNotIn("&PORE_PRESSURE_PARAMETERS", changed)

    def test_rejects_missing_namelist(self):
        with self.assertRaises(ValueError):
            without_pore_pressure("&NEWRUN_PARAMETERS\n PORE_PRESSURE_FLAG=T,\n/\n")

    def test_rejects_missing_flag(self):
        with self.assertRaises(ValueError):
            without_pore_pressure("&PORE_PRESSURE_PARAMETERS\n /\n")


if __name__ == "__main__":
    unittest.main()
