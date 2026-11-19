#!/usr/bin/env python3
"""Unit checks for the isolated FINH no-pore-pressure input transformation."""

import unittest

from compare_finh_pore_toggle import select_variants, without_pore_pressure


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

    def test_uses_both_variants_when_input_is_original(self):
        original = (
            "&NEWRUN_PARAMETERS\n PORE_PRESSURE_FLAG = T,\n /\n"
            "&PORE_PRESSURE_PARAMETERS\n PORE_PRESSURE%F_INHIBIT_MODE='OFF',\n /\n"
        )
        variants = select_variants(original)
        self.assertEqual([name for name, _ in variants], ["original", "no_pore_pressure"])
        self.assertIn("PORE_PRESSURE_FLAG = F,", variants[1][1])

    def test_accepts_manually_disabled_fixture_without_git(self):
        disabled = "&NEWRUN_PARAMETERS\n PORE_PRESSURE_FLAG = F,\n /\n"
        self.assertEqual(select_variants(disabled), (("no_pore_pressure", disabled),))

    def test_rejects_inconsistent_disabled_fixture(self):
        invalid = (
            "&NEWRUN_PARAMETERS\n PORE_PRESSURE_FLAG = F,\n /\n"
            "&PORE_PRESSURE_PARAMETERS\n /\n"
        )
        with self.assertRaises(ValueError):
            select_variants(invalid)

    def test_rejects_missing_namelist(self):
        with self.assertRaises(ValueError):
            without_pore_pressure("&NEWRUN_PARAMETERS\n PORE_PRESSURE_FLAG=T,\n/\n")

    def test_rejects_missing_flag(self):
        with self.assertRaises(ValueError):
            without_pore_pressure("&PORE_PRESSURE_PARAMETERS\n /\n")


if __name__ == "__main__":
    unittest.main()
