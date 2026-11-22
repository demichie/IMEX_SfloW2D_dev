#!/usr/bin/env python3
"""Solver-independent tests of the extended two-layer evidence checks."""

import csv
from pathlib import Path
import tempfile
import unittest

from extended_validation import (
    BUDGET_COLUMNS, RUN, budgets, cells, compare_snapshots, set_input,
)

def fixed_record(values):
    return "".join(f"{value:20.12E}" for value in values) + "\n"


class ExtendedEtnaValidationTest(unittest.TestCase):
    def test_binary_restart_input_does_not_keep_l2_text_path(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            set_input(root, 0.01, 0.005, restart_binary=True)
            contents = (root / "IMEX_SfloW2D.inp").read_text()
            self.assertIn("&RESTART_L1_PARAMETERS", contents)
            self.assertIn("N_RESTART_FILES=1,", contents)
            self.assertIn("RESTART_FILES='restart.bin',", contents)
            self.assertNotIn("&RESTART_L2_PARAMETERS", contents)
            self.assertNotIn("RESTART_FILES_L2=", contents)
            self.assertNotIn("N_RESTART_FILES_L2=", contents)
            self.assertNotIn("RESTART_FILES='OUTPUT/", contents)
            self.assertNotIn("OUTPUT_OVERWRITE=T", contents)
            self.assertNotIn("etna_L2_0000.q_2d", contents)
            self.assertIn("N_LAYERS = 2", contents)
            self.assertIn("OUTPUT_RUNOUT_FLAG=F,", contents)

    def test_text_initialization_keeps_both_layer_restart_groups(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            set_input(root, 0.01, 0.005)
            contents = (root / "IMEX_SfloW2D.inp").read_text()
            self.assertIn("&RESTART_L1_PARAMETERS", contents)
            self.assertIn("RESTART_FILES='./DEM/init_volume.asc',", contents)
            self.assertIn("&RESTART_L2_PARAMETERS", contents)
            self.assertIn("RESTART_FILES_L2='NONE',", contents)
            self.assertNotIn("RESTART_FILES='restart.bin',", contents)

    def test_three_budget_records_are_consistent(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "OUTPUT" / RUN
            (root / "metadata").mkdir(parents=True)
            (root / "fields").mkdir()
            (root / "metadata" / "output_manifest.csv").write_text(
                "snapshot_index,time_seconds\n0,0\n1,0.005\n2,0.01\n"
            )
            for n, time in enumerate((0.0, 0.005, 0.01)):
                values = [
                    time, 10.0 - n, float(n), 10.0,
                    0.2 - 0.02 * n, 0.02 * n, 0.2,
                    5.0 - 0.5 * n, 0.5 * n, 0.0, 5.0,
                ]
                with (root / "fields" / f"{RUN}_budget_{n:04d}.csv").open("w", newline="") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(BUDGET_COLUMNS)
                    writer.writerow(values)
            result = budgets(Path(directory))
            self.assertEqual(len(result), 3)
            self.assertEqual(result[2]["mass_L2"], 2.0)
            with (root / "fields" / f"{RUN}_budget_0002.csv").open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(BUDGET_COLUMNS)
                writer.writerow([0.01, 8.0, 2.0, 11.0, 0.16, 0.04, 0.2, 4.0, 1.0, 0.0, 5.0])
            with self.assertRaisesRegex(AssertionError, "mass_total"):
                budgets(Path(directory))

    def test_snapshot_comparison_rejects_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            roots = [Path(directory) / name for name in ("a", "b")]
            for root in roots:
                folder = root / "OUTPUT" / RUN / "snapshots" / "0001_t_0"
                folder.mkdir(parents=True)
                for layer in (1, 2):
                    (folder / f"{RUN}_L{layer}_0001.q_2d").write_text(
                        fixed_record((0.0, 0.0, 12.0, 0.0, 0.0, 0.0, 5.0)) + "\n",
                        encoding="ascii"
                    )
            compare_snapshots(roots[0], roots[1], 1)
            mismatch = roots[1] / "OUTPUT" / RUN / "snapshots" / "0001_t_0" / f"{RUN}_L2_0001.q_2d"
            mismatch.write_text(
                fixed_record((0.0, 0.0, 12.1, 0.0, 0.0, 0.0, 5.0)) + "\n",
                encoding="ascii"
            )
            with self.assertRaisesRegex(AssertionError, "Snapshot index"):
                compare_snapshots(roots[0], roots[1], 1)
            mismatch.write_text(
                fixed_record((0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 5.0)),
                encoding="ascii"
            )
            with self.assertRaisesRegex(AssertionError, "Negative"):
                list(cells(mismatch))


if __name__ == "__main__":
    unittest.main()
