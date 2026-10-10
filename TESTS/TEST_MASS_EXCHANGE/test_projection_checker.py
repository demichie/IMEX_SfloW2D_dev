"""Negative controls for N8-B evidence completeness and frozen error limits."""
import contextlib
import io
from pathlib import Path
import tempfile
import unittest

from check_projection_stress import SIZES, STEPS, check


class ProjectionCheckerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="imex-n8-b-checker-")
        self.addCleanup(self.temporary.cleanup)
        self.path = Path(self.temporary.name)/"projection.log"
        (self.path.parent/"projection_snapshot.bin").write_bytes(b"checker fixture only")
        self.lines = ["Actual OpenMP team: 1"]
        for dimension, nx, ny in SIZES:
            for name, steps in STEPS.items():
                l1 = 1/nx**2 if name == "smooth" else 1/nx
                for iteration in range(1, steps+1):
                    self.lines.append(f"N8B_RECORD {name} {dimension} {nx} {ny} {iteration} "
                                      f"0 0 0 {l1} 0.001 0")

    def evaluate(self, lines):
        self.path.write_text("\n".join(lines)+"\n")
        with contextlib.redirect_stdout(io.StringIO()):
            return check(self.path, "strict", 1)

    def test_complete_evidence(self):
        result = self.evaluate(self.lines)
        self.assertEqual(len(result["records"]), 588)
        self.assertTrue(result["unequal_area_algebra"]["wrong_weights_rejected"])

    def test_missing_case_or_wrong_team(self):
        for lines in (self.lines[:-1], ["Actual OpenMP team: 4"]+self.lines[1:]):
            with self.assertRaises(AssertionError):
                self.evaluate(lines)

    def test_nonfinite_or_wrong_volume_or_nonlocal_mismatch(self):
        for index, value in ((6, "NaN"), (7, "1"), (11, "0.001")):
            lines = self.lines.copy()
            fields = lines[1].split(); fields[index] = value; lines[1] = " ".join(fields)
            with self.assertRaises(AssertionError):
                self.evaluate(lines)

    def test_missing_l1_refinement(self):
        lines = self.lines.copy()
        for i, line in enumerate(lines):
            if line.startswith("N8B_RECORD smooth "):
                fields = line.split(); fields[9] = "0.01"; lines[i] = " ".join(fields)
        with self.assertRaises(AssertionError):
            self.evaluate(lines)


if __name__ == "__main__":
    unittest.main()
