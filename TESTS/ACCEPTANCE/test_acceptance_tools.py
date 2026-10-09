"""Focused checks of audit parsing and evidence classification, not model physics."""

import json
from pathlib import Path
import tempfile
import unittest

from capture_solver import digest, fortran_float
from run_acceptance import compare_threads, execute, overlay_worktree


class AuditToolsTests(unittest.TestCase):
    def test_fortran_exponents(self):
        self.assertEqual(fortran_float("0.987774513049-133"), 0.987774513049e-133)
        self.assertEqual(fortran_float(b"-0.322725514122-157"), -0.322725514122e-157)
        self.assertEqual(fortran_float("2.5D+03"), 2500.0)
        self.assertEqual(fortran_float("0.0"), 0.0)
        with self.assertRaises(ValueError):
            fortran_float("1-3")

    @staticmethod
    def pair():
        return [
            {"profile": "strict", "test": "case", "outer_threads": threads,
             "status": "pass", "calls": [{"actual_threads": threads,
                 "files": {"case_0001.q_2d": {"sha256": "same"}}}]}
            for threads in (1, 4)]

    def test_actual_thread_evidence_required(self):
        values = self.pair()
        self.assertEqual(compare_threads(values)[0]["status"], "pass")
        values[1]["calls"][0]["actual_threads"] = 1
        self.assertEqual(compare_threads(values)[0]["status"], "fail")

    def test_empty_or_different_outputs_do_not_pass(self):
        values = self.pair()
        values[1]["calls"][0]["files"]["case_0001.q_2d"]["sha256"] = "different"
        self.assertEqual(compare_threads(values)[0]["status"], "fail")
        for value in values:
            value["calls"][0]["files"] = {}
        self.assertEqual(compare_threads(values)[0]["status"], "fail")

    def test_failed_trajectory_is_not_equivalence(self):
        values = self.pair()
        values[1]["status"] = "fail"
        self.assertEqual(compare_threads(values)[0]["status"], "not_completed")

    def test_command_log_and_nonzero_exit(self):
        import os
        import sys
        with tempfile.TemporaryDirectory(prefix="imex-audit-test-") as temporary:
            directory = Path(temporary)
            log = directory / "test.log"
            result = execute([sys.executable, "-c", "print('evidence'); raise SystemExit(7)"],
                             directory, log, dict(os.environ), 10)
            self.assertEqual(result["returncode"], 7)
            self.assertEqual(log.read_text().strip(), "evidence")
            self.assertEqual(result["log_sha256"], digest(log))

    def test_plan_links(self):
        tools = Path(__file__).resolve().parent
        plan = json.loads((tools / "acceptance_plan.json").read_text())
        ids = [criterion["id"] for criterion in plan["criteria"]]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertTrue(all(criterion["gap"] for criterion in plan["criteria"]))
        names = plan["unit_tests"] + [test["name"] for test in plan["solver_tests"]]
        for name in names:
            self.assertTrue((tools.parents[1] / "TESTS" / name / "run_test.sh").is_file())

    def test_candidate_overlay_excludes_ignored_and_tracks_removals(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory(prefix="imex-audit-overlay-") as temporary:
            root = Path(temporary)
            repo, export = root / "repo", root / "export"
            repo.mkdir()
            export.mkdir()
            (repo / "source.f90").write_text("candidate source\n")
            (repo / "new_test.f90").write_text("new test\n")
            (repo / "ignored.o").write_bytes(b"stale object")
            (export / "removed.f90").write_text("old source\n")
            answers = [b"source.f90\0removed.f90\0", b"new_test.f90\0",
                       b"candidate patch", " M source.f90\n D removed.f90\n?? new_test.f90\n"]
            with patch("run_acceptance.subprocess.check_output", side_effect=answers):
                metadata = overlay_worktree(repo, export)
            self.assertEqual((export / "source.f90").read_text(), "candidate source\n")
            self.assertTrue((export / "new_test.f90").is_file())
            self.assertFalse((export / "removed.f90").exists())
            self.assertFalse((export / "ignored.o").exists())
            self.assertTrue(metadata["base_revision_only"])
            self.assertEqual(metadata["patch_sha256"], digest(root / "candidate.patch"))


if __name__ == "__main__":
    unittest.main()
