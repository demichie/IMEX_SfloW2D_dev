"""Focused checks of audit parsing and evidence classification, not model physics."""

import json
import math
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

    def test_closure_inventory_and_action_mapping(self):
        """Keep the original inventories complete and satisfied contracts free of new blockers."""
        tools = Path(__file__).resolve().parent
        execution = json.loads((tools / "acceptance_plan.json").read_text())
        closure = json.loads((tools / execution["closure_plan"]).read_text())
        remaining = closure["remaining_actions"]
        completed = closure.get("completed_actions", [])
        actions = {item["id"]: item for item in remaining + completed}
        self.assertEqual(len(actions), len(remaining) + len(completed))
        remaining_ids = {item["id"] for item in remaining}
        self.assertTrue(all(item["status"] == "satisfied" for item in completed))
        self.assertEqual(sum(a["kind"] == "test" for a in actions.values()), 6)
        self.assertEqual(sum(a["kind"] == "decision" for a in actions.values()), 2)
        self.assertEqual(closure["milestone_status"]["N9"], "open")
        records = []
        for group, prefix, count in (("mass_exchange", "I", 13), ("hydrodynamics", "H", 14)):
            items = closure["required_tests"][group]
            self.assertEqual([i["id"] for i in items], [f"{prefix}{n:02}" for n in range(1, count + 1)])
            records.extend(items)
        criteria = {item["id"]: item for item in closure["criteria"]}
        expected = {f"N7-{n:02}" for n in range(1, 10)} | {f"N8-{n:02}" for n in range(1, 6)}
        self.assertEqual(set(criteria), expected)
        self.assertEqual(len(criteria), len(closure["criteria"]))
        for milestone, group in (("N7", "hydrodynamics"), ("N8", "mass_exchange")):
            contracts_closed = all(c["closure_status"] == "satisfied" for c in criteria.values()
                                   if c["id"].startswith(milestone + "-"))
            cases_closed = all(c["status"] == "satisfied" for c in closure["required_tests"][group]
                               if not (milestone == "N7" and c["id"] == "H14"))
            self.assertEqual(closure["milestone_status"][milestone],
                             "satisfied" if contracts_closed and cases_closed else "open")
        records.extend(criteria.values())
        for item in records:
            status = item.get("closure_status", item.get("status"))
            self.assertIn(status, ("satisfied", "partial", "missing", "decision"))
            if status == "satisfied":
                self.assertEqual(item["remaining_actions"], [], item["id"])
            else:
                self.assertTrue(item["remaining_actions"], item["id"])
            self.assertTrue(set(item["remaining_actions"]) <= remaining_ids, item["id"])
        for item in execution["criteria"]:
            if item["id"] in criteria:
                for field in ("closure_status", "remaining_actions"):
                    self.assertEqual(item[field], criteria[item["id"]][field])
        required_ids = {item["id"] for group in closure["required_tests"].values() for item in group}
        for action in actions.values():
            self.assertTrue(set(action["criteria"]) <= set(criteria))
            self.assertTrue(set(action["normative_items"]) <= required_ids)
            self.assertTrue(action["exit_conditions"])
            for criterion in action["criteria"]:
                if action["id"] in remaining_ids:
                    self.assertIn(action["id"], criteria[criterion]["remaining_actions"])
                else:
                    self.assertNotIn(action["id"], criteria[criterion]["remaining_actions"])

    def test_closure_evidence_provenance(self):
        """Validate archived evidence hashes and links without claiming a fresh solver execution."""
        tools = Path(__file__).resolve().parent
        repo = tools.parents[1]
        closure = json.loads((tools / "n7_n8_closure_plan.json").read_text())
        self.assertTrue((repo / closure["report"]).is_file())
        artifacts = {}
        for name, record in closure["evidence_artifacts"].items():
            path = repo / record["path"]
            self.assertEqual(digest(path), record["sha256"])
            artifacts[name] = json.loads(path.read_text())
        for criterion in closure["criteria"]:
            self.assertTrue(set(criterion["evidence"]) <= set(artifacts))
        latest = artifacts["latest"]
        self.assertEqual(latest["summary"]["failed"], 0)
        self.assertEqual(latest["summary"]["passed"], 36)
        self.assertFalse(latest["summary"]["n9_accepted"])
        self.assertEqual(len(latest["thread_comparisons"]), 6)
        self.assertTrue(all(p["status"] == "pass" for p in latest["thread_comparisons"]))
        n8 = artifacts["n8"]
        records = n8["mass_exchange_production_evidence"]["records"]
        self.assertEqual({(r["profile"], r["actual_threads"]) for r in records},
                         {("strict", 1), ("strict", 4), ("optimized", 1), ("optimized", 4)})
        self.assertEqual([len(r["diagnostics"]) for r in records], [12] * 4)
        for filename in ("src/equation_terms_2d.f90", "src/mass_exchange_2d.f90",
                         "src/geometry_2d.f90", "TESTS/TEST_MASS_EXCHANGE/test_mass_exchange.f90"):
            self.assertEqual(n8["source_and_fixture_sha256"][filename],
                             latest["source_and_fixture_sha256"][filename])
        for group in closure["required_tests"].values():
            for item in group:
                for test in item["evidence_tests"]:
                    self.assertTrue((repo / "TESTS" / test / "run_test.sh").is_file())

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

    def test_n8_a_completion_evidence(self):
        """Close only the bounded cell-law package and preserve historical payloads and limits."""
        repo = Path(__file__).resolve().parents[2]
        closure = json.loads((repo / "TESTS/ACCEPTANCE/n7_n8_closure_plan.json").read_text())
        evidence = json.loads((repo / closure["evidence_artifacts"]["n8_a"]["path"]).read_text())
        latest = json.loads((repo / closure["evidence_artifacts"]["latest"]["path"]).read_text())
        historical = json.loads((repo / closure["evidence_artifacts"]["n8"]["path"]).read_text())
        self.assertEqual([a["id"] for a in closure["completed_actions"]], ["N8-A"])
        self.assertEqual(closure["completed_actions"][0]["evidence"], ["n8_a"])
        self.assertEqual(next(c for c in closure["criteria"] if c["id"] == "N8-01")["closure_status"],
                         "satisfied")
        self.assertFalse(evidence["summary"]["solver_changes"])
        self.assertFalse(evidence["summary"]["full_36_run_audit_repeated"])
        for name, value in evidence["source_and_fixture_sha256"].items():
            if name.startswith("src/"):
                self.assertEqual(value, latest["source_and_fixture_sha256"][name])
        records = evidence["mass_exchange_production_evidence"]["records"]
        self.assertEqual({(r["profile"], r["actual_threads"]) for r in records},
                         {("strict", 1), ("strict", 4), ("optimized", 1), ("optimized", 4)})
        self.assertEqual(len(records), 4)
        for record in records:
            old = next(r for r in historical["mass_exchange_production_evidence"]["records"]
                       if (r["profile"], r["actual_threads"]) == (record["profile"], record["actual_threads"]))
            self.assertEqual(len(record["diagnostics"]), 18)
            cases = [d["case"] for d in record["diagnostics"]]
            self.assertEqual(cases[:12], [d["case"] for d in old["diagnostics"]])
            self.assertEqual(cases[12:], ["gas_packing", "gas_reserve", "gas_exhausted", "gas_inflow",
                                          "gas_combined", "gas_masked"])
            self.assertEqual(record["flat_cases"], cases)
            for name, value in old["sha256"].items():
                self.assertEqual(record["sha256"][name], value)
            for diagnostic in record["diagnostics"]:
                self.assertEqual(len(diagnostic["production"]), 5)
                self.assertEqual(len(diagnostic["reference"]), 5)
                for actual, reference in zip(diagnostic["production"], diagnostic["reference"]):
                    self.assertTrue(math.isfinite(actual) and math.isfinite(reference))
                    self.assertLessEqual(abs(actual-reference), 2e-12*max(1.0, abs(reference)))
        for profile in ("strict", "optimized"):
            pair = [r for r in records if r["profile"] == profile]
            self.assertEqual(pair[0]["sha256"], pair[1]["sha256"])
        self.assertEqual({c["case"] for c in evidence["negative_controls"]},
                         {"missing_flat_assertion", "missing_gas_case", "wrong_actual_team",
                          "corrupt_volume_diagnostic", "nonfinite_reference"})
        self.assertTrue(all(c["status"] == "rejected" for c in evidence["negative_controls"]))
        for milestone in ("N7", "N8", "N9"):
            self.assertEqual(evidence["summary"][milestone], "open")


if __name__ == "__main__":
    unittest.main()
