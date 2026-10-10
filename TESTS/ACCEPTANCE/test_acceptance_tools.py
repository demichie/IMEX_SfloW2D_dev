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
        action = next(a for a in closure["completed_actions"] if a["id"] == "N8-A")
        self.assertEqual(action["evidence"], ["n8_a"])
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

    def test_n8_b_completion_evidence(self):
        """Close N8 only after refinement, repeated geometry, area approval and full regression."""
        repo = Path(__file__).resolve().parents[2]
        closure = json.loads((repo / "TESTS/ACCEPTANCE/n7_n8_closure_plan.json").read_text())
        evidence = json.loads((repo / closure["evidence_artifacts"]["n8_b"]["path"]).read_text())
        full = json.loads((repo / closure["evidence_artifacts"]["n8_b_full"]["path"]).read_text())
        self.assertEqual(closure["milestone_status"], {"N7": "open", "N8": "satisfied", "N9": "open"})
        completed = {a["id"] for a in closure["completed_actions"]}
        self.assertTrue({"N8-A", "N8-B", "D-N8-AREA"} <= completed)
        self.assertEqual(evidence["area_decision"]["status"], "accepted_by_user")
        self.assertTrue(evidence["area_decision"]["algebra"]["wrong_weights_rejected"])
        self.assertEqual(len(evidence["area_decision"]["algebra"]["records"]), 3)
        self.assertEqual(full["summary"]["passed"], 36)
        self.assertEqual(full["summary"]["failed"], 0)
        self.assertFalse(full["summary"]["n9_accepted"])
        self.assertEqual(set(full["production_changes"]), {"src/geometry_2d.f90"})
        self.assertEqual(len(full["comparison_to_baseline"]["cases"]), 16)
        self.assertTrue(all(c["status"] == "unchanged" for c in full["comparison_to_baseline"]["cases"]))
        records = evidence["records"]
        self.assertEqual({(r["profile"], r["actual_threads"]) for r in records},
                         {("strict", 1), ("strict", 4), ("optimized", 1), ("optimized", 4)})
        self.assertEqual(len(records), 4)
        self.assertEqual(sum(r["record_count"] for r in records), 2352)
        for profile in ("strict", "optimized"):
            pair = sorted((r for r in records if r["profile"] == profile), key=lambda r: r["actual_threads"])
            self.assertEqual(pair[0]["snapshot_sha256"], pair[1]["snapshot_sha256"])
            self.assertEqual(pair[1]["identical_metrics_reference"], profile+":1")
            self.assertEqual(len(pair[0]["groups"]), 24)
            self.assertEqual(sum(g["updates"] for g in pair[0]["groups"]), 588)
            for series in pair[0]["refinements"]:
                self.assertEqual(len(series["L1"]), 3)
                self.assertTrue(all(0 < ratio < 0.8 for ratio in series["L1_ratios"]))
            for group in pair[0]["groups"]:
                self.assertLessEqual(group["max_volume_balance_error"], 2048*2.220446049250313e-16)
                self.assertTrue(all(math.isfinite(v) for v in group["metrics_min"]+group["metrics_max"]))
        for filename, value in evidence["source_and_fixture_sha256"].items():
            if filename.startswith("src/"):
                self.assertEqual(value, full["source_and_fixture_sha256"][filename])

    def test_stationary_composition_correction_evidence(self):
        """Keep the verified stationary subset distinct from contact and CFL acceptance."""
        repo = Path(__file__).resolve().parents[2]
        closure = json.loads((repo / "TESTS/ACCEPTANCE/n7_n8_closure_plan.json").read_text())
        record = closure["evidence_artifacts"]["hydrostatic_roundoff"]
        evidence = json.loads((repo / record["path"]).read_text())
        self.assertEqual(digest(repo / record["path"]), record["sha256"])
        self.assertEqual(closure["milestone_status"], {"N7": "open", "N8": "satisfied", "N9": "open"})
        self.assertEqual({a["id"] for a in closure["remaining_actions"]},
                         {"N7-B", "N7-C", "N7-D", "D-N7-CFL"})
        cases = {c["id"]: c for c in closure["required_tests"]["hydrodynamics"]}
        self.assertEqual(cases["H04"]["status"], "satisfied")
        self.assertEqual(cases["H04"]["remaining_actions"], [])
        self.assertEqual(cases["H04"]["evidence_tests"], ["TEST_HYDROSTATIC_ROUNDOFF"])
        self.assertEqual(cases["H05"]["status"], "missing")
        self.assertEqual(cases["H05"]["remaining_actions"], ["N7-B"])
        self.assertEqual(evidence["summary"]["N7_B"], "open")
        self.assertEqual(evidence["production_changes"], ["src/hyperbolic_2d.f90"])
        contract = evidence["guard_contract"]
        self.assertEqual(contract["correction"]["roundoff_multiplier"], 64)
        self.assertEqual(contract["controls"]["velocities"], [1e-20, 1e-12])
        self.assertEqual(contract["controls"]["transverse_cells"], 5)
        frozen = repo / "TESTS/TEST_N7_COMPOSITION/contract.json"
        self.assertEqual(digest(frozen), evidence["equilibrium_contract_sha256"])
        for profile in evidence["profiles"].values():
            self.assertEqual(profile["equilibrium_count"], 90)
            self.assertEqual(profile["control_count"], 36)
            for metric in ("maximum_equilibrium_error", "maximum_force_scaled_residual",
                           "maximum_component_budget", "maximum_boundary_mass_flux"):
                self.assertEqual(profile[metric], 0)
            self.assertGreater(profile["minimum_physical_nonequilibrium_speed"], 1e-10)
        full = evidence["full_audit"]
        self.assertEqual(full["summary"],
                         {"test_runs": 38, "failed": 0, "passed": 38, "n9_accepted": False})
        self.assertEqual(len(full["tests"]), 38)
        self.assertTrue(all(t["status"] == "pass" for t in full["tests"]))
        self.assertEqual(len(full["thread_comparisons"]), 6)
        self.assertTrue(all(t["status"] == "pass" for t in full["thread_comparisons"]))
        comparisons = evidence["comparison_to_parent_full_audit"]
        self.assertEqual(len(comparisons), 16)
        self.assertTrue(all(t["effective_inputs_identical"] for t in comparisons))
        self.assertEqual(sum(t["status"] == "identical" for t in comparisons), 10)
        self.assertEqual({t["test"] for t in comparisons if t["status"] == "changed"},
                         {"TEST_PCCU_GATE_H", "TEST_PCCU_INCLINED_EXCAVATION"})

    def test_n7_a_completion_evidence(self):
        """Close the frozen N7-A cases, retaining output changes and remaining gates explicitly."""
        repo = Path(__file__).resolve().parents[2]
        closure = json.loads((repo / "TESTS/ACCEPTANCE/n7_n8_closure_plan.json").read_text())
        evidence = json.loads((repo / closure["evidence_artifacts"]["n7_a"]["path"]).read_text())
        spatial = json.loads((repo / closure["evidence_artifacts"]["spatial"]["path"]).read_text())
        action = next(a for a in closure["completed_actions"] if a["id"] == "N7-A")
        self.assertEqual(action["status"], "satisfied")
        self.assertEqual(action["evidence"], ["n7_a"])
        self.assertEqual(closure["milestone_status"], {"N7": "open", "N8": "satisfied", "N9": "open"})
        self.assertEqual({a["id"] for a in closure["remaining_actions"]},
                         {"N7-B", "N7-C", "N7-D", "D-N7-CFL"})
        self.assertEqual(next(c for c in closure["criteria"] if c["id"] == "N7-01")["closure_status"],
                         "satisfied")
        for item in closure["required_tests"]["hydrodynamics"]:
            if item["id"] in {"H01", "H02", "H03", "H07", "H11"}:
                self.assertEqual(item["status"], "satisfied")
                self.assertEqual(item["remaining_actions"], [])
                self.assertIn("TEST_N7_EQUILIBRIUM", item["evidence_tests"])
        summary = evidence["summary"]
        self.assertEqual(summary["N7_A"], "satisfied")
        self.assertEqual(summary["production_solver_runs"], 408)
        self.assertEqual(summary["bitwise_thread_pairs"], 204)
        self.assertEqual(summary["observed_steps"], 27360)
        self.assertEqual(len(evidence["case_inventory"]), 102)
        self.assertEqual(len(set(evidence["case_inventory"])), 102)
        contract = json.loads((repo / "TESTS/TEST_N7_EQUILIBRIUM/contract.json").read_text())
        self.assertEqual(evidence["contract"], contract)
        self.assertEqual(contract["stages"], [2, 3, 4])
        self.assertEqual(contract["roundoff_epsilon_multiplier"], 32768)
        roundoff = 32768 * 2.220446049250313e-16
        records = evidence["records"]
        self.assertEqual({r["profile"] for r in records}, {"strict", "optimized"})
        self.assertEqual(len(records), 2)
        self.assertEqual(sum(r["solver_runs"] for r in records), 408)
        self.assertEqual(2*sum(g["one_team_steps"] for r in records for g in r["groups"].values()),
                         27360)
        for record in records:
            self.assertEqual(record["actual_threads"], [1, 4])
            self.assertEqual(record["case_count"], 102)
            self.assertEqual({name: g["cases"] for name, g in record["groups"].items()},
                             {"lake": 72, "advection": 9, "ritter": 9, "excavation": 12})
            for name, group in record["groups"].items():
                diagnostics = group["admissibility_envelope"]
                for state in ("known", "solved", "raw_final"):
                    self.assertGreaterEqual(diagnostics[f"minimum_{state}_h"]["minimum"], -roundoff)
                self.assertEqual(diagnostics["failed_local_solves"]["maximum"], 0)
                self.assertLessEqual(diagnostics["maximum_relative_mass_drift"]["maximum"], roundoff)
                self.assertLessEqual(diagnostics["maximum_dt_over_CFL"]["maximum"], 1+roundoff)
                for interval in diagnostics.values():
                    self.assertTrue(all(math.isfinite(v) for v in interval.values()))
                if name == "lake":
                    for field in ("maximum_equilibrium_error", "maximum_production_residual",
                                  "maximum_reference_residual"):
                        self.assertLessEqual(group[field], roundoff)
                    self.assertEqual(diagnostics["maximum_final_repair"]["maximum"], 0)
                else:
                    self.assertLessEqual(group["maximum_field_error_over_frozen_limit"], 1)
                if name == "excavation":
                    self.assertEqual(group["maximum_relative_uphill_mass"], 0)
            self.assertEqual(len(record["refinements"]), 6)
            for series in record["refinements"]:
                self.assertEqual(len(series["L1"]), 3)
                limit = contract[series["case"]]["maximum_refinement_ratio"]
                self.assertTrue(all(0 < ratio < limit for ratio in series["ratios"]))
        full = evidence["full_audit"]
        self.assertEqual(full["summary"],
                         {"test_runs": 37, "failed": 0, "passed": 37, "n9_accepted": False})
        self.assertEqual(len(full["tests"]), 37)
        self.assertTrue(all(t["status"] == "pass" for t in full["tests"]))
        self.assertEqual(len(full["thread_comparisons"]), 6)
        self.assertTrue(all(p["status"] == "pass" for p in full["thread_comparisons"]))
        self.assertEqual(set(evidence["production_changes"]),
                         {"src/hp_reconstruction_2d.f90", "src/state_conversion_2d.f90"})
        for name in ("hp_pccu_1d_core.py", "hp_pccu_2d_core.py"):
            path = "TESTS/TEST_SPATIAL_OPERATOR/reference/"+name
            self.assertEqual(evidence["source_and_fixture_sha256"][path],
                             spatial["source_and_fixture_sha256"][path])
        comparisons = evidence["comparison_to_N8_B"]["cases"]
        self.assertEqual(len(comparisons), 16)
        self.assertEqual(sum(c["status"] == "unchanged" for c in comparisons), 9)
        self.assertEqual(sum(c["status"] == "changed" for c in comparisons), 7)
        strict = [c for c in comparisons if c["profile"] == "strict_debug"]
        self.assertEqual(len(strict), 8)
        self.assertTrue(all(c["status"] == "unchanged" for c in strict))
        replays = evidence["optimized_field_characterization"]["records"]
        self.assertEqual(len(replays), 4)
        for replay in replays:
            self.assertEqual(len(replay["snapshots"]), replay["compared_snapshots"])
            self.assertGreater(replay["compared_snapshots"], 0)
            self.assertTrue(math.isfinite(replay["maximum_old_new_Linf"]))
            for snapshot in replay["snapshots"]:
                width = snapshot["shape"][1]-2
                for field in ("old_new_component_Linf", "old_new_component_mean_L1",
                              "old_optimized_to_strict_Linf", "new_optimized_to_strict_Linf"):
                    self.assertEqual(len(snapshot[field]), width)
                    self.assertTrue(all(math.isfinite(v) and v >= 0 for v in snapshot[field]))


if __name__ == "__main__":
    unittest.main()
