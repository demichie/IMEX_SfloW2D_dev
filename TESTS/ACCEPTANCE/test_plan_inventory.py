#!/usr/bin/env python3
"""Validate that the current acceptance plan names executable, existing suites."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PLAN = ROOT / "TESTS/ACCEPTANCE/acceptance_plan.json"


def main():
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    unit = plan["unit_tests"]
    solver = [item["name"] for item in plan["solver_tests"]]
    active = unit + solver
    errors = []
    for name in active:
        if not (ROOT / "TESTS" / name / "run_test.sh").is_file():
            errors.append(f"{name}: run_test.sh missing")
    if len(active) != len(set(active)):
        errors.append("duplicate active suite in acceptance plan")
    for criterion in plan["criteria"]:
        missing = set(criterion["tests"]) - set(active)
        if missing:
            errors.append(f'{criterion["id"]}: criteria reference inactive suites: {sorted(missing)}')
        historical = criterion.get("historical_tests", [])
        if set(historical) & set(active):
            errors.append(f'{criterion["id"]}: historical suite also declared active')
        pending = criterion.get("closure_status") == "historical_only_pending_current_revalidation"
        if bool(historical) != pending:
            errors.append(f'{criterion["id"]}: historical coverage status mismatch')
        if pending and criterion.get("coverage") != "historical_evidence_pending_current_revalidation":
            errors.append(f'{criterion["id"]}: historical coverage not identified')
    for error in errors:
        print("FAIL:", error)
    if errors:
        return 1
    print(f"PASS: {len(unit)} unit, {len(solver)} solver suites; "
          f"{sum(bool(c.get('historical_tests')) for c in plan['criteria'])} N7 criteria awaiting revalidation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
