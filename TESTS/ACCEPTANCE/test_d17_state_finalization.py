#!/usr/bin/env python3
"""Static D19 regression guard for Lie state admissibility and atomic publication."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]


def load(path):
    return (ROOT / path).read_text(encoding="utf-8")


def body(text, name):
    # Use anchored lines: END SUBROUTINE may omit the routine name.
    opening = re.search(
        rf"(?im)^[ \t]*subroutine[ \t]+{re.escape(name)}[ \t]*\(",
        text)
    assert opening, f"missing routine {name}"
    closing = re.search(
        rf"(?im)^[ \t]*end[ \t]+subroutine(?:[ \t]+{re.escape(name)})?[ \t]*$",
        text[opening.end():])
    assert closing, f"missing END SUBROUTINE for {name}"
    return text[opening.end():opening.end() + closing.start()]


def ordered(code, patterns):
    offsets = []
    for pattern in patterns:
        match = re.search(pattern, code, re.I | re.S)
        assert match, f"missing invariant: {pattern}"
        offsets.append(match.start())
    assert offsets == sorted(offsets), "Lie publication order is incorrect"


def check():
    state = load("src/state/state_2d.f90")
    time = body(load("src/numerics/time_integration_2d.f90"), "advance_lie")
    executor = body(load("src/numerics/lie_source_executor_2d.f90"),
                    "execute_lie_source_groups")
    conservative = re.search(
        r"(?is)logical\s+function\s+conservative_mass_admissible\b.*?"
        r"end\s+function\s+conservative_mass_admissible", state)
    assert conservative, "missing conservative mass admissibility guard"
    assert "idx_add_gas_first:idx_add_gas_last" in conservative.group(0)
    assert "idx_solidEqn_first:idx_solidEqn_last" in conservative.group(0)

    ordered(time, [
        r"this%q_before\s*=\s*state%q",
        r"call\s+spatial_operator%evaluate",
        r"this%q_trial\(1:n_eqns,j,k\)\s*=",
        r"conservative_mass_admissible\(this%q_trial",
        r"call\s+state%snapshot_substrate",
        r"call\s+begin_geometry_transaction",
        r"call\s+execute_lie_source_groups",
        r"if\s*\(source_step_failed\)",
        r"state%q\s*=\s*this%q_trial",
        r"state%qp\s*=\s*this%qp_trial",
        r"call\s+commit_geometry_transaction",
    ])
    assert "CALL state%restore_substrate(substrate_before)" in time
    assert "CALL rollback_geometry_transaction" in time
    assert "CALL initialize_dry_physical_cell(this%qp_before(:,j,k))" in time
    assert "CALL initialize_dry_physical_cell(this%qp_trial(:,j,k))" in time
    ordered(executor, [
        r"q_group\s*=\s*q",
        r"qp_group\s*=\s*qp",
        r"cell_loop:\s*do",
        r"if\s*\(step_failed\)\s*then",
        r"if\s*\(with_exchange\s*\.and\.\s*topo_change_flag\)",
        r"call\s+assign_topography_vertices",
        r"state%deposit\(j,k,1:n_solid\)\s*=",
        r"call\s+state%refresh_substrate_totals",
        r"q\s*=\s*q_group",
        r"qp\s*=\s*qp_group",
    ])
    assert "CALL initialize_dry_physical_cell(qp_candidate)" in executor
    assert "CALL initialize_dry_physical_cell(qp_group(:,j,k))" in executor
    print("PASS: D19 Lie post-spatial guards, group publication and rollback")


if __name__ == "__main__":
    check()
