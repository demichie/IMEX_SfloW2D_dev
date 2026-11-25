#!/usr/bin/env python3
"""Static D19/D20 guard for admissible cell candidates and atomic Lie publication."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]


def load(path):
    return (ROOT / path).read_text(encoding="utf-8")


def body(text, name):
    opening = re.search(
        rf"(?im)^[ \t]*(?:pure[ \t]+)?subroutine[ \t]+{re.escape(name)}[ \t]*\(",
        text)
    assert opening, f"missing routine {name}"
    closing = re.search(
        rf"(?im)^[ \t]*end[ \t]+subroutine(?:[ \t]+{re.escape(name)})?[ \t]*$",
        text[opening.end():])
    assert closing, f"missing END SUBROUTINE for {name}"
    return text[opening.end():opening.end() + closing.start()]


def ordered(code, patterns):
    positions = []
    for pattern in patterns:
        match = re.search(pattern, code, re.I | re.S)
        assert match, f"missing invariant: {pattern}"
        positions.append(match.start())
    assert positions == sorted(positions), "Lie execution or publication order is incorrect"


def check():
    composition = load("src/state/state_composition.f90")
    integrator_source = load("src/numerics/time_integration_2d.f90")
    time = body(integrator_source, "advance_lie")
    spatial_cell = body(integrator_source, "update_spatial_cell")
    executor = body(load("src/numerics/lie_source_executor_2d.f90"),
                    "execute_lie_source_groups")

    conservative = re.search(
        r"(?is)logical\s+function\s+conservative_mass_admissible\b.*?"
        r"end\s+function\s+conservative_mass_admissible", composition)
    assert conservative, "missing conservative mass admissibility guard"
    assert "idx_add_gas_first:idx_add_gas_last" in conservative.group(0)
    assert "idx_solidEqn_first:idx_solidEqn_last" in conservative.group(0)
    assert "SUM(q_cell(idx_add_gas_first:idx_add_gas_last))" in conservative.group(0)

    # Spatial candidates are validated before entering the Lie source sweep.
    assert re.search(
        r"conservative_mass_admissible\(q_candidate\(1:n_eqns\)\)",
        spatial_cell, re.I)
    assert "CALL initialize_dry_physical_cell(qp_candidate)" in spatial_cell
    assert "qp_pre_all=>state%qp" in time
    assert "qp_pre_all=>this%qp_before" in time
    assert "this%cell_status(l)=0" in time

    # There is no redundant full-domain conservative snapshot or retry rollback.
    for obsolete in ("this%q_before", "this%q_trial=state%q",
                     "this%qp_trial=state%qp", "state%q=this%q_trial",
                     "state%qp=this%qp_trial", "rollback_geometry_transaction",
                     "restore_substrate"):
        assert obsolete not in time, f"retired transaction path: {obsolete}"

    # The fatal failure check must precede authoritative publication.
    ordered(time, [
        r"call\s+execute_lie_source_groups\(",
        r"if\s*\(source_step_failed\)\s*call\s+fatal_error",
        r"state%q\(:,:,j,k\)\s*=\s*this%q_trial\(:,:,j,k\)",
        r"state%qp\(:,:,j,k\)\s*=\s*this%qp_trial\(:,:,j,k\)",
        r"call\s+commit_geometry_transaction",
    ])
    assert "IF (needs_geometry_update) CALL begin_geometry_transaction" in time

    # Each OpenMP worklist cell owns the candidates of its configured layers.
    for fragment in ("q_layers=q(:,:,j,k)", "qp_layers=qp(:,:,j,k)",
                     "qpentry=qp_original(:,layer,j,k)",
                     "q_local=q_candidate", "qp_local=qp_candidate"):
        assert fragment.lower() in executor.lower(), f"missing cell-local Lie contract: {fragment}"
    ordered(executor, [
        r"cell_loop:\s*do",
        r"group_loop:\s*do",
        r"if\s*\(cell_error==0\)",
        r"q\(:,:,j,k\)\s*=\s*q_layers",
        r"qp\(:,:,j,k\)\s*=\s*qp_layers",
        r"if\s*\(first_failure/=huge\(0\)\)",
        r"call\s+apply_local_bed_rate",
        r"call\s+state%refresh_substrate_totals",
    ])
    assert "CALL initialize_dry_physical_cell(qp_candidate)" in executor
    print("PASS: D19/D20 layer-local Lie guards, fatal rejection and atomic publication")


if __name__ == "__main__":
    check()
