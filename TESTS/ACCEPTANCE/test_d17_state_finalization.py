#!/usr/bin/env python3
"""D17.19-10a source-level regression guard for post-RK finalization ordering.

This is deliberately a structural test, not a numerical equivalence test.
The observer must see raw conservative assembly, before any repairs.
"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
state = (ROOT / "src/state/state_2d.f90").read_text()
time = (ROOT / "src/numerics/time_integration_2d.f90").read_text()
step = (ROOT / "src/numerics/step_source_executor_2d.f90").read_text()


def body(text, name):
    pattern = rf"(?is)\bsubroutine\s+{name}\s*\(.*?\n(.*?)\bend\s+subroutine\s+{name}\b"
    match = re.search(pattern, text)
    assert match, f"missing routine {name}"
    return match.group(1)


def ordered(text, patterns):
    locations = []
    for expression in patterns:
        match = re.search(expression, text, re.I | re.S)
        assert match, f"missing invariant: {expression}"
        locations.append(match.start())
    assert locations == sorted(locations), "finalization steps are out of order"


def check():
    final = body(state, "finalize_imex_cell_state")
    ordered(final, [
        r"ISNAN\(q_cell\)",
        r"call\s+clip_roundoff_negative_component_masses\(q_cell\)",
        r"negative_thickness_check:",
        r"negative_alpha_check:",
        r"call\s+qc_to_qp\(q_cell\s*,\s*qp_cell",
        r"call\s+project_dry_cell_momenta",
        r"qp_cell\(4\)\s*\.LT\.\s*273\.0_wp",
        r"if\s*\(\s*SUM\(q_cell\(5:4\+n_solid\)\)",
        r"IF\s*\(\s*nodata\s*\)",
    ])
    for term in ("-1.0E-7_wp", "1.0E-10_wp",
                 "qc_to_qp(q_previous", "solid_excess_roundoff"):
        assert term.lower() in final.lower(), f"missing legacy behavior: {term}"

    rk = body(time, "imex_RK_solver")
    ordered(rk, [
        r"q\(1:n_vars,j,k\)\s*=\s*q_old_cell\s*-\s*dt\*residual_cell",
        r"CALL\s+observer\(3,\s*0,\s*j,\s*k,\s*q\(:,j,k\),\s*0\)",
        r"CALL\s+finalize_imex_cell_state\(q\(:,j,k\)",
    ])
    assert "CALL initialize_dry_physical_cell(this%qp_rk(:,j,k))" in rk
    assert "CALL initialize_dry_physical_cell(this%qp_rk(1:n_vars+2,j,k))" in rk

    dry = body(state, "initialize_dry_physical_cell")
    ordered(dry, [r"qp_cell\s*=\s*0\.0_wp", r"qp_cell\(4\)\s*=\s*T_ambient"])
    assert step.count("CALL initialize_dry_physical_cell(") == 2
    print("PASS: D17 post-RK observer, repair order, dry halo and STEP guards")


if __name__ == "__main__":
    check()
