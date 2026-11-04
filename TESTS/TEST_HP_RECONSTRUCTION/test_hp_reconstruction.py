#!/usr/bin/env python3
"""Guard removal of legacy HP weights and preserve the selected-depth production path.

Physical regression is exercised by TEST_SPATIAL_OPERATOR and TEST_PCCU_GATE_H.
"""
from pathlib import Path
root = Path(__file__).resolve().parents[2]
r = (root / "src/spatial/reconstruction_2d.f90").read_text()
f = (root / "src/spatial/hyperbolic_2d.f90").read_text()
b = (root / "src/Makefile.am").read_text()
assert not (root / "src/spatial/hp_reconstruction_2d.f90").exists()
for s in (r, f, b):
    assert "hp_reconstruction_2d" not in s
for old in ("reconstruct_hp_line", "evaluate_hp_cell", "commit_hp_cell",
            "hp_blended", "hp_scratch", "hp_dynamic_residual_threshold"):
    assert old not in r, old
for expected in ("discharge_bc_to_velocity", "enforce_boundary_discharge",
                 "froude_h_blend_weight", "reconstruct_kt_eta_traces"):
    assert expected in r, expected

# The regular 2-D interior must use contiguous component stencils. The
# boundary path retains the scalar, condition-aware reconstruction.
interior = r.split("IF ( regular_interior ) THEN", 1)[1].split(
    "vars_loop:DO i=1,n_vars", 1
)[0]
for component in ("center_state", "west_state", "east_state",
                  "south_state", "north_state"):
    assert component in interior, component
assert "DO CONCURRENT (i=1:n_vars)" in interior
assert "reconstruction_cell_value" not in interior.split("DO CONCURRENT (i=1:n_vars)", 1)[1], (
    "interior limiter still fetches scalar neighbours repeatedly"
)
assert "CALL limit([west_state(i),center_state(i),east_state(i)]" in interior
assert "CALL limit([south_state(i),center_state(i),north_state(i)]" in interior
geometry = (root / "src/geometry/geometry_2d.f90").read_text()
for pure_declaration in ("PURE SUBROUTINE limit(", "PURE REAL(wp) FUNCTION minmod(",
                         "PURE REAL(wp) function maxmod("):
    assert pure_declaration in geometry, pure_declaration

# Depth slopes in the hydrostatic reconstruction are required only by the
# minmod and Froude-blended transport-depth selections.
assert r.count("kt_trace_reconstruction=='minmod' .OR.") >= 2
for expected in ("CASE ('h_x')", "CASE ('minmod')", "CASE ('froude_blend')",
                 "qp_adv_left_use(1)=h_eta_left",
                 "CALL qp_to_qc(qp_adv_left_use"):
    assert expected in f, expected
print("PASS: selected-depth primitive reconstruction with discharge boundary conditions")
