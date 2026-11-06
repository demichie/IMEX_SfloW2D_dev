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

# The first face-depth estimate may be dry even when the accepted HP depth
# is wet. Zero provisional velocity must not be treated as an invalid inlet.
# At the final HP face a dry Neumann boundary has zero transport; a genuinely
# prescribed nonzero Dirichlet discharge remains an explicit input error.
first_pass = r.split("SUBROUTINE discharge_bc_to_velocity(", 1)[1].split(
    "END SUBROUTINE discharge_bc_to_velocity", 1
)[0]
final_pass = r.split("SUBROUTINE enforce_boundary_discharge(", 1)[1].split(
    "END SUBROUTINE enforce_boundary_discharge", 1
)[0]
assert "Nonzero boundary discharge at dry first-pass reconstruction" not in first_pass
assert "velocity_bc%value=0.0_wp" in first_pass
assert "discharge_bc%flag==0" in final_pass
assert "Nonzero Dirichlet discharge at a dry HP boundary" in final_pass
assert "trace(idx)=0.0_wp" in final_pass

# The regular 2-D interior must use contiguous component stencils. The
# boundary path retains the scalar, condition-aware reconstruction.
interior = r.split("IF ( regular_interior ) THEN", 1)[1].split(
    "vars_loop:DO i=1,n_vars", 1
)[0]
for component in ("center_state", "west_state", "east_state",
                  "south_state", "north_state"):
    assert component in interior, component
assert "uniform_generalized_minmod" in interior
assert "DO CONCURRENT (i=1:n_vars)" in interior
assert "qrec_prime_x=minmod(" in interior
assert "qrec_prime_y=minmod(" in interior
assert "reconstruction_cell_value" not in interior.split("DO CONCURRENT (i=1:n_vars)", 1)[1], (
    "interior limiter still fetches scalar neighbours repeatedly"
)
assert "CALL limit([west_state(i),center_state(i),east_state(i)]" in interior
assert "CALL limit([south_state(i),center_state(i),north_state(i)]" in interior
geometry = (root / "src/geometry/geometry_2d.f90").read_text()
for pure_declaration in ("PURE SUBROUTINE limit(", "PURE ELEMENTAL REAL(wp) FUNCTION minmod(",
                         "PURE REAL(wp) function maxmod("):
    assert pure_declaration in geometry, pure_declaration

# Depth slopes in the hydrostatic reconstruction are required only by the
# minmod and Froude-blended transport-depth selections.
assert r.count("kt_trace_reconstruction=='minmod' .OR.") >= 2
# Trace-mode dispatch lives in the reconstruction workspace. The
# selected-depth conservative conversion lives in the PCCU interface.
for expected in ("CASE ('h_x')", "CASE ('minmod')", "CASE ('froude_blend')"):
    assert expected in r, expected
for expected in ("qp_adv_left_use(1)=h_eta_left",
                 "CALL primitive_to_conservative(qp_adv_left_use"):
    assert expected in f, expected
print("PASS: selected-depth primitive reconstruction with discharge boundary conditions")
