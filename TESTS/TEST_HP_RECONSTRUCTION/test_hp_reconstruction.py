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
for expected in ("trace(2) = trace(1)*trace(idx_u)",
                 "trace(3) = trace(1)*trace(idx_v)",
                 "froude_h_blend_weight", "reconstruct_kt_eta_traces"):
    assert expected in r, expected
for expected in ("CASE ('h_x')", "CASE ('minmod')", "CASE ('froude_blend')",
                 "qp_adv_left_use(2)=h_eta_left*qp_adv_left(idx_u)",
                 "qp_adv_left_use(3)=h_eta_left*qp_adv_left(idx_v)"):
    assert expected in f, expected
print("PASS: only selected-depth velocity/momentum reconstruction remains")
