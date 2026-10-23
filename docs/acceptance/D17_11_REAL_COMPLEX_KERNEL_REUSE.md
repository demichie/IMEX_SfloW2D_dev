# D17.11 — REAL/COMPLEX kernel reuse

## Scope

The lahar carrier-fluid viscosity correlation had separate REAL and COMPLEX
implementations of the same piecewise temperature law. The COMPLEX
implementation is now the single formula; the existing REAL generic overload
wraps it and returns its real component. Public call signatures and the
generic `eval_lahar_fluid_viscosity` interface remain unchanged.

No `.inc` file was needed for this single correlation. The D17.10 implicit
evaluator already computes through a shared COMPLEX path and projects the REAL
rate at its REAL interface; the lava cooling and pore-pressure diffusion laws
likewise do not maintain a separate REAL formula to consolidate here.

## Local verification

`TESTS/TEST_RHEOLOGY_CELL_LOCAL/run_test.sh` passes in strict and optimized
profiles, with 255/255 rheology fixtures per profile. The source-module
regression additionally compares REAL and COMPLEX viscosity values in both
branches of the piecewise temperature law and checks the complex-step
temperature derivative against a centered finite difference.

`TESTS/TEST_IMEX_STAGES/run_test.sh` also passes its full strict and optimized
stage, exact-flow, mixed-group, CFL, boundary and temporal-order regressions.

GNU Fortran strict and optimized results are local validation only. This
report does not claim a complete D17.20 acceptance matrix or cross-platform
acceptance.
