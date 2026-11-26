# D17.11 — REAL/COMPLEX kernel reuse

## Scope

The lahar carrier-fluid viscosity correlation originally used a COMPLEX-only
implementation whose REAL overload converted inputs to COMPLEX and projected
the REAL part. This avoided duplicate formulas but imposed COMPLEX
transcendental arithmetic on REAL mass-exchange and rheology evaluations.

The two overloads now instantiate a single
`src/physics/rheology_laws/rheology_lahar_viscosity.inc` expression block.
The REAL overload uses exclusively REAL intermediates, while the COMPLEX
overload retains complex-step sensitivity. Temperature branching, constants,
algebraic expression order, public signatures and generic
`eval_lahar_fluid_viscosity` remain unchanged.

The separate `rheology_lahar_drag.inc` continues to share the momentum
drag expression across REAL and COMPLEX source kernels.

## Verification status

The original D17.11 implementation passed strict and optimized
`TESTS/TEST_RHEOLOGY_CELL_LOCAL/run_test.sh` regressions, including
255/255 rheology fixtures per profile, as well as the
`TESTS/TEST_IMEX_STAGES/run_test.sh` suite. Those results describe the
earlier implementation, not automatic validation of this replacement.

For this change, the source-module regression covers matching REAL/COMPLEX
viscosity values in both temperature branches and complex-step derivatives
with respect to both temperature and solid fraction. An isolated GNU Fortran
harness also confirmed matching COMPLEX values against the previous
formula and matching REAL values in optimized and debug builds.
Full-repository and cross-platform regressions still require CI execution.
