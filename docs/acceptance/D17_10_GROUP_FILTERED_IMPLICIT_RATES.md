# D17.10 — Group-filtered semi-implicit and implicit rates

## Scope

The semi-implicit and implicit equation-term evaluators now accept a Newton
group and its selected solver descriptor. They include only terms registered
in both descriptors. Existing calls without descriptors retain the legacy
whole-source behavior for compatibility; the production Newton stage passes
the descriptors through rate evaluation, residual assembly, complex-step
Jacobian evaluation, and line search.

The semi-implicit REAL path selects rheology terms by role. The implicit REAL
and COMPLEX paths select the matching rheology, the complete radiative/
convective cooling bundle, or pore-pressure dissipation. Incomplete cooling
bundles and mismatched group/solver descriptors fail explicitly.

The Newton stage coordinator executes each group's registered solvers in
execution order. Each Newton solve sees only its selected rate terms and
equations; a legacy coupled block remains one solver over its full equation
set. Semi-implicit momentum contributions are applied only to the current
solver's equations before that solve.

## Local verification

- `TESTS/TEST_RHEOLOGY_CELL_LOCAL/run_test.sh`: strict and optimized profiles
  pass; 255/255 rheology fixtures per profile. The grouped-stage regression
  checks term isolation for rheology, cooling, and pore pressure, including
  REAL and COMPLEX implicit rates and the semi-implicit REAL rate. It also
  executes three ordered Newton solvers in one group and checks the thermal
  and pore-pressure residuals against the final sequential state. The existing
  joint Newton regression continues to check the coupled legacy block.
- `TESTS/TEST_SOURCE_REGISTRY/run_test.sh`: strict and optimized profiles pass
  for 864 legacy configurations, 14,857 assertions, and 28 rejected invalid
  configurations per profile.
- `TESTS/TEST_IMEX_STAGES/run_test.sh`: full strict and optimized stage,
  boundary, CFL, exact-flow, mixed-group, and temporal-order checks pass.

These are local GNU Fortran results. No Intel ifx, Linux, CI, or cross-platform
acceptance is claimed. The full D17.20 acceptance matrix remains open.
