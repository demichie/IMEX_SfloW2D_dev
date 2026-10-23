# D17-C analytic stage bridge

Status: **implemented; local GNU Fortran validation passed**. This phase
connects the validated D17 stage groups to the existing analytic physical
kernels. It does not add native input syntax or route explicit, boundary, or
post-RK mass-exchange evaluation through the registry.

## Implementation

- `equation_terms_2d` now exposes `eval_analytic_terms` and a read-only cell
  context for bed slopes, curvature, gravity, and the frozen stochastic value.
  The bridge accepts the physical interval, preserves the full conservative
  state, and updates only the velocity, temperature, or pore-pressure target.
- `nonlinear_solver_2d` consumes the normalized stage-group order and analytic
  solver descriptors instead of dispatching analytic groups from the legacy
  equation partition. Direct calls to rheology, thermal, and pore-pressure
  analytic laws have moved behind the bridge. The friction-momentum projection
  is also accessed through `equation_terms_2d`.
- The simulation-owned registry is checked and copied once at the start of an
  IMEX advance, before cell loops. Newton groups still pass their complete
  equation block to the existing joint Newton solve. Semi-implicit
  initialization, post-Newton projection, rollback, status codes, analytic
  velocity → temperature → pore-pressure order, and the effective IM increment
  remain unchanged.
- Production entry points, standalone test runners, and the source dependency
  graph now include the registry metadata required by stage dispatch.

## Local validation

GNU Fortran 15.2.0 on macOS ARM64:

| Check | Result |
|---|---|
| Rheology/cell-local suite, strict and optimized | 255/255 fixtures pass per profile; analytic velocity, cooling, pressure, and grouped-stage checks pass |
| IMEX stage integration, strict and optimized | Exact-flow and mixed ANALYTIC/NEWTON stage traces pass for `N_RK=2/3/4`; group-order and reconstruction-boundary checks pass; invalid stage-count cases are rejected as expected |
| D17 source registry, strict and optimized | 864 legacy configurations, 14,857 assertions, and 28 rejected invalid configurations per profile |
| Direct `equation_terms_2d` runners | State conversion, initialization, workspace lifecycle, stochastic evolution, curvature, and Froude regressions pass |
| Production entry point | `IMEX_SfloW2D.f90` compiles against the updated integration API |
| Whitespace validation | `git -c core.whitespace=cr-at-eol diff --check` (preserves the existing CRLF source file) |

The IMEX test harness's Python comparison step requires Python 3.10 or newer
and NumPy. The Python 3.11 Conda environment was used. Initial attempts with the
workspace's Python 3.9 environment and a system Python lacking NumPy did not
run the comparison step; no repository or dependency files were changed for
those environment issues.

No Intel ifx, Linux, CI, or portability acceptance is claimed. The subsequent
numbered implementation of the explicit stage interface is recorded in
[D17.9 explicit stage routing](D17_9_EXPLICIT_STAGE_ROUTING.md), and the
group-filtered semi-implicit/implicit rate interfaces in
[D17.10 group-filtered implicit rates](D17_10_GROUP_FILTERED_IMPLICIT_RATES.md).
Multiple Newton-solver execution, remaining stage/step routing, native input
and full acceptance remain open.
