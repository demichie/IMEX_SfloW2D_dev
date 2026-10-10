# IMEX stage admissibility and temporal order audit

The covered production IMEX states remain admissible and reproduce the
independent stage calculations for `N_RK=2,3,4`. Temporal refinement confirms
order one for two stages and order two for three/four stages, including a
coupled explicit-transport/implicit-drag problem. However, `N_RK=1` skips its
intended explicit update. That defect is reproduced, not corrected or accepted.
The combined CFL policy also remains a scientific decision; N9 stays open.

## Changes and observation contract

[`time_integration_2d.f90`](../../src/time_integration_2d.f90) has an optional
read-only `observer` argument on `advance`. Its three events expose the known
conservative stage before primitive conversion, the solved stage, and the
raw final assembly before clipping or component repair. A solved-stage status
distinguishes no local solve, convergence and failure. Calls occur inside the
existing OpenMP loops, so callbacks must be thread-safe and must not mutate
configuration or states. The test writes only each callback's unique cell.

Normal production callers omit the callback. There are no new diagnostic
allocations in the integrator, and no tableau, arithmetic update, nonlinear
tolerance, CFL formula or input default is changed. The `n_RK` declaration
comment now correctly identifies a stage count rather than an order.

## Stage comparisons and admissibility

[`TEST_IMEX_STAGES`](../../TESTS/TEST_IMEX_STAGES/run_test.sh) compiles actual
production modules in strict/debug and optimized profiles, with floating-point
traps and OpenMP in both. Its assertion driver remains strictly compiled.
Eight observer-off comparisons per profile cover resting and dynamic states
with every supported stage count; the returned conservative/primitive fields
must be byte-identical with and without observation.

Four stage fixtures cover rough-bed equilibrium, smooth dynamics, rough-bed
drainage and compact wet support. They use the same pinned, unmodified Python
spatial cores as the
[spatial-operator audit](SPATIAL_OPERATOR_2026_10_09.md). Each fixture runs with
slope/curvature both disabled or both enabled, four stage counts, directional
minimum or harmonic timestep, two profiles, and actual teams of one/four
threads. This gives 128 thread-paired configurations, or 256 executions.
Worksets are built by the real domain routine and depend on the stage count.
Any independent reference state outside the active workset must remain zero.

The independent tableaux assemble reference states from the pinned spatial
RHS, rather than reading private production coefficients. Every conservative
component is compared separately: the much larger thermal-energy scale cannot
hide momentum errors. The fixed comparison limit is
`32768*epsilon(float64)*max(1,max(abs(reference_component)))`.
The largest normalized discrepancy is `9.36e-13`.

Raw known, solved and final states are checked for finite values, nonnegative
thickness and component/carrier masses, and admissible temperature. All
observed final repair increments are zero and no local Newton solve reports
failure in these fixtures. Closed-collar mass conservation and lake-at-rest
continuation are also checked. Seven deliberately corrupted observations or
admissibility statistics must fail the comparator in each profile.

The closure is liquid with one zero solid component. These tests are not a
positivity theorem, a variable-composition transport gate, a characterization
of every rheology, or a validation of explicitly time-dependent inlet sources.

## Temporal refinement and implicit solves

The autonomous order fixtures have an 8 by 7 rectangular grid, first-order
spatial reconstruction, fixed spatial resolution and final time `0.32 s`.
Step counts are 8, 16, 32, 64 and 128, all within the initial production timestep
bound. Keeping the spatial operator fixed isolates temporal convergence.

Thermal transport uses a constant-density moving liquid and varying
temperature. Its independent oracle is the exact matrix exponential of the
linear central-upwind thermal operator, including Neumann boundary rows.
A resting temperature perturbation is not used as an order oracle: HP's
equilibrium detection intentionally suppresses scalar fluxes at rest.

The drag oracle is the existing quadratic-friction law, model 6, with
`h=1 m`, `friction_factor=1` and initial velocity `(0.6,0.8) m/s`. Its exact
velocity is the initial vector divided by `1+t`. The mixed oracle evolves the
independent thermal operator with this analytic velocity using RK4; agreement
between 2048 and 4096 oracle steps is checked before comparison.

| Stage count | Thermal transport order | Quadratic drag order | Mixed order |
| --- | --- | --- | --- |
| 2 | 1.0015 | 0.9996 | 1.0010 |
| 3 | 2.0044 | 2.0008 | 2.0039 |
| 4 | 2.0067 | 2.0005 | 2.0087 |

These are the finest-pair slopes; strict and optimized values agree at the
displayed precision. All eighteen order series pass the unchanged target
window of expected order plus/minus `0.25` for the last two refinement pairs.
Temperature errors use a `1 K` scale, momentum errors a
`rho_l*(1 m/s)` scale, and the mixed test takes the maximum of the two.
Independent first/second-order and mixed colored-tree tableau conditions are
checked as well.

Nonzero implicit drag stages are additionally checked against their quadratic
equation and independent algebraic root, including known/final momentum
assembly. Production Newton convergence permits either its residual criterion
or a normalized correction below `tol_rel=1e-5`; `tol_abs` is also `1e-5`.
Some captured stages converge through the correction criterion without meeting
the residual-only test. This is recorded explicitly and is not mislabeled as
a machine-precision implicit solve. Neither tolerance is tightened by the test.

## Reproduced single-stage defect

The `N_RK=1` branch stores an explicit diagonal coefficient of one, with
explicit final weight one and zero implicit tableau/weights. Stage assembly
uses only earlier-stage coefficients, so its single stage is the initial
state. The final stiff-accuracy test nevertheless compares the stored diagonal
against the final weight, declares a match and returns that unchanged stage
instead of applying the evaluated explicit residual.

The dynamic and thermal tests reproduce the no-op with both compiler profiles
and both thread counts. The intended forward-Euler result is retained alongside
the observed result in
[`results_imex_stages_2026_10_10.json`](results_imex_stages_2026_10_10.json).
The default suite checks this known behavior as characterization; it does not
declare `N_RK=1` scientifically valid. The separate
`--require-design-order` gate fails with an assertion naming the defect.

The next numerical correction should remove this erroneous assembly selection
and update the regression to require Euler advancement. Whether a one-stage
option should integrate implicit physics or reject such configurations is a
separate policy choice: its existing implicit weights are zero.

## CFL interpretation

At the same initial state, the retained directional-minimum step and the
reference harmonic step are different parameterizations. Stage-dependent
wave bounds are evaluated independently for every observed reference stage.
Only stages whose explicit RHS actually participates should be interpreted
as explicit CFL checks; unused final stages are recorded separately.
The defective one-stage branch discards its computed RHS, so it is excluded
from these participating-stage statistics.

For the retained step, the largest ratio to a directional minimum in an
explicitly used stage is `1.194`, in the drainage case with three stages.
Its largest combined Courant number over used stages is about `0.4473`.
The corresponding maxima across harmonic-policy cases are `0.6677` and
`0.2532`. Thus even the harmonic step can exceed its initial combined target
`0.24` after stage velocities change. Neither policy produced an inadmissible
state in this finite matrix, but neither result proves general admissibility.

No timestep is changed midway through a tableau. A future admissibility/retry
policy must reject and recompute the whole step if needed, and should be tested
with the actual IMEX method rather than inheriting an SSPRK2 proof. Changing
the CFL formula alone would not resolve all these questions.

## Evidence and reproduction

The candidate is based on `84ef3018`, in the actual repository
`/Users/demichie/Codes/GIT/IMEX_SfloW2D_dev`. The JSON artifact retains the full
acceptance manifest, stage/order records, input/output fingerprints, the
design-gate failure, and comparisons with the preceding baseline. Observed
hashes are not approved golden references.

From the repository root:

```sh
sh TESTS/TEST_IMEX_STAGES/run_test.sh
```

Use `KEEP_TEST_WORKDIR=1` to retain generated fixtures, executables and payloads.
The test requires gfortran, LAPACK, Python 3 with NumPy, `nf-config` and
`nc-config`. With a retained executable, the comparator's `--section order`
isolates temporal tests; adding `--require-design-order` reproduces the expected
scientific-gate failure. Run it from a fresh temporary directory.

Full audit logs remain in `/tmp/imex-stages-acceptance-20261010-final` and the
separate design-gate log in `/tmp/imex-stages-design-gate-20261010.log`.
These are temporary, not permanent raw-output archives. The
[acceptance plan](../../TESTS/ACCEPTANCE/acceptance_plan.json) keeps the
single-stage defect, CFL policy, nonautonomous-source and mixed-composition
requirements open before N9 or the layer migration.

The complete audit passes 20 unit/documentation scripts and 16 full-solver
runs, 36 in total; seven audit-tool self-tests and all six Gate-H/restart
cross-thread pairs also pass. All 16 full-solver cases retain identical
canonical q/checkpoint fingerprints against the preceding spatial-operator
baseline. These successful regression/characterization checks coexist with
the failed design-order gate and do not override it.

For a fresh full audit of the uncommitted candidate:

```sh
python3 TESTS/ACCEPTANCE/run_acceptance.py \
  --working-tree --revision HEAD \
  --output /tmp/imex-stages-repeat \
  --netcdf /opt/homebrew/opt/netcdf-fortran --jobs 3
```

After committing, use the resulting SHA without `--working-tree`.
