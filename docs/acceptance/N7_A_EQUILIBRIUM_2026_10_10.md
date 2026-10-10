# N7-A: equilibrium and one-dimensional dynamics

N7-A is satisfied in its frozen fixtures. The focused suite passes 408 production
solver runs: 102 cases per strict/optimized profile, each with actual one/four
thread teams. All 204 complete payload pairs match bitwise. Across these runs,
27360 time steps are observed. The full regression passes 37/37 tests and all
six Gate-H/deterministic-restart/stochastic-restart thread comparisons.

Two required regressions exposed production defects and justified bounded fixes.
Optimized dynamic outputs consequently change; this is not an output-neutral
cleanup. N7 as a whole and N9 scientific baseline approval remain open. N8 stays
satisfied. The next bounded package is N7-B, not layer storage or new physics.

The [machine-readable evidence](results_n7_a_2026_10_10.json) pins production and
numerical fixture hashes, the frozen contract, case inventory, execution
fingerprints, raw-evidence hashes, all audit results and baseline comparisons.
The [test README](../../TESTS/TEST_N7_EQUILIBRIUM/README.md) explains execution and
the [closure checklist](N7_N8_CLOSURE_2026_10_10.md) maps original requirements.
The base revision is `026061fe`; candidate hashes, not that base SHA alone,
identify the executed changes. Final documentation/closure metadata and
acceptance-tool assertions are added after numerical execution.

## Frozen cases and independent oracles

All cases use `N_RK=2,3,4`. The new script reuses the unchanged raw-stage driver
from `TEST_IMEX_STAGES`. The checksum-pinned clean 1D/2D spatial kernels remain
unchanged; independent adapters use the current IMEX assembly, not the prototype's
SSPRK2 stepper or retry. No CFL, stage-count, physical-law or global compiler-policy
change is introduced, and D-N7-CFL is not resolved by this execution.

| Family | Cases per profile/team | Grid and time | Oracle and original requirement |
| --- | --- | --- | --- |
| Lake | 72 | 32 cells in 1D; 32 x 24 in 2D; spacing 0.25; t=0.02 | Smooth, continuous one-cell and fully wet parabolic beds; all four slope/curvature pairs; component-wise unchanged state, raw stages and force-scaled residual. H03/H11. |
| Smooth advection | 9 | 64/128/256 cells over length 16; t=0.4 | Uniform moving liquid and Gaussian temperature pulse; independent linear scalar central-upwind operator plus analytic translating Gaussian. H01. |
| Ritter | 9 | 40/80/160 cells over length 20; t=0.2 | Exact cell integrals of the dry-bed rarefaction, including ambient-density gravity correction; separate clean-core finite-time field/stage oracle. H02. |
| Excavation | 12 | 60 cells over length 30; t=0.01 | Ten-degree slope with two depth-2 walls spanning one cell each; all four flags; finite-time reference, total mass, zero uphill transfer and positive downslope acceleration. H07. |

Smooth advection is passive thermal transport, not a composition-contact gate.
Its limiter-zero refinement test demonstrates decreasing spatial error, not
second-order spatial convergence. The already established IMEX temporal-order
tests are retained. H04/H05 variable-composition trajectories belong to N7-B.

The [contract](../../TESTS/TEST_N7_EQUILIBRIUM/contract.json) was fixed before the
first execution and was not changed after failures. With epsilon for float64,
raw thickness positivity, relative mass budget, uphill mass and lake equilibrium
use `32768*epsilon`. Finite-time field errors use
`32768*epsilon*steps` with separate initial mass, momentum, thermal and solid
scales; the thermal-energy scale never hides momentum or mass errors.
Advection requires L1 <= 2 K, Linf <= 8 K; Ritter requires mean thickness error
<= 0.12 and mean momentum/c0 error <= 0.15. Both require successive L1 ratios
below 0.9 on three refinements.

Known, solved and raw-final minimum mass, raw temperature/component fractions,
local solve status, used/computed CFL and mass budgets are checked every step.
Complete final and last-step raw arrays are compared separately. The checker
self-tests include deliberately corrupted fields and diagnostics.

## Raw states and existing dry cleanup

The rest-transport guard can leave a pressure impulse in a cell with exactly
zero mass. At the first coarse Ritter RK2 step the raw dry momentum is about
24.49557 in conservative units, before production canonicalizes that dry state.
The independent reference has the same impulse. It is not a raw negative-mass
failure and cannot be treated as physical velocity at zero mass.

The initial checker incorrectly prohibited every final repair. It now compares
the raw reference states before sanitization and applies the existing explicit
dry mapping only for mass below epsilon. Every step's cleanup magnitude must
match that independently predicted mapping. No negative-mass clipping is used
to demonstrate positivity and no acceptance tolerance was relaxed.

## Required production corrections

In `mixt_var`, refined Ritter generated a positive subnormal speed squared.
`MIN(cap, buoyancy_depth/speed_squared)` overflowed before the cap was applied;
the strict build stopped at that quotient. The existing 1e15 upper cap is now
tested before division. Ordinary states and the exact-rest convention Ri=0
retain their previous formula. A focused conversion regression exercises all
three branches with floating-point traps against both production profiles.

In HP momentum reconstruction, optimized fused multiply-add retained a tiny
discharge when an admissible face velocity bound was exactly zero. That residue
disabled the downstream exact-rest scalar-flux guard and created spurious mass
transport. The all-flag excavation test failed its pre-existing frozen field
limit: scaled mass error 2.3049e-8 versus 7.2760e-10.

Only the two face mean-momentum products are now stored in thread-private local
`VOLATILE` scalars before the final addition/subtraction. This preserves exact
cancellation without changing bounds, formulas, OpenMP or all-module compiler
flags. As a diagnostic, compiling only the old HP object without FMA removed
the discrepancy. The new exact-zero-bound unit test fails against the old
optimized object and passes with the bounded fix. Both signs are tested.

## Observed results

The maximum relative mass drift is 9.095e-16; uphill mass in the excavation is
exactly zero. Lake equilibrium error is at most 1.209e-16 and the inferred
production force-scaled residual at most 1.390e-16. All raw mass minima are
nonnegative in these cases and all local implicit solves succeed.

Advection L1 errors decrease from approximately 0.383 K to 0.217 K to 0.117 K;
successive ratios are approximately 0.568 and 0.536. Ritter thickness L1 errors
range from 0.00557--0.00621 on the coarsest grid to 0.00221--0.00265 on the finest;
all ratios are below 0.69. The largest component-scaled reference discrepancy is
approximately 1.74e-11, within the pre-frozen step-scaled allowance. Strict and
optimized executions independently pass these limits.

## Changed optimized baseline: characterization, not approval

Against N8-B, all eight strict full-solver fingerprints and the optimized lake
remain unchanged: 9/16 cases. Seven optimized dynamic case fingerprints change:
the two-thread configurations of Gate H, deterministic restart and stochastic
restart, plus the embedded excavation. Within the candidate revision all
restart and one/four-thread comparisons still pass.

An additional one-thread replay retains complete historical/candidate/strict
fields using unchanged physical inputs. These are descriptive measurements,
not newly fitted acceptance thresholds. For each field column the old/new
Linf difference is divided by `max(1, max(abs(old_column)))`; the maximum below
is over retained snapshots and columns, not a local relative velocity error.
Output columns include auxiliary/deposit states when the fixture writes them.

| Full-test family | Snapshots compared | Maximum normalized old/new Linf |
| --- | --- | --- |
| Deterministic restart | 27 | 9.180e-3 |
| Gate H | 12 | 1.582e-2 |
| Inclined excavation | 18 | 2.961e-12 |
| Stochastic restart | 22 | 1.171e-2 |

The candidate moves closer to strict in the deterministic/Gate-H aggregate
norms and removes the excavation discrepancy. This is not uniform cross-profile
identity: stochastic auxiliary output already has an old optimized/strict
normalized difference of approximately 1.323, which remains approximately
1.323 after this correction. The archived per-column values preserve that
distinction; no explanation or cross-profile scientific acceptance is inferred
from the restart test alone. Direct strict/optimized equality is not claimed.

Long EXAMPLE_ETNA and EXAMPLE_LAVA_IMO characterization should be repeated before
adopting candidate optimized outputs as a scientific baseline. The previous long
runs do not cover these two new production changes. Remaining N7-B/C/D gates,
D-N7-CFL and N9 approval are explicitly left open.
