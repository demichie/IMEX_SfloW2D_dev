# N7 variable composition equilibrium and transport diagnostics

The 90-case equilibrium subset now passes after a bounded correction in
`hyperbolic_2d`. N7-B remains open because the moving-contact reference needs a
better isolated fixture. This package does not certify the full planned
equilibrium/transport inventory. The CFL and clean reference kernels are
unchanged. See the [historical diagnosis](../../docs/acceptance/N7_B_DIAGNOSIS_2026_10_10.md)
and [correction report](../../docs/acceptance/HYDROSTATIC_ROUNDOFF_2026_10_10.md).

## Run the verified equilibrium subset

From the repository root:

```sh
KEEP_TEST_WORKDIR=1 sh TESTS/TEST_N7_COMPOSITION/run_test.sh --equilibrium-only
```

This checks all original 90 equilibria in both profiles, with unchanged states
and tolerances, and writes `equilibrium_evidence.json`. The separate
`TEST_HYDROSTATIC_ROUNDOFF/run_test.sh` also contrasts 36 slow-motion, tangential
and real-pressure-force controls against the checksum-pinned old hyperbolic
module, using identical other current objects. Both scripts retain actual
one/four-thread binary comparisons.

Without an option, `run_test.sh` attempts the original full inventory; its
contact comparison is still a draft and must not be treated as a valid solver
acceptance gate. It exits nonzero on the mismatch and stores `failure.json`.
The pre-correction equilibrium failure is archived in the historical diagnosis,
not relabeled as a failure of the corrected equilibrium subset.

## Run the partial investigation

```sh
KEEP_TEST_WORKDIR=1 sh TESTS/TEST_N7_COMPOSITION/run_test.sh --diagnose
```

This mode builds both strict and historical optimized production modules and
runs 36 diagnostic cases per profile with actual one/four-thread teams. It
records failing criteria in `diagnosis.json` without stopping at the first
failure. A successful diagnostic command only means the investigation ran;
it does not mean the numerical gates passed.

The diagnostic subset comprises both 1D equilibrium families with slope and
curvature disabled: three closures, one/two-step RK2 prefixes and 20-step
RK2/3/4 trajectories. Six moving-contact probes cover both signs and all three
closures at 80 cells with RK3. The 2D/flag matrix, other contact refinements and
rotation are outside this diagnostic subset. Complete `result.bin` and
`composition.bin` payloads must match bitwise between thread counts.

## Frozen construction and criteria

[`contract.json`](contract.json) fixes the numerical cases and limits before
execution. Its 162-case inventory per profile consists of 90 equilibria and
72 contacts. Only the equilibrium subset is verified. The general acceptance
plan runs it through `TEST_HYDROSTATIC_ROUNDOFF`; it does not run the draft
contact comparison as a passing gate.

The three current closures use two solids with distinct densities and heat
capacities: liquid carrier; ambient-air carrier with one added gas; and the
gas-containing mixture with explicit liquid. Independent specific-volume and
heat-capacity sums construct and decode states without clipping. Each transported
component and the residual liquid/air carrier has its own inventory budget.

For the Q1-bed equilibrium, the varying fractions lie in the constant-density
specific-volume nullspace at uniform temperature. Constant free surface is
valid because Gamma is then constant. The second, flat-bed family deliberately
has varying Gamma and depth, with `Gamma*h**2` constant. Constant free surface
alone is not its equilibrium condition. That family uses piecewise-constant
traces so the endpoint pressure identity applies directly.

Equilibrium, raw face fractions, TVD bounds and component budgets retain the
fixed `32768*float64_epsilon` limits. Moving-contact field comparisons retain
that limit multiplied by the step count, plus the frozen analytic-error and
refinement gates. No criterion is relaxed after a failure.

The moving contact uses an independently integrated scalar CU/TVD reference
with the current explicit IMEX weights, not the one-solid prototype as a gas
reference. Its prescribed velocity tapers to zero away from the pulse for
closed boundaries. The current comparison interval is insufficiently isolated
from that taper at roundoff-level precision; this is a fixture/reference defect,
not evidence of a production contact defect. A revised reference contract must
preserve this failed construction and be justified before a new run.

## Observation format

The shared `TEST_IMEX_STAGES` driver has an optional sixth argument selecting
the closure. Omission or zero preserves its original pure-liquid layout and
binary payload. Nonzero modes append no data to `result.bin`; instead,
`composition.bin` separately stores all-step component budgets and extrema,
last-stage raw cell-face traces, conservative interface states and the spatial
term. The spatial term has the production `+div(flux)` sign and is subtracted
by the IMEX assembly.

Known, solved and raw-final states are observed before final repair. A separate
test-only workspace reconstructs each solved stage, checks fractions before
thermodynamic closure and independently verifies the final conservative faces.
It does not change the production spatial operator or per-cell source laws.
Five checker unit tests exercise equilibrium identities, each component,
analytic signed transport, deliberate diagnostic corruption and the fixed
roundoff/reference contract.

N7-B, D-N7-CFL and N9 remain unresolved by these diagnostics.
