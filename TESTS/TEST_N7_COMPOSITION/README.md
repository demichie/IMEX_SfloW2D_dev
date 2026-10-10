# N7 variable composition equilibrium and transport regression

The full acceptance inventory checks 90 equilibria and 72 signed moving
contacts per compiler profile, with actual one/four-thread teams. It uses
independent thermodynamic identities, analytic translated cell averages and a
scalar CU/TVD reference with the current IMEX weights. CFL and the clean reference
kernels are unchanged. The complete matrix passes in the
[geometry and transport completion](../../docs/acceptance/N7_B_GEOMETRY_AND_TRANSPORT_2026_10_10.md).
The original failed contact and its diagnostic evidence
remain archived separately; they are not relabeled as passing tests.

## Run the full acceptance inventory

From the repository root:

```sh
KEEP_TEST_WORKDIR=1 sh TESTS/TEST_N7_COMPOSITION/run_test.sh
```

The script requires all 162 cases in each profile before writing `evidence.json`.
Any failure exits nonzero and writes `failure.json`; a partial diagnostic or
equilibrium-only execution cannot close N7-B. The general acceptance plan runs
this full gate as well as the separate stationary-roundoff controls.

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

The original equilibrium construction and limits are unchanged. See the
[historical diagnosis](../../docs/acceptance/N7_B_DIAGNOSIS_2026_10_10.md) and
[stationary correction](../../docs/acceptance/HYDROSTATIC_ROUNDOFF_2026_10_10.md)
for the preceding failure and its bounded production fix.

## Run the partial investigation

```sh
KEEP_TEST_WORKDIR=1 sh TESTS/TEST_N7_COMPOSITION/run_test.sh --diagnose
```

This mode retains the original contact construction in `contract.json`, not
the revised acceptance contact. It builds strict and historical optimized modules and
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
72 contacts. The current contact construction is separately fixed in
[`contact_contract_v3.json`](contact_contract_v3.json); the original contract
is not overwritten, so the historical evidence and equilibrium hashes remain
valid. The failed intermediate construction remains in
[`contact_contract_v2.json`](contact_contract_v2.json).

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

The revised moving contact translates the same pulse and the complete 16 m
comparison interval into a 480 m domain. Longitudinal grids have 480/960/1920
cells along x and 960 cells along y, preserving the original spacings 1/0.5/0.25 m,
time 0.1 s, timestep rule, limiter and all acceptance limits. The velocity
transitions occupy 110--115 m and 365--370 m. No comparison cells are removed.

Before launching each solver, an index-space check requires both the comparison
interval and physical boundaries to lie outside a conservative dependency cone
of `4*n_RK*steps` cells. Fully wet flat-bed reconstruction and fluxes have local
stencils; the prescribed timestep removes global timestep dependence. The
initially quiet boundary buffers must remain quiet throughout the test, which
still checks zero boundary flux and budgets over the entire domain. This is a
finite-time isolated transport fixture, not a general boundary-condition test.

The v2 attempt isolated the central interval but retained transitions too near
the physical boundaries; its measured nonzero boundary flux failed the unchanged
criterion. V3 extends the quiet boundary buffers as well. Each revised
construction was frozen before its execution; no tolerance was enlarged.

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
Six checker unit tests exercise equilibrium identities, each component,
analytic signed transport, deliberate diagnostic corruption and the fixed
roundoff/reference contract, plus preservation of resolution/limits and rejection
of an insufficient dependency-cone separation.

Only a full passing inventory establishes N7-B within its defined fixtures.
D-N7-CFL and N9 remain separate decisions.
