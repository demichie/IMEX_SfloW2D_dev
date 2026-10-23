# D17.13 — Legacy-compatible default configuration

## Scope

The configured-source bridge now has regression coverage for the compiled
legacy IMEX defaults:

- velocity group 3, energy group 1 and pore-pressure group 2;
- ANALYTIC, ANALYTIC and NEWTON methods for those configured group IDs;
- an inactive-physics case, which retains the zero-rate implicit momentum
  block, reserves explicit stage group 4 and registers no inactive `step`
  terms;
- an active lava/energy/pore-pressure case, which registers energy cooling,
  pore-pressure dissipation and the complete lava rheology bundle in their
  expected analytic/Newton groups and execution order.

The tests exercise the actual `initialize_configured_sources` bridge and
registry validation. The source-term inventory, activations, equation
footprints, representations and atomic mass-exchange batch remain defined by
the [D17 source contract](D17_SOURCE_CONTRACT_2026_10_23.md) and its
[machine-readable manifest](D17_SOURCE_TERMS_MANIFEST_2026_10_23.json).

The grouped-stage runtime regression also builds the registry through
`initialize_configured_sources` with the legacy defaults and runs the ordered
analytic energy, analytic pore-pressure and Newton velocity stage solves.
Existing exact-flow oracles verify the energy/pressure targets; the Newton
momentum residual and effective stage increment are checked against the
assembled source rates. The separate joint-Newton regression continues to
cover explicitly configured shared Newton groups.

## Local verification

`TESTS/TEST_SOURCE_REGISTRY/run_test.sh` passes in strict and optimized
profiles. Each profile covers 864 legacy configurations, 14,873 assertions
and rejects 28 invalid configurations. The added checks verify the configured
default schedule.

`TESTS/TEST_RHEOLOGY_CELL_LOCAL/run_test.sh` passes in strict and optimized
profiles, with 255/255 fixtures per profile, including configured-default
runtime solves. `TESTS/TEST_IMEX_STAGES/run_test.sh` passes in strict and
optimized profiles, including stage trajectories for N_RK=2/3/4, mixed
analytic/Newton groups, CFL and boundary regressions.

These tests do not claim end-to-end numerical equivalence for every source
activation or complete D17.20 acceptance.

No new cross-platform or CI acceptance is claimed.
