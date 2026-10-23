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

## Local verification

`TESTS/TEST_SOURCE_REGISTRY/run_test.sh` passes in strict and optimized
profiles. Each profile covers 864 legacy configurations, 14,873 assertions
and rejects 28 invalid configurations. The added checks verify the configured
default schedule; they do not claim end-to-end numerical equivalence for every
source activation or complete D17.20 acceptance.

No new cross-platform or CI acceptance is claimed.
