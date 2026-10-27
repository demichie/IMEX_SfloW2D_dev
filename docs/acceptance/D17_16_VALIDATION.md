# D17.16 / D18 — Physical-source input and registry validation

## Current input contract

Physical sources are configured through their owning input namelists.
`METHOD`, `CADENCE`, and `GROUP` (as the corresponding physical-source fields)
determine the schedule. `source_configuration_2d%initialize_configured_sources`
constructs both the equation partition and the validated source registry.
There is no separate `SOURCE_GROUPS` input or numerical IMEX group assignment.

Methods are `EXPLICIT`, `ANALYTIC`, and `NEWTON`, subject to each source's
catalogued capabilities. Cadences are `STAGE` and `STEP`. For the current
implementation, physical velocity, thermal, and pore-pressure sources with
`CADENCE=STEP` are explicitly rejected; the enabled post-RK mass-exchange
batch uses `EXPLICIT/STEP`. Stage groups are ordered by method class and
group identity. Explicit source groups and implicit equation blocks are
validated before integration.

The `source_input_2d` reader rejects the removed `SOURCE_GROUPS` namelist.
The numerical input reader validates physical method, cadence, and group
settings. It does not translate the previously supported IMEX group fields.
Input templates and generators must use only the physical-source settings.

## Regression coverage

`TESTS/TEST_SOURCE_REGISTRY/run_test.sh` exercises the source catalog,
physical triplet combinations and invalid registry metadata in strict and
optimized compilation profiles. `TESTS/TEST_IMEX_STAGES/run_test.sh`
exercises actual input parsing, the production timestep-event sequence,
and analytic/Newton group trajectories. `TESTS/TEST_MASS_EXCHANGE/run_test.sh`
exercises atomic post-RK exchange, budgets and checkpoint continuation.

The original D17.16 acceptance report described a previous input contract.
Its former examples are not valid input for this branch. This page specifies
the current D18 behavior; previous CI evidence should be interpreted only
at the commits where it was collected.

## Remaining work

Complete the D18 `STEP` execution modes that are not yet implemented and
extend numerical acceptance across the new physical-source combinations.
Passing CI does not itself establish equivalence for untested combinations.
