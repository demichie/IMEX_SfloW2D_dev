# Source registry and D18 physical configuration

Run the registry tests from the repository root with GNU Fortran:

```sh
D17_WORK_DIR="$PWD/build-local/d17-registry/focused" sh TESTS/TEST_SOURCE_REGISTRY/run_test.sh
```

The runner checks the physical source catalog, source activation, validation of group
membership and capabilities, and 864 combinations of physical METHOD/GROUP
settings for velocity, cooling and pore pressure. The registry is constructed
through `initialize_configured_sources`, with CADENCE determining STAGE and
post-RK STEP scheduling. Invalid physical source entries are rejected.

`source_catalog_2d` declares the canonical physical sources and their supported
interfaces. `source_configuration_2d` builds the stage and step groups and
the equation partition from physical scheduling choices. The simulation owns
and finalizes the resulting `source_registry_type`. No alternative input
registry initializer or SOURCE_GROUPS namelist is supported.

The tests exercise both analytic and Newton target groups and the atomic
post-RK mass-exchange batch. The strict and optimized GNU profiles retain
logs and compiled objects under D17_WORK_DIR (or a temporary directory).
New physical source families should be added to the canonical catalog and
covered by relevant numerical tests; do not introduce numeric physics IDs
or alternative source configuration paths.

## Paired interlayer sources (block 4)

The same runner builds and executes `test_interlayer_mass_exchange`,
`test_interlayer_lie_transaction`, and
`test_interlayer_momentum_exchange` with strict and optimized profiles.

The cell-local mass-exchange test covers the four layer-presence cases,
two simultaneous gross solid streams, both carrier modes, distinct
settling/elutriation limit factors, independent upward-carrier limiting
in `ALLOW_CONCENTRATED`, and different donor/recipient carrier densities.
It checks both final packing constraints, per-class donor caps,
antisymmetric mass, momentum and thermal energy, and additional gas
species transport.

The transaction test exercises the actual Lie-group entry point with
group-zero snapshots and later-group candidates, verifies atomic rejection,
then applies analytic Chezy drag to the accepted mass-exchange candidates.
It also checks that optional Darcy and OU conservative components remain
local to their respective layers. These tests do not enable the
production `N_LAYERS=2` gate. OpenMP integrated coverage remains in
`TEST_SPATIAL_OPERATOR`.
