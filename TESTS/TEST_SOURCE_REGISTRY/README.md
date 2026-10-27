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
