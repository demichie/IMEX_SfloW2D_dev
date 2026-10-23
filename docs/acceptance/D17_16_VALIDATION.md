# D17.16 — Input and source-schedule validation

## Implemented

The legacy input path continues to validate `IMEX_*` group IDs (1–3) and the
supported methods (`ANALYTIC` or `NEWTON`). Sentinel values distinguish
explicitly supplied legacy fields from defaults, so native and explicitly
configured legacy schedules cannot be combined accidentally.

The native input is a separate optional `SOURCE_GROUPS` namelist. For example,
the velocity block remains present even when no physical source rates are
active:

```fortran
&SOURCE_GROUPS
  source_group_count=1,
  source_group_id(1)=1,
  source_group_order(1)=1,
  source_group_cadence(1)='STAGE',
  source_group_method(1)='NEWTON',
  source_group_term_count(1)=0,
  source_group_solver_count(1)=1,
  source_solver_order(1,1)=1,
  source_solver_target_count(1,1)=1,
  source_solver_targets(1,1,1)='VELOCITY'
/
```

Supported cadences are `STAGE` and `STEP`; methods are `EXPLICIT`, `ANALYTIC`
and `NEWTON`; solver targets are `VELOCITY`, `TEMPERATURE` and
`PORE_PRESSURE`. `STEP` groups must be explicit. The named term list must
include every physical term activated by the other input options exactly once;
omissions do not disable physics. Duplicate or unknown terms, missing implicit
targets, duplicate group/order identities, invalid target layouts, malformed
or repeated namelists, and conflicts with explicit legacy `IMEX_*` settings
are rejected.

Input storage has fixed catalog-bounded capacities: at most
`SOURCE_PHYSICAL_TERM_COUNT` groups and terms per group, three solvers per
group and three targets per solver. The native builder constructs the equation
partition and registry only after validating the complete active-term
inventory.

## Local evidence

`TESTS/TEST_SOURCE_REGISTRY/run_test.sh` passed in strict and optimized
profiles. Each profile covered 864 legacy configurations, 14,893 assertions,
28 invalid registry configurations and parser cases for a valid namelist,
absence, malformed input, duplicates and a native/legacy conflict.

`TESTS/TEST_MASS_EXCHANGE/run_test.sh` passed in strict and optimized profiles.
It compiles the production module order (including `source_configuration_2d`,
`simulation_2d` and `inpout_2d`) and exercises mass-exchange and checkpoint
continuation.

## Remaining acceptance gap

The parser cases exercise the production namelist-reader helper using scratch
input units; the complete `read_param` startup flow is not yet run by a
dedicated input fixture. Native schedule equivalence across physical model
combinations and cross-platform/CI acceptance remain open. D17.16 is therefore
**partially complete**.
