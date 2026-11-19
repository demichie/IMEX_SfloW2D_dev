# Repository documentation

## Purpose and role in IMEX-SfloW2D

This index distinguishes **current implementation instructions** from
**historical acceptance and migration records**. Dated evidence is kept for
reproducibility; it is not a second implementation roadmap or a source of
obsolete compatibility requirements.

## Current engineering guides

- [Development roadmap](IMEX_SfloW2D_Development_Roadmap.md): consolidated priorities and development history.
- [Two-layer architecture and roadmap](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md): authoritative current two-layer design and outstanding gate.
- [Two-layer code audit](TWO_LAYER_CODE_AUDIT.md): architecture and validation findings at its stated baseline.
- [State architecture](STATE_ARCHITECTURE_REDESIGN_ROADMAP.md): state ownership and scientific invariants; dated baseline, check current source before acting.
- [Lie timestep transaction](LIE_TIMESTEP_TRANSACTION_ROADMAP.md): current step-acceptance and failure contract.
- [Fortran style and Codee](FORTRAN_STYLE.md): mandatory code and Doxygen conventions.
- [Adding a rheology model](ADDING_A_RHEOLOGY_MODEL.md): extension procedure and tests.
- [NetCDF output](NETCDF_OUTPUT.md): physical-field output schema and binary restart relationship.
- [Lie geometry plan](LIE_CELL_OUTER_GEOMETRY_FREEZE_PLAN.md): rationale and implementation record; status statements have their own date/baseline.

## Source and example documentation

- [Source directory map](../src/physics/README.md) and the individual
  `src/*/README.md` guides document every source directory, including the
  constitutive equations of physical-law families.
- [Examples](../EXAMPLES/) use `README.md` as the single run guide in
  each example. The
  [README migration audit](../EXAMPLES/README_SECOND_PASS_AUDIT.md)
  retains the comparison evidence for the removed plain-text READMEs.
- [Examples input audit](EXAMPLES_INPUT_COMPATIBILITY_AUDIT.md) records
  a dated input-validation step; re-run its checker before relying on its
  historical conclusions.

## Historical evidence (retained intentionally)

- [Acceptance records](acceptance/) contain dated test results, manifests,
  reference inputs and supporting logs. These are **not** current
  compatibility formats or prescribed scientific configurations.
- [Archived roadmaps](archive/roadmaps/) record completed or superseded
  proposals. Do not implement them instead of the current roadmap.
- [Debug warning cleanup](DEBUG_WARNING_CLEANUP_2026_11_06.md) is a
  snapshot of a particular compiler-warning baseline, not a perpetual
  warning count.

For any new change, inspect the current branch, declarations, tests and
GitHub Actions results before trusting a dated document's status.
