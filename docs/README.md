# Repository documentation

This directory contains **current developer contracts** and a separately
identified **historical acceptance record**. The repository's source code
and CI for the exact active commit remain the source of truth.

## Current documentation

- [Development roadmap](IMEX_SfloW2D_Development_Roadmap.md): current TL13 priorities.
- [Two-layer architecture](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md):
  authoritative design, physical assumptions and remaining gates.
- [Two-layer closure audit](TWO_LAYER_CODE_AUDIT.md): technical snapshot;
  historical status claims must be checked against current code.
- [State architecture](STATE_ARCHITECTURE_REDESIGN_ROADMAP.md): ownership,
  conversion and regression contract.
- [Lie timestep transaction](LIE_TIMESTEP_TRANSACTION_ROADMAP.md):
  accepted-step publication and fatal rejection policy.
- [Fortran style](FORTRAN_STYLE.md): mandatory formatting and Doxygen rules.
- [Adding a rheology model](ADDING_A_RHEOLOGY_MODEL.md): extension guide
  with illustrative, **not registered**, physics.
- [NetCDF output](NETCDF_OUTPUT.md): named physical visualization fields;
  binary checkpoints hold the full restart state.
- [Visualization with ParaView and QGIS](website/getting-started/visualization.md):
  menu-by-menu instructions for scientific plots and geographic textures.

## Other authoritative sources

- [Physics/source documentation](../src/physics/README.md) and
  [state documentation](../src/state/README.md) describe owners and APIs.
- Each example's `EXAMPLES/*/README.md` is the canonical run guide.
- [Website sources](website/) are built by the MkDocs publishing workflow.
- [Dated acceptance reports](acceptance/README.md) and their pinned JSON
  fixtures are historical evidence; some are still **required by the
  executable acceptance tests** in `TESTS/ACCEPTANCE/`.

## Retired documentation and reproducibility

Superseded roadmap archives, one-off cleanup plans and the obsolete
large state-inventory manifest are no longer maintained in the active
tree. They remain accessible through the immutable
[pre-cleanup Git snapshot](https://github.com/demichie/IMEX_SfloW2D_dev/tree/d91843a0faf8c8491cf50c240cdd50a7cd006407/docs).
Test-consumed acceptance JSON fixtures remain local; links in dated
reports also point to their original pinned revision for reproducibility.
Old binary/restart/NetCDF versions in historical reports do **not**
imply ongoing compatibility requirements.

Do not use a dated migration plan to reinstate obsolete interfaces.
Check the current branch, source, tests and CI results for the exact SHA.
