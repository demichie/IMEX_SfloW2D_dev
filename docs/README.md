# Repository documentation

This directory contains the **current scientific and developer documentation**
for IMEX-SfloW2D. The source, tests and CI results for the exact revision
remain authoritative for implemented behavior.

## Current specifications and guides

| Document | Purpose |
| --- | --- |
| [Development roadmap](IMEX_SfloW2D_Development_Roadmap.md) | Integrated functionality and outstanding scientific/engineering work |
| [Two-layer architecture](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md) | Authoritative ownership, spatial coupling and layer source contracts |
| [Two-layer hyperbolicity](TWO_LAYER_HYPERBOLICITY.md) | Density/shear characteristic limits and conditional CFL arguments |
| [Interlayer physical closures](TL13B_C_PHYSICAL_CLOSURES.md) | Source budgets, accepted assumptions and calibration gaps |
| [Initial-state input and restart](RESTART_INPUT.md) | Layer initialization and binary checkpoint continuation |
| [NetCDF output](NETCDF_OUTPUT.md) | Physical visualization fields, units and layer names |
| [Adding a rheology model](ADDING_A_RHEOLOGY_MODEL.md) | Source extension guide (illustrative laws are not registered) |
| [Fortran style](FORTRAN_STYLE.md) | Mandatory source and Doxygen conventions |

Implementation-level contracts are maintained **beside the code**:

- [State and conversions](../src/state/README.md)
- [Finite-volume spatial methods](../src/spatial/README.md)
- [Geometry and active-domain OpenMP worksets](../src/geometry/README.md)
- [Time integration and Lie source transactions](../src/numerics/README.md)
- [Physical laws and source treatments](../src/physics/README.md)
- [I/O and restart implementation](../src/io/README.md)

## User documentation

The MkDocs website is a **user manual**: document available physics,
simulation setup, numerical limitations, outputs and visualization there.
Keep implementation roadmaps, dated acceptance evidence, CI/CD
instructions and developer audits in the repository, not the public
website navigation. Public explanations of two-layer equations and
physical assumptions must link to other website pages rather than
requiring readers to consult internal GitHub development documents.


- [Website sources](website/) and [published MkDocs](https://demichie.github.io/IMEX_SfloW2D_dev/)
- [Doxygen entry page](doxygen/mainpage.md) and [published API](https://demichie.github.io/IMEX_SfloW2D_dev/doxygen/)
- [Visualization](website/getting-started/visualization.md) and [Docker](website/getting-started/docker.md); publication workflows are described in [Actions](../.github/ACTIONS.md)
- [Runnable examples](../EXAMPLES/), each documented by its own README

## Historical evidence

The [acceptance directory](acceptance/README.md) retains dated scientific
reports and machine-readable evidence. **Some JSON fixtures are consumed by
active acceptance tests**. Preserve them until their dependencies are
explicitly migrated; dated results apply only to their recorded revisions.

Superseded state-redesign, Lie-transaction and active-domain development
roadmaps, and the historical TL13 implementation diary, are no longer
maintained as current specifications. Their still-relevant invariants are
included in the module guides above. Original texts remain accessible in
[the fixed pre-cleanup revision](https://github.com/demichie/IMEX_SfloW2D_dev/tree/cccdee3ef8d54b7efea92ace97d4beabc0c52ac5/docs).
An older [migration-document archive](https://github.com/demichie/IMEX_SfloW2D_dev/tree/d91843a0faf8c8491cf50c240cdd50a7cd006407/docs)
is also available.

Historical records do not authorize restoring retired interfaces, earlier
restart formats or alternate production implementations.
