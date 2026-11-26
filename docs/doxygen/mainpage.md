# IMEX-SfloW2D Fortran API

IMEX-SfloW2D is a depth-averaged finite-volume solver for geophysical
and volcanic flows. This reference is generated from the Fortran source
of the same revision as the public documentation website.

## Entry points and data ownership

| Location | Responsibility |
| --- | --- |
| `src/IMEX_SfloW2D.f90` | Main executable and top-level run loop |
| `src/core/simulation_2d.f90` | Simulation lifecycle and canonical per-layer workspaces |
| `src/core/model_layout_2d.f90` | Resolved layer count, state dimensions and indices |
| `src/state/state_fields.f90` | Authoritative conservative/primitive storage |
| `src/state/state_conversions.f90` | State conversion and admissibility |
| `src/spatial/spatial_operator_2d.f90` | Spatial residual, reconstruction and coupled pressure |
| `src/numerics/time_integration_2d.f90` | Accepted timestep and source/spatial coordination |
| `src/numerics/source_catalog_2d.f90` | Canonical source identities, roles and methods |
| `src/numerics/source_registry_2d.f90` | Setup-time source-group validation |
| `src/numerics/lie_source_executor_2d.f90` | Per-cell source transactions and thread-owned scratch |
| `src/physics/` | Cell-local constitutive, thermal, pressure and exchange laws |

The solver supports both `N_LAYERS=1` and `N_LAYERS=2` without a special
activation flag. Both use one canonical `q/qp(component,layer,x,y)` state.
Per-layer physical parameters are independently configured; the terrain,
mesh and substrate inventory are shared.

### Uniform physical-treatment interfaces

Physical laws retain the standard explicit, semi-implicit, implicit
REAL/COMPLEX and analytic entry points, even where a contribution is
identically zero or a method reports that it is unsupported.
These intentional interfaces should not be removed as unused code.

## Scientific limits

The coupled two-layer gravitational characteristic spectrum is not
unconditionally real. Density inversion and some interlayer shear states
can produce complex eigenvalues; a repeated real spectrum may lack a
complete eigenbasis. The implemented total-depth CFL envelope does not
certify hyperbolicity. A successful run is not a proof of physical validity.

- [Architecture](https://github.com/demichie/IMEX_SfloW2D_dev/blob/master/docs/TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md)
- [Two-layer hyperbolicity](https://github.com/demichie/IMEX_SfloW2D_dev/blob/master/docs/TWO_LAYER_HYPERBOLICITY.md)
- [Physical closure assumptions](https://github.com/demichie/IMEX_SfloW2D_dev/blob/master/docs/TL13B_C_PHYSICAL_CLOSURES.md)
- [Fortran conventions](https://github.com/demichie/IMEX_SfloW2D_dev/blob/master/docs/FORTRAN_STYLE.md)
- [MkDocs user guide](https://demichie.github.io/IMEX_SfloW2D_dev/)

## Browse the generated reference

Use the **Modules**, **Files** and **Search** sections of the Doxygen
navigation to inspect individual procedures, arguments and call relations.
The implementation and tests for the selected source revision remain
authoritative; this reference is not a substitute for numerical validation.
