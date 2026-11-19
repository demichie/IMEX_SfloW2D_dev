# core: Solver-wide configuration, model layout, execution lifecycle, diagnostics and initialization

## Purpose and role in IMEX-SfloW2D

Solver-wide configuration, model layout, execution lifecycle, diagnostics and initialization. It owns the simulation context and global control flow, not constitutive formulas.

## Files in this directory

| File | Responsibility |
| --- | --- |
| [`complexify.f90`](./complexify.f90) | REAL/COMPLEX differentiation support. |
| [`constitutive_parameters_2d.f90`](./constitutive_parameters_2d.f90) | Shared material and ambient parameters. |
| [`diagnostics_2d.f90`](./diagnostics_2d.f90) | Error reporting and fatal diagnostics. |
| [`init_2d.f90`](./init_2d.f90) | Initial-condition construction and release. |
| [`layer_physics_2d.f90`](./layer_physics_2d.f90) | Per-layer physical configuration selection. |
| [`model_layout_2d.f90`](./model_layout_2d.f90) | Resolved component counts and indices. |
| [`parameters_2d.f90`](./parameters_2d.f90) | Working precision, global run settings and equation-index vocabulary. |
| [`runtime_2d.f90`](./runtime_2d.f90) | Simulation time and step bookkeeping. |
| [`simulation_2d.f90`](./simulation_2d.f90) | Simulation context, initialization, advance and finalization. |
| [`solid_transport_parameters_2d.f90`](./solid_transport_parameters_2d.f90) | Solid-species material and packing parameters. |

## Mathematical and architectural contracts

The simulation object aggregates state, per-layer geometry, source registry and workspaces. Model layout is resolved before allocation; do not duplicate its component indices or move the numerical kernels here.

## Interfaces, invariants and checks

- Respect [the project Fortran style](../../docs/FORTRAN_STYLE.md) and the official Codee configuration.
- Preserve physics, conservative/primitive ownership and REAL/COMPLEX behavior when modifying any implementation.
- Validate affected optimized/debug builds, source or spatial regression tests and OpenMP behavior as appropriate.
