# numerics: Time integration, ordered Lie physical-source groups, source capability validation and nonlinear solves; physics is evaluated by the physics modules

## Purpose and role in IMEX-SfloW2D

Time integration, ordered Lie physical-source groups, source capability validation and nonlinear solves; physics is evaluated by the physics modules.

## Files in this directory

| File | Responsibility |
| --- | --- |
| [`lapack_interfaces_2d.f90`](./lapack_interfaces_2d.f90) | Explicit interfaces to linear algebra routines. |
| [`layer_contacts_2d.f90`](./layer_contacts_2d.f90) | Two-layer contact and exposure decisions. |
| [`lie_source_executor_2d.f90`](./lie_source_executor_2d.f90) | Run ordered Lie groups and cell-local trial transactions; hold worker scratch and exchange-free source masks. |
| [`nonlinear_solver_2d.f90`](./nonlinear_solver_2d.f90) | Newton solvers and complex-step Jacobian machinery. |
| [`source_catalog_2d.f90`](./source_catalog_2d.f90) | Single immutable set of schedulable source identities and available methods. |
| [`source_configuration_2d.f90`](./source_configuration_2d.f90) | Translate configured physics into ordered source groups. |
| [`source_input_2d.f90`](./source_input_2d.f90) | Source-related input/configuration bindings. |
| [`source_registry_2d.f90`](./source_registry_2d.f90) | Validate source groups during setup; the Lie executor owns the operational copy. |
| [`time_integration_2d.f90`](./time_integration_2d.f90) | Time integration and accepted-step control. |

## Mathematical and architectural contracts

Group 0 combines explicit spatial/source activity; positive groups represent subsequent Lie operations. A model's supported EXPLICIT, SEMI_IMPLICIT, NEWTON or ANALYTIC pathways are catalog properties, not assumptions. Neither numerical scheduling nor mass-transfer law implementations should be duplicated here.

## Interfaces, invariants and checks

- Respect [the project Fortran style](../../docs/FORTRAN_STYLE.md) and the official Codee configuration.
- Preserve physics, conservative/primitive ownership and REAL/COMPLEX behavior when modifying any implementation.
- Validate affected optimized/debug builds, source or spatial regression tests and OpenMP behavior as appropriate.
