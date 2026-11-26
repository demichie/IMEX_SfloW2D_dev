# Numerical integration and Lie source execution

`src/numerics/` owns time stepping, ordered physical-source groups,
nonlinear solves and source capability validation. Constitutive laws
remain in `src/physics/`, not duplicate implementations here.

## Modules

| File | Responsibility |
| --- | --- |
| [`time_integration_2d.f90`](time_integration_2d.f90) | Trial conservative/primitive fields and accepted-step publication |
| [`lie_source_executor_2d.f90`](lie_source_executor_2d.f90) | Ordered cell-local Lie groups and OpenMP worker scratch |
| [`source_catalog_2d.f90`](source_catalog_2d.f90) | Source identities and supported numerical treatments |
| [`source_configuration_2d.f90`](source_configuration_2d.f90) | Convert physical settings to source groups |
| [`source_registry_2d.f90`](source_registry_2d.f90) | Validate groups at simulation initialization |
| [`source_input_2d.f90`](source_input_2d.f90) | Source configuration input bindings |
| [`nonlinear_solver_2d.f90`](nonlinear_solver_2d.f90) | Newton solves and complex-step Jacobians |
| [`lapack_interfaces_2d.f90`](lapack_interfaces_2d.f90) | Linear algebra interfaces |
| [`layer_contacts_2d.f90`](layer_contacts_2d.f90) | Lower/upper layer exposure and contact |

## Timestep and failure contracts

- `time_integration_2d` owns `q_trial` and `qp_trial`. **Authoritative
  fluid `q/qp` are published only after the spatial stage and every
  scheduled Lie group succeeds.**
- Group zero evaluates rates on the pre-spatial state `q^n` and applies
  accepted increments to the post-spatial candidate `q^*`.
  Positive-numbered groups execute in the configured order.
- The catalog validates each law's genuine `EXPLICIT`,
  `SEMI_IMPLICIT`, `NEWTON` or `ANALYTIC` capability. Uniform
  per-law treatment interfaces remain, including intentional zero
  contributions or explicit unsupported-treatment errors.
- Invalid spatial updates, source groups or state conversions are
  **fatal**. There is no recoverable global substrate rollback or
  timestep retry; shared substrate inventories may change during
  the cell sweep. Geometry-local deferred publication is separate.
- One- and two-layer configurations share this integration path.
  Interlayer transfers and drag obey the pair-wide source contracts in
  [the architecture specification](../../docs/TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md).

## OpenMP and diagnostics

Each worklist cell owns its candidate and visits its physical groups in
order. The executor uses `SCHEDULE(DYNAMIC,16)` for heterogeneous cell
work without changing a cell's source order. Failed-cell selection uses
the earliest worklist index with synchronized witness capture.

The failure-only `D19_SUBSTRATE` diagnostic reports the first rejected
substrate candidate and its inventories/rates. Optional
`IMEX_PROFILE_LIE=1` measures approximate worker cell-compute balance;
it is not a measurement of whole-solver parallel efficiency.

## Verification

Preserve numerical evaluation order, REAL/COMPLEX equivalence,
conservation, dry states and disjoint OpenMP writes. Use the
`TEST_SPATIAL_OPERATOR`, `TEST_MASS_EXCHANGE`,
`TEST_WORKSPACE_LIFECYCLE` and source registry tests as appropriate;
verify optimized and strict/debug builds plus representative
single-/two-layer restarts and end-to-end OpenMP runs.

See [the state implementation](../state/README.md),
[physical source laws](../physics/README.md) and
[manually dispatched validation workflows](../../.github/ACTIONS.md).
