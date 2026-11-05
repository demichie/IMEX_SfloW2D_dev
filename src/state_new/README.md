# state_new — canonical solver state

This directory contains the sole production state implementation.

- `state_new_fields.f90` owns grid fields (`q`, `qp`), substrate inventories,
  accumulated diagnostics and substrate transaction snapshots.
- `state_new_types.f90` owns the REAL/COMPLEX cell-local physical state types
  and reusable per-worker workspaces.
- `state_new_composition.f90`, `state_new_mixture.f90`,
  `state_new_conversions.f90` and `state_new_physical.f90`
  implement canonical cell-local transformations and evaluations.
- `state_new.f90` exports these modules as one public facade.

Production code must not import or link `state_2d` or
`state_conversion_2d`; the historical `src/state/` implementation was
removed. The REAL/COMPLEX shared `.inc` kernels remain in this directory.

This change intentionally does not alter the Lie splitting, settling,
source rates or global time-step rollback rules.
