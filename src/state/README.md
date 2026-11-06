# state — canonical solver state

This directory contains the sole production state implementation.

- `state_fields.f90` owns grid fields (`q`, `qp`), substrate inventories,
  accumulated diagnostics and substrate transaction snapshots.
- `state_types.f90` owns the REAL/COMPLEX cell-local physical state types,
  initialized once and reusable directly per worker.
- `state_composition.f90`, `state_mixture.f90`,
  `state_conversions.f90` and `state_physical.f90`
  implement canonical cell-local transformations and evaluations.
- `state.f90` exports these modules as one public facade.

Production code must not import or link `state_2d` or
`state_conversion_2d`; the historical `src/state/` implementation was
removed. The REAL/COMPLEX shared `.inc` kernels live in `src/state/inc/`.
The name `kinematics.inc` is retained because it computes velocities as well
as curvature acceleration.

This change intentionally does not alter the Lie splitting, settling,
source rates or global time-step rollback rules.

`evaluate_mixture_from_volume` is restricted to initialization/source input
where constituent volume fractions are prescribed. It first obtains component
masses per unit volume from the canonical gas equation of state and solid
densities, then calls `evaluate_mixture` for all final thermodynamic outputs.
Normal cell operations start from mass fractions in `qp` and do not use this
inverse conversion.

Dry-state thermodynamics for REAL and COMPLEX builders and refresh are
shared by `inc/dry_thermodynamics.inc`. The wet-state mixture, heat-capacity
validation and reduced-gravity assembly shared by qp builders and refresh is
implemented once in `inc/wet_thermodynamics.inc` (REAL and COMPLEX).

The former single-field workspace containers and their forwarding initializers
were removed; callers initialize the physical state directly. The conservative
mass-admissibility function is exported by `state_composition`, not `state_fields`.
