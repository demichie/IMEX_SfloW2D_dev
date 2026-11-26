# state/inc: Shared cell-local mathematical kernels included in REAL and COMPLEX implementations, not separately compiled Fortran modules

## Purpose and role in IMEX-SfloW2D

Shared cell-local mathematical kernels included in REAL and COMPLEX implementations, not separately compiled Fortran modules.

## Files in this directory

| File | Responsibility |
| --- | --- |
| [`conservative_from_primitive.inc`](./conservative_from_primitive.inc) | Primitive-to-conservative cell conversion. |
| [`density_from_mass.inc`](./density_from_mass.inc) | Composition-derived mixture density. |
| [`dry_thermodynamics.inc`](./dry_thermodynamics.inc) | Thermodynamic conventions for dry states. |
| [`gas_specific_volume.inc`](./gas_specific_volume.inc) | Gas-phase specific volume. |
| [`kinematics.inc`](./kinematics.inc) | Cell-local velocity and speed invariants. |
| [`mixture_from_mass.inc`](./mixture_from_mass.inc) | Mixture composition from component masses. |
| [`mixture_heat_capacity.inc`](./mixture_heat_capacity.inc) | Mixture heat capacity. |
| [`reduced_gravity.inc`](./reduced_gravity.inc) | Reduced-gravity evaluation. |
| [`richardson.inc`](./richardson.inc) | Local Richardson-number calculation. |
| [`temperature_from_energy.inc`](./temperature_from_energy.inc) | Temperature from conservative thermal content. |
| [`velocity_from_momentum.inc`](./velocity_from_momentum.inc) | Velocity from conservative momentum. |

## Mathematical and architectural contracts

These includes share substantive formulas across typed REAL and COMPLEX caller contexts. They must not allocate a second canonical state or read cell indices. With `M=rho_m h`, the common conversions include `u=m_x/M`, `v=m_y/M` and `T=Q_T/(M c_{p,m})` for admissible wet states; dry-cell branches follow `state_fields`.

## Interfaces, invariants and checks

- Respect [the project Fortran style](../../../docs/FORTRAN_STYLE.md) and the official Codee configuration.
- Preserve physics, conservative/primitive ownership and REAL/COMPLEX behavior when modifying any implementation.
- Validate affected optimized/debug builds, source or spatial regression tests and OpenMP behavior as appropriate.
