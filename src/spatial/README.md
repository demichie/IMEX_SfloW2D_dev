# spatial: Finite-volume spatial discretization: reconstruction, numerical fluxes, nonconservative contributions, CFL estimates and explicit transport

## Purpose and role in IMEX-SfloW2D

Finite-volume spatial discretization: reconstruction, numerical fluxes, nonconservative contributions, CFL estimates and explicit transport.

## Files in this directory

| File | Responsibility |
| --- | --- |
| [`hyperbolic_2d.f90`](./hyperbolic_2d.f90) | Hyperbolic flux and characteristic evaluations. |
| [`nonconservative_2d.f90`](./nonconservative_2d.f90) | Nonconservative spatial contributions. |
| [`pccu_2d.f90`](./pccu_2d.f90) | PCCU numerical spatial-flux machinery. |
| [`reconstruction_2d.f90`](./reconstruction_2d.f90) | Cell-to-face reconstruction, boundary traces and admissibility. |
| [`spatial_operator_2d.f90`](./spatial_operator_2d.f90) | Model-facing spatial residual/CFL API and workset orchestration. |

## Mathematical and architectural contracts

In schematic form the semi-discrete balance is `dQ_{j,k}/dt = -(F_{j+1/2,k}-F_{j-1/2,k})/dx -(G_{j,k+1/2}-G_{j,k-1/2})/dy + spatial/source terms`. Actual path contributions, boundary fluxes and reconstruction constraints are implemented here. The spatial operator owns the transport residual, not physical constitutive parameters.

## Interfaces, invariants and checks

- Respect [the project Fortran style](../../docs/FORTRAN_STYLE.md) and the official Codee configuration.
- Preserve physics, conservative/primitive ownership and REAL/COMPLEX behavior when modifying any implementation.
- Validate affected optimized/debug builds, source or spatial regression tests and OpenMP behavior as appropriate.
