# geometry: Computational geometry, topography and active-domain worksets consumed by finite-volume operators

## Purpose and role in IMEX-SfloW2D

Computational geometry, topography and active-domain worksets consumed by finite-volume operators.

## Files in this directory

| File | Responsibility |
| --- | --- |
| [`domain_2d.f90`](./domain_2d.f90) | Wet/active cells, interfaces and compact work lists for numerical kernels. |
| [`geometry_2d.f90`](./geometry_2d.f90) | Grid coordinates, terrain geometry, slopes, curvature and updates following bed changes. |

## Mathematical and architectural contracts

The geometric coefficients are supplied to local physical kernels. Cell traversal and wet-cell worksets remain outside cell-local constitutive laws.

## Interfaces, invariants and checks

- Respect [the project Fortran style](../../docs/FORTRAN_STYLE.md) and the official Codee configuration.
- Preserve physics, conservative/primitive ownership and REAL/COMPLEX behavior when modifying any implementation.
- Validate affected optimized/debug builds, source or spatial regression tests and OpenMP behavior as appropriate.
