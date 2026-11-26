# stochastic: Stochastic Ornstein-Uhlenbeck forcing and reproducible random-number generation used by configured physical laws

## Purpose and role in IMEX-SfloW2D

Stochastic Ornstein-Uhlenbeck forcing and reproducible random-number generation used by configured physical laws.

## Files in this directory

| File | Responsibility |
| --- | --- |
| [`stochastic_module.f90`](./stochastic_module.f90) | OU stochastic state and effective rheological multiplier/perturbation. |
| [`stochastic_random_2d.f90`](./stochastic_random_2d.f90) | Random-number generation utilities and streams. |

## Mathematical and architectural contracts

Distinguish the transported stochastic variable from the effective value used by rheology. Generate random increments outside parallel per-cell source kernels; preserve OpenMP reproducibility and accepted-step semantics.

## Interfaces, invariants and checks

- Respect [the project Fortran style](../../docs/FORTRAN_STYLE.md) and the official Codee configuration.
- Preserve physics, conservative/primitive ownership and REAL/COMPLEX behavior when modifying any implementation.
- Validate affected optimized/debug builds, source or spatial regression tests and OpenMP behavior as appropriate.
