# substrate

## Purpose and role in IMEX-SfloW2D

Shared terrain/substrate configuration for both layers; not a second mobile material layer.

## Files in this directory

| File | Responsibility |
| --- | --- |
| [`substrate_parameters_2d.f90`](./substrate_parameters_2d.f90) | Canonical shared terrain/material configuration |

## Physical equations

Erosion and deposition update substrate inventories. Schematically, with eroded solid volume rate $E_i$ and deposited rate $D_i$, a bed-species storage $B_{s,i}$ satisfies

$$
\dot B_{s,i}=D_i-E_i.
$$

This is a balance illustration: substrate parameters do **not** integrate it. Actual porosity conversion, inventory and topographic updates are defined in state, geometry and mass-exchange source code.

## Architectural and numerical contracts

The substrate is canonical and shared. Its NAMELIST-accessible erodible fractions are kept separate from derived types with ALLOCATABLE components where compiler input handling requires it.

Numerical thresholds, dry-state branches, regularization, and available methods are governed by the linked implementation. The equations above describe implemented physical terms; they are not a replacement for checking each kernel. See the [parent physics guide](../README.md) and [Fortran style](../../../docs/FORTRAN_STYLE.md).
