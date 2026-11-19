# io: Parse input, manage output paths, save fields and restart records; connect disk formats to canonical runtime/state representations

## Purpose and role in IMEX-SfloW2D

Parse input, manage output paths, save fields and restart records; connect disk formats to canonical runtime/state representations.

## Files in this directory

| File | Responsibility |
| --- | --- |
| [`inpout_2d.f90`](./inpout_2d.f90) | Input NAMELIST parsing, output schedules, NetCDF, binary restart and diagnostic output. |
| [`output_dirs.c`](./output_dirs.c) | Filesystem directory helper used by output paths. |
| [`output_paths.f90`](./output_paths.f90) | Output filename/path construction and directory naming. |

## Mathematical and architectural contracts

Input NAMELIST aliases for single-layer use are part of the intended input model; they are not historical physics objects. Restart and NetCDF schemas follow the current development format rather than preserving obsolete alternatives.

## Interfaces, invariants and checks

- Respect [the project Fortran style](../../docs/FORTRAN_STYLE.md) and the official Codee configuration.
- Preserve physics, conservative/primitive ownership and REAL/COMPLEX behavior when modifying any implementation.
- Validate affected optimized/debug builds, source or spatial regression tests and OpenMP behavior as appropriate.
