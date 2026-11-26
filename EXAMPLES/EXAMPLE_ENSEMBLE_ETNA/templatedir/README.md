# Etna ensemble realization template

## Purpose and role in IMEX-SfloW2D

This directory is a **template**, not an independent simulation example.
`EXAMPLE_ENSEMBLE_ETNA/create_inputfiles.py` copies it into separate ensemble
member directories and replaces sampled `ENSEMBLE_*` input placeholders.
The parent [ensemble README](../README.md) is the authoritative guide to
preparing, executing and aggregating the members.

## Files in this directory

| File | Role |
| --- | --- |
| `IMEX_SfloW2D.inp` | Parameter template, including deliberate ensemble placeholders. |
| `input_ellipsoid.py` | Per-member initial ellipsoid/excavation configuration. |
| `create_input_ellipsoid.py` | Per-member initial-volume and DEM preprocessing. |
| `run.sh` | Per-member preprocessing followed by the solver. |
| `plot_overlay.py` | Optional plotting utility. |

## Shared terrain and execution

Before creating members, prepare the Etna DEM once according to the parent
guide and establish the intended symbolic link in `templatedir/DEM/`.
Do not commit another copy of the large raster here. The executable symlink
is prepared in the template prior to member generation; do not run the
placeholder-bearing input as if it were a completed realization.
