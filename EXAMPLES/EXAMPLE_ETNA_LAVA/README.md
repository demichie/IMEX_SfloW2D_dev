---
title: "Etna Lava Flow with Two Bottom Fissures"
summary: "One-hour lava-flow simulation on original Etna DEM with temperature-dependent rheology."
category: "Volcanic flows"
tags:
  - "Etna"
  - "lava"
  - "fissures"
  - "thermal rheology"
---

# Etna Lava Flow with Two Bottom Fissures

## Overview

A one-hour lava-flow simulation using two bottom fissural sources and temperature-dependent rheology model 3. This example replaces the previous `EXAMPLE_LAVA_IMO` selection for the N9 lava validation case. Unlike the Etna avalanche example, it uses the original topography rather than excavating an initial volume.

## Main parameters

| Parameter | Value |
| --- | --- |
| Grid | 250 × 150 cells |
| Cell size | 10 m |
| Duration / output cadence | 3600 s / 60 s |
| Liquid density | 2600 kg/m³ |
| Fissure temperature | 1538 K |
| Number of fissures | 2 |
| Rheology | Model 3, temperature-dependent |
| `N_RK` / CFL policy | 2 / `directional_min` (default) |
| `KT_TRACE_RECONSTRUCTION` | `h_x` |
| Slope correction / curvature | Disabled |

## How to run

From the repository root, run:

```sh
cd EXAMPLES/EXAMPLE_ETNA_LAVA
python3 prepare.py
OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE ../../bin/IMEX_SfloW2D
```

`prepare.py` extracts `Etna2014_crop.asc` from `../EXAMPLE_ETNA/DEM/Etna2014_crop.zip`, verifying both archive and raster SHA-256 checksums. Only the Python standard library is needed. For a separate clean run, call `python3 EXAMPLES/EXAMPLE_ETNA_LAVA/prepare.py --directory /absolute/path/to/new-run` from the repository root, then run the executable by absolute path from the prepared directory. **Do not run `create_input_ellipsoid.py` for this lava case**: it changes the initial topography and volume.

## Expected output

The solver writes under `OUTPUT/exampleEtna_lava/`: the time-series NetCDF is in `fields/exampleEtna_lava.nc`, ESRI snapshots are in indexed directories under `snapshots/`, terrain exports are in `terrain/`, runout diagnostics are in `diagnostics/`, and the binary checkpoint is in `restart/restart.bin`. Conservative text output is disabled in the supplied input. Open `OUTPUT/exampleEtna_lava/fields/exampleEtna_lava.nc` in ParaView; generated DEM and output files are ignored by Git but not deleted.

## Validation evidence and limits

The current input explicitly selects `h_x` for lava trace reconstruction. The retained validation runs predate this setting; `fixture.json` records both the original and current input checksums.

The original README points to [the ETNA_LAVA validation report](../../docs/acceptance/ETNA_LAVA_EXAMPLE_2026_10_11.md) and to [`fixture.json`](fixture.json), which retains input/DEM checksums and relevant parameter values. The 3600 s runs were checked under strict/optimized profiles and actual one/four-thread teams, with bitwise-equal decoded fields/checkpoints within each profile. The lava front reaches the domain boundary, so a full mass budget must include boundary fluxes; the retained N9 closed-domain mass check ends at the last output before contact. Scientific N9 approval is documented separately from successful execution of this example.

## Documentation verification note

**N9 scientific signoff.** The older README was authored before N9 signoff. The independent N9 v5 acceptance documentation now records approval of its pinned scope. The included validation notes are retained and must not be interpreted as approval of changed input data.

## Additional original README details

The DEM source raster is **3000 × 2000 pixels at 1 m** and the computational grid lies within its extent. The preparation script verifies SHA-256 checksums, and accepts an existing DEM only if identical; it does not overwrite differing existing inputs/DEMs. The solver imports the original terrain to the computational grid. Files generated during the example run are ignored by Git but not deleted. The documented N9 mass check uses saved output only up to 3180 s, the last snapshot before domain-boundary contact; the complete field and restart tests continue to 3600 s.
