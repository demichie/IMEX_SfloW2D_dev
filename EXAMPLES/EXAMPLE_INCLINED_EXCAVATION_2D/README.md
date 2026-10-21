---
title: "Inclined Excavation Wet–Dry HP-PCCU Test"
summary: "Two-dimensional well-balanced wet–dry test with Q1 bed geometry and excavated plateau."
category: "Numerical verification"
tags:
  - "HP-PCCU"
  - "wet-dry"
  - "well-balanced"
  - "topography"
---

# Inclined Excavation Wet–Dry HP-PCCU Test

## Overview

This numerical-verification example tests spurious transport across wet–dry interfaces on an inclined bed and the hydrostatic well-balanced reconstruction in HP-PCCU. It combines the particle-free setup of the BUMP subcritical case with the two-dimensional restart and boundary layout of EXAMPLE_2D. The original planar bed is `B0(x)=10−tan(10°)x`; a rectangular excavation of depth 2 m produces a wet interior and continuous one-cell Q1 ramps.

## Main parameters

| Parameter | Value |
| --- | --- |
| x range | 0–30 m |
| y range | −10 to +10 m |
| Default grid | 120 × 80 |
| dx = dy | 0.25 m |
| Excavation | 10 ≤ x < 20 m, −5 ≤ y < 5 m |
| Initial thickness | 2 m inside; `h=B0−B` on ramps |
| CFL | 0.24 |
| Limiter / θ | Generalized minmod (3) / 1.3 |
| `N_RK` | 2 |
| Slope/curvature corrections | Disabled |

## Dependencies

The input generator requires Python 3 and NumPy. Matplotlib is required only when the `--plot` option is enabled.

## How to run

From this directory run `python3 create_example.py 120 --plot`; then run the compiled solver with the generated input. For the original verification case, preserve the documented solver flags and initial-state conditions.

## Expected output

`create_example.py` produces `topography_dem.asc`, `inclined_excavation_2D_0000.q_2d`, `IMEX_SfloW2D.inp`, `initial_condition_reference.npz`, and (with `--plot`) `initial_thickness.png` and `initial_centerline.png`.

## Geometry and Q1 consistency

The DEM grid is aligned with computational **vertices**, and its interior values match the authoritative `B_vertex` values. Cell-center elevations and shared face elevations are all derived from this single continuous Q1 topographic field. No independent discontinuous bed values are assigned to the two sides of a shared internal face.

The excavation is rectangular with one-cell linear ramps and bilinear corner patches. On the interior plateau, `h=2 m`; on the ramps, `h=B0−B`; outside it, `h=0`. Consequently `η=B+h=B0` over the wet region. The slope descends toward positive x; initially `u=v=0`, and gravity must accelerate fluid eastward, not laterally.

The problem uses `N_SOLID=0`, `N_ADD_GAS=0`, `LIQUID_FLAG=T`, `GAS_FLAG=F`, `RHEOLOGY_FLAG=F`, thermal energy as the fourth equation, `ENTRAINMENT_FLAG=F`, `LOSS_RATE=0`, `RECONSTR_COEFF=1`, `SLOPE_CORRECTION_FLAG=F` and `CURVATURE_TERM_FLAG=F`. A future variant with slope correction enabled is noted by the source README, but is not part of this default setup.

## Expected transport behavior and diagnostics

- **West/upslope wall** (`x=10 m`): leakage for `x<10 m` is spurious before any physically driven overtopping.
- **North and south walls** (`y=±5 m`): premature transport to the lateral dry regions is spurious; the initial free-surface gradient has no y component.
- **East/downslope wall** (`x=20 m`): acceleration within the excavation precedes accumulation and overtopping of the eastern crest.

The Q1 shared-face bed and positivity-limited `η−B` reconstruction produce zero reconstructed wet/dry endpoint thickness initially at every outer crest. The scheme is intended to produce zero initial PCCU mass flux without a special discontinuous-step classifier.

The original README proposes the following integral measures: `M_uphill` for `x<10 m`; `M_downhill` for `x≥20 m`; `M_south` for `y<−5 m` and `10≤x<20`; `M_north` for `y≥5 m` and `10≤x<20`. Before physically driven overtopping, uphill and lateral masses should remain zero up to roundoff, while downhill mass grows only after eastern overtopping.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 21 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.
