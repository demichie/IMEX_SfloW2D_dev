---
title: "One-Dimensional Dam Break with Mass Exchange"
summary: "One-dimensional dam break over flat terrain with deposition, entrainment, and friction."
category: "Numerical benchmarks"
tags:
  - "dam break"
  - "1D"
  - "entrainment"
---

# One-Dimensional Dam Break with Mass Exchange

## Overview

A 1-D dam-break problem on flat terrain, including deposition, entrainment, and basal friction.

## Main parameters

| Parameter | Value |
| --- | --- |
| Geometry | 1-D, flat bed |
| Particle volume fraction | Generator argument |
| Temperature | Generator argument (K) |
| x-resolution | Generator argument |

## Dependencies

Requires Python 3 and NumPy. Matplotlib is required only when the plotting flag is `true`.

## Source integration

The table follows ascending **Group** numbers (the order of Lie source updates after spatial transport). **Scope** identifies the layer or layer pair affected; **Source (input setting)** reproduces the relevant input parameter names and values; **Method** is the configured numerical treatment; **Status** distinguishes enabled sources from settings that do not activate a process; **Role** explains the physical contribution. An entry marked *Configured only* is not an executed source group. Group 0 denotes the initial explicit source stage, not a time-integration stage count.

| Group | Scope | Source (input setting) | Method | Status | Role |
| ---: | --- | --- | --- | --- | --- |
| 0 | Flow | `GEOMETRY_CURVATURE_METHOD='EXPLICIT'` | `EXPLICIT` | Configured only | Curvature correction is not enabled. |
| 1 | Flow | `RHEOLOGY%RHEOLOGY_MODEL='quadratic'` | `NEWTON` | Active | Quadratic basal resistance modifies momentum. |
| 2 | Flow | `THERMAL%THERMAL_METHOD='NEWTON'` | `NEWTON` | Configured only | Cooling is not enabled. |
| 4 | Flow | `MASS_EXCHANGE_L1%DEPOSITION%SETTLING_FLAG=T`, `MASS_EXCHANGE_L1%ENTRAINMENT%FLAG=T` | `EXPLICIT` | Active | Settling and atmospheric entrainment exchange material with the flow. |

## How to run

From this example directory, first generate the inputs, then run the solver.

### 1. Generate the input files

The generator takes **four positional arguments**, in this order:

1. Number of computational cells along x.
2. Initial solid-particle volume fraction (dimensionless, 0–1).
3. Initial flow temperature (K).
4. Plot the initial solution (`true` or `false`).

The script prepares the input for the 1-D flat-bed dam-break with deposition, entrainment and friction.

**Generator command:**

```bash
./create_example.py 400 0.5 300 false
```

### 2. Run the solver

**Solver command:**

```bash
../../src/IMEX_SfloW2D
```

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` may be inspected with ParaView.

## Documentation provenance

This Markdown guide supersedes the earlier plain-text instructions. Differences in historical commands or parameters are documented here and in [the examples documentation audit](../README_SECOND_PASS_AUDIT.md).

## Topography files

Run `python create_example.py` from this example directory. The generator creates
`DEM/topography_dem.asc` (and the directory if necessary). The generated
`IMEX_SfloW2D.inp` points to that path; the DEM is not written to the
example root. Keep numerical outputs outside `DEM/`.
