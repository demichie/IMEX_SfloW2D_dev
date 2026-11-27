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

The table lists the **active physical source processes** in ascending **Group** order, complementing the Physical processes description above. **Scope** identifies the affected layer or layer pair; **Source** gives the exact input parameter name (and the selected rheology model where relevant), without repeating physical coefficients; **Method** is the integration treatment; **Role** briefly explains the physical effect. Different processes may share a group.

| Group | Scope | Source | Method | Role |
| ---: | --- | --- | --- | --- |
| 1 | Flow | `RHEOLOGY%RHEOLOGY_MODEL='quadratic'` | `NEWTON` | Quadratic basal resistance modifies momentum. |
| 4 | Flow | `MASS_EXCHANGE_L1%DEPOSITION%SETTLING_FLAG` | `EXPLICIT` | Settles suspended solid particles onto the substrate. |
| 4 | Flow | `MASS_EXCHANGE_L1%ENTRAINMENT%FLAG` | `EXPLICIT` | Entrains ambient gas into the flowing mixture. |


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
