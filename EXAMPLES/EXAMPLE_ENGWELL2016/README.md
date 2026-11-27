---
title: "Radial Flow Benchmark Inspired by Engwell et al. (2016)"
summary: "Radial supercritical source with particle sedimentation and air entrainment."
category: "Volcanic flows"
tags:
  - "radial source"
  - "coignimbrite"
  - "supercritical"
---

# Radial Flow Benchmark Inspired by Engwell et al. (2016)

## Overview

A supercritical radial inlet over flat terrain, intended to reproduce the behavior of Figure 3 in Engwell et al. (2016). Sedimentation and atmospheric-air entrainment are enabled; a simple friction law is used. The source must be supercritical (Ri < 1) to avoid characteristics entering the prescribed boundary from the computational domain.

## Main parameters

| Parameter | Value |
| --- | --- |
| Source radius | Generator argument |
| Initial thickness | Generator argument |
| Initial temperature | Generator argument (K) |
| Richardson number | Generator argument, < 1 |
| Solid mass fraction (`XS_SOURCE`) | Sixth generator argument (`xs_init`) |

## Dependencies

The input generator requires Python 3 and NumPy.

## Source integration

The table lists the **active physical source processes** in ascending **Group** order, complementing the Physical processes description above. **Scope** identifies the affected layer or layer pair; **Source** gives the exact input parameter name (and the selected rheology model where relevant), without repeating physical coefficients; **Method** is the integration treatment; **Role** briefly explains the physical effect. Different processes may share a group.

| Group | Scope | Source | Method | Role |
| ---: | --- | --- | --- | --- |
| 0 | Flow | `RADIAL_SOURCE_FLAG` | `EXPLICIT` | Radial cell injection supplies the supercritical inflow. |
| 1 | Flow | `RHEOLOGY%RHEOLOGY_MODEL='quadratic'` | `NEWTON` | Quadratic basal resistance modifies momentum. |
| 4 | Flow | `MASS_EXCHANGE%DEPOSITION%SETTLING_FLAG` | `EXPLICIT` | Settles suspended solid particles onto the substrate. |
| 4 | Flow | `MASS_EXCHANGE%ENTRAINMENT%FLAG` | `EXPLICIT` | Entrains ambient gas into the flowing mixture. |


## How to run

From this example directory, first generate the inputs and then run the solver.

### 1. Generate the input files

The generator expects **6 positional arguments**, in this exact order:

1. Number of cells along x (positive integer).
2. Source radius (m; positive).
3. Source thickness (m; positive).
4. Initial source temperature (K; positive).
5. Initial Richardson number (dimensionless; below 1 for the supercritical inlet).
6. Initial solid **mass fraction** (dimensionless, between 0 and 1).

The sixth value is written to `XS_SOURCE`, a **mass fraction**; the older README also calls it a volume fraction, but that statement conflicts with the current generator. The supercritical regime (`Ri<1`) prevents characteristics from propagating upstream across the prescribed source.

**Generator command:**

```bash
python create_example.py 200 2000 2000 900 0.9 0.8
```

### 2. Run the solver

**Solver command:**

```bash
../../src/IMEX_SfloW2D
```

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` can be viewed in ParaView.

## References

1. Bursik, M. I., and Woods, A. W. (1996). The dynamics and thermodynamics of large ash flows. *Bulletin of Volcanology*, 58, 175–193.
2. Engwell, S. L., de' Michieli Vitturi, M., Esposti Ongaro, T., and Neri, A. (2016). Insights into the formation and dynamics of coignimbrite plumes from one-dimensional models. *Journal of Geophysical Research: Solid Earth*, 121(6), 4211–4231.

## Documentation provenance

This Markdown guide supersedes the earlier plain-text instructions. Historical discrepancies are recorded here and in [the examples documentation audit](../README_SECOND_PASS_AUDIT.md).

## Topography files

Run `python create_example.py` from this example directory. The generator creates
`DEM/topography_dem.asc` and the `DEM/` directory when needed. The input
references that location; do not place generated terrain in the example root.
