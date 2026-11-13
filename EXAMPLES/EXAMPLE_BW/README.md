---
title: "One-Dimensional Thermally Evolving Flow"
summary: "1-D ash-flow-inspired test of thermal transport, friction, deposition and entrainment."
category: "Volcanic flows"
tags:
  - "thermal energy"
  - "ash flow"
  - "1D"
---

# One-Dimensional Thermally Evolving Flow

## Overview

A 1-D flow down a gentle slope with friction, sedimentation and air entrainment. The thermal-energy transport equation is retained. Deposition does not change the topography when `SUBSTRATE%TOPO_CHANGE_FLAG=F`; the original README recommends changing this flag to study its effect. The parameters are similar to Bursik and Woods (1996).

## Main parameters

| Parameter | Value |
| --- | --- |
| Geometry | 1-D gentle slope |
| Thermal-energy transport | Enabled |
| Bed adjustment by deposition | Disabled by default |
| x-resolution, solids fraction, temperature | Generator arguments |

## Dependencies

Requires Python 3 and NumPy. Matplotlib is required only when the plotting flag is `true`.

## How to run

From this example directory, first generate the inputs, then run the solver.

### 1. Generate the input files

The generator takes **four positional arguments**, in this order:

1. Number of computational cells along x.
2. Initial solid-particle volume fraction (dimensionless, 0–1).
3. Initial flow temperature (K).
4. Plot the initial solution (`true` or `false`).

The example includes friction, settling/deposition and entrainment, and transports thermal content.

**Generator command:**

```bash
./create_example.py 400 0.01 900 false
```

### 2. Run the solver

**Solver command:**

```bash
../../bin/IMEX_SfloW2D
```

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc`, viewable in ParaView.

## References

1. Bursik, M. I., and Woods, A. W. (1996). The dynamics and thermodynamics of large ash flows. *Bulletin of Volcanology*, 58, 175–193.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 21 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.

## Topography files

Run `python create_example.py` from this example directory. The generator creates
`DEM/topography_dem.asc` and the `DEM/` directory when needed. The input
references that location; do not place generated terrain in the example root.
