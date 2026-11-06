---
title: "Supercritical Flow over a Bump"
summary: "Supercritical shallow-water flow with weak bump response."
category: "Analytical benchmarks"
tags:
  - "SWASHES"
  - "bump"
  - "shallow water"
---

# Supercritical Flow over a Bump

## Overview

Supercritical shallow-water flow with weak bump response. The flow uses a single phase with nominal density 1000 kg/m³, no gas, no entrainment, erosion, sedimentation or friction (`RHEOLOGY_FLAG=F`). The Richardson number is very small and the bump causes only a modest change in thickness.

## Main parameters

| Parameter | Value |
| --- | --- |
| Regime | Ri < 1 (typically Ri < 10^-2) |
| Density | 1000 kg/m³ (as stated in the original README) |
| Friction | Disabled |
| Boundaries | Left: hu = 10 m²/s and h = 1 m. Right: zero-gradient outflow. |

## Dependencies

Requires Python 3 and NumPy. Matplotlib is required only when the plotting flag is `true`.

## How to run

From this example directory, first prepare the inputs and then run the solver.

### 1. Generate the input files

The generator accepts **2 positional arguments**, in this order:

1. Number of computational cells along x (positive integer).
2. Whether to plot the initial solution (`true` or `false`).

The Python script sets up the specified SWASHES bump-flow test using a single constant-density phase, with neither friction nor mass-exchange processes.

**Generator command:**

```bash
./create_example.py 200 false
```

### 2. Run the solver

**Solver command:**

```bash
../../../bin/IMEX_SfloW2D
```

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` is viewable in ParaView.

## References

1. Delestre, O., Lucas, C., Ksinant, P. A., Darboux, F., Laguerre, C., Vo, T. N. T., et al. (2013). SWASHES: a compilation of shallow water analytic solutions for hydraulic and environmental studies. *International Journal for Numerical Methods in Fluids*, 72(3), 269–300.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 21 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.
