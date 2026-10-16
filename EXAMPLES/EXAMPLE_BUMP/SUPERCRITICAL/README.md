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

Generate inputs with `./create_example.py 200 false` (x-cell count and plotting flag) and run `../../../bin/IMEX_SfloW2D` from this case directory.

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` is viewable in ParaView.

## References

1. Delestre, O., Lucas, C., Ksinant, P. A., Darboux, F., Laguerre, C., Vo, T. N. T., et al. (2013). SWASHES: a compilation of shallow water analytic solutions for hydraulic and environmental studies. *International Journal for Numerical Methods in Fluids*, 72(3), 269–300.

## Legacy documentation and preservation

The original [`README.txt`](README.txt) in this example directory remains the archival source for the complete historical narrative, parameter discussion, scientific references and commands. It has not been deleted or rewritten. This standardized README must not replace it as the sole provenance source until a full section-by-section textual comparison has been accepted.
