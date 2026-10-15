---
title: "Transcritical Bump Flow without Shock"
summary: "Transcritical transition from subcritical to supercritical flow without a shock."
category: "Analytical benchmarks"
tags:
  - "SWASHES"
  - "bump"
  - "shallow water"
---

# Transcritical Bump Flow without Shock

## Overview

Transcritical transition from subcritical to supercritical flow without a shock. The flow uses a single phase with nominal density 1000 kg/m³, no gas, no entrainment, erosion, sedimentation or friction (`RHEOLOGY_FLAG=F`). The original documentation cites SWASHES test 3.1.4.

## Main parameters

| Parameter | Value |
| --- | --- |
| Regime | Ri > 1 to Ri < 1 |
| Density | 1000 kg/m³ (as stated in the original README) |
| Friction | Disabled |
| Boundaries | Left: prescribed specific discharge. Right: zero-gradient outflow. |

## How to run

Generate inputs with `./create_example.py 200 false` (x-cell count and plotting flag) and run `../../../bin/IMEX_SfloW2D` from this case directory.

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` is viewable in ParaView.

## References

1. Delestre, O., Lucas, C., Ksinant, P. A., Darboux, F., Laguerre, C., Vo, T. N. T., et al. (2013). SWASHES: a compilation of shallow water analytic solutions for hydraulic and environmental studies. *International Journal for Numerical Methods in Fluids*, 72(3), 269–300.

## Legacy documentation and preservation

The original [`README.txt`](README.txt) in this example directory remains the archival source for the complete historical narrative, parameter discussion, scientific references and commands. It has not been deleted or rewritten. This standardized README must not replace it as the sole provenance source until a full section-by-section textual comparison has been accepted.
