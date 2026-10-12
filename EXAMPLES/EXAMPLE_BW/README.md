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

A 1-D flow down a gentle slope with friction, sedimentation and air entrainment. The thermal-energy transport equation is retained. Deposition does not change the topography when `TOPO_CHANGE_FLAG=F`; the original README recommends changing this flag to study its effect. The parameters are similar to Bursik and Woods (1996).

## Main parameters

| Parameter | Value |
| --- | --- |
| Geometry | 1-D gentle slope |
| Thermal-energy transport | Enabled |
| Bed adjustment by deposition | Disabled by default |
| x-resolution, solids fraction, temperature | Generator arguments |

## How to run

Generate with `./create_example.py 400 0.01 900 false` (cells, solids volume fraction, temperature, plotting flag), then execute `../../bin/IMEX_SfloW2D`.

## Expected output

Native NetCDF (`<run_name>.nc`), viewable in ParaView.

## References

1. Bursik, M. I., and Woods, A. W. (1996). The dynamics and thermodynamics of large ash flows. *Bulletin of Volcanology*, 58, 175–193.

## Legacy documentation and preservation

The original [`README.txt`](README.txt) in this example directory remains the archival source for the complete historical narrative, parameter discussion, scientific references and commands. It has not been deleted or rewritten. This standardized README must not replace it as the sole provenance source until a full section-by-section textual comparison has been accepted.
