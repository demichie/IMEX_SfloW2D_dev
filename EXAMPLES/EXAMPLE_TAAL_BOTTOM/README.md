---
title: "Bottom Elliptical PDC Source at Taal"
summary: "Gas–particle flow with pore-pressure evolution entering vertically from an elliptical bottom source."
category: "Volcanic flows"
tags:
  - "Taal"
  - "bottom source"
  - "pore pressure"
  - "pulses"
---

# Bottom Elliptical PDC Source at Taal

## Overview

A dilute gas–particle flow over Taal terrain, with an elliptical bottom-source area and vertical injection velocity. A pore-pressure equation is included. The elliptical source is parameterized by its major and minor semi-axes and its orientation.

## Main parameters

| Parameter | Value |
| --- | --- |
| Grid | 500 × 700, cells of 50 m |
| Duration | 50 s |
| Source center | (282345, 1548428) |
| Ellipse radii | 1500 m and 500 m |
| Orientation | 135 degrees |
| Source velocity / temperature | 10 m/s / 350 K |
| `XS_SOURCE` | 0.976 |
| Pore pressure | Enabled, fraction 0.5 |
| `TIME_PARAM` | Period 30 s, active 20 s, ramp 0 s, phase offset 150 s |

## How to run

Create an executable link with `ln -s ../../bin/IMEX_SfloW2D .`, then execute `./IMEX_SfloW2D` from this example directory.

## Expected output

The committed input enables NetCDF and runout diagnostics but disables ESRI raster output (`OUTPUT_ESRI_FLAG=F`).

## Elliptical source and pulses (preserved)

```fortran
&RADIAL_SOURCE_PARAMETERS
 X_SOURCE=282345.0D0,
 Y_SOURCE=1548428.0D0,
 R_SOURCE=1500.0D0,
 R2_SOURCE=500.0D0,
 ANGLE_SOURCE=135.0,
 VEL_SOURCE=10.0D0,
 T_SOURCE=350.0D0,
 XS_SOURCE=0.976D0,
 TIME_PARAM=30.D0,20.D0,0.D0,150.D0,
/
```

For the bottom radial source, `TIME_PARAM(1:3)` define the cycle period (30 s), injection duration (20 s) and ramp duration (0 s). Unlike the lateral source, `TIME_PARAM(4)` is a phase offset used in the periodic window, not an absolute cut-off or total source duration. The tracked value (150 s) is exactly five periods, so it causes no net phase shift: the source is active during 0–20 s and 30–50 s, until the simulation ends. Removing the timing line gives a steady source. The legacy README labels the fourth value a total duration; under the current bottom-source implementation its single-pulse example (`300,20,0,150`) instead shifts the pulse to 150–170 s, if the simulation runs that long. The physical interpretation of `XS_SOURCE` must not be conflated with `ALPHAS_SOURCE` without reading the model equation.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 24 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.
