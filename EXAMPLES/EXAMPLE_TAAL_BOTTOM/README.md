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
| `TIME_PARAM` | (30, 20, 0, 150) s |

## How to run

Create an executable link with `ln -s ../../bin/IMEX_SfloW2D .`, then execute `./IMEX_SfloW2D` from this example directory.

## Expected output

The committed input enables NetCDF and runout diagnostics but disables ESRI raster output (`OUTPUT_ESRI_FLAG=F`).


Native NetCDF and runout output are enabled; the current input disables ESRI rasters.

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

`TIME_PARAM` defines the cycle period (30 s), injection duration (20 s), ramp duration (0 s) and nominal total source duration (150 s). Removing it gives a steady source. The original README proposes `TIME_PARAM=300.D0,20.D0,0.D0,150.D0` for a single 20 s pulse. Note that the simulation itself ends at 50 s in the checked input, so the full nominal source duration is not simulated. The physical interpretation of `XS_SOURCE` must not be conflated with `ALPHAS_SOURCE` without reading the model equation.

## Legacy documentation and preservation

The original [`README.txt`](README.txt) in this example directory remains the archival source for the complete historical narrative, parameter discussion, scientific references and commands. It has not been deleted or rewritten. This standardized README must not replace it as the sole provenance source until a full section-by-section textual comparison has been accepted.

## Documentation verification note

**Input-file output flags.** The tracked input enables NetCDF and runout output, with `OUTPUT_ESRI_FLAG=F`; generated ESRI files should not be promised as default outputs.
