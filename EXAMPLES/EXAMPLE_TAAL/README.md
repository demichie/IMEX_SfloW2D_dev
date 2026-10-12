---
title: "Pulsed Radial PDC Source at Taal"
summary: "Supercritical radial dilute gas–particle mixture spreading over Taal Volcano topography."
category: "Volcanic flows"
tags:
  - "Taal"
  - "radial source"
  - "pulses"
  - "DEM"
---

# Pulsed Radial PDC Source at Taal

## Overview

A supercritical radial source injects a dilute gas–particle mixture onto the topography of the Taal volcanic complex, Philippines, through repeated pulses until the source is stopped.

## Main parameters

| Parameter | Value |
| --- | --- |
| Grid | 500 × 700 |
| Cell size | 50 m |
| Duration | 200 s |
| Source coordinates | (282345, 1548428) |
| Source thickness / radius | 40 m / 250 m |
| Source velocity / temperature | 50 m/s / 350 K |
| Source solids volume fraction | 0.02 |
| `TIME_PARAM` | (30, 20, 4, 150) s |

## How to run

Create the local executable link with `ln -s ../../bin/IMEX_SfloW2D .` and run `./IMEX_SfloW2D` from the example directory. The input and `topography_dem.asc` are already tracked.

## Expected output

The committed input uses `OUTPUT_NETCDF_FLAG=T` and `OUTPUT_ESRI_FLAG=F`; NetCDF is the configured field output.


Native NetCDF output is enabled (`OUTPUT_NETCDF_FLAG=T`), while `OUTPUT_ESRI_FLAG=F` for the current input.

## Pulse semantics (preserved)

```fortran
&RADIAL_SOURCE_PARAMETERS
 X_SOURCE=282345.0D0,
 Y_SOURCE=1548428.0D0,
 H_SOURCE=40.0D0,
 R_SOURCE=250.0D0,
 VEL_SOURCE=50.0D0,
 T_SOURCE=350.0D0,
 ALPHAS_SOURCE=0.02D0,
 TIME_PARAM=30.D0,20.D0,4.0D0,150.D0,
/
```

The four timing parameters mean: time between cycles (30 s); injection duration per cycle (20 s, followed by 10 s off); ramp time for increasing/decreasing source intensity (4 s), scaling height and velocity to preserve Richardson number; and total source activation duration (150 s). For steady injection remove the timing line. The legacy README says to make the period longer than the source duration for a single pulse and gives `TIME_PARAM=100.D0,20.D0,0.D0,200.D0`; it claims this gives one 20 s pulse, but the parameter relationship should be checked because the fourth value exceeds the first.

## Legacy documentation and preservation

The original [`README.txt`](README.txt) in this example directory remains the archival source for the complete historical narrative, parameter discussion, scientific references and commands. It has not been deleted or rewritten. This standardized README must not replace it as the sole provenance source until a full section-by-section textual comparison has been accepted.

## Documentation verification note

**Input-file output flags.** The tracked input enables NetCDF but sets `OUTPUT_ESRI_FLAG=F`; therefore ESRI raster outputs should not be promised as enabled by default.
