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

The committed input enables native NetCDF output and disables ESRI raster output (`OUTPUT_ESRI_FLAG=F`).

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

For this lateral radial source, the four timing parameters mean: period between cycles (30 s); injection duration per cycle (20 s, followed by 10 s off); ramp time for increasing/decreasing source intensity (4 s), scaling height and velocity to preserve Richardson number; and absolute source cut-off time (150 s). Thus the tracked input produces five pulses before the source stops at 150 s; the simulation continues to 200 s. For steady injection remove the timing line. The legacy README says to make the period longer than the source duration for a single pulse, but its example `TIME_PARAM=100.D0,20.D0,0.D0,200.D0` has a 100 s period and a 200 s cut-off; under the current implementation it produces pulses at 0–20 s and 100–120 s, not a single pulse.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 21 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.
