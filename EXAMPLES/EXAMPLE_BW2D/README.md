---
title: "Two-Dimensional Pulsed Bottom-Source Flow"
summary: "Thermally evolving, radially spreading flow supplied from a circular bottom source."
category: "Volcanic flows"
tags:
  - "thermal energy"
  - "bottom source"
  - "pulses"
---

# Two-Dimensional Pulsed Bottom-Source Flow

## Overview

A 2-D radial flow with basal friction, sedimentation and entrainment above sloping terrain, based on Bursik and Woods (1996). In contrast to a conventional radial inlet, mixture enters vertically from the bottom of a circular region without horizontal momentum; sedimentation and entrainment are disabled within the source area. At steady state the outward mass flux matches the bottom input flux. Retained thermal-energy transport is enabled; topography is fixed by `TOPO_CHANGE_FLAG=F`.

## Main parameters

| Parameter | Value |
| --- | --- |
| Domain | 60 km × 60 km in documented invocation |
| Grid | 300 × 300 cells in documented invocation |
| Source radius | 2000 m |
| Slope | 2 degrees |
| Source temperature | 330 K |
| Source velocity | 100 m/s |
| Solid volume fraction | 8 × 10^-5 |
| Particle size | 2.5 × 10^-4 m |

## Dependencies

The input generator requires Python 3 and NumPy.

## How to run

Run `python create_example.py 300 30000.0 2000.0 2.0 330.0 100.0 0.00008 2.5e-4`, then `../../bin/IMEX_SfloW2D`. The eight arguments are number of x cells, half-domain size (m), source radius (m), slope (degrees), source temperature (K), velocity (m/s), solid volume fraction and particle size (m). The original README incorrectly states six arguments.

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc`, viewable in ParaView.

## Source pulse timing in the original documentation

The original README contains two descriptions of the four `TIME_PARAM` settings. An earlier description uses a 30 s period and 20 s duration; the later, explicit generated namelist states:

```fortran
&RADIAL_SOURCE_PARAMETERS
 X_SOURCE=0.0D0,
 Y_SOURCE=0.0D0,
 R_SOURCE=2000.0D0,
 VEL_SOURCE=100.0D0,
 T_SOURCE=330.0D0,
 ALPHAS_SOURCE=8e-05,
 TIME_PARAM=600.D0, 100.D0, 50.0D0, 50.D0,
/
```

The generated template sets `BOTTOM_RADIAL_SOURCE_FLAG=T`, so `TIME_PARAM(4)` is a phase offset, not a cut-off. With the explicit values, the 50 s phase equals the ramp duration: the source starts at its peak, declines linearly to zero over the first 50 s, and stays off for the remainder of the 600 s period. The earlier 30 s/20 s prose does not describe this generated configuration.

For steady-state injection, remove the `TIME_PARAM` line from the radial-source namelist.

## References

1. Bursik, M. I., and Woods, A. W. (1996). The dynamics and thermodynamics of large ash flows. *Bulletin of Volcanology*, 58, 175–193.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 24 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.

## Documentation verification note

**Source timing.** The documented invocation supplies **eight** arguments. The template-generated `TIME_PARAM = 600, 100, 50, 50` describes a 600 s period, a 100 s active window, a 50 s linear ramp and a 50 s phase offset. The contradictory 30/20 s prose is historical, not the active configuration.
