---
title: "Supercritical Flow over a Slope Discontinuity"
summary: "Supercritical incoming flow traversing a slope, flat region and bed discontinuity."
category: "Numerical benchmarks"
tags:
  - "supercritical"
  - "topography"
  - "advection"
---

# Supercritical Flow over a Slope Discontinuity

## Overview

Supercritical inflow (`Ri<1`) enters from the left and passes an initial slope, a flat region, and a topographic discontinuity. Sedimentation and entrainment are absent, basal friction is disabled (`RHEOLOGY_FLAG=F`), and the retained thermal-energy equation is advected by the hyperbolic operator.

## Main parameters

| Parameter | Value |
| --- | --- |
| Flow regime | Supercritical, Ri < 1 |
| Friction | Disabled |
| Sedimentation / entrainment | Disabled |
| Thermal energy | Hyperbolic advection |

## Dependencies

Requires Python 3 and NumPy. Matplotlib is required only when the plotting flag is `true`.

## How to run

Generate with `./create_example.py 400 1.0 300 false` (x cells, solids volume fraction, temperature in K, plot flag), then run `../../bin/IMEX_SfloW2D`.

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` is viewable in ParaView.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 24 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.
