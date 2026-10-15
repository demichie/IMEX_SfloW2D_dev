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

## How to run

Generate with `./create_example.py 400 1.0 300 false` (x cells, solids volume fraction, temperature in K, plot flag), then run `../../bin/IMEX_SfloW2D`.

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` is viewable in ParaView.

## Legacy documentation and preservation

The original [`README.txt`](README.txt) in this example directory remains the archival source for the complete historical narrative, parameter discussion, scientific references and commands. It has not been deleted or rewritten. This standardized README must not replace it as the sole provenance source until a full section-by-section textual comparison has been accepted.
