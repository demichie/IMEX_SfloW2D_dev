---
title: "Dilute PDC Collapse at Campi Flegrei"
summary: "Collapse of a dilute gas–particle mixture over real Campi Flegrei topography."
category: "Volcanic flows"
tags:
  - "Campi Flegrei"
  - "PDC"
  - "collapse"
---

# Dilute PDC Collapse at Campi Flegrei

## Overview

Dilute gas–particle mixture collapsing over the Campi Flegrei area. Initial conditions are similar to those studied by Ongaro, Orsucci and Cornolti.

## Main parameters

| Parameter | Value |
| --- | --- |
| Grid | 350 × 260 cells |
| Cell size | 50 m |
| Simulation duration | 150 s |
| Collapse radius | 500 m |
| Collapse height | 400 m |
| Initial solid volume fraction | 0.015 |
| Topography | `DEM/topo.zip` |

## Source integration

The table follows ascending **Group** numbers (the order of Lie source updates after spatial transport). **Scope** identifies the layer or layer pair affected; **Method** is the configured numerical treatment; **Status** distinguishes enabled sources from settings that do not activate a process; **Role** explains the physical contribution. An entry marked *Configured only* is not an executed source group. Group 0 denotes the initial explicit source stage, not a time-integration stage count.

| Group | Scope | Method | Status | Role |
| ---: | --- | --- | --- | --- |
| 0 | Flow | `EXPLICIT` | Active | Terrain-curvature source contributes to momentum. |
| 1 | Flow | `NEWTON` | Configured only | Quadratic rheology is disabled. |
| 2 | Flow | `NEWTON` | Configured only | Cooling is not enabled. |
| 4 | Flow | `EXPLICIT` | Active | Particle settling removes suspended solids into deposits. |

## How to run

From this example directory, execute the following commands in order:

**Extract the topography archive:**

```bash
unzip DEM/topo.zip -d DEM
```

The input uses `DEM/topography_dem.asc`; unpack the archive only into `DEM/`.

**Run the solver:**

```bash
../../src/IMEX_SfloW2D
```


## Expected output

The committed input enables native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` and disables ESRI raster output (`OUTPUT_ESRI_FLAG=F`). The older README's `.asc` description does not match the current input; ESRI output requires enabling that flag.

## References

1. Esposti Ongaro, T., Orsucci, S., and Cornolti, F. A fast, calibrated model for pyroclastic density currents kinematics and hazard. *Journal of Volcanology and Geothermal Research*, 327, 257–272.

## Documentation provenance

This Markdown guide supersedes the earlier plain-text instructions. Differences in historical commands or parameters are documented here and in [the examples documentation audit](../README_SECOND_PASS_AUDIT.md).
