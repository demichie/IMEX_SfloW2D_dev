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
| Topography | `topo.zip` |

## How to run

First run `unzip topo.zip` in this directory, then `../../bin/IMEX_SfloW2D`.

## Expected output

The committed input sets `OUTPUT_ESRI_FLAG=F` and `OUTPUT_NETCDF_FLAG=T`. The older README describes ESRI `.asc` files; that description reflects a possible configuration rather than the default output flags currently tracked.


Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` is enabled by the input. The older README also mentions ESRI `.asc` files, although `OUTPUT_ESRI_FLAG=F` in the current tracked input.

## References

1. Esposti Ongaro, T., Orsucci, S., and Cornolti, F. A fast, calibrated model for pyroclastic density currents kinematics and hazard. *Journal of Volcanology and Geothermal Research*, 327, 257–272.

## Legacy documentation and preservation

The original [`README.txt`](README.txt) in this example directory remains the archival source for the complete historical narrative, parameter discussion, scientific references and commands. It has not been deleted or rewritten. This standardized README must not replace it as the sole provenance source until a full section-by-section textual comparison has been accepted.

## Documentation verification note

**Input-file output flags.** The tracked input enables NetCDF but sets `OUTPUT_ESRI_FLAG=F`; ESRI rasters are not enabled in the committed input even though the old README describes them.
