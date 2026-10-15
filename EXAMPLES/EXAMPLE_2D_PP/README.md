---
title: "Granular Avalanche with Pore Pressure"
summary: "Inclined granular avalanche with transported pore pressure and Voellmy-Salm-based rheology."
category: "Granular flows"
tags:
  - "pore pressure"
  - "avalanche"
  - "Voellmy-Salm"
---

# Granular Avalanche with Pore Pressure

## Overview

An avalanche initially held within a hemispherical shell slides down a 35-degree inclined plane into a horizontal runout plane, as in the geometric setup of Wang et al. (2004). Sedimentation and air entrainment are disabled. A pore-pressure transport equation is solved (Gueugneau et al., 2017, Equation 12), and pore pressure affects the Voellmy-Salm-based friction model.

## Main parameters

| Parameter | Value |
| --- | --- |
| Slope | 35 degrees |
| Initial geometry | Hemispherical shell |
| Pore pressure | Transported |
| Sedimentation / air entrainment | Disabled |
| Solid volume fraction, temperature, resolution | Preprocessing arguments |

## How to run

Run `python create_example.py 100 0.5 300 false`, then `../../bin/IMEX_SfloW2D`. Script arguments: x cells, solid volume fraction, temperature (K), plot flag. NumPy and Basemap are listed by the legacy README.

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` is viewable in ParaView.

## References

1. Wang, Y., Hutter, K., and Pudasaini, S. P. (2004). The Savage-Hutter theory: A system of partial differential equations for avalanche flows of snow, debris, and mud. *ZAMM*, 84(8), 507–527. DOI: 10.1002/zamm.200310123.
2. Gueugneau, V., Kelfoun, K., Roche, O., and Chupin, L. (2017). Effects of pore pressure in pyroclastic flows: numerical simulation and experimental validation. *Geophysical Research Letters*, 44(5), 2194–2202.

## Legacy documentation and preservation

The original [`README.txt`](README.txt) in this example directory remains the archival source for the complete historical narrative, parameter discussion, scientific references and commands. It has not been deleted or rewritten. This standardized README must not replace it as the sole provenance source until a full section-by-section textual comparison has been accepted.
