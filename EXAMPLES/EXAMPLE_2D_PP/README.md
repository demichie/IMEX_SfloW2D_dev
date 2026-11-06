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

From this example directory, first generate the inputs, then run the solver.

### 1. Generate the input files

The generator takes **four positional arguments**, in this order:

1. Number of computational cells along x.
2. Initial solid-particle volume fraction (dimensionless, 0–1).
3. Initial flow temperature (K).
4. Plot the initial solution (`true` or `false`).

Requires Python 3 and NumPy; Matplotlib is needed only for plotting. The historical Basemap dependency is obsolete.

**Generator command:**

```bash
python create_example.py 100 0.5 300 false
```

### 2. Run the solver

**Solver command:**

```bash
../../bin/IMEX_SfloW2D
```

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` is viewable in ParaView.

## References

1. Wang, Y., Hutter, K., and Pudasaini, S. P. (2004). The Savage-Hutter theory: A system of partial differential equations for avalanche flows of snow, debris, and mud. *ZAMM*, 84(8), 507–527. DOI: 10.1002/zamm.200310123.
2. Gueugneau, V., Kelfoun, K., Roche, O., and Chupin, L. (2017). Effects of pore pressure in pyroclastic flows: numerical simulation and experimental validation. *Geophysical Research Letters*, 44(5), 2194–2202.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 21 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.
