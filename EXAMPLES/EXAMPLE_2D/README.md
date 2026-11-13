---
title: "Granular Avalanche on an Inclined Plane"
summary: "Finite granular avalanche descending a 35-degree slope onto a horizontal runout plane."
category: "Granular flows"
tags:
  - "avalanche"
  - "inclined plane"
  - "entrainment"
---

# Granular Avalanche on an Inclined Plane

## Overview

A finite granular mass initially confined in a hemispherical shell is released onto a plane inclined at 35 degrees, connected smoothly to a horizontal runout surface. The configuration resembles Wang et al. (2004), Example 4.1, but uses a different rheology, following Bursik and Woods (1996). Particle sedimentation and atmospheric-air entrainment are represented.

## Main parameters

| Parameter | Value |
| --- | --- |
| Slope | 35 degrees |
| Initial geometry | Hemispherical shell |
| Solid volume fraction | User-specified |
| Temperature | User-specified, K |
| Grid resolution | Set by preprocessing argument |

## How to run

From this example directory, first generate the inputs, then run the solver.

### 1. Generate the input files

The generator takes **four positional arguments**, in this order:

1. Number of computational cells along x.
2. Initial solid-particle volume fraction (dimensionless, 0–1).
3. Initial flow temperature (K).
4. Plot the initial solution (`true` or `false`).

The tracked generator is `create_example.py`, not the historical `create_example4.py`. Requires NumPy; Matplotlib is used when plotting. Basemap is no longer used.

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

The generated input enables both native NetCDF and ESRI raster output. NetCDF is written to `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` and can be viewed in ParaView; the ESRI rasters can be opened in a GIS.

The initial conditions and topography resemble Example 4.1 from Wang et al. (2004), but the rheological model is not the same as the original reference. The original README calls the script `create_example4.py` in its description, although the actual script and command are `create_example.py`.

## References

1. Wang, Y., Hutter, K., and Pudasaini, S. P. (2004). The Savage-Hutter theory: A system of partial differential equations for avalanche flows of snow, debris, and mud. *ZAMM*, 84(8), 507–527. DOI: 10.1002/zamm.200310123.
2. Bursik, M. I., and Woods, A. W. (1996). The dynamics and thermodynamics of large ash flows. *Bulletin of Volcanology*, 58, 175–193.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 21 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.

## Documentation verification note

**Generator command.** The tracked Python generator is `create_example.py`; the legacy mention of `create_example4.py` is a stale filename.

## Topography files

Run `python create_example.py` from this example directory. The generator creates
`DEM/topography_dem.asc` (and the directory if necessary). The generated
`IMEX_SfloW2D.inp` points to that path; the DEM is not written to the
example root. Keep numerical outputs outside `DEM/`.
