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

From this directory, generate the input files with `python create_example.py 100 0.5 300 false`, then run `../../bin/IMEX_SfloW2D`. The generator requires NumPy; the legacy README also lists Basemap. Arguments: number of x cells, initial solid volume fraction, temperature (K), plot flag.

## Expected output

The simulation produces native NetCDF output (`<run_name>.nc`), viewable in ParaView.

The initial conditions and topography resemble Example 4.1 from Wang et al. (2004), but the rheological model is not the same as the original reference. The original README calls the script `create_example4.py` in its description, although the actual script and command are `create_example.py`.

## References

1. Wang, Y., Hutter, K., and Pudasaini, S. P. (2004). The Savage-Hutter theory: A system of partial differential equations for avalanche flows of snow, debris, and mud. *ZAMM*, 84(8), 507–527. DOI: 10.1002/zamm.200310123.
2. Bursik, M. I., and Woods, A. W. (1996). The dynamics and thermodynamics of large ash flows. *Bulletin of Volcanology*, 58, 175–193.

## Legacy documentation and preservation

The original [`README.txt`](README.txt) in this example directory remains the archival source for the complete historical narrative, parameter discussion, scientific references and commands. It has not been deleted or rewritten. This standardized README must not replace it as the sole provenance source until a full section-by-section textual comparison has been accepted.

## Documentation verification note

**Generator command.** The tracked Python generator is `create_example.py`; the legacy mention of `create_example4.py` is a stale filename.
