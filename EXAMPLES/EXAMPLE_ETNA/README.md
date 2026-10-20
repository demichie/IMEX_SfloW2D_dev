---
title: "Etna Pyroclastic Avalanche with Voellmy–Salm Rheology"
summary: "Avalanche on real Etna terrain with an initial volume defined by ellipsoid–topography intersection."
category: "Volcanic flows"
tags:
  - "Etna"
  - "Voellmy-Salm"
  - "DEM"
  - "avalanche"
---

# Etna Pyroclastic Avalanche with Voellmy–Salm Rheology

## Overview

A pyroclastic avalanche is simulated over real Mount Etna terrain with Voellmy–Salm rheology. The initial mass occupies a niche obtained by intersecting the terrain with an ellipsoid. Two user-selected points in the x–y plane define the lowest and highest positions of the volume and two ellipsoid semi-axes; a third semi-axis is supplied. The preprocessing script searches the niche concavity to obtain a desired volume, creating both the initial thickness and the modified terrain.

## Main parameters

| Parameter | Value |
| --- | --- |
| Rheology | Voellmy–Salm |
| Topography | `DEM/Etna2014_crop.zip` |
| Initial shape | Ellipsoidal excavation |
| Example target volume | 250000 m³ (configurable) |
| Example semi-width | 150 m (configurable) |

## Dependencies

Requires Python 3, NumPy, and pandas for preprocessing. Matplotlib is required only for the optional raster plot.

## How to run

From this directory, extract the DEM archive into `DEM/`, edit `input_ellipsoid.py` and run `python create_input_ellipsoid.py`. Then execute `../../bin/IMEX_SfloW2D`. For cleaning, the tracked script is `./cleanFolder.sh` (the original README incorrectly wrote `./cleanfolder.sh`). For raster plotting, create a symbolic link with `ln -s ../../UTILS/plot_overlay.py .` and run `python plot_overlay.py`.

## Expected output

Native NetCDF is enabled, and the legacy README describes ESRI ASCII `.asc` products and their visualization using GIS software or `plot_overlay.py`.

## Original ellipsoid configuration example

The values explicitly supplied in the original README are:

```python
DEM_folder = './DEM/'
DEM_file = 'Etna2014_crop.asc'

x1 = 500300
y1 = 4177600
x2 = 500430
y2 = 4177500
semi_width = 150.0
vol = 250000
```

These are documentation example values, not a replacement for checking `input_ellipsoid.py` and the current model input. The method uses two points to establish the longitudinal orientation and extents, user-chosen third semi-axis, and an automatic concavity search to match the specified volume.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 24 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.

## Documentation verification note

**Cleanup command.** The tracked cleanup utility is `cleanFolder.sh` (case sensitive); the old `cleanfolder.sh` spelling is not the tracked filename.
