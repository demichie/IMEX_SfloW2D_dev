---
title: "Etna Pyroclastic Avalanche Ensemble"
summary: "Ensemble workflow for sampling uncertain source and rheology parameters and generating probabilistic hazard maps."
category: "Ensemble simulations"
tags:
  - "Etna"
  - "uncertainty"
  - "hazard maps"
  - "ensemble"
---

# Etna Pyroclastic Avalanche Ensemble

## Overview

This example sets up multiple pyroclastic-avalanche simulations at Etna, samples uncertain initial and rheological parameters (including two locations with different weights), runs a collection of realizations, and aggregates outputs into probabilistic thickness and dynamic-pressure-related maps. Its directory hierarchy minimizes repeated copies of DEMs by using symbolic links.

## Main parameters

| Parameter | Value |
| --- | --- |
| Site | Mount Etna |
| Input DEM | `DEM/Etna2014_crop.zip` |
| Parameter sampling | `create_ensemble.py` |
| Case generation | `create_inputfiles.py` |
| Parallel launch | `launch_jobs.py`, configurable `max_processes` |
| Hazard aggregation | `create_prob_maps.py`, `plot_ensemble.py` |

## How to run

From `EXAMPLE_ENSEMBLE_ETNA`, follow the sequence below.

1. Link `../../../bin/IMEX_SfloW2D` as `templatedir/IMEX_SfloW2D`.
2. Unzip `DEM/Etna2014_crop.zip`, create `templatedir/DEM/`, and link its extracted raster into that folder. The original README's commands are reproduced in the workflow notes below.
3. Execute `python create_ensemble.py` to generate the CSV of sampled uncertain parameters. Sampling ranges and sample count are configured in this script.
4. Execute `python create_inputfiles.py` to create `ensemble.xxxxx` directories by copying `templatedir` and replacing placeholders of the form `ENSEMBLE_var1`.
5. Execute `python launch_jobs.py` to run the `run.sh` script in each member directory. Configure `max_processes` in that script.
6. Execute `python create_prob_maps.py` to aggregate raster products, optionally grouping realizations with pandas `groupby`.
7. Execute `python plot_ensemble.py` to make PNG maps.

## Expected output

The workflow produces per-realization solver outputs, probability rasters (`.asc`) for maximum thickness and dynamic-pressure/thickness thresholds, grouped directories `result.xxx`, and PNG summary maps.

## Original directory layout and setup commands

The legacy README recommends maintaining the directory organization:

```text
EXAMPLE_ENSEMBLE_ETNA/
├── DEM/
└── templatedir/
    └── DEM/
```

Its executable-link and DEM preparation instructions are:

```sh
cd templatedir
ln -s ../../../bin/IMEX_SfloW2D .
cd ..
cd DEM
unzip Etna2014_crop.zip
cd ../templatedir
mkdir DEM
cd DEM
ln -s ../../DEM/Etna2014_crop.asc .
cd ../..
```

The original README specifies substitution tokens in the form `ENSEMBLE_xyz` and gives the following parameter example (retained here because it explains the scientific interpretation of the sampling):

```fortran
&RESTART_PARAMETERS
 RESTART_FILE="./DEM/init_volume.asc",
 T_INIT=373.D0,
 T_AMBIENT=300.D0,
 SED_VOL_PERC=50.0,
 u_init=ENSEMBLE_vx,
 v_init=ENSEMBLE_vy,
/
&RHEOLOGY_PARAMETERS
 RHEOLOGY_MODEL=1,
 MU=ENSEMBLE_mu,
 XI=ENSEMBLE_xi
/
```

`templatedir/run.sh` orchestrates the ellipsoidal-volume preprocessing, solver execution and plotting. Postprocessing can be omitted by commenting out the corresponding lines *before* copying the template folder. The original documentation skips step 7 in its numbering and contains the misspellings `ENSAMBLE` and `create_ensamble.py`; these are documentation errors, not names to adopt as new conventions.

## Legacy documentation and preservation

The original [`README.txt`](README.txt) in this example directory remains the archival source for the complete historical narrative, parameter discussion, scientific references and commands. It has not been deleted or rewritten. This standardized README must not replace it as the sole provenance source until a full section-by-section textual comparison has been accepted.

## Documentation verification note

**Run order and placeholders.** The placeholder prefix used by `create_inputfiles.py` is `ENSEMBLE_` (not the typo `ENSAMBLE_`). The generation sequence retains symbolic links, template copies, per-member `run.sh`, raster aggregation and probability plotting.
