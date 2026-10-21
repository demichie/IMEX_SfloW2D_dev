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

## Dependencies

Requires Python 3 and the IMEX_SfloW2D executable. NumPy and pandas are used for sampling, input generation and aggregation; Matplotlib is used for plots. The optional `create_hazard_maps.py` additionally requires netCDF4 and GDAL's Python bindings (`osgeo`).

## How to run

From `EXAMPLE_ENSEMBLE_ETNA`, follow the sequence below.

1. Link `../../../bin/IMEX_SfloW2D` as `templatedir/IMEX_SfloW2D`.
2. Unzip `DEM/Etna2014_crop.zip`, create `templatedir/DEM/`, and link its extracted raster into that folder. The original README's commands are reproduced in the workflow notes below.
3. Execute `python create_ensemble.py` to generate the CSV of sampled uncertain parameters. Sampling ranges and sample count are configured in this script.
4. Execute `python create_inputfiles.py` to create `ensemble.xxxxx` directories by copying `templatedir` and replacing placeholders of the form `ENSEMBLE_var1`.
5. Execute `python launch_jobs.py` to run the `run.sh` script in each member directory. Configure `max_processes` in that script. The launcher waits for every member and exits with an error if any run fails.
6. Execute `python create_prob_maps.py` to aggregate raster products, optionally grouping realizations with pandas `groupby`.
7. Execute `python plot_ensemble.py` to make PNG maps.
8. Optionally execute `python create_hazard_maps.py` to calculate exceedance maps from maximum-thickness rasters.

## Expected output

Each realization writes solver output beneath `ensemble.xxxxx/OUTPUT/<RUN_NAME>/`: snapshots in indexed directories under `snapshots/`, maximum-thickness products and runout tables in `diagnostics/`, vulnerability rasters in `vulnerability/`, terrain exports in `terrain/`, and input backups/metadata in `metadata/`. The stock template sets `OUTPUT_CONS_FLAG=F`, so it does not write NetCDF fields unless that flag is enabled. The aggregation scripts read the per-realization outputs and produce probability rasters (`.asc`) in `result.xxx`; `plot_ensemble.py` makes PNG summary maps, and the optional hazard-map script writes to `exceedance_maps/`.

The template currently sets `RUN_NAME="exampleEtna"`, so the default output root for each member is `ensemble.xxxxx/OUTPUT/exampleEtna/`.

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
 RHEOLOGY_MODEL = 'voellmy',
 MU=ENSEMBLE_mu,
 XI=ENSEMBLE_xi
/
```

`templatedir/run.sh` orchestrates the ellipsoidal-volume preprocessing and solver execution. Its optional per-run plotting command is commented out by default; enable it *before* copying the template folder to run the overlay postprocessor for each realization. The original documentation skips step 7 in its numbering and contains the misspellings `ENSAMBLE` and `create_ensamble.py`; these are documentation errors, not names to adopt as new conventions.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 21 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.

## Documentation verification note

**Run order and placeholders.** `create_inputfiles.py` substitutes sampled values for matching `ENSEMBLE_` tokens (not the typo `ENSAMBLE_`). The `ENSEMBLE_vx` and `ENSEMBLE_vy` tokens are intentionally resolved later by `create_input_ellipsoid.py`, invoked by each member's `run.sh`; they are expected to remain after `create_inputfiles.py` alone. The generation sequence retains symbolic links, template copies, per-member runs, raster aggregation and probability plotting.
