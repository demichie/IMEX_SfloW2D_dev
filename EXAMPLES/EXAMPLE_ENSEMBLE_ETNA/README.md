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

Perform these operations in order from the
`EXAMPLES/EXAMPLE_ENSEMBLE_ETNA` directory. The shared DEM and
template executable link must exist before generating any ensemble
member.

### 1. Prepare the shared DEM and executable link

The suggested directory organization is:

```text
EXAMPLE_ENSEMBLE_ETNA/
├── DEM/
└── templatedir/
    └── DEM/
```

From `EXAMPLES/EXAMPLE_ENSEMBLE_ETNA`, create the executable link and unpack/link the shared DEM **before** generating any member directories:

```bash
cd templatedir
```

```bash
ln -s ../../../bin/IMEX_SfloW2D .
```

```bash
cd ..
```

```bash
cd DEM
```

```bash
unzip Etna2014_crop.zip
```

```bash
cd ../templatedir
```

```bash
mkdir DEM
```

```bash
cd DEM
```

```bash
ln -s ../../DEM/Etna2014_crop.asc .
```

```bash
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

### 2. Generate ensemble members and run simulations

The number of samples, uncertain parameters and their ranges are
defined in `create_ensemble.py`. Each parameter has a unique keyword
(for example `var1`); this example also combines two possible source
locations with different weights. The template contains placeholders
such as `ENSEMBLE_var1` which are substituted when creating member
directories. The `ENSEMBLE_vx` and `ENSEMBLE_vy` placeholders are
resolved later by member preprocessing.

**Step 2a — Sample uncertain parameters.** Sampling ranges and sample count are set
in the script.

```bash
python create_ensemble.py
```

**Step 2b — Create realization directories.** This copies `templatedir` into
`ensemble.xxxxx` directories and replaces the relevant `ENSEMBLE_var1`
placeholders.

```bash
python create_inputfiles.py
```

**Step 2c — Launch simulations.** Configure `max_processes` in the launcher. It waits
for every realization and reports any failed run.

```bash
python launch_jobs.py
```

**Step 2d — Aggregate probability rasters.** Grouping realizations with pandas
`groupby` is optional.

```bash
python create_prob_maps.py
```

**Step 2e — Plot the ensemble results.**

```bash
python plot_ensemble.py
```

**Step 2f — Optionally compute exceedance maps** from maximum-thickness rasters.

```bash
python create_hazard_maps.py
```


## Expected output

Each realization writes solver output beneath `ensemble.xxxxx/OUTPUT/<RUN_NAME>/`: snapshots in indexed directories under `snapshots/`, maximum-thickness products and runout tables in `diagnostics/`, vulnerability rasters in `vulnerability/`, terrain exports in `terrain/`, and input backups/metadata in `metadata/`. The stock template sets `OUTPUT_CONS_FLAG=F`, so it does not write NetCDF fields unless that flag is enabled. The aggregation scripts read the per-realization outputs and produce probability rasters (`.asc`) in `result.xxx`; `plot_ensemble.py` makes PNG summary maps, and the optional hazard-map script writes to `exceedance_maps/`.

The template currently sets `RUN_NAME="exampleEtna"`, so the default output root for each member is `ensemble.xxxxx/OUTPUT/exampleEtna/`.


## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 21 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.

## Documentation verification note

**Run order and placeholders.** `create_inputfiles.py` substitutes sampled values for matching `ENSEMBLE_` tokens (not the typo `ENSAMBLE_`). The `ENSEMBLE_vx` and `ENSEMBLE_vy` tokens are intentionally resolved later by `create_input_ellipsoid.py`, invoked by each member's `run.sh`; they are expected to remain after `create_inputfiles.py` alone. The generation sequence retains symbolic links, template copies, per-member runs, raster aggregation and probability plotting.
