---
title: "Fluidized Granular Dam Break with Dynamic Permeability"
summary: "Two-dimensional dam-break channel experiment with pore pressure, compaction and regularized μ(I) rheology."
category: "Granular flows"
tags:
  - "dam break"
  - "pore pressure"
  - "mu(I)"
  - "compaction"
---

# Fluidized Granular Dam Break with Dynamic Permeability

## Overview

This long-channel dam break is designed to examine simultaneous compaction and gas diffusion in thin fluidized granular flows. It uses a regularized μ(I) rheology (Barker and Gray, 2017), pore-pressure loss, and a selectable inhibition function for the transition from unhindered to inhibited gas escape as solids packing approaches its maximum.

## Main parameters

| Parameter | Value |
| --- | --- |
| Rheology model | 11 (regularized μ(I)) |
| MU_S, MU_2 | 0.55, 0.73 |
| MUI_INF | 0.05 (0 would recover unregularized form) |
| I_0 | 0.279 |
| Hydraulic permeability | Template value 1.10 × 10^-11; unused with dynamic permeability enabled |
| GAS_LOSS_FLAG | T |
| F_INHIBIT_MODE | DYNAMIC |
| DYNAMIC_PERMEABILITY_FLAG | T |
| Momentum boundaries | Template `FLAG=1` is outflow-only in the current solver; it does not impose impermeable walls. |

## Dependencies

The setup generator requires Python 3 and NumPy; Matplotlib is required only when its plot flag is `true`. The centerline postprocessor requires NumPy, Matplotlib, and either netCDF4 or h5py to read the NetCDF file. MP4 output additionally requires `imageio[ffmpeg]`, `imageio-ffmpeg`, and system `ffmpeg`; without FFmpeg the script falls back to GIF.

## Physical-model namelists in the legacy example

The original documentation provides the following regularized
`mu(I)` coefficients (Barker and Gray, 2017):

```fortran
&RHEOLOGY_PARAMETERS
 RHEOLOGY_MODEL = 'mui',
 MU_S = 0.55D0,
 MU_2 = 0.73D0,
 MUI_INF = 0.05D0,
 I_0 = 0.279D0,
/
```

The value `MUI_INF=0.05` regularizes the model; the legacy text notes
that setting `MUI_INF=0` recovers the unregularized form.

The documented pore-pressure and permeability settings are:

```fortran
&PORE_PRESSURE_PARAMETERS
 hydraulic_permeability = 1.10D-11,
 GAS_LOSS_FLAG = T,
 ALPHA_TRANS = 0.2,
 F_INHIBIT_MODE = 'DYNAMIC',
 DYNAMIC_PERMEABILITY_FLAG = T,
/
```

These blocks reproduce the historical settings for reference; consult
the tracked input/template for the active run.


## How to run

Run these steps from the example directory, in the order shown.

### 1. Generate the dam-break inputs

The generator `create_dambreak.py` takes **ten positional arguments**:

1. **NX** — Number of computational cells along the channel (x direction).
2. **NY** — Number of computational cells across the channel (y direction).
3. **X_LENGTH** — Total channel length (m).
4. **Y_WIDTH** — Total channel width (m).
5. **BED_LENGTH** — Initial dam extent along x (m).
6. **BED_HEIGHT** — Initial dam height (m).
7. **ALFAS** — Initial solid-particle volume fraction (dimensionless, 0–1).
8. **T** — Initial/ambient temperature (K).
9. **PORE_PRESSURE_FRACT** — Initial pore-pressure fraction relative to hydrostatic pressure (0–1).
10. **PLOT_FLAG** — Whether to create a preview PNG (`true` or `false`).

The first command writes `IMEX_SfloW2D.inp`, `topography_dem.asc`,
the initial conservative state, and (when requested) the preview figure.
The legacy README mentions `create_dambreak_static.py`; the tracked
script to use is `create_dambreak.py`.

**Generator command:**

```bash
python create_dambreak.py 300 15 3.0 0.15 0.20 0.40 0.58 300 1.0 true
```

### 2. Link the compiled solver (first run only)

**Create the executable link:**

```bash
ln -s ../../bin/IMEX_SfloW2D .
```

### 3. Run the simulation

**Solver command:**

```bash
./IMEX_SfloW2D IMEX_SfloW2D.inp
```

The simulation writes `OUTPUT/dambreak2D/fields/dambreak2D.nc`,
viewable with ParaView, and snapshots under the run output directory.

### 4. Post-process the channel centerline

The postprocessor extracts the flow height at the specified fraction
across the channel: `0` is the lower side, `0.5` the centerline,
and `1` the upper side. `--animate` also creates an MP4 movie
(or GIF when FFmpeg is unavailable).

**Postprocessing command:**

```bash
python height_summary_nc.py OUTPUT/dambreak2D/fields/dambreak2D.nc --center_fraction 0.5 --animate
```

For MP4 output, install the optional Python packages:

```bash
pip install -U "imageio[ffmpeg]" imageio-ffmpeg
```

On Ubuntu/Debian, the historical README also recommends installing
the system FFmpeg executable:

```bash
sudo apt-get update && sudo apt-get install -y ffmpeg
```

### 5. Optional cleanup

Preview which generated files would be removed:

```bash
python clean_folder.py
```

Only after checking the dry-run output, apply the cleanup:

```bash
python clean_folder.py --apply
```


## Expected output

The generator writes `IMEX_SfloW2D.inp`, `topography_dem.asc`, an initial conservative-state file and optional preview. The run writes `OUTPUT/dambreak2D/fields/dambreak2D.nc`; `height_summary_nc.py` produces `frames/dambreak2D_hcenter_####.png` and an MP4 (or GIF fallback). For MP4, install `imageio[ffmpeg]`, `imageio-ffmpeg`, and system `ffmpeg`.

## Pore-pressure and permeability options

`F_INHIBIT_MODE` selects the inhibition factor for gas escape as solid fraction approaches maximum packing:

- `OFF`: inhibition is always 1, even above maximum packing.
- `STEP`: discontinuous switch at maximum packing.
- `STATIC`: smooth transition independent of height, starting at user-provided `ALPHA_TRANS`.
- `DYNAMIC`: height-dependent smooth transition derived from flow height. The parser still requires `ALPHA_TRANS > 0` even though this branch uses the height-derived transition; the template sets `ALPHA_TRANS=0.2`.

The original README warns that `GAS_LOSS_FLAG=T` is required to avoid nonphysical pressure decrease without gas loss. When `DYNAMIC_PERMEABILITY_FLAG=T`, the prescribed permeability is ignored and a local value is derived from solids fraction and Sauter diameter. The tracked template sets `hydraulic_permeability=1.10D-11`; the legacy text's `5.85e-12` is stale, and neither value controls this dynamic-permeability run.

## Provenance and documentation notes

The original README attributes the compaction/permeability formulation to Claudia Elijas-Parra and the channel setup to ECP Breard. It also uses historical names `create_dambreak_static.py` and `create_example_Roche_static.py`. Its claim of impermeable side walls and a no-flux right boundary conflicts with the tracked template: momentum boundaries use `FLAG=1`, which the current solver treats as outflow-only.

## References

Barker and Gray (2017), cited in the original README for regularized μ(I). Complete bibliographic details require verification.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 21 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.

## Documentation verification note

**Generator filename.** The existing setup program is `create_dambreak.py` (the names `create_dambreak_static.py` and `create_example_Roche_static.py` are historical references in the source/README). The supplied CLI has 10 positional arguments. Any alternative names in embedded docstrings are not working paths.
