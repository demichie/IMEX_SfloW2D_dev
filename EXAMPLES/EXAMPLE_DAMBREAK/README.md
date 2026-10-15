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
| Hydraulic permeability | 1.10 × 10^-11 (may be overridden dynamically) |
| GAS_LOSS_FLAG | T |
| F_INHIBIT_MODE | DYNAMIC |
| DYNAMIC_PERMEABILITY_FLAG | T |

## How to run

Generate using `python create_dambreak.py 300 15 3.0 0.15 0.20 0.40 0.58 300 1.0 true`. Arguments: NX, NY, channel length (m), channel width (m), dam length (m), dam height (m), solids volume fraction, ambient temperature (K), hydrostatic pore-pressure fraction, plot flag. From the example directory, create a solver link (`ln -s ../../bin/IMEX_SfloW2D .`) and run `./IMEX_SfloW2D IMEX_SfloW2D.inp`. For centerline postprocessing, use `python height_summary_nc.py OUTPUT/dambreak2D/fields/dambreak2D.nc --center_fraction 0.5 --animate`. Clean with `python clean_folder.py` (dry run), followed by `python clean_folder.py --apply`.

## Expected output

The generator writes `IMEX_SfloW2D.inp`, `topography_dem.asc`, an initial conservative-state file and optional preview. The run writes `OUTPUT/dambreak2D/fields/dambreak2D.nc`; `height_summary_nc.py` produces `frames/dambreak2D_hcenter_####.png` and an MP4 (or GIF fallback). For MP4, install `imageio[ffmpeg]`, `imageio-ffmpeg`, and system `ffmpeg`.

## Pore-pressure and permeability options

`F_INHIBIT_MODE` selects the inhibition factor for gas escape as solid fraction approaches maximum packing:

- `OFF`: inhibition is always 1, even above maximum packing.
- `STEP`: discontinuous switch at maximum packing.
- `STATIC`: smooth transition independent of height, starting at user-provided `ALPHA_TRANS`.
- `DYNAMIC`: height-dependent smooth transition, no manually set `ALPHA_TRANS` needed.

The original README warns that `GAS_LOSS_FLAG=T` is required to avoid nonphysical pressure decrease without gas loss. When `DYNAMIC_PERMEABILITY_FLAG=T`, the prescribed permeability is ignored and a local value is derived from solids fraction and Sauter diameter. The original text mentions both `1.10D-11` and `5.85e-12` as supplied values; this inconsistency is retained here as an audit note for verification.

## Provenance and documentation notes

The original README attributes the compaction/permeability formulation to Claudia Elijas-Parra and the channel setup to ECP Breard. It also uses historical names `create_dambreak_static.py` and `create_example_Roche_static.py` and describes a two-dimensional long rectangular channel with impermeable side walls and a no-flux right boundary.

## References

Barker and Gray (2017), cited in the original README for regularized μ(I). Complete bibliographic details require verification.

## Legacy documentation and preservation

The original [`README.txt`](README.txt) in this example directory remains the archival source for the complete historical narrative, parameter discussion, scientific references and commands. It has not been deleted or rewritten. This standardized README must not replace it as the sole provenance source until a full section-by-section textual comparison has been accepted.

## Documentation verification note

**Generator filename.** The existing setup program is `create_dambreak.py` (the names `create_dambreak_static.py` and `create_example_Roche_static.py` are historical references in the source/README). The supplied CLI has 10 positional arguments. Any alternative names in embedded docstrings are not working paths.
