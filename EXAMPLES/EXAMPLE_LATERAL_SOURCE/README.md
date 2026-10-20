---
title: "Supercritical Lateral Source"
summary: "Lateral supercritical mass-flux source on flat terrain with sedimentation and entrainment."
category: "Source-driven flows"
tags:
  - "lateral source"
  - "supercritical"
  - "entrainment"
---

# Supercritical Lateral Source

## Overview

Supercritical lateral inflow of a gas–particle mixture on flat terrain. Sedimentation and air entrainment are modeled using a simple friction law. The source must be supercritical (`Ri<1`) so that no characteristic returns upstream into the prescribed inflow.

## Main parameters

| Parameter | Value |
| --- | --- |
| Source extent | 250 m in the documented invocation |
| Source mass-flow rate | 10,000,000 kg/s in the documented invocation |
| Temperature | 300 K in the documented invocation (ambient temperature) |
| Initial Richardson number | Generator argument, < 1 |
| Solid mass fraction (`XS_SOURCE`) | 0.01 in the documented invocation |

## Dependencies

The input generator requires Python 3 and NumPy.

## How to run

Run `python create_example.py 200 250 10000000 300 0.5 0.01` and then `../../bin/IMEX_SfloW2D`. The six arguments are x cells, source extent (m), total source mass-flow rate (kg/s), temperature (K), Richardson number and solid mass fraction. The generator puts the third argument in `MFR_SOURCE`; the solver derives the source thickness from that rate. The sixth argument is written to `XS_SOURCE`, a mass fraction. The legacy README's "initial thickness" description for the third argument and "kg/s" units for the last argument are errors.

## Expected output

The documented 300 K invocation uses the same temperature as the 300 K ambient air. With `XS_SOURCE=0.01`, the solids make the source mixture slightly denser than ambient: approximately 1.188 kg/m³ versus 1.176 kg/m³, giving positive reduced gravity of about 0.098 m/s². The run reaches 100 s, writes nonempty NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc`, and produces nonzero mass. Air entrainment is enabled, so the final mixture mass can exceed the imposed `MFR_SOURCE * T_END`.

The earlier zero-mass result at 300 K was caused by the `froude_blend` reconstruction replacing the prescribed boundary-source depth with a zero reconstructed depth. The solver now preserves that prescribed trace when computing wave speeds, advective fluxes and the hydrostatic path. The validated 300 K run reported final mixture mass of about `1.331e9 kg` and solid mass of about `9.764e6 kg`.

The former 900 K invocation is not retained as the active setup: it gives negative source reduced gravity and is rejected by the solver when a Richardson number is prescribed.

### Original explanatory note

The legacy README says the example uses a simple friction model "as in [1]", but no corresponding bibliography entry is supplied. Do not infer a specific reference until it has been checked.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record, including its historical 900 K invocation. This guide uses the approved 300 K setup instead; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.

## Documentation verification note

**Parameter interpretation.** The generator maps its third argument to `MFR_SOURCE`; the solver reports this rate in kg/s and derives `H_SOURCE` from it. The final argument maps to `XS_SOURCE`, which is a source mass fraction. For the documented pressure, ambient temperature and 1% solids, the neutral-density threshold is approximately 303.03 K; 300 K is selected to match ambient temperature and yields positive reduced gravity. The former 900 K invocation yields negative reduced gravity and is invalid for the prescribed positive Richardson-number setup.
