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
| Temperature | Generator argument (K) |
| Initial Richardson number | Generator argument, < 1 |
| Solid mass fraction (`XS_SOURCE`) | 0.01 in the documented invocation |

## Dependencies

The input generator requires Python 3 and NumPy.

## How to run

Run `python create_example.py 200 250 10000000 900 0.5 0.01` and then `../../bin/IMEX_SfloW2D`. The six arguments are x cells, source extent (m), total source mass-flow rate (kg/s), temperature (K), Richardson number and solid mass fraction. The generator puts the third argument in `MFR_SOURCE`; the solver derives the source thickness from that rate. The sixth argument is written to `XS_SOURCE`, a mass fraction. The legacy README's "initial thickness" description for the third argument and "kg/s" units for the last argument are errors.

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` is enabled; open it using ParaView. File creation alone does not establish a successful source run. With the documented 900 K invocation, the source mixture density is below ambient density, giving negative reduced gravity. The mass-flow/Richardson-number setup would evaluate `sqrt(RI_SOURCE / reduced_gravity)` while deriving source thickness. The solver now rejects this combination before that calculation, prints an explicit error and exits unsuccessfully; the 900 K setup remains physically invalid.

The earlier zero-mass result at 300 K was caused by the `froude_blend` reconstruction replacing the prescribed boundary-source depth with a zero reconstructed depth. The solver now preserves that prescribed trace when computing wave speeds, advective fluxes and the hydrostatic path. With the same diagnostic parameters as the documented command except `T_SOURCE=300 K`, the run reached 100 s, wrote a nonempty NetCDF file and produced nonzero mass. The final reported mixture mass was about `1.331e9 kg`; air entrainment is enabled, so total mixture mass can exceed the imposed `MFR_SOURCE * T_END` of `1.0e9 kg`. The solid mass was about `9.764e6 kg`.

This diagnostic does not validate the documented 900 K configuration or establish that 300 K is a scientifically approved replacement. The 900 K input remains invalid because it gives negative source reduced gravity and non-finite source properties.

### Original explanatory note

The legacy README says the example uses a simple friction model "as in [1]", but no corresponding bibliography entry is supplied. Do not infer a specific reference until it has been checked.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 24 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.

## Documentation verification note

**Parameter interpretation.** The generator maps its third argument to `MFR_SOURCE`; the solver reports this rate in kg/s and derives `H_SOURCE` from it. The final argument maps to `XS_SOURCE`, which is a source mass fraction. The documented 900 K input yields negative source reduced gravity. The separate 300 K diagnostic now injects mass after the boundary-trace fix, but is not an approved physical replacement.
