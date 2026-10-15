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
| Source extent | Generator argument (m) |
| Source mass-flow rate | Third generator argument; units to verify |
| Temperature | Generator argument (K) |
| Initial Richardson number | Generator argument, < 1 |
| Solid fraction | Generator argument; basis to verify |

## How to run

Run `python create_example.py 200 250 10000000 900 0.5 0.01` and then `../../bin/IMEX_SfloW2D`. The script accepts x cells, source extent, source mass flow rate, temperature, Richardson number and solids fraction. The original README also labels the third argument "initial thickness" and gives incompatible units for the last argument; those contradictions require second-pass verification.

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` is enabled; open it using ParaView.

### Original explanatory note

The legacy README says the example uses a simple friction model "as in [1]", but no corresponding bibliography entry is supplied. Do not infer a specific reference until it has been checked.

## Legacy documentation and preservation

The original [`README.txt`](README.txt) in this example directory remains the archival source for the complete historical narrative, parameter discussion, scientific references and commands. It has not been deleted or rewritten. This standardized README must not replace it as the sole provenance source until a full section-by-section textual comparison has been accepted.

## Documentation verification note

**Review mass-flow-rate parameter.** The third `create_example.py` argument is the **source mass-flow rate**, not initial thickness. The legacy README lists incompatible descriptions and units. Use the source mass-flow-rate designation; the precise SI unit and how the generator converts it to inflow state should be confirmed against the full script.
