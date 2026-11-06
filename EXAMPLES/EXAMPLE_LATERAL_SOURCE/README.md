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
| Temperature | 900 K in the documented invocation |
| Initial Richardson number | Generator argument, < 1 |
| Solid volume fraction (`ALPHAS_SOURCE`) | 0.01 in the documented invocation |

## Dependencies

The input generator requires Python 3 and NumPy.

## How to run

From this example directory, first generate the inputs and then run the solver.

### 1. Generate the input files

The generator expects **6 positional arguments**, in this exact order:

1. Number of cells along x (positive integer).
2. Lateral source extent (m; positive).
3. Total source **mass-flow rate** (kg/s; positive), written to `MFR_SOURCE`. It is not a source thickness.
4. Source temperature (K; positive).
5. Initial Richardson number (dimensionless; below 1 for a supercritical source).
6. Initial solid **volume fraction** (dimensionless, 0–1), written to `ALPHAS_SOURCE`.

The legacy README inconsistently describes argument 3 as an initial thickness and argument 6 as a mass fraction. The active generator instead writes `MFR_SOURCE` and `ALPHAS_SOURCE`. `XS_SOURCE`, if used, would denote a solid mass fraction.

**Generator command:**

```bash
python create_example.py 200 250 10000000 900 0.5 0.01
```

### 2. Run the solver

**Solver command:**

```bash
../../bin/IMEX_SfloW2D
```

## Expected output

For the documented pressure, solid density and ambient conditions, the gas density at 900 K is approximately 0.392 kg/m³. A 1% solid volume fraction (`ALPHAS_SOURCE=0.01`, `RHO_S=2000 kg/m³`) gives a source mixture density of approximately 20.39 kg/m³, compared with ambient air at approximately 1.176 kg/m³. The resulting reduced gravity is positive, approximately 9.24 m/s². The earlier negative-gravity calculation used `XS_SOURCE=0.01`, incorrectly treating the intended volume fraction as a mass fraction.

The solver distinguishes the two composition inputs: `ALPHAS_SOURCE` is the solid volume fraction, while `XS_SOURCE` is the solid mass fraction. The active example now retains 900 K and uses the intended 1% volume fraction. A run with the documented parameters reached 100 s, wrote a nonempty NetCDF file and reported final mixture and solid masses of about `2.287e8 kg` and `2.025e8 kg`, respectively. Air entrainment is enabled, so the final mixture mass can exceed the imposed `MFR_SOURCE * T_END`. The solver also preserves the prescribed boundary-source depth when computing wave speeds, advective fluxes and the hydrostatic path.

The former 900 K invocation with `XS_SOURCE=0.01` is a different, mass-basis mixture and is correctly rejected when its reduced gravity is non-positive and a Richardson number is prescribed.

### Original explanatory note

The legacy README says the example uses a simple friction model "as in [1]", but no corresponding bibliography entry is supplied. Do not infer a specific reference until it has been checked.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record, including its historical description of the concentration as a mass fraction. The active guide corrects the composition to 1% by volume and records the distinction between `ALPHAS_SOURCE` and `XS_SOURCE`; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.

## Documentation verification note

**Parameter interpretation.** The generator maps its third argument to `MFR_SOURCE`; the solver reports this rate in kg/s and derives `H_SOURCE` from it. The final argument maps to `ALPHAS_SOURCE`, a solid volume fraction. At 900 K, 1% solid volume yields a mixture density of approximately 20.39 kg/m³ and positive reduced gravity of approximately 9.24 m/s². By contrast, `XS_SOURCE=0.01` means 1% solids by mass and gives negative reduced gravity at 900 K; that is not the intended example composition.
