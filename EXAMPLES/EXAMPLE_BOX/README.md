---
title: "One-Dimensional Dam Break with Mass Exchange"
summary: "One-dimensional dam break over flat terrain with deposition, entrainment, and friction."
category: "Numerical benchmarks"
tags:
  - "dam break"
  - "1D"
  - "entrainment"
---

# One-Dimensional Dam Break with Mass Exchange

## Overview

A 1-D dam-break problem on flat terrain, including deposition, entrainment, and basal friction.

## Main parameters

| Parameter | Value |
| --- | --- |
| Geometry | 1-D, flat bed |
| Particle volume fraction | Generator argument |
| Temperature | Generator argument (K) |
| x-resolution | Generator argument |

## Dependencies

Requires Python 3 and NumPy. Matplotlib is required only when the plotting flag is `true`.

## How to run

From this example directory, first generate the inputs, then run the solver.

### 1. Generate the input files

The generator takes **four positional arguments**, in this order:

1. Number of computational cells along x.
2. Initial solid-particle volume fraction (dimensionless, 0–1).
3. Initial flow temperature (K).
4. Plot the initial solution (`true` or `false`).

The script prepares the input for the 1-D flat-bed dam-break with deposition, entrainment and friction.

**Generator command:**

```bash
./create_example.py 400 0.5 300 false
```

### 2. Run the solver

**Solver command:**

```bash
../../bin/IMEX_SfloW2D
```

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` may be inspected with ParaView.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 21 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.
