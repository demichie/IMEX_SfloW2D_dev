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

Run `./create_example.py 400 0.5 300 false` to prepare the simulation, followed by `../../bin/IMEX_SfloW2D`. Arguments: x-cell count, particle volume fraction, temperature (K), plotting flag.

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` may be inspected with ParaView.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 24 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.
