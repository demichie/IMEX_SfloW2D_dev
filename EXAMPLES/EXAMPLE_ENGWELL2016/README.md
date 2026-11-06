---
title: "Radial Flow Benchmark Inspired by Engwell et al. (2016)"
summary: "Radial supercritical source with particle sedimentation and air entrainment."
category: "Volcanic flows"
tags:
  - "radial source"
  - "coignimbrite"
  - "supercritical"
---

# Radial Flow Benchmark Inspired by Engwell et al. (2016)

## Overview

A supercritical radial inlet over flat terrain, intended to reproduce the behavior of Figure 3 in Engwell et al. (2016). Sedimentation and atmospheric-air entrainment are enabled; a simple friction law is used. The source must be supercritical (Ri < 1) to avoid characteristics entering the prescribed boundary from the computational domain.

## Main parameters

| Parameter | Value |
| --- | --- |
| Source radius | Generator argument |
| Initial thickness | Generator argument |
| Initial temperature | Generator argument (K) |
| Richardson number | Generator argument, < 1 |
| Solid mass fraction (`XS_SOURCE`) | Sixth generator argument (`xs_init`) |

## Dependencies

The input generator requires Python 3 and NumPy.

## How to run

From this example directory, first generate the inputs and then run the solver.

### 1. Generate the input files

The generator expects **6 positional arguments**, in this exact order:

1. Number of cells along x (positive integer).
2. Source radius (m; positive).
3. Source thickness (m; positive).
4. Initial source temperature (K; positive).
5. Initial Richardson number (dimensionless; below 1 for the supercritical inlet).
6. Initial solid **mass fraction** (dimensionless, between 0 and 1).

The sixth value is written to `XS_SOURCE`, a **mass fraction**; the older README also calls it a volume fraction, but that statement conflicts with the current generator. The supercritical regime (`Ri<1`) prevents characteristics from propagating upstream across the prescribed source.

**Generator command:**

```bash
python create_example.py 200 2000 2000 900 0.9 0.8
```

### 2. Run the solver

**Solver command:**

```bash
../../bin/IMEX_SfloW2D
```

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` can be viewed in ParaView.

## References

1. Bursik, M. I., and Woods, A. W. (1996). The dynamics and thermodynamics of large ash flows. *Bulletin of Volcanology*, 58, 175–193.
2. Engwell, S. L., de' Michieli Vitturi, M., Esposti Ongaro, T., and Neri, A. (2016). Insights into the formation and dynamics of coignimbrite plumes from one-dimensional models. *Journal of Geophysical Research: Solid Earth*, 121(6), 4211–4231.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 21 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.
