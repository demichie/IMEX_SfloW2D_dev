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

Run `python create_example.py 200 2000 2000 900 0.9 0.8` and then `../../bin/IMEX_SfloW2D`. The script expects six arguments: x cells, source radius, source thickness, temperature, Richardson number and solid mass fraction. The sixth argument is substituted into `XS_SOURCE`, which the solver reads as a mass fraction; the legacy description calling it a volume fraction is incorrect.

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` can be viewed in ParaView.

## References

1. Bursik, M. I., and Woods, A. W. (1996). The dynamics and thermodynamics of large ash flows. *Bulletin of Volcanology*, 58, 175–193.
2. Engwell, S. L., de' Michieli Vitturi, M., Esposti Ongaro, T., and Neri, A. (2016). Insights into the formation and dynamics of coignimbrite plumes from one-dimensional models. *Journal of Geophysical Research: Solid Earth*, 121(6), 4211–4231.

## Legacy documentation and preservation

The paired [`README.txt`](README.txt) is retained as a verbatim archival record. This guide was compared section by section with the legacy file on 21 October 2026; verified differences in commands, dependencies and references are recorded in the Examples second-pass audit.
