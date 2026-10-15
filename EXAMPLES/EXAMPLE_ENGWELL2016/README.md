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
| Solid mass fraction (`xs_init`) | Sixth generator argument; verify transformation in solver setup |

## How to run

Run `python create_example.py 200 2000 2000 900 0.9 0.8` and then `../../bin/IMEX_SfloW2D`. The script expects six arguments: x cells, source radius, thickness, temperature, Richardson number and solids fraction. The legacy README inconsistently describes the last quantity as both volume and mass fraction; confirm against code before making a definitive claim.

## Expected output

Native NetCDF output at `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` can be viewed in ParaView.

## References

1. Bursik, M. I., and Woods, A. W. (1996). The dynamics and thermodynamics of large ash flows. *Bulletin of Volcanology*, 58, 175–193.
2. Engwell, S. L., de' Michieli Vitturi, M., Esposti Ongaro, T., and Neri, A. (2016). Insights into the formation and dynamics of coignimbrite plumes from one-dimensional models. *Journal of Geophysical Research: Solid Earth*, 121(6), 4211–4231.

## Legacy documentation and preservation

The original [`README.txt`](README.txt) in this example directory remains the archival source for the complete historical narrative, parameter discussion, scientific references and commands. It has not been deleted or rewritten. This standardized README must not replace it as the sole provenance source until a full section-by-section textual comparison has been accepted.

## Documentation verification note

**Review solid fraction definition.** The last generator argument is named `xs_init` in the script and is therefore documented here as a source solid **mass fraction**, rather than an unspecified volume fraction. The README originally uses both terms; the definition must be checked against the source-to-conservative-variable conversion before claiming the formula.
