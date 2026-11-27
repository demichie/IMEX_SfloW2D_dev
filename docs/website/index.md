# IMEX-SfloW2D

!!! danger "Development version — not for production use"
    **IMEX-SfloW2D is under active development and is not intended for production use.** Its numerical methods, physical models, interfaces, and results may change. Do not rely on this development version for operational hazard assessment, safety-critical decisions, or other production applications. Validate the software and its results independently before research use.

**Finite-volume modelling of volcanic and other geophysical flows.**

IMEX-SfloW2D is an open-source Fortran model for two-dimensional, depth-averaged
flows over complex terrain. It can simulate a single flowing mixture or two
interacting layers, with configurable material properties, basal resistance,
transport, deposition and other physical processes.

[Get started](getting-started/index.md){ .md-button .md-button--primary }
[Explore examples](examples/index.md){ .md-button }
[Understand the model](model/index.md){ .md-button }

## What can you simulate?

- **Single-layer flows:** gravity-driven granular avalanches, sediment-rich
  lahars and other configured mixtures, including lava applications.
- **Two-layer flows:** a lower dense mixture interacting with an upper flow
  through pressure coupling, particle transfer and interlayer drag.
- **Terrain-dependent dynamics:** flow over a supplied DEM, with optional
  erosion, deposition and terrain evolution according to the selected laws.
- **Transport and thermal behaviour:** mixture composition, optional pore
  pressure, thermal exchange and other configured physical sources.

The applicable equations and source laws depend on the case configuration.
An example is a starting point, not evidence that its coefficients are
calibrated for a particular event.

## Using the documentation

| Your goal | Start with |
| --- | --- |
| Install the solver and run your first case | [Getting started](getting-started/index.md) |
| Choose a simulation to adapt | [Example catalogue](examples/index.md) |
| Understand the equations and numerical treatment | [Physical and numerical model](model/index.md) |
| Compare rheological and environmental processes | [Rheological models](model/rheology.md) and [erosion and deposition](model/mass-exchange.md) |
| Configure two interacting layers | [Two-layer model](model/two-layer.md) |
| Read and visualize NetCDF results | [Output files](getting-started/output-files.md) and [visualization](getting-started/visualization.md) |
| Cite the software | [References and citation](references/index.md) |

!!! warning "Interpretation of two-layer results"
    The two-layer pressure system is not hyperbolic for every combination of
    layer density and relative velocity. A completed simulation does not
    establish that its local characteristic speeds are real or that its
    physical closures are calibrated. Read the
    [characteristic limits](model/hyperbolicity.md) before interpreting
    results from strongly sheared or density-inverted flows.

For the implementation itself, see the
[Fortran API reference](https://demichie.github.io/IMEX_SfloW2D_dev/doxygen/).
Source code, license and development guidance are available in the
[GitHub repository](https://github.com/demichie/IMEX_SfloW2D_dev).
