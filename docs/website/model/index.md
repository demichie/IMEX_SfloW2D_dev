# Physical and numerical model

IMEX-SfloW2D solves depth-averaged conservation equations for mass, momentum and, where configured, additional transported mixture variables and thermal effects. It uses finite-volume spatial discretization together with implicit-explicit and analytic source treatments.

## Architecture

| Component | Responsibility |
| --- | --- |
| `src/state/` | Canonical conservative and primitive states, conversions and admissibility |
| `src/spatial/` | Reconstruction, numerical fluxes, topographic and coupled pressure terms |
| `src/physics/` | Cell-local rheology, thermal, pore-pressure, injection and mass-exchange laws |
| `src/numerics/` | Time integration, Lie source ordering, implicit solves and accepted-state publication |
| `src/core/` | Model configuration, layer descriptors, lifecycle and solver diagnostics |

A two-layer calculation shares terrain and geometry but preserves independently configured layer states and source laws. It does not run two separate single-layer simulations.

- [Spatial and temporal integration](numerics.md)
- [Two-layer coupling and scientific limits](two-layer.md)
- [Current architecture specification](https://github.com/demichie/IMEX_SfloW2D_dev/blob/master/docs/TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md)
- [Source-code reference](https://demichie.github.io/IMEX_SfloW2D_dev/doxygen/)

For physically specific examples, use the [case catalogue](../examples/index.md). Dated acceptance reports are historical evidence tied to their recorded revisions, not unconditional validation of every flow regime.
