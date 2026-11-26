# Two-layer formulation

IMEX-SfloW2D can evolve two interacting, depth-averaged mixtures using
the same executable as a single-layer simulation.

Set `N_LAYERS=2` in `IMEX_SfloW2D.inp`. For a single layer, use
`N_LAYERS=1`. There is no separate two-layer executable or
activation environment variable.

## Layer geometry

The lower layer is **L1**, with depth $h_1$; the upper layer is
**L2**, with depth $h_2$. Both share the terrain elevation $B$.

The interface between layers lies at $B+h_1$. The upper free
surface is

$$
\eta_{\mathrm{top}}=B+h_1+h_2.
$$

The model uses one computational grid and shared terrain/substrate
inventories. Each layer has its own conservative and primitive
variables, composition and configurable physical properties.

## Layer configurations

| Aspect | L1 (lower) | L2 (upper) |
| --- | --- | --- |
| Flow state | Depth, velocity, mixture and optional variables | Separate depth, velocity, mixture and optional variables |
| Rheology and cooling | Layer-specific parameters | Layer-specific parameters |
| Initial condition | Own initial raster or complete state | Own initial raster or complete state |
| Coupling | Receives pressure and exchanges from L2 | Receives pressure and exchanges from L1 |
| Substrate deposition | Supported for L1 | L2 settling transfers solids to L1, not directly to the substrate |

Layer contact with the terrain or atmosphere depends on whether
the other layer is present. For example, an isolated wet upper
layer can contact terrain where the configured law allows it.

## Input groups

For a two-layer run, supply the appropriate L1 and L2 groups,
such as `INITIAL_L1_PARAMETERS` and `INITIAL_L2_PARAMETERS`.
The shared `INTERLAYER_PARAMETERS` group configures the
pair-wide transfer and drag laws.

For `N_LAYERS=1`, the unsuffixed forms such as
`INITIAL_PARAMETERS` and `RHEOLOGY_PARAMETERS` are intentional,
supported names.

Two-layer primitive initialization currently requires the
layer release times to be zero; independent delayed layer
releases are not supported. See
[initial conditions and restart](../getting-started/initial-conditions.md).

## What couples the layers?

**Pressure coupling** is part of the finite-volume spatial equations:
both layers affect horizontal pressure forces. The exact terms and
terrain metric are described on the
[coupled pressure page](coupling.md).

**Mass exchange** may transfer solids and associated carrier
between the layers through settling and elutriation. **Interlayer drag**
relaxes relative motion while conserving combined momentum at fixed
layer masses. Their assumptions are explained under
[interlayer processes](interlayer-processes.md).

## Scientific applicability

!!! warning "Density and shear restrictions"
    The coupled gravitational system does not always have real
    characteristic speeds. Inverted density ordering or sufficiently
    strong interlayer shear can produce complex eigenvalues; even equal
    density can lead to a defective repeated characteristic.
    Successful execution or satisfaction of a CFL condition does not
    establish general physical validity.

Read the [characteristic analysis](hyperbolicity.md), particularly
when layer densities are close, inverted, or relative velocities
are large. Interlayer transfer laws also require scenario-specific
calibration.

## Try a simulation

- [Two-layer 1D dam-break](../examples/EXAMPLE_DAMBREAK_TWO_LAYER_1D/index.md):
  begin with a lower layer and an initially dry upper layer.
- [Two-layer Etna collapse](../examples/EXAMPLE_ETNA_TWO_LAYER/index.md):
  explore lower-to-upper particle transfer and a dilute upper current.

These are model demonstrations, not universal calibration benchmarks.
