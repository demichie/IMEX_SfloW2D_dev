# Two-layer formulation

Set `N_LAYERS=2` in `IMEX_SfloW2D.inp` to use the integrated coupled solver. No environment activation flag is required. `N_LAYERS=1` remains the standard single-layer configuration.

## Shared and per-layer quantities

- **Shared:** terrain elevation, grid and geometry, active-cell domain, material/species catalogue and substrate inventories.
- **Per layer:** conservative and primitive states, rheology, thermal/pore-pressure configurations, and supported external source laws.
- **Interlayer:** mass-transfer and momentum-drag laws use the single pair-wide `INTERLAYER_PARAMETERS` configuration.

The four local wetness cases (both wet, L1 only, L2 only, both dry) determine which layer touches terrain or atmosphere. Interlayer drag requires both layers present.

## Physical coupling

With terrain `B` and layer depths `h1` and `h2`, the two free surfaces are `eta1=B+h1` and `eta2=B+h1+h2`. Coupled pressure enters the spatial operator, not a separate Lie source. Interlayer mass transfer and drag are handled by coupled source transactions that enforce donor availability and update both layers consistently.

## Important applicability limits

!!! warning "Coupled characteristics are not always real"
    Inverted density ordering (`rho2>rho1`) can produce complex gravitational eigenvalues even at zero relative velocity. Sufficiently strong shear can produce complex eigenvalues even for stable density ordering. Equal-density wet/wet states may have repeated real eigenvalues without a complete eigenbasis.

The current total-depth CFL speed envelope is not a general hyperbolicity test. Reducing the timestep cannot make a non-hyperbolic state hyperbolic. The documented spectrum is an analysis of a frozen smooth gravitational subsystem, not a global certificate for the full discrete solver.

**Read before interpreting results:**

- [Scientific hyperbolicity analysis](https://github.com/demichie/IMEX_SfloW2D_dev/blob/master/docs/TWO_LAYER_HYPERBOLICITY.md)
- [Architecture and conservation contracts](https://github.com/demichie/IMEX_SfloW2D_dev/blob/master/docs/TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md)
- [Physical closure assumptions](https://github.com/demichie/IMEX_SfloW2D_dev/blob/master/docs/TL13B_C_PHYSICAL_CLOSURES.md)

## Runnable examples

- [Etna granular collapse and upper dilute current](../examples/EXAMPLE_ETNA_TWO_LAYER/index.md)
- [1D dam-break with initially dry upper layer](../examples/EXAMPLE_DAMBREAK_TWO_LAYER_1D/index.md)

These examples demonstrate supported workflows; they are not empirical validation for every density, shear or exchange regime.
