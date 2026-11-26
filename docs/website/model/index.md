# Physical and numerical model

IMEX-SfloW2D evolves depth-integrated mixtures over a two-dimensional
Cartesian grid and a supplied terrain surface. The model combines
finite-volume transport with constitutive laws and source terms selected
in the simulation input.

## Flow variables

The principal physical quantities are flow depth, horizontal velocity,
mixture density and composition. Depending on the configuration, the
simulation may also transport thermal content, multiple solid or gas
species, excess pore pressure and a stochastic variable.

The model advances conservative mass, momentum and thermal/constituent
quantities and converts them to physical variables such as depth,
velocity and temperature for constitutive laws and output. Dry cells
have no flowing material.

## Physical processes

The available source laws describe, where enabled:

- Resistance to motion through the configured rheology, including
  frictional or velocity-dependent formulations.
- Exchange with the substrate or atmosphere, including erosion,
  deposition, entrainment and carrier loss when supported by the
  selected flow model.
- Thermal cooling and other temperature-dependent effects.
- Optional pore-pressure and stochastic processes.
- Sources of new material, such as configured injections.

Not every law or time-integration treatment is valid for every
material. Choose a model from an existing example and verify that
its physical assumptions fit your application.

## Single layer and two layers

With `N_LAYERS=1`, one depth-averaged flow exchanges material and momentum
with its surroundings according to its configured laws.

With `N_LAYERS=2`, two flowing layers share a computational terrain and
are coupled by hydrostatic pressure and, when enabled, interlayer
mass and momentum exchange. Layer configurations may differ.
The two-layer equations have important restrictions related to density
ordering and relative velocity.

## Read more

- [Numerical methods](numerics.md) — reconstruction, fluxes, time stepping
  and timestep selection.
- [Two-layer formulation](two-layer.md) — layer variables and input.
- [Coupled pressure equations](coupling.md) — the implemented pressure
  coupling between layers.
- [Interlayer processes](interlayer-processes.md) — settling, elutriation,
  carrier transfer, drag and energy accounting.
- [Characteristic limits](hyperbolicity.md) — when the coupled
  gravitational wave speeds can become complex.

To connect the equations to a concrete configuration, use the
[example catalogue](../examples/index.md).
