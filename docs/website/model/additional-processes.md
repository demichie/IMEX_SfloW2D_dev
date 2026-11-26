# Pore pressure, stochastic forcing and material injection

In addition to rheology, environmental exchange and heat transfer,
IMEX-SfloW2D provides optional physical processes used in
particular simulation configurations.

## Pore-pressure evolution

A Darcy-based formulation describes an optional pore-pressure
component and associated carrier-gas evolution. Where the selected
rheology supports the coupling, pore pressure reduces the effective
normal stress and changes basal resistance.

The Darcy source includes a supported analytic update over a finite
time increment. Its exact permeability and gas/thermal closures
depend on the configured material parameters.

## Stochastic variation

The model can evolve an Ornstein–Uhlenbeck-type stochastic variable
and use a transformed fluctuation in selected rheological properties.
The transported stochastic variable and the effective value applied
to friction are distinct. Stochastic forcing must be configured and
interpreted as a model assumption, not as an automatic representation
of all physical variability.

## External material injections

Material can enter through configured circular or fissural bottom
sources, radial sources, or boundary injections. Source descriptions
include injection geometry, activation time and properties of the
injected mixture.

Bottom and radial injections supply cell-source material; boundary
injections impose boundary traces and should not be interpreted as
the same cell-volume operation.

## Choosing compatible processes

Different source laws support different explicit, semi-implicit,
implicit or analytic treatments. The solver validates these
choices rather than silently substituting another method.
Source groups are applied in the configured order, as described
under [numerical methods](numerics.md).

The [example catalogue](../examples/index.md) illustrates
combinations relevant to particular flow problems.
