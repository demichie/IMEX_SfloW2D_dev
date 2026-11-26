# Lava, lahar and other geophysical-flow applications

IMEX-SfloW2D is a configurable flow solver, not a single
material-specific equation. Applications differ through their
rheology, composition, source processes, initial conditions and
boundary forcing.

## Lava flows

Lava configurations may combine temperature-dependent viscous
and yield-stress resistance with thermal cooling and lava
deposition. These processes play different roles:

- **Rheology** determines how resistance changes with temperature,
  depth and motion.
- **Cooling** changes the thermal content of moving material.
- **Lava deposition** removes mobile lava and can increase bed elevation.

A lava simulation does not automatically activate every available
thermal or deposition law. Verify the selected processes in its input.

## Lahars and sediment-rich flows

Lahar configurations may combine concentration- and
temperature-sensitive resistance with solids transport,
erosion, deposition and carrier exchange. The moving mixture
composition can evolve, so fixed rheological parameters alone
do not describe every possible run.

For quantitative interpretation, distinguish the assumed rheology
from sediment availability, settling, erosion and flow-boundary
conditions.

## Granular avalanches and other mixtures

For granular flows, alternatives include Voellmy, Froude-dependent,
Lucas or granular µ(I)-type resistance, subject to each law's
assumptions. Additional pore-pressure or stochastic processes
may be useful in specifically configured cases.

## A model configuration, not just an application label

The terms *lava*, *lahar* and *granular avalanche* describe
applications; they do not uniquely define a complete parameter
set. Identify the actual selected [rheology](rheology.md),
[mass-exchange processes](mass-exchange.md),
[thermal processes](thermal.md) and
[numerical treatment](numerics.md).

Browse the [examples](../examples/index.md) to see those
choices assembled into runnable cases.
