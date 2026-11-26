# Getting started

You can run IMEX-SfloW2D by compiling the Fortran executable or using
a prebuilt Docker image. Each example includes its own simulation
settings, terrain prerequisites and preparation commands.

## Run your first simulation

1. [Build from source](installation.md) or [use Docker](docker.md).
2. Choose a case from the [examples catalogue](../examples/index.md).
3. Prepare the required DEM and initial-state files as described by
   that example. Preparation scripts generate inputs; they do not
   run the solver unless their instructions explicitly say so.
4. Open the example directory containing `IMEX_SfloW2D.inp` and
   run the solver there.
5. Inspect the [output directories](output-files.md) and use
   [ParaView or QGIS](visualization.md) or the
   [Python plotting tools](postprocessing.md) to examine the results.

## Before you run

- Check that the case's terrain and input files are available.
- Verify that `N_LAYERS` matches your intended single- or
  two-layer calculation.
- Read the [initial conditions and restart guide](initial-conditions.md)
  if changing the initial state or continuing a run.
- Check the [two-layer applicability limits](../model/hyperbolicity.md)
  when using interacting layers.

The examples are provided to demonstrate how to configure and use
the model. Quantitative predictions require appropriate physical
parameters, validation and numerical convergence checks.
