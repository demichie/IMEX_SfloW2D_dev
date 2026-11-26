# Getting started

Choose an installation method and then follow the README of a specific example.

1. [Build the solver from source](installation.md) using gfortran, NetCDF-Fortran and LAPACK/BLAS, or use the [Docker container](docker.md) when a published image is available.
2. Open the [examples catalogue](../examples/index.md). Each case documents its input prerequisites; some require copying and extracting a shared DEM.
3. Prepare input rasters or initial conditions with the selected example's script. Such scripts do **not** run the solver unless explicitly stated.
4. Run the executable from the case directory, where `IMEX_SfloW2D.inp` is located.
5. Consult [output files](output-files.md) and [visualization](visualization.md) to inspect NetCDF, diagnostics and restart checkpoints.

For two-layer cases, consult the [coupled-model limits](../model/two-layer.md) before interpreting results.
