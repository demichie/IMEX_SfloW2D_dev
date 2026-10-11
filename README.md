# Depth-averaged gas-particles model

Shallow water model for multiphase flow (gas+particles) with density of gas
temperature-dependent.

To compile the code you need a Fortran compiler and the NetCDF library for Fortran. You can install both with anaconda, by creating an anaconda environment and activating it:

> conda create -n fortran_env conda-forge::gfortran_linux-64  sysroot_linux-64  make  conda-forge::netcdf-fortran  conda-forge::liblapack  conda-forge::libblas
> 
> conda activate fortran_env

or, on a osx computer:

> conda create -n fortran_env conda-forge::gfortran make conda-forge::netcdf-fortran conda-forge::liblapack conda-forge::libblas
> 
> conda activate fortran_env


To compile:

> autoreconf -i

Then on linux (replace USERNAME with you user account):

> ./configure --with-netcdf=/home/USERNAME/anaconda3/envs/fortran_env

On OSX (replace USERNAME with you user account):

> ./configure --with-netcdf=/USERS/USERNAME/anaconda3/envs/fortran_env

To compile the code with OpenMP add the following flag in src/Makefile:
1) with gfortran: -fopenmp
2) with intel: -qopenmp

> make

> make install

The executable is copied in the bin folder.

Several examples can be found in the EXAMPLES folder.

The [Etna lava example](EXAMPLES/EXAMPLE_ETNA_LAVA/README.md) uses two bottom
fissures and the original Etna DEM before avalanche preprocessing. Together
with `EXAMPLE_ETNA`, it is selected for the N9 long-case validation; it replaces
`EXAMPLE_LAVA_IMO` in that selection without removing the existing IMO example.

## CFL timestep policies

`NUMERIC_PARAMETERS` accepts `CFL_POLICY='directional_min'` (the default,
including when omitted) or `CFL_POLICY='combined'`. Both use the same final
HP-PCCU face speeds. With nonnegative directional speed maxima `lambda_x/y`,
the spatial suggestion is respectively:

- `min(MAX_DT, CFL*dx/lambda_x, CFL*dy/lambda_y)`, ignoring zero-speed directions;
- `min(MAX_DT, CFL/(lambda_x/dx + lambda_y/dy))`, using only active directions.

If no direction propagates, the suggestion is `MAX_DT`. Both policies agree
on a one-dimensional x/y grid. The existing admissible CFL range remains
`0 <= CFL <= 0.5` in 1D and `0 <= CFL <= 0.25` in 2D. Precisely `CFL=-1`
is also accepted: it returns `MAX_DT` without evaluating reconstruction/speeds.
Other negative values and unknown policy names are rejected.

The accepted timestep still obeys output/runout/probe events, the final time
and the existing growth limiter. Neither policy is a general stability proof
for IMEX stages; the supported stage counts remain `N_RK=2,3,4`. The selected
policy is written to the effective input backup (`RUN_NAME.bak`).

## Runtime diagnostics

`VERBOSE_LEVEL` controls how much diagnostic information is printed and never
pauses an execution. Interactive debug pauses can be enabled independently in
the `RUN_PARAMETERS` namelist:

```fortran
INTERACTIVE_DEBUG_FLAG = T,
```

The flag defaults to `F` and should remain disabled for batch and production
runs. Interactive pauses requested inside an active OpenMP parallel region are
skipped to avoid unsafe concurrent access to standard input. Invalid input and
nonphysical solver states terminate with a nonzero error status after printing
their diagnostic context.

## Vertical closure status

The legacy velocity/concentration vertical-profile implementation has been
removed. The current solver uses the standard depth-averaged equations, and
the input options `VERTICAL_PROFILES_FLAG` and
`VERTICAL_PROFILES_PARAMETERS` are no longer accepted.

Vertical-structure effects will be reintroduced through a separate closure
interface that maps conservative state and local context to coefficients used
by the equation terms. This will allow analytical/quadrature and data-driven
backends to share the same solver-facing contract.

## Bottom fissural sources

Bottom fissural sources are enabled with `BOTTOM_FISSURAL_SOURCE_FLAG = T` and
`N_FISSURES` in `NEWRUN_PARAMETERS`. They cannot be combined with radial or
lateral source flags. Each fissure is a finite-width rectangular strip defined
by two centerline endpoints and a width. Its cell coverage is computed by exact
polygon/cell intersection and kept separate from the other fissures.

`FISSURAL_SOURCE_PARAMETERS` accepts endpoint coordinates, widths, a linear
injection velocity and temperature per fissure, plus four time parameters per
fissure: period, active duration, ramp duration and phase. If the time
parameters are omitted, the source is active at constant strength for the run.
All fissures use the same source composition, specified with the existing
`ALPHAS_SOURCE`/`XS_SOURCE`, `ALPHAG_SOURCE`/`XG_SOURCE`, and liquid fraction
fields in the fissural namelist. Volumetric-flow-rate input is not implemented.

```fortran
&NEWRUN_PARAMETERS
	BOTTOM_FISSURAL_SOURCE_FLAG = T,
	N_FISSURES = 1,
	...
/
&FISSURAL_SOURCE_PARAMETERS
	X_FISSURES_END_POINTS = 100.0, 300.0,
	Y_FISSURES_END_POINTS = 200.0, 200.0,
	WIDTH_FISSURES = 20.0,
	LINEAR_VEL_FISSURES = 1.0,
	T_FISSURES = 1200.0,
	TIME_PARAM_FISSURES = 60.0, 45.0, 5.0, 0.0,
	ALPHAS_SOURCE = 0.02,
/
```
