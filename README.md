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
