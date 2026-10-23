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

Generated Autotools files such as `configure` and `Makefile.in` are not tracked
in Git. Regenerate them from a checkout with:

> autoreconf -fi

The build uses the Fortran compiler selected by `FC` and the NetCDF Fortran
library. For a Conda environment, pass its prefix as `--with-netcdf`:

> ./configure --with-netcdf=/path/to/fortran_env

The historical optimized profile is the default. To enable runtime checks and
floating-point traps:

> ./configure --with-netcdf=/path/to/fortran_env --enable-strict-debug

`FCFLAGS` supplied by the environment is respected. Non-GNU toolchains can
override `OPENMP_FCFLAGS`; `NETCDF_LIBS` and `LAPACK_LIBS` can likewise be
provided when their library names or link dependencies differ from the
defaults.

> make

> make check

> make install

The executable is copied in the bin folder.

`make check` runs the legacy example integration cases in temporary
directories, leaving the source fixtures untouched. The CI workflow also runs
binary/stochastic restart and HP-PCCU Gate-H checks with actual one- and
four-thread OpenMP teams in both build profiles.

The legacy test runner can optionally run the focused `TEST_DYNAMIC_FINH`
fixture with `TESTS/run_tests.sh /path/to/IMEX_SfloW2D TEST_DYNAMIC_FINH`.
The Intel ifx portability workflow exposes the same focused diagnostic through
the `diagnose_dynamic_fin_h` manual-dispatch input and uploads its solver log
and elapsed-time record. The focused diagnostic run is limited to five minutes
so a stalled test still produces diagnostic artifacts promptly.

Several examples can be found in the EXAMPLES folder.

The [Etna lava example](EXAMPLES/EXAMPLE_ETNA_LAVA/README.md) uses two bottom
fissures and the original Etna DEM before avalanche preprocessing. Together
with `EXAMPLE_ETNA`, it is selected for the N9 long-case validation; it replaces
`EXAMPLE_LAVA_IMO` in that selection without removing the existing IMO example.

## Release downloads

GitHub automatically provides **Source code (zip)** and **Source code
(tar.gz)** snapshots for tags. These are Git snapshots, not the prepared
Autotools source distribution; they may not contain generated files such as
`configure`. When available, use the explicitly named
`IMEX_SfloW2D-<version>.tar.gz` and
`IMEX_SfloW2D-<version>-examples.tar.gz` assets attached to a validated
[release](https://github.com/demichie/IMEX_SfloW2D_dev/releases).

## Per-run output directories

Solver-generated files are grouped below `OUTPUT/<RUN_NAME>/` by default.
`OUTPUT_BASE_DIRECTORY` in `RUN_PARAMETERS` selects a different parent
directory, and `OUTPUT_OVERWRITE=T` explicitly allows reusing an existing run
directory. Reuse does not remove files left by an earlier run. Without the
overwrite flag, an existing run directory is rejected before solver outputs
are opened.

The run directory contains `metadata/` (effective input backup and snapshot
manifest), `fields/` (NetCDF), `snapshots/` (conservative and raster snapshots,
each in an index/time-stamped directory), `diagnostics/`, `probes/`,
`vulnerability/`, `terrain/`, and `restart/` (the written `restart.bin`).
The manifest records each snapshot index and its exact time. Input files such
as topography and externally supplied restart checkpoints remain at their
configured input paths. For an in-place binary continuation, the basename
`restart.bin` resolves to `OUTPUT/<RUN_NAME>/restart/restart.bin` when that
checkpoint exists; otherwise the configured input path is used.

`RUN_NAME` accepts ASCII letters, digits, `_` and `-`. Output path components
are created without following symlinks; parent-directory (`..`) components in
the configured output base are rejected. Before reusing an existing run
directory, its contents are checked for symlinks, hard-linked files and
non-regular entries.

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
for IMEX stages; the supported stage counts remain `N_RK=2,3,4`. The selected policy is written to the effective input backup
(`metadata/RUN_NAME.bak` within the per-run output directory).

## Runout and binary continuation

The binary checkpoint format is now **v5**. It preserves the original runout
origin and initial extents, previous barycentre, diagnostic sampling history,
maximum runout and steady-stop state together with the physical state.
Earlier binary formats, including v4, are rejected: they cannot provide an
exact diagnostic continuation. ASCII initial conditions remain supported and
are not treated as a binary continuation merely because `RESTART=T`.

Runout is written only to `RUN_NAME_runout.csv`; the accompanying barycentre
and directional extents remain in `RUN_NAME_mass_center.csv`. Their existing
five-column schemas and formatting are unchanged. `runout.txt` is neither read
nor written. Existing legacy files are not automatically deleted.

For a binary continuation, both existing CSVs are checked before either is
opened for append. Headers, finite ordered records and the exact last row must
match the restored diagnostic history. The last sample need not occur at the
checkpoint time. No header or checkpoint-time row is duplicated. A missing
CSV starts a new segment with a header; its earlier rows cannot be recovered
from a checkpoint.

Incompatible, incomplete or later CSV histories cause an error without
truncation or overwriting. Use matching histories or a new execution directory
to start a segment. Keep the runout origin configuration, `DT_RUNOUT`,
`ALPHAS_THRESHOLD`, `EPS_STOP` and `OUTPUT_RUNOUT_FLAG` unchanged for an exact
continuation. Changes to the simulation end time are allowed; an already
requested steady stop retains its original deadline.

This policy applies to the runout and mass-centre CSVs. Probe CSV handling is
unchanged by this package. See the
[restart verification](docs/acceptance/N9_RESTART_V5_2026_10_11.md) for the
verified cases and remaining acceptance scope.

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
