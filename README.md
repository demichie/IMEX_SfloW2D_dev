# IMEX-SfloW2D

IMEX-SfloW2D is a Fortran finite-volume, depth-averaged model for geophysical
and volcanic flows, using high-resolution spatial discretization and
implicit-explicit (IMEX) time integration. Supported single-layer carrier
configurations include gas, liquid and lava, with model-specific physics.

The solver supports both `N_LAYERS=1` and `N_LAYERS=2` directly
through the input file, without an additional activation variable.
The two-layer formulation has documented scientific limitations,
particularly for density inversion and interlayer shear; see
[two-layer hyperbolicity](docs/TWO_LAYER_HYPERBOLICITY.md) and the
[two-layer architecture](docs/TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md).
The existing single-layer numerical path remains unchanged.

## Build

Requirements: a Fortran compiler, Autoconf/Automake, NetCDF-Fortran,
LAPACK/BLAS and Make. On Ubuntu, equivalent packages include `gfortran`,
`autoconf`, `automake`, `libnetcdff-dev`, `liblapack-dev` and
`libblas-dev`. Use libraries built for your selected Fortran compiler.

```sh
autoreconf -fi
./configure --with-netcdf=/usr FC=gfortran
make -j2
make check
```

For a strict/debug build, configure with `--enable-strict-debug`.
The default build is optimized. Configure-time overrides include `FCFLAGS`,
`OPENMP_FCFLAGS`, `NETCDF_LIBS` and `LAPACK_LIBS`. For Conda or another
custom NetCDF-Fortran prefix, set `--with-netcdf=/path/to/prefix`.
Generated Autotools files are not tracked in Git.

`make install` installs the compiled executable. Focused numerical, physics,
restart and OpenMP checks are available under `TESTS/`; see
[contributing](CONTRIBUTING.md) and the [manual CI workflows](.github/README.md).
GitHub Actions in this repository run only on manual dispatch.

## Documentation and examples

- [Developer documentation](docs/README.md)
- [Documentation website](https://demichie.github.io/IMEX_SfloW2D_dev/)
- [Example configurations](EXAMPLES/) (each runnable example has its own README)
- [Model architecture and implementation status](docs/TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md)
- [NetCDF fields and binary restart](docs/NETCDF_OUTPUT.md)
- [Citation metadata](CITATION.cff) and [GPLv3 license](COPYING)

## Releases

GitHub automatically supplies source snapshots for tags; these are not
Autotools-generated distribution archives and may omit `configure`. When
available, download the separately prepared source and examples archives from
the [releases page](https://github.com/demichie/IMEX_SfloW2D_dev/releases).

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

Binary restart is a **single development-only checkpoint schema**. It stores
the canonical conservative and auxiliary state together with the runtime and
runout sampling history required for continuation. Checkpoints written by
superseded development implementations are not accepted; regenerate them
instead of maintaining a versioned reader.

Runout diagnostics are written to `RUN_NAME_runout.csv` and
`RUN_NAME_mass_center.csv`. On continuation, their existing headers and
sampling history must be consistent with the restored checkpoint; incompatible
histories are rejected without truncation. A missing CSV begins a new output
segment and cannot reconstruct earlier diagnostic rows.

Keep the runout origin and relevant diagnostic scheduling parameters unchanged
for a continuous restart. The dated
[restart acceptance report](docs/acceptance/N9_RESTART_V5_2026_10_11.md)
records historical validation evidence, **not** an additional supported
checkpoint format. Consult [the current output guide](docs/NETCDF_OUTPUT.md)
for physical-field and restart-file responsibilities.

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

The former velocity/concentration vertical-profile implementation has been
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
