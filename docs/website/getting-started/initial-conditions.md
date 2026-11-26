# Initial conditions and restart

IMEX-SfloW2D reads its simulation settings from
`IMEX_SfloW2D.inp` in the current case directory. This page explains
how to supply initial flowing material and how to resume a simulation.

The examples provide complete, working input files; start from one
that matches the intended physics.

## Starting from a raster or conservative snapshot

The `RESTART` flag enables initialization from files. **An input raster
or a `.q_2d` snapshot may specify an initial state rather than continue
a previous simulation**.

For a one-layer case, use exactly one of
`RESTART_PARAMETERS` or `RESTART_L1_PARAMETERS` and the matching
`INITIAL_PARAMETERS` or `INITIAL_L1_PARAMETERS` group.
For example:

```fortran
&RESTART_PARAMETERS
 N_RESTART_FILES=2,
 RESTART_FILES='first.asc','second.asc',
 RELEASE_TIME=0.D0,10.D0,
/
```

This configuration supplies two layer-thickness rasters and a
release schedule. Required files must exist at their specified
paths relative to the working directory.

## Two-layer initial states

A two-layer case requires separate L1 and L2 groups, even if one
layer starts dry:

```fortran
&RESTART_L1_PARAMETERS
 N_RESTART_FILES=2,
 RESTART_FILES='first_l1.asc','second_l1.asc',
 RELEASE_TIME=0.D0,0.D0,
/
&RESTART_L2_PARAMETERS
 N_RESTART_FILES_L2=1,
 RESTART_FILES_L2='NONE',
 RELEASE_TIME_L2=0.D0,
/
```

Here the upper layer is initially dry. Material properties,
composition, temperature and velocity are configured using
`INITIAL_L1_PARAMETERS` and `INITIAL_L2_PARAMETERS`.

For an initial state supplied by `.q_2d` conservative snapshots,
provide exactly one complete snapshot for each layer, with the
same output index. Do not mix `.q_2d` and `.asc` inputs in
one primitive initialization.

!!! warning "Two-layer delayed release"
    Nonzero L1 or L2 release times are not supported for two-layer
    primitive initialization. Set both release schedules to zero.

## Continuing from a binary checkpoint

The solver writes a complete simulation checkpoint to
`OUTPUT/<RUN_NAME>/restart/restart.bin`.

A binary `.bin` checkpoint contains the complete state for
**all configured layers**, together with the information required
to continue the run. Supply it through the L1/unsuffixed restart
group, without adding `RESTART_L2_PARAMETERS`.

When resuming in the same output directory,
`RESTART_FILES='restart.bin'` selects that run's existing
checkpoint. A checkpoint at an explicitly specified other path
remains an external input.

A binary checkpoint is different from a visualization NetCDF
file. **NetCDF output cannot be used as a complete solver restart.**

## Practical checks

- Check the case's `N_LAYERS` and initial-state group names.
- Prepare all referenced DEM and initial-state files before running.
- For two layers, check that both initial fields use the same grid.
- Keep the full input, version and outputs of the original run when
  restarting or comparing numerical results.

For worked examples, see [Getting started](index.md) and
[Two-layer formulation](../model/two-layer.md).
