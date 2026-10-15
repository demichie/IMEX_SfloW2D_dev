# Output files and formats

The solver writes generated files below a directory dedicated to each run:

```text
OUTPUT/<RUN_NAME>/
├── diagnostics/    maximum-thickness raster and runout CSV files
├── fields/         time-series NetCDF file
├── metadata/       effective input backup and snapshot manifest
├── probes/         probe CSV files
├── restart/        binary checkpoint, restart.bin
├── snapshots/      conservative and raster snapshots
├── terrain/        exported terrain rasters
└── vulnerability/ vulnerability tables and rasters
```

`OUTPUT_BASE_DIRECTORY` in `RUN_PARAMETERS` changes the parent directory;
the default is `OUTPUT`. It may be relative or absolute, but its path cannot
contain a `..` component or symlink. `RUN_NAME` may contain ASCII letters,
digits, underscores and hyphens.
Existing run directories are checked for symlinks, hard-linked files and
non-regular entries before they may be reused.

`OUTPUT_OVERWRITE` defaults to `.FALSE.`. A run fails if its output directory
already exists unless `OUTPUT_OVERWRITE=T` is explicitly selected. Reusing a
directory does not delete stale files from earlier runs, so use a new
`RUN_NAME` when a clean output set is required.

Each solver snapshot is written below `snapshots/` in a folder named with its
zero-based index and a readable encoding of the scheduled time. The exact
snapshot index and time in seconds are recorded in
`metadata/output_manifest.csv`. The full time series remains a single file,
`fields/<RUN_NAME>.nc`.

The effective input, including initialized defaults, is saved as
`metadata/<RUN_NAME>.bak`. Runout and mass-centre histories are written to
`diagnostics/`; probe histories are kept separately in `probes/`.

The solver writes binary checkpoints to the fixed path
`restart/restart.bin`. When continuing a run in place and
`RESTART_FILES="restart.bin"` is selected, an existing checkpoint in that
run's output directory is resolved automatically. An explicitly configured
checkpoint elsewhere remains an input and is not moved. In-place binary
continuation reuses the existing run directory without requiring
`OUTPUT_OVERWRITE=T`; other existing run directories require explicit
overwrite permission.

Input rasters and externally prepared files remain at their configured paths.
The solver does not relocate topography, erodibility, or externally supplied
restart files.

Post-processing scripts are separate from solver output management. The
`UTILS/plot_overlay.py`, `plot_runout.py` and `plot_probes.py` utilities save
figures beside the input file they process; pass files from the appropriate
run subdirectory to keep those figures with the run. `flip_texture.py` writes
its derived images in the current working directory.
