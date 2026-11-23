# Vesuvius LAHAR — 600 s single-layer example

This example runs the supplied Vesuvius LAHAR configuration with the
production single-layer solver. It uses a reduced, grid-aligned computational
domain to avoid processing terrain far from the flow. Run all commands from
`EXAMPLES/EXAMPLE_LAHAR`.

## Input and computational domain

The supplied `IMEX_SfloW2D.inp` is retained, including its physical and
numerical parameters. Only the obsolete introductory migration-status comment
was updated.

| Setting | Value |
| --- | --- |
| Run name | `VESUVIO_001_001` |
| Model | Single layer (`N_LAYERS=1`), six solid species |
| Lower-left corner | `X0=440995 m`, `Y0=4515991 m` |
| Grid | `340 x 280` cells, `50 m` resolution |
| Domain extent | `[440995, 457995] x [4515991, 4529991] m` |
| Time | `T_START=0 s`, `T_END=600 s` |
| Output interval | `600 s` |
| Time stepping | `CFL=0.24`, `MAX_DT=5 s` |
| Rheology | `lahar` |

The input enables settling, erosion, and liquid-carrier loss, and reads
the initial thickness through the existing ASCII restart/initial-state path.
The terrain modification flag is disabled. Preserve the original input
parameters when using this example as a performance or numerical reference.

## Prepare the DEM rasters

The repository stores each large ASCII raster as an individual ZIP archive
in `DEM/`. The three archives must be added separately; they are **not**
included in the initial example commit.

After the ZIP files have been uploaded, extract them:

```bash
unzip -o DEM/ves_asciiNEW.zip -d DEM
unzip -o DEM/VESUVIO_001_001_erodible.zip -d DEM
unzip -o DEM/masked_init_vol_thick_001_001.zip -d DEM
```

The input expects these exact paths:

- `DEM/ves_asciiNEW.asc` — terrain elevation;
- `DEM/VESUVIO_001_001_erodible.asc` — erodible-bed thickness;
- `DEM/masked_init_vol_thick_001_001.asc` — initial LAHAR thickness.

See [DEM/README.md](DEM/README.md) for archive contents, naming, and
verification. The shared `EXAMPLES/.gitignore` excludes extracted
`DEM/*.asc` rasters, so only the three ZIP files should be committed.

## Run

With the three rasters extracted and the solver already built:

```bash
../../src/IMEX_SfloW2D
```

Output is written below `OUTPUT/VESUVIO_001_001/`, including NetCDF
fields and the final binary restart. No preprocessing script or preliminary
single-layer run is required for the supplied ASCII initial thickness.

For an explicit OpenMP thread count, for example:

```bash
OMP_NUM_THREADS=4 ../../src/IMEX_SfloW2D
```

## Local reference and validation status

The supplied, grid-aligned input was run locally through `t=600 s`.
The user-reported log contains **898 timesteps**, **6.250 s** of wall time
for iterations, and **7.217 s** total elapsed wall time. These times are
machine-dependent observations, **not** regression thresholds or evidence
of a particular OpenMP speedup.

The flow was reported to remain within the reduced domain. A cellwise
comparison of conservative fields against the larger-domain reference and
a repository CI run with all three raster archives are **not yet recorded**.
Do not treat the example as an automated scientific regression until those
checks and reference data have been added.
