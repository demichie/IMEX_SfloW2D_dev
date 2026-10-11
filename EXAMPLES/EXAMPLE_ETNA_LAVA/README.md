# Etna lava example

This one-hour lava simulation uses two bottom fissures and temperature-dependent
rheology 3. It replaces `EXAMPLE_LAVA_IMO` as the lava case selected for N9,
alongside the existing `EXAMPLE_ETNA` avalanche case. Scientific-baseline
approval remains a separate N9 decision.

The supplied input is preserved byte for byte: 250 by 150 cells at 10 m,
`T_END=3600 s`, outputs every 60 s, liquid density 2600 kg/m3 and fissure
temperature 1538 K. Slope correction and curvature are OFF, as supplied.
`CFL_POLICY` is omitted and therefore selects `directional_min`; `N_RK=2`.
Parameter values and input/DEM checksums are recorded in [fixture.json](fixture.json).

## Prepare the original DEM

From the repository root:

```sh
cd EXAMPLES/EXAMPLE_ETNA_LAVA
python3 prepare.py
```

The preparation script extracts `Etna2014_crop.asc` from the existing
`../EXAMPLE_ETNA/DEM/Etna2014_crop.zip` into `./DEM/`. It verifies both archive
and raster SHA-256 checksums and accepts an existing DEM only if identical.
The source raster is 3000 by 2000 pixels at 1 m; the simulation grid lies
inside its extent. Only the Python standard library is required.

Do not run `create_input_ellipsoid.py` for this lava example. That avalanche
preprocessing excavates a niche and generates a different bed and initial
volume; this case starts from the original topography with fissural injection.
The solver still performs its normal topography import onto the simulation grid.

## Run the simulation

After building the solver in this repository, from the example directory:

```sh
OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE ../../bin/IMEX_SfloW2D
```

The input enables native NetCDF, ESRI thickness/temperature rasters and runout
diagnostics. Open `exampleEtna_lava.nc` directly in ParaView. Conservative text
output is disabled in the supplied input; the solver writes `restart.bin`.
Generated DEM and run outputs are ignored by Git, but are not deleted.

For an isolated validation directory, without running inside the repository:

```sh
python3 EXAMPLES/EXAMPLE_ETNA_LAVA/prepare.py --directory /absolute/path/to/new-run
```

This copies the unchanged input and original DEM into the selected directory.
Existing differing input/DEM files are never overwritten. Run the compiled
executable there by its absolute path.

## Validation evidence

The [validation report](../../docs/acceptance/ETNA_LAVA_EXAMPLE_2026_10_11.md)
records full 3600 s runs in strict/optimized builds with actual one/four-thread
teams and bitwise-equal decoded fields/checkpoints within each profile.
The flow reaches the domain boundary: the N9 mass budget must include numerical
boundary fluxes rather than treating the domain as closed. N9 approval remains
separate from successful example execution.
