# Post-processing utilities

These optional Python utilities operate on **existing solver outputs**. They
do not participate in the numerical solver or its build.

Install NumPy and Matplotlib for plotting. From the repository root:

```sh
python3 UTILS/plot_overlay.py EXAMPLES/EXAMPLE_ETNA/OUTPUT/exampleEtna
python3 UTILS/plot_runout.py EXAMPLES/EXAMPLE_ETNA/OUTPUT/exampleEtna/diagnostics/exampleEtna_runout.csv
python3 UTILS/plot_probes.py PATH/TO/OUTPUT/RUN_NAME/probes/RUN_NAME_prb_01.csv --rate
```

Replace example paths with the actual run name and output directory;
`OUTPUT_BASE_DIRECTORY` may change its location. Graphics are written next
to the corresponding source raster or CSV. Plotting existing PNG/PDF files
does not modify solver outputs.

- [ParaView and NetCDF](README_ParaView.md)
- [Georeferenced textures](README_texture.md)

`plot_overlay.py` plots ESRI ASCII rasters using the exported computational
terrain (`terrain/dem_esri.asc`) and the exact snapshot time in
`metadata/output_manifest.csv`. Select `--kind snapshots`, `diagnostics`,
or `vulnerability` to limit outputs. Snapshot rasters require
`OUTPUT_ESRI_FLAG=T`; the NetCDF-only workflow should use ParaView instead.
The utility does not construct combined vulnerability tables.

`plot_probes.py` reads the solver's three-line CSV header (coordinates,
variable names, units). `--rate` adds finite-difference plots.
`plot_runout.py` reads the runout CSV header and data columns.

Run focused script checks with:

```sh
python3 -m unittest discover -s UTILS -p 'test_*.py'
```
