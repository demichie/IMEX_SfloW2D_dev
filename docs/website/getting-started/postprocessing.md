# Post-processing with Python

IMEX-SfloW2D provides optional Python plotting utilities for
simulation outputs. They run **after** a simulation; none of
these tools advances the numerical solver.

Install Python 3 with NumPy and Matplotlib, then run the tools
from the repository root. Replace the example paths and run
names with your own output directory.

## Plot ESRI ASCII results

```bash
python3 UTILS/plot_overlay.py EXAMPLES/EXAMPLE_ETNA/OUTPUT/exampleEtna
```

The script uses the saved terrain and snapshot manifest to
associate the raster data with the correct simulation time.
Use `--kind snapshots`, `diagnostics` or `vulnerability`
to limit which rasters are plotted.

Snapshot rasters must have been enabled using
`OUTPUT_ESRI_FLAG=T`. For NetCDF-only runs, use
[ParaView](visualization.md).

## Plot runout histories

```bash
python3 UTILS/plot_runout.py \
  EXAMPLES/EXAMPLE_ETNA/OUTPUT/exampleEtna/diagnostics/exampleEtna_runout.csv
```

This operates on a saved runout CSV, where output has been
configured. Figures are written beside the selected input.

## Plot probe histories

```bash
python3 UTILS/plot_probes.py \
  PATH/TO/OUTPUT/RUN_NAME/probes/RUN_NAME_prb_01.csv --rate
```

The script reads the probe coordinates, variable names and
units from the CSV header. The optional `--rate` argument
also produces finite-difference rate plots.

The utilities create figures alongside their inputs; they
do not modify the saved physical state.

For 3D terrain, flow surfaces, velocity arrows or satellite
overlays, follow [Visualizing results](visualization.md).
