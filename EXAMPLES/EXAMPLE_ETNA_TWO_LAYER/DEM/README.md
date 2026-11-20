# Etna DEM input

From `EXAMPLES/EXAMPLE_ETNA_TWO_LAYER`, copy the **only shared input file**
into this directory:

```bash
cp ../EXAMPLE_ETNA/DEM/Etna2014_crop.zip DEM/
```

The example-local `prepare.py` extracts `Etna2014_crop.asc` and runs the
example-local `create_input_ellipsoid.py` and `input_ellipsoid.py` to produce
`new_dem.asc`, `init_volume.asc`, and `ellipsoid.asc`. Do not commit the
copied ZIP or generated rasters. All other scripts and inputs are stored
inside this example.
