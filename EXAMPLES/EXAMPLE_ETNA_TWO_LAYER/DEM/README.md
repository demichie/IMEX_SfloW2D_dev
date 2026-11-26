# Etna DEM input

From `EXAMPLES/EXAMPLE_ETNA_TWO_LAYER`, copy **only** the shared archive:

```bash
cp ../EXAMPLE_ETNA/DEM/Etna2014_crop.zip DEM/
python3 prepare.py
```

`prepare.py` uses the local ellipsoid scripts, extracts
`Etna2014_crop.asc` and produces `new_dem.asc`, `init_volume.asc`
and `ellipsoid.asc`. It does **not** execute the Fortran solver.
The two-layer executable reads `init_volume.asc` directly as the
initial L1 thickness. L2 starts dry via `RESTART_FILES='...','NONE'`.
