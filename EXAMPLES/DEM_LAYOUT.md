# Example terrain layout and ownership

Each example with a DEM has its own `DEM/` directory. Topography generators
write `DEM/topography_dem.asc` directly; their input templates specify the
same path with `TOPOGRAPHY_FILE` in `NEWRUN_PARAMETERS`. Keep checkpoint,
result and solver output files outside `DEM/`.

## Large and shared DEMs

| Owner | Committed archive | Reusing example | Preparation |
| --- | --- | --- | --- |
| `EXAMPLE_TAAL` | `DEM/topography_dem.zip` | `EXAMPLE_TAAL_BOTTOM` | Copy/unzip into its own `DEM/` as described in both README files. |
| `EXAMPLE_ETNA` | `DEM/Etna2014_crop.zip` | `EXAMPLE_ENSEMBLE_ETNA`, `EXAMPLE_ETNA_LAVA` | Copy/unzip or run the lava preparation script. |
| `EXAMPLE_CF` | `DEM/topo.zip` | — | Unzip in the owning `DEM/`. |

Generated ASC files and extracted shared archives are intentionally untracked.
The source raster is immutable: moving or extracting it must not change
its numeric values, shape, extent, resolution or NoData conventions. The
`sectors.asc` raster in the Etna ensemble is a categorical sector map, not
the topographic DEM, and remains at its original location.

The example root README and the README inside each `DEM/` describe the
required local preparation. Scripts should create `DEM/` automatically.
