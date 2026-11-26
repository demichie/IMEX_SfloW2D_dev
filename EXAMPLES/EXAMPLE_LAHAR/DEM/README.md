# DEM inputs for EXAMPLE_LAHAR

This directory contains the three cropped ESRI ASCII rasters required by
`../IMEX_SfloW2D.inp`. Large rasters are committed as separate ZIP files,
**one ZIP per ASC**, rather than as uncompressed `.asc` files.

| ZIP to upload | ASC expected by the solver | Meaning |
| --- | --- | --- |
| `ves_asciiNEW.zip` | `ves_asciiNEW.asc` | Terrain elevation |
| `VESUVIO_001_001_erodible.zip` | `VESUVIO_001_001_erodible.asc` | Erodible-bed thickness |
| `masked_init_vol_thick_001_001.zip` | `masked_init_vol_thick_001_001.asc` | Initial LAHAR thickness |

Each ZIP archive contains its corresponding `.asc` raster at the archive
root, without an additional `DEM/` directory.

From `EXAMPLES/EXAMPLE_LAHAR`, test and extract the archives:

```bash
unzip -t DEM/ves_asciiNEW.zip
unzip -t DEM/VESUVIO_001_001_erodible.zip
unzip -t DEM/masked_init_vol_thick_001_001.zip

unzip -o DEM/ves_asciiNEW.zip -d DEM
unzip -o DEM/VESUVIO_001_001_erodible.zip -d DEM
unzip -o DEM/masked_init_vol_thick_001_001.zip -d DEM

test -f DEM/ves_asciiNEW.asc
test -f DEM/VESUVIO_001_001_erodible.asc
test -f DEM/masked_init_vol_thick_001_001.asc
```

The simulation grid is `340 x 280` cells at `50 m`, with origin
`(440995, 4515991)`. The raster crop should cover this domain and
preserve original elevation/thickness values, NODATA values, grid
registration, and the original DEM resolution. Cropping the raster
does not change the solver grid or the physical input parameters.

The repository-wide `EXAMPLES/.gitignore` excludes extracted
`DEM/*.asc` files. Commit only the three ZIP archives listed above
and this README. The solver reads the extracted ASCII files, not ZIPs.
