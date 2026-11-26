# Canonical Taal DEM

`topography_dem.zip` is the only tracked copy of the Taal terrain
(`topography_dem.asc`, stored with lossless ZIP compression).

From `EXAMPLES/EXAMPLE_TAAL`:

```bash
unzip DEM/topography_dem.zip -d DEM
```

The solver input selects `DEM/topography_dem.asc`. The uncompressed raster
is generated locally; do not commit it. `EXAMPLE_TAAL_BOTTOM` must obtain
its terrain from this directory rather than keeping a duplicate archive.
