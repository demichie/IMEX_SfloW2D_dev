# Shared Taal DEM (not duplicated in Git)

The canonical raster is in `../../EXAMPLE_TAAL/DEM/`. Prepare this example
from `EXAMPLES/EXAMPLE_TAAL_BOTTOM` with either option:

**Copy the ZIP and extract it locally:**

```bash
cp ../EXAMPLE_TAAL/DEM/topography_dem.zip DEM/
unzip DEM/topography_dem.zip -d DEM
```

**Or copy the already extracted ASC:**

```bash
cp ../EXAMPLE_TAAL/DEM/topography_dem.asc DEM/
```

The second option requires extracting the parent archive first. The solver
expects `DEM/topography_dem.asc`. Neither the local extracted ASC nor the
copied ZIP should be committed here; the authoritative archive stays in
`EXAMPLE_TAAL/DEM`.
