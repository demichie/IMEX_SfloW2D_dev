# Etna DEM shared from the canonical example

The tracked archive exists only at `../../EXAMPLE_ETNA/DEM/Etna2014_crop.zip`.
From `EXAMPLE_ENSEMBLE_ETNA`:

```bash
cp ../EXAMPLE_ETNA/DEM/Etna2014_crop.zip DEM/
unzip DEM/Etna2014_crop.zip -d DEM
```

Alternatively, extract the canonical Etna archive and copy
`../EXAMPLE_ETNA/DEM/Etna2014_crop.asc` here. The ensemble template
uses a symlink to `../../DEM/Etna2014_crop.asc` from `templatedir/DEM`.
Do not commit a duplicated archive or uncompressed DEM in this example.
