# Georeferenced terrain texture with QGIS and ParaView

A satellite/orthophoto texture can be overlaid on a terrain surface
provided its coordinate reference system (CRS), geographic extent and
orientation match the solver DEM. The numerical grid itself does not
supply a CRS: determine it from the source DEM metadata and the original
data provider; do not assume a UTM zone.

1. In QGIS, load the source DEM from the example's `DEM/` directory.
   Assign its **actual** CRS (assigning a wrong CRS does not reproject
   coordinates).
2. Add a licensed orthophoto/satellite imagery source using a supported
   QGIS raster or XYZ connection. Ensure that its terms permit export
   and reuse.
3. Reproject the imagery into the DEM's CRS if required.
4. Use the DEM's geographic extent when clipping/exporting the image;
   preserve pixel orientation, aspect ratio, and any associated
   georeferencing information.
5. In ParaView, load the solver's NetCDF terrain (`b`), create a
   **Warp By Scalar** terrain, and apply **Texture Map to Plane**
   with its plane aligned to the DEM's x/y domain.
6. Load the exported image as the texture and verify correspondence at
   several known geographic features before using it in analysis.

The computational terrain exported at
`OUTPUT/<RUN_NAME>/terrain/dem_esri.asc` can help confirm the plotting
extent. Avoid trial-and-error flips or transposes: correct the CRS,
raster geotransform or mapping-plane orientation at the source.
