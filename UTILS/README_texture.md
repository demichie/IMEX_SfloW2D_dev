# QGIS and ParaView: satellite texture for the simulation terrain

This guide preserves the original menu-by-menu workflow while updating
paths and clarifying the coordinate-system checks. QGIS, QuickMapServices
and ParaView versions may show slightly different menu labels.

## 1. Identify the coordinate reference system

The solver stores Cartesian `x`/`y` coordinates but does **not**
define a geographic CRS. Read the source DEM metadata or its data
provider's documentation to determine the correct CRS (for example,
a particular WGS 84 / UTM zone only **if that is the source CRS**).
Do not infer the zone from an approximate location or change the CRS
label to make rasters appear aligned.

The simulation grid export is usually
`OUTPUT/<RUN_NAME>/terrain/dem_esri.asc`. The original terrain
can also be under the example's `DEM/` folder. Use the computational
export when you need **exactly the simulated footprint**.

## 2. Load the DEM into QGIS

1. Open **QGIS**.
2. Choose **Layer > Add Layer > Add Raster Layer** (or press
   **Ctrl+Shift+R**, where supported).
3. In the raster **Source** selector, browse to `dem_esri.asc`
   or the prepared source DEM. Click **Add**, then **Close**.
4. Check the CRS in the lower-right status bar. Click that indicator
   or choose **Project > Properties > CRS**.
5. Search for and select the DEM's verified CRS (e.g. the correct
   **WGS 84 / UTM zone**), then **Apply > OK**. If QGIS asks for the
   layer's CRS separately, use the layer context menu
   **Layer CRS > Set Layer CRS** to assign its *known* CRS.
   This does **not** reproject the raster.
6. Right-click the DEM in **Layers** and select **Zoom to Layer**.
   In **Layer Properties > Information**, confirm extent and pixel size.

## 3. Add satellite or orthophoto imagery

### Option A: QuickMapServices (original workflow)

1. Choose **Plugins > Manage and Install Plugins**.
2. Search for **QuickMapServices** and install/enable it.
3. If the installed version exposes them, choose
   **Web > QuickMapServices > Settings > More services >
   Get contributed pack** to add community map providers.
4. Open **Web > QuickMapServices**, select an appropriate imagery
   provider and add its satellite/orthophoto layer. The old guide
   specifically referred to **Google > Google.cn Satellite**;
   this exact entry is not guaranteed to exist in current provider
   lists and is not required.
5. Check that the imagery source permits your intended use and export.

### Option B: QGIS XYZ Tiles / raster source

1. In the **Browser** panel, expand **XYZ Tiles**.
2. Choose an existing permitted imagery connection, or right-click
   **XYZ Tiles > New Connection** and enter your provider's URL
   and settings.
3. Add the layer to the map. QGIS can reproject it on screen for
   display, but exports still require careful CRS and extent checks.

## 4. Export an image covering the simulation domain

The original workflow used **Convert map to raster**, then clipped
the result and exported an image. Preserve those options:

1. With the DEM and imagery loaded, open
   **Processing > Toolbox**.
2. Search for **Convert map to raster** in the toolbox's
   raster tools (the exact group name can vary).
3. For **Minimum extent to render** / **Extent**, use the
   **...** menu and choose **Calculate from Layer** or
   **Use Layer Extent**, selecting `dem_esri`.
4. Set **Map units per pixel** to the target resolution, for example
   `10` if the DEM CRS is in metres and that resolution is suitable.
   Choose an output raster and click **Run**.
5. If necessary, select **Raster > Extraction >
   Clip Raster by Extent**; for **Clipping extent** choose
   **Calculate from Layer / Use Layer Extent** with the DEM,
   then **Run**. Some installations expose this tool through
   **Processing > Toolbox > GDAL > Raster extraction**.
6. For a plain image suitable for a ParaView texture, choose
   **Project > Import/Export > Export Map to Image**.
   In the export dialog set the extent from the DEM layer,
   choose dimensions or DPI, and **Save**.
7. Open the resulting image and verify that its four corners,
   width/height ratio and north-up orientation match the DEM.
   Preserve the CRS/georeference separately for reproducibility.

**Important:** a screenshot or map export may include map decorations,
resampling or a different extent. For quantitative matching, prefer a
georeferenced raster processing/export workflow; use the plain image
only as a visualization texture.

## 5. Apply the texture in ParaView

1. Start ParaView and choose **File > Open**.
2. Open `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc`,
   uncheck **Spherical Coordinates** if the reader offers it,
   and click **Apply**.
3. In **Pipeline Browser**, select the original NetCDF input.
   Choose **Filters > Alphabetical > Warp By Scalar**,
   select `b` under **Scalars**, use vertical normal
   `0, 0, 1`, and click **Apply**.
4. Switch the render-view camera from 2D to 3D if needed.
5. In **Pipeline Browser**, select the terrain
   **Warp By Scalar**, then choose
   **Filters > Alphabetical > Texture Map to Plane**
   (or use **Filters > Search**). Click **Apply**.
6. In **Properties** for **Texture Map to Plane**, enable
   **Advanced Properties** using the gear/wheel icon if the
   plane controls are hidden.
7. Set the texture plane's origin and two axes to span the
   DEM's `x`/`y` extent; do not rely on automatic orientation
   if the image is rotated or mirrored.
8. In the **Properties** / **Display** section, locate
   **Texture** (often under **Miscellaneous** in older ParaView).
   Choose **Load...** and select the image exported from QGIS.
9. If the texture is hidden by scalar colours, set
   **Coloring > Solid Color** or disable scalar colouring
   for this display. Use **Reset Camera** if needed.
10. Compare at least three identifiable geographic features
    (including ones away from the domain centre) to confirm
    translation, rotation, scale and orientation.

If the overlay is mirrored, **first** check the CRS,
raster export's top/bottom row orientation and texture mapping
plane. Flipping an image manually should be a documented,
last-resort correction rather than a trial-and-error procedure.

See [the ParaView guide](README_ParaView.md) for colour maps,
wet-cell thresholding, velocities, deposits and time animation.
