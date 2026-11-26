# Visualizing IMEX-SfloW2D results with ParaView and QGIS

This guide brings together the complete step-by-step ParaView workflow and
the QGIS-to-ParaView satellite-texture workflow. It preserves the explicit
menu paths and filter controls from the original utility guides. Menu names
may differ across versions of ParaView, QGIS and QuickMapServices.

For plotting ESRI ASCII rasters or diagnostic CSV series with the repository's
Python scripts, see the [post-processing utilities](https://github.com/demichie/IMEX_SfloW2D_dev/tree/master/UTILS).

## ParaView: simulation fields and 3D visualization

This guide follows the original ParaView workflow (terrain, flow surface,
threshold, colours, opacity, velocity arrows, deposits and time annotation).
It updates the field names to the current solver output. Menu labels can
differ between ParaView versions; use **Filters > Search** (or the
**Alphabetical** submenu) when a named filter is not shown.

### 1. Open the NetCDF simulation

1. Start ParaView and choose **File > Open**.
2. Browse to `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc`. Select the file
   and click **OK**.
3. In **Properties** on the left, choose the appropriate NetCDF reader
   when prompted. If **Spherical Coordinates** is offered, uncheck it
   for a Cartesian model grid, then click **Apply**.
4. In **Pipeline Browser**, select the loaded dataset. Check the
   **Information** tab for the domain dimensions and available arrays.
   Use the time selector on the top toolbar to choose a snapshot.
5. If the view is flat, select the **2D/3D** camera control (or the
   3D interaction mode) above the render view.

Current single-layer variables include:

| Variable | Meaning | Units |
| --- | --- | --- |
| `b` | Terrain elevation | m |
| `h` | Flow depth | m |
| `eta_top` | Upper free-surface elevation | m |
| `T` | Flow temperature | K |
| `u`, `v` | Horizontal velocity components | m/s |
| `rho_m` | Bulk mixture density | kg/m³ |
| `solid_fraction` | Total solid volume fraction, when configured | 1 |
| `deposit_tot` | Total deposited thickness | m |
| `hMax` | Maximum thickness | m |
| `pDynMax` | Maximum dynamic pressure | Pa |

See [the full NetCDF reference](https://github.com/demichie/IMEX_SfloW2D_dev/blob/master/docs/NETCDF_OUTPUT.md). For two-layer simulations, layer-specific fields use `_l1` and `_l2`;
`N_LAYERS=2` runs through the standard integrated solver. Interpret the coupled
fields within the [two-layer hyperbolicity limits](../model/two-layer.md).

### 2. Build a 3D terrain surface

1. In **Pipeline Browser**, select the original `.nc` dataset, not
   another filter.
2. Choose **Filters > Alphabetical > Warp By Scalar**, or use
   **Filters > Search** and search for `Warp By Scalar`.
3. In **Properties**, set **Scalars** to `b`. Set **Normal** to
   `0, 0, 1` (or equivalent vertical displacement) and
   **Scale Factor** to `1` for metre coordinates and elevations.
4. Click **Apply**. Use **Reset Camera** if the terrain is off screen;
   select 3D interaction mode if necessary.
5. In the colouring drop-down at the top of the render view (often
   displaying **Solid Color**), choose `b`.
6. Click **Edit Color Map** to open **Color Map Editor**. Use the
   **Choose Preset** control (folder/heart icon in some versions)
   to select a colour palette. Close the editor when finished.

If `b` is stored as **Cell Data** and **Warp By Scalar** cannot select it,
apply **Filters > Alphabetical > Cell Data to Point Data** to the input
first, then apply the warp.

### 3. Show the flow surface above the terrain

1. Re-select the **original NetCDF input** in **Pipeline Browser**.
2. Apply **Filters > Alphabetical > Warp By Scalar** again.
3. Set **Scalars** to `eta_top`, **Normal** to `0, 0, 1`,
   **Scale Factor** to `1`, then click **Apply**.
4. Select this second **Warp By Scalar** filter in **Pipeline Browser**
   and choose **Filters > Alphabetical > Threshold**.
5. In **Properties**, set **Scalars** to `h`. Set the lower
   threshold just above zero (for example `0.001` m if appropriate
   to the simulation) and retain the upper range. Click **Apply**.
   This avoids showing the dry area; unlike the old procedure,
   **do not threshold temperature against ambient temperature**.
6. In the colouring drop-down choose `rho_m`, `solid_fraction`,
   `T` or another meaningful physical field. If a field is not
   available for the configuration, choose another.
7. Open **Edit Color Map**, choose the desired preset, and, for
   transparency, enable **opacity mapping for surfaces** if supported.
   In **Opacity Transfer Function** adjust opacity points. You may
   also use **Properties > Display > Opacity** for uniform transparency.

Keep the terrain and flow pipelines both visible using the eye icons
in **Pipeline Browser**. A small vertical offset may help distinguish
coincident surfaces visually, but must not be mistaken for physical depth.

### 4. Add velocity arrows

1. Select the flow **Threshold** result in **Pipeline Browser**.
2. Choose **Filters > Alphabetical > Calculator**.
3. In **Properties**, set **Result Array Name** to `Vel`, enter
   `u*iHat + v*jHat` in the expression field, and click **Apply**.
4. Select **Calculator**, then choose
   **Filters > Alphabetical > Glyph**.
5. In **Properties**, select **Glyph Type > Arrow**,
   **Orientation Array > Vel**, and (optionally) **Scale Array > Vel**.
   Adjust **Scale Factor** to the model length scale and click **Apply**.
6. If too many or no arrows appear, adjust **Glyph Mode** or
   **Masking** (e.g. **Every Nth Point** or a maximum number of
   samples, depending on version).
7. In **Coloring**, choose `Vel` or **GlyphVector** when available.

The legacy UI used a **Scaling > vector** selector; newer versions may
show **Orientation Array**, **Scale Array**, and **Glyph Mode** instead.

### 5. Display deposited material

1. Select the original NetCDF input and apply
   **Filters > Alphabetical > Warp By Scalar** with `b` again.
2. In **Coloring**, choose `deposit_tot`.
3. Select this new warp and apply **Filters > Alphabetical > Threshold**.
4. Select `deposit_tot` as the threshold variable, choose a positive
   minimum thickness appropriate to the run, and click **Apply**.
5. Adjust the colour preset or the **Color Map Editor** as desired.

### 6. Add simulation time and save the view

1. Select an active pipeline object and choose
   **Filters > Alphabetical > Annotate Time** (called
   **Annotate Time Filter** in some versions).
2. Set **Format** to `Time: %5.2f s` when this option is present,
   and click **Apply**. Position the text in the render view.
3. Use the top time controls to move between output frames.
4. Save a still image via **File > Save Screenshot**; to save a
   movie choose **File > Save Animation**.
5. To reuse the complete pipeline in ParaView, choose
   **File > Save State**, and reopen it with **File > Load State**.

Before interpreting results, verify the grid orientation, the displayed
time, the field units and that vertical warp uses the intended geometry.

---

## QGIS and ParaView: georeferenced satellite textures

This guide preserves the original menu-by-menu workflow while updating
paths and clarifying the coordinate-system checks. QGIS, QuickMapServices
and ParaView versions may show slightly different menu labels.

### 1. Identify the coordinate reference system

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

### 2. Load the DEM into QGIS

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

### 3. Add satellite or orthophoto imagery

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

### 4. Export an image covering the simulation domain

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

### 5. Apply the texture in ParaView

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

For colour maps, wet-cell thresholds, velocity glyphs, deposited material and
animation, follow the [ParaView section](#paraview-simulation-fields-and-3d-visualization)
above.
