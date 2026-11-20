# ParaView: step-by-step visualization of IMEX-SfloW2D output

This guide follows the original ParaView workflow (terrain, flow surface,
threshold, colours, opacity, velocity arrows, deposits and time annotation).
It updates the field names to the current solver output. Menu labels can
differ between ParaView versions; use **Filters > Search** (or the
**Alphabetical** submenu) when a named filter is not shown.

## 1. Open the NetCDF simulation

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

See [the full NetCDF reference](../docs/NETCDF_OUTPUT.md). For prepared
two-layer fields, layer-specific names use `_l1` and `_l2`;
production two-layer execution remains disabled.

## 2. Build a 3D terrain surface

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

## 3. Show the flow surface above the terrain

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

## 4. Add velocity arrows

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

## 5. Display deposited material

1. Select the original NetCDF input and apply
   **Filters > Alphabetical > Warp By Scalar** with `b` again.
2. In **Coloring**, choose `deposit_tot`.
3. Select this new warp and apply **Filters > Alphabetical > Threshold**.
4. Select `deposit_tot` as the threshold variable, choose a positive
   minimum thickness appropriate to the run, and click **Apply**.
5. Adjust the colour preset or the **Color Map Editor** as desired.

## 6. Add simulation time and save the view

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
