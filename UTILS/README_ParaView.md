# Visualizing IMEX-SfloW2D NetCDF in ParaView

Open `OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc` in ParaView. In readers
that offer **Spherical Coordinates**, disable it for Cartesian model grids.
The solver writes `x`, `y` and `time` coordinates; inspect time steps
with the animation controls.

Current single-layer physical variables include `h` (depth, m),
`u` and `v` (horizontal velocity, m/s), `T` (temperature, K),
`rho_m` (mixture density, kg/m³), `solid_fraction` (when solids are
configured), and `mass` (kg/m²). Shared variables include `b`
(terrain elevation, m), `eta_top` (upper surface elevation, m),
`deposit_tot`, `hMax` and `pDynMax`. See
[NetCDF field documentation](../docs/NETCDF_OUTPUT.md) for the full
authoritative list and units. The prepared two-layer layout appends
`_l1`/`_l2` to layer-specific fields; two-layer execution is not yet
enabled for production.

Suggested visualization:

1. Colour the terrain grid by `b`. To create a relief surface, apply
   **Warp By Scalar** using `b` as the elevation scalar, selecting an
   appropriate scale factor.
2. Create a second surface with **Warp By Scalar** using `eta_top`.
   Apply **Threshold** on `h` to display only wet cells; prefer a
   physically meaningful positive depth threshold rather than
   thresholding temperature.
3. Colour the flow by `rho_m`, `T`, `h` or `solid_fraction`,
   depending on the configured physics. Use the colour and opacity
   controls for the desired view.
4. For velocity glyphs, use **Calculator** to define
   `u*iHat + v*jHat` and then **Glyph**. Adjust glyph density to avoid
   masking the surface.
5. For deposition, colour the terrain by `deposit_tot` and threshold
   positive thickness when desired.

ParaView menu labels and reader choices vary by installed version. Verify
that arrays are read as a Cartesian grid and that the selected time is
correct before interpreting the visualization.
