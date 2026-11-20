# NetCDF physical fields for visualization

Set `OUTPUT_NETCDF_FLAG=T` under `RUN_PARAMETERS` to write
`OUTPUT/<run_name>/fields/<run_name>.nc`.

The output is a single schema of scalar physical fields with Fortran dimensions
`(x,y,time)`, readable as arrays in ParaView. The internal conservative `q`
and primitive `qp` arrays are **not** copied as multidimensional NetCDF
tensors. They remain authoritative in memory and in the binary checkpoint.

For one layer the field names are unsuffixed. With `N_LAYERS=2`, layer-specific
fields have suffix `_l1` or `_l2`; both layers share the same field convention.

| Variable (single layer) | Meaning | Units |
| --- | --- | --- |
| `h`, `u`, `v` | Depth and horizontal velocity components | m, m/s, m/s |
| `T` | Temperature | K |
| `mass` | Layer mass per unit area | kg/m² |
| `rho_m` | Bulk density, computed as mass/depth in wet cells; zero in dry cells | kg/m³ |
| `solid_fraction` | Total solid volume fraction (when solids are configured) | 1 |
| `pore_pressure` | Excess pore pressure (only when its state component is configured) | Pa |
| `stochastic` | Transported stochastic state (when configured) | 1 |

For multiple solid classes, `solid_fraction_class1`, `solid_fraction_class2`,
etc. expose individual volume fractions. Extra gas classes use
`gas_fraction_class1`, etc. Shared solid inventories are
`deposit_class1`, `erosion_class1`, `erodible_class1`, etc.
Each is a separate `(x,y,time)` field; no extra `component`, `layer`
or `solid_class` dimensions occur in visualization variables.

The shared fields `b`, `h_total`, `eta_top`, `mass_total`,
`deposit_tot`, `hMax`, `pDynMax`, `modVelMax` retain their
previous meanings and units. Coordinates are `x`, `y`, and `time`.

The binary checkpoint remains the **only** restart state file. It preserves the
complete authoritative conservative arrays and all restart metadata; NetCDF is
for visualization and analysis only. Restart continuation appends records to
an existing file that already uses the current visualization schema; there is
no compatibility path for pre-migration NetCDF histories. Regenerate those
outputs with the development solver.

The dam-break `height_summary_nc.py` reads the physical `h` variable.

## Verification

The structural two-layer `TEST_SPATIAL_OPERATOR` test checks scalar fields
in both layers, absent obsolete tensors, solid inventories and appended time
records. Run the relevant actions and a single-layer dam-break example to
verify the full output chain in the current environment.

## Degenerate spatial axes (1D simulations)

When `COMP_CELLS_X=1` or `COMP_CELLS_Y=1`, native NetCDF visualization
uses two samples on that axis. Every physical field is duplicated along the
missing dimension, without changing the computational cells, state, budgets,
or binary restart. The second visualization coordinate is offset by 4% of
the physical span along the other axis, matching the v2 Python converter;
a 1-by-1 computational grid uses its cell spacing as a fallback. NetCDF
history continuation validates the *visualization* dimensions, not the
computational dimensions. Regenerate files created with the previous
single-sample NetCDF schema before appending.
