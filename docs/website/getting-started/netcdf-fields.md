# NetCDF fields and units

Enable `OUTPUT_NETCDF_FLAG=T` in `RUN_PARAMETERS` to write

```text
OUTPUT/<RUN_NAME>/fields/<RUN_NAME>.nc
```

The file is intended for plotting, scientific analysis and
visualization. Its physical fields use coordinates `x`, `y` and
`time`, with one scalar field per variable.

## Flow variables

For a single layer, fields have the following names. For two
layers, **layer-specific fields** receive suffix `_l1` or `_l2`.

| Variable | Meaning | Units |
| --- | --- | --- |
| `h` | Flow depth | m |
| `u`, `v` | Horizontal velocities | m/s |
| `T` | Mixture temperature | K |
| `mass` | Mobile layer mass per unit area | kg/m² |
| `rho_m` | Bulk mixture density (zero in dry cells) | kg/m³ |
| `solid_fraction` | Total solid **volume** fraction, when present | 1 |
| `solid_mass_fraction` | Total solid **mass** fraction, when present | 1 |
| `pore_pressure` | Excess pore pressure, when configured | Pa |
| `stochastic` | Transported stochastic variable, when configured | 1 |

For multiple solid classes, `solid_fraction_class1`,
`solid_fraction_class2` and so on hold individual volume
fractions. Their mass-fraction counterparts begin with
`solid_mass_fraction_class`.

For additional gas species, `gas_fraction_class1` and
`gas_mass_fraction_class1` distinguish gas volume and mass
fractions. Fields appear only when the corresponding
components are configured.

## Shared terrain and diagnostic fields

| Variable | Meaning |
| --- | --- |
| `b` | Terrain elevation [m] |
| `h_total` | Total flowing depth [m] |
| `eta_top` | Uppermost free-surface elevation [m] |
| `mass_total` | Total mobile mixture mass per unit area [kg/m²] |
| `deposit_tot` | Cumulative deposited material thickness [m] |
| `hMax` | Maximum flowing thickness [m] |
| `pDynMax` | Maximum dynamic pressure [Pa] |
| `modVelMax` | Maximum velocity magnitude [m/s] |

When several solid classes are present, separate substrate
inventory fields include `deposit_class1`, `erosion_class1`
and `erodible_class1`, with one field per configured class.

## One-dimensional simulations

A simulation with one computational cell in `x` or `y` is
exported with **two visualization samples** along that axis.
Values are duplicated only for display; the solver does not
add computational cells.

## Which file should you use?

- **NetCDF:** visualization and quantitative analysis of
  physical fields at output times.
- **`snapshots/`:** additional raster or conservative snapshots
  when their output options are enabled.
- **`restart/restart.bin`:** complete binary continuation state
  for restarting the same solver.

A NetCDF file is **not** a substitute for the binary checkpoint.
See [Output directories](output-files.md),
[Initial conditions and restart](initial-conditions.md)
and the [ParaView guide](visualization.md).
