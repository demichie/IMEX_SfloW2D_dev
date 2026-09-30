# Native NetCDF output

Set `OUTPUT_NETCDF_FLAG = T` in `RUN_PARAMETERS` to write the scientific
solution directly to `<run_name>.nc`. NetCDF is the only multidimensional
physical-solution format produced by the solver.

The conservative `.q_2d` files remain available through `OUTPUT_CONS_FLAG`.
They have a distinct purpose: restart and conservative-state diagnostics.
ESRI rasters, probes, runout, mass-centre and other specialized outputs are
also unchanged.

## Legacy `.p_2d` field audit

The following matrix records the migration of every field in the former
formatted physical output. Coordinates are one-dimensional NetCDF coordinate
variables; all physical fields have dimensions `(x, y, time)`.

| Former `.p_2d` field | Native NetCDF variable |
| --- | --- |
| `x`, `y` | `x`, `y` |
| flow thickness | `h` |
| bed and free surface | `b`, `eta` |
| velocity components | `u`, `v`, `w` |
| solid fractions | `solid_frac_XX` |
| additional-gas fractions | `add_gas_frac_XX` |
| liquid fraction | `alphal` |
| temperature | `T` |
| mixture density | `rhomix` |
| reduced gravity | `red grav` |
| deposit per solid class | `deposit__XX` |
| erosion per solid class | `erosion_XX` |
| total available erodible thickness | `erodible` |
| basal shear velocity | `shear_velocity` |
| Richardson number | `Ri` |
| Rouse number per solid class | `Rouse_XX` |
| stochastic variables | `Zs_XX` |
| absolute pore pressures | `PorePres_XX` |
| effective friction coefficient | `mu eff` |
| maximum flow thickness | `hMax` |
| maximum dynamic pressure | `pDynMax` |
| maximum velocity magnitude | `modVelMax` |
| granular inertial number | `inertial_number` |

`XX` is the one-based, zero-padded class or variable index. Fields that are
not active for a selected physical model are retained in the schema and
written as zero, matching the former fixed-column output behavior.

The free-surface elevation is `eta = b + h`. The vertical velocity `w` is
computed at cell centers as `u * db/dx + v * db/dy` from the horizontal
velocity and the bed-slope fields.

NetCDF files written with the previous schema, where `w` stored free-surface
elevation, cannot be appended to with the current schema. Use a new output
filename when continuing a run from one of those files.
