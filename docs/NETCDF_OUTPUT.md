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
| dispersed-liquid fraction (removed) | not emitted; no dispersed-liquid conservative component in the current state |
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

## Prepared two-layer output schema

The production `N_LAYERS=2` gate remains closed until the functional and
scientific acceptance gates. The I/O layer is prepared and exercised by
standalone two-layer fixtures.

For two-layer runs, **all existing NetCDF variables keep their historical
L1 meaning** and the existing `(x,y,time)` dimensions. In addition, the
same file contains:

| Variable | Dimensions | Meaning |
| --- | --- | --- |
| `q_layers` | `(component,layer,x,y,time)` | All conservative variables in the canonical solver order |
| `qp_layers` | `(component,layer,x,y,time)` | All primitive variables in the canonical solver order |
| `eta_L2` | `(x,y,time)` | Upper free surface, `B+h1+h2` |
| `h_total` | `(x,y,time)` | Total fluid depth, `h1+h2` |
| `mass_L1` | `(x,y,time)` | L1 conserved mass per horizontal area, `q(1,1,:,:)` |
| `mass_L2` | `(x,y,time)` | L2 conserved mass per horizontal area, `q(1,2,:,:)` |
| `mass_total` | `(x,y,time)` | L1+L2 conserved mass per horizontal area |

The `component` dimension is `n_vars`; `layer` has length 2 and
stores L1 before L2. Both tensors preserve inactive optional components
as zero. The bed, erosion/deposition and substrate fields are shared,
not duplicated per layer. Restart appends require the same two-layer
schema; a single-layer NetCDF file cannot be silently reused as a
two-layer file or vice versa.

If `OUTPUT_CONS_FLAG=T`, the original `<run>_<index>.q_2d`
still denotes L1 and `<run>_L2_<index>.q_2d` contains L2, with the
same component ordering and snapshot index. The extra ESRI rasters are
`<run>_L2_<index>.asc`, `<run>_T_L2_<index>.asc`,
`<run>_total_depth_<index>.asc` and `<run>_eta_L2_<index>.asc`.
The original DEM is still written only once. In addition, every
two-layer output produces `fields/<run>_budget_<index>.csv`, recording
L1/L2/total mass and volume and the total solid inventory including
the shared deposited and erodible substrate.

## Prepared two-layer binary restart

Single-layer checkpoints retain **format v7, byte for byte**. Two-layer
checkpoints use format **v8**, with an explicit layer count in the
dimension header and both conservative layer arrays. The shared terrain,
substrate, output schedule and hazard histories are written once.
The L1 OU field and intrinsic random stream retain their previous
positions and meanings. V8 appends the L2 OU field and the independent
L2 RNG stream. The reader validates file version, layer count, component
layout, physical flags and random-state length before using the data;
it reconstructs derived topographic geometry rather than saving a
second geometry representation.

For a two-layer **text cold-start**, provide exactly two `.q_2d`
files using `RESTART_FILES(1)` and `RESTART_FILES(2)` and set
`N_RESTART_FILES=2`. The files must use the same output index and
common equation layout, with disabled optional components zero.
Legacy single-layer `.asc` and `.q_2d` inputs are unaffected.
A pair of thickness-only `.asc` files is **not** an accepted
two-layer initialization: it cannot independently describe all
physical and optional layer components.

The integrated spatial-operator regression writes and reads a v8
checkpoint, compares the independent L2 RNG and OU state, reads back
both NetCDF tensors and the total-mass field, and checks separate
L1/L2 text initialization. Existing single-layer restart regressions
continue to verify the v7 byte-level contract.
