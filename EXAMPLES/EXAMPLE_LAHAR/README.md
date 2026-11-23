# Vesuvius LAHAR — single-layer simulation

## Overview

This example simulates the propagation of a **lahar** (a water–sediment
mixture) over Vesuvius terrain using the two-dimensional, depth-averaged
single-layer model (`N_LAYERS=1`). The flow is initialized from a spatially
distributed thickness raster and contains **six solid grain-size classes**
suspended in a liquid carrier. The model accounts for momentum resistance,
sediment erosion and deposition, and loss of the liquid carrier.

The simulation covers **600 s** of physical time. The computational domain
is a grid-aligned portion of the original terrain, chosen to encompass the
flow while avoiding unnecessary calculation outside the area of interest.

## Domain and initial conditions

| Parameter | Value |
| --- | --- |
| Run name | `VESUVIO_001_001` |
| Lower-left corner | `X0=440995 m`, `Y0=4515991 m` |
| Grid | `340 x 280` cells at `50 m` |
| Extent | `440995–457995 m` (X), `4515991–4529991 m` (Y) |
| Initial velocity | `u=v=0 m/s` |
| Initial mixture temperature | `373 K` |
| Initial solid volume fraction | `51%` in total, distributed among six classes |
| Simulation interval | `0–600 s` |
| Output interval | `600 s` |
| CFL and maximum timestep | `0.24` and `5 s` |

All four sides use the homogeneous boundary values specified in
`IMEX_SfloW2D.inp`. The initial flow geometry is given by a raster rather
than a prescribed geometric collapse or an active injection source.

## Input rasters

Three ESRI ASCII (`.asc`) rasters are read from `DEM/`:

| Raster | Purpose |
| --- | --- |
| `ves_asciiNEW.asc` | Ground-surface elevation used to construct the computational topography and slopes. |
| `VESUVIO_001_001_erodible.asc` | Initial thickness of the erodible substrate available to the sediment-erosion process. |
| `masked_init_vol_thick_001_001.asc` | Initial lahar-flow thickness, combined with the initial temperature, velocity and solid composition specified in the input. |

The input enables `SUBTRACT_INIT_FLAG=T`: the initial flow thickness is
subtracted during terrain/erodible-bed initialization, as defined by the
solver's initialization procedure. The model then uses the resulting
topography and substrate configuration.

The rasters are distributed as three separate ZIP archives; their original
`.asc` filenames must be preserved on extraction. See
[DEM/README.md](DEM/README.md) for the archive names and details.

## Physical processes

**Mixture and rheology.** The flow comprises a liquid carrier
(`RHO_L=1000 kg/m³`) and six solid classes, with particle diameters
from **16 mm to 0.015625 mm** and densities from **1400 to 2200 kg/m³**.
The selected `lahar` rheology combines concentration- and
temperature-dependent carrier viscosity, viscous and turbulent resistance,
and a concentration-dependent yield-stress contribution. Its coefficients
are specified in `RHEOLOGY_PARAMETERS`.

**Erosion.** Bed erosion is enabled with coefficient `1e-4`. Eroded
material is drawn from the finite erodible inventory, whose composition is
given by `ERODIBLE_FRACT` and whose initial thickness is read from
`VESUVIO_001_001_erodible.asc`. The configured substrate porosity is
`0.21`.

**Deposition.** Sediment settling/deposition is enabled with
`MASS_EXCHANGE%DEPOSITION%SETTLING_FLAG=T`. The process transfers
suspended solid material to the substrate, accounting for the six
particle classes.

**Liquid-carrier loss.** The prescribed carrier-liquid withdrawal rate is
`1e-8 m/s`. Erosion, deposition and liquid loss are handled in the
explicit mass-exchange source group (`METHOD='EXPLICIT'`, `GROUP=4`).

**Disabled processes.** Atmospheric entrainment is disabled, and there is
no radial injection or collapsing-volume source. Dynamic topographic
evolution is disabled (`SUBSTRATE%TOPO_CHANGE_FLAG=F`), although
sediment inventories and deposition/erosion diagnostics are still tracked.
No additional gas phase or second flowing layer is used.

## Run the example

From the repository root:

```bash
cd EXAMPLES/EXAMPLE_LAHAR
```

Extract the three DEM archives into the same `DEM/` directory:

```bash
unzip -o DEM/ves_asciiNEW.zip -d DEM
unzip -o DEM/VESUVIO_001_001_erodible.zip -d DEM
unzip -o DEM/masked_init_vol_thick_001_001.zip -d DEM
```

With the IMEX-SfloW2D executable built at `../../src/IMEX_SfloW2D`,
start the simulation:

```bash
../../src/IMEX_SfloW2D
```

The solver uses OpenMP by default. To specify the thread count explicitly,
for example:

```bash
OMP_NUM_THREADS=4 ../../src/IMEX_SfloW2D
```

## Output files

The simulation writes to `OUTPUT/VESUVIO_001_001/`. With the output flags
in the supplied input file, the principal files are:

| Path relative to `OUTPUT/VESUVIO_001_001/` | Contents |
| --- | --- |
| `fields/VESUVIO_001_001.nc` | NetCDF fields over time: flow thickness, velocities, temperature, mass, mixture density, solid fractions by class, terrain, deposition, erosion, remaining erodible sediment, and maximum-flow diagnostics. |
| `restart/restart.bin` | Binary checkpoint used to continue the simulation. |
| `metadata/VESUVIO_001_001.bak` | Backup of the run parameters written by the solver. |
| `metadata/output_manifest.csv` | Output snapshot indices and their physical times. |

`OUTPUT_NETCDF_FLAG=T` enables the NetCDF output, with
`DT_OUTPUT=600 s`. The current configuration has
`OUTPUT_CONS_FLAG=F`, `OUTPUT_ESRI_FLAG=F` and
`OUTPUT_RUNOUT_FLAG=F`, so it does **not** write separate conservative
`.q_2d` snapshots, ESRI field rasters, or runout CSV files.
