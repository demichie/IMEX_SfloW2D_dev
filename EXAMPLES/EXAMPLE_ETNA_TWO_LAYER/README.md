# Etna two-layer — standalone experimental example

Run every command from **EXAMPLES/EXAMPLE_ETNA_TWO_LAYER**. All scripts and
inputs needed for preprocessing are inside this directory. The **only file**
copied from `EXAMPLE_ETNA` is the shared DEM ZIP; it is not tracked twice.

This is an experimental TL13 two-layer example, **not a scientifically
validated prediction**. Running `N_LAYERS=2` still requires
`IMEX_TWO_LAYER_VALIDATION=1`.

## 1. Copy the shared DEM once

After entering this directory from the repository root:

```bash
cd EXAMPLES/EXAMPLE_ETNA_TWO_LAYER
cp ../EXAMPLE_ETNA/DEM/Etna2014_crop.zip DEM/\nunzip DEM/Etna2014_crop.zip -d DEM/
```

Requirements for preprocessing: Python 3, NumPy and pandas.
The compiled solver is expected at `../../src/IMEX_SfloW2D`.

The DEM extraction is a manual prerequisite, not part of `prepare.py`.\nCheck that `DEM/Etna2014_crop.asc` exists before proceeding.\n\n## 2. Prepare only the topography and initial thickness

```bash
python3 prepare.py
```

This **never starts the Fortran program**. It runs the local
`create_input_ellipsoid.py` with parameters from `input_ellipsoid.py`
and produces:

- `DEM/Etna2014_crop.asc`: original topography (extracted manually beforehand).
- `DEM/new_dem.asc`: modified terrain.
- `DEM/init_volume.asc`: initial L1 thickness raster.
- `DEM/ellipsoid.asc`: excavation geometry for inspection.

The ellipsoid uses the original Etna center `(500430, 4177500)` and target
volume `100000 m³`. If changing the DEM or excavation parameters, delete
the generated rasters before rerunning preprocessing. No other example
scripts or input files are used.

## 3. Run the two-layer simulation

No separate conservative-state initialization step is needed.
The Fortran solver reads the primitive initial conditions directly from
`IMEX_SfloW2D.inp` and converts them using its **existing** mixture and
conservative-state conversion logic.

```bash
IMEX_TWO_LAYER_VALIDATION=1 ../../src/IMEX_SfloW2D
```

The current Etna example starts with:

| Primitive initial state | L1 | L2 |
|---|---|---|
| Thickness | `DEM/init_volume.asc` | `NONE` (dry) |
| Temperature | 373 K | Unused while dry |
| Horizontal velocity | 0 m/s | 0 m/s |
| Solids by volume | 50% | 0% |

The input uses independent `RESTART_L1_PARAMETERS` and
`RESTART_L2_PARAMETERS` groups. Each layer has its own file list and count:
`RESTART_FILES` belongs only to L1, while `RESTART_FILES_L2` belongs only
to L2. Use an ESRI `.asc` thickness raster, or a sole `'NONE'` entry for
a dry layer. Independent primitives are defined by `INITIAL_L1_PARAMETERS`
and `INITIAL_L2_PARAMETERS`. To start only L2, set L1's `RESTART_FILES='NONE'`
and L2's `RESTART_FILES_L2='./DEM/l2_thickness.asc'`. Both layers can also
start wet. Only same-time initial rasters are supported in the coupled model;
nonzero per-layer release times require future activation logic and are
rejected explicitly. The percentages refer to the *common* solid species.
See [restart input](../../docs/RESTART_INPUT.md) for the complete contract.

These primitive inputs are for **new initial states**, not continuation
from a previous run. Conservative `.q_2d` snapshots and binary restart
files remain available via the existing restart mechanisms.

The example requests `T_END=200 s`, `DT_OUTPUT=5 s`, and `CFL=0.24`.
A shorter exploratory run can set `T_END=10.0` and `DT_OUTPUT=1.0`
but neither time horizon establishes scientific validation. Outputs are
written to `OUTPUT/exampleEtna_two_layer/`.

L1 models the dense granular avalanche with Voellmy friction (`mu=0.3`,
`xi=300`). L2 uses the quadratic rheology with the layer-specific basal
`RHEOLOGY_L2%FRICTION_FACTOR=0.02` and Newton treatment (Lie group 1).
This is the same quadratic-friction coefficient used in `EXAMPLE_TAAL`.

L2 begins strictly dry and can form as L1→L2 elutriation develops. The
`LINEAR_USTAR` elutriation settings are lifted solid volume fraction
`0.15`, efficiency `0.01`, and critical Rouse number `100`
(explicit Lie group 4). Settling and independent substrate mass exchange
in L2 are disabled.

The explicit `ELUTRIATION_SOLID_VOLUME_FRACTION=0.15` keeps the lifted
parcel composition prescribed. Omitting this field instead selects the
current local L1 solid volume fraction after the spatial update, and its
gas/solid co-lift ratio varies from cell to cell. This is an alternative
physical assumption, not a calibration of the Etna case.

Interlayer momentum drag is enabled using the Chézy slip-velocity law,
with dimensionless `CHEZY_DRAG_COEFF=0.03` and a separate analytic Lie
group 5. This coefficient is an exploratory choice used in the
interlayer-drag numerical tests, **not** a calibration for Etna.
The Chézy drag acts only where both layers are thicker than 1 mm; it
exchanges opposite momentum impulses, conserving the pair's total momentum.
Its coefficient is **independent** of L2's quadratic basal
`FRICTION_FACTOR`.

## 4. Optional checks — commands from this directory

Unit checks (no solver or DEM required):

```bash
python3 -m unittest discover -s . -p 'test_*.py'
```

Short controlled integration smoke:

```bash
python3 run_validation.py \
  --solver ../../src/IMEX_SfloW2D \
  --workdir ./validation-smoke
```

Extended conservation/restart/OpenMP checks:

```bash
python3 extended_validation.py \
  --solver ../../src/IMEX_SfloW2D \
  --workdir ./validation-extended
```

The checkers run from this example and create isolated, new or empty
work directories. They copy the already extracted local DEM for preprocessing but **do
not** invoke a preliminary single-layer solver to construct initial
snapshots. The two-layer executable alone constructs its conservative
initial state when each test begins.

The short test checks snapshots, ESRI/NetCDF output, budgets and binary
restart presence. The extended test compares continuous vs restarted
trajectories, 1/2/4 OpenMP threads and L2 formation over a 2 s
trajectory. These are **functional** checks only; TL13b coupled wave
and physical validation is still pending.

## Source files

`prepare.py` generates only DEM rasters. `IMEX_SfloW2D.inp` defines the
two-layer run, including its primitive initial conditions. The local
`input_ellipsoid.py` and `create_input_ellipsoid.py` generate the
original Etna excavation without using another example's scripts.
`run_validation.py`, `extended_validation.py`, and `test_*.py`
provide optional regression evidence.

The `.asc` rasters and `OUTPUT/` are generated and untracked. The
only copied source data is `DEM/Etna2014_crop.zip`.

## Positivity diagnosis of L2

On a rejected spatial update, the solver prints `KP_POSITIVITY` records
after `LIE_SPATIAL_FAILURE` when `N_LAYERS=2`. They report the
cell-owned W/E/S/N transport depths, their four-face mean (equation 3.31
of Kurganov–Petrova), conservative face masses obtained from the same
selected traces as the PCCU mass flux, and the characteristic speed
bounds and CFL values. The diagnostic also reconstructs the eight mass
terms analogous to equation 3.32, distinguishing the four cell-owned
terms from the four incoming neighbor terms. `conservation_offset`
reports the difference between the published cell mass and the mean
of the reconstructed face masses; a nonzero offset means the eight-term
positivity argument does not apply without that extra contribution.

The eight-term decomposition assumes the unmodified central-upwind
mass flux. Dry-donor suppression and a degenerate-speed fallback can
alter it; compare the reported decomposition with the actual
`LIE_SPATIAL_FAILURE` mass update and `LIE_FACE_FAILURE` fluxes.
This is read-only failure evidence, not a positivity fix.

## L2 face-base and reconstruction audit

The failure-only `KP_BASE` diagnostic uses the completed L1 reconstruction as
the effective L2 base. It prints, in W/E/S/N order, each cell-owned base,
the opposite reconstructed base at the same interface, their jump
(`opposite - own`), and the difference between the cell-owned base
average and `B_cent + h1_cell`. These are distinct questions:
a jump between adjacent cells is allowed by the current reconstruction
contract, while the cell-owned four-face average must still reproduce
the cell-center effective base.

The `KP_TRACE` records separate first-pass primitive reconstructions
(`raw_h`, `raw_mass`, `raw_density`) from the selected PCCU transport
traces (`KP_POSITIVITY h_W,E,S,N`, `mass_W,E,S,N`, and
`selected_density`). Both mass sets use the canonical
`primitive_to_conservative` conversion. Compare their differences
from the published conservative cell mass to determine whether a
large reconstruction offset is already present before the L2
hydrostatic/transport selection. The test
`TESTS/TEST_SPATIAL_OPERATOR` also checks the effective base's
cell-owned average identity and independently retained face jumps.

These diagnostics do not prove that the two-layer pressure system
satisfies the Kurganov-Petrova theorem; they test the necessary mass
reconstruction and numerical speed hypotheses without modifying the
numerical method or mass-exchange sources.
