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

L1 reproduces the original Etna dense avalanche (Voellmy `mu=0.3`,
`xi=300`). L2 begins strictly dry, and L1→L2 elutriation may create
it dynamically after shear develops. The exploratory `LINEAR_USTAR`
coefficients are lifted solid fraction `0.15`, efficiency `0.01`,
Rouse threshold `100`. L2 settling and Chezy drag are disabled in this
example to isolate formation of the upper current; those processes remain
available in the two-layer implementation.

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
