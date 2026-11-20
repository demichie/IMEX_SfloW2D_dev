# Etna two-layer — standalone experimental example

This example runs from **its own directory**. It contains its own
topography-preparation scripts, preprocessing parameters, single-layer seed
input and two-layer input. **The only file to copy from another example is
the shared Etna DEM ZIP.** No scripts are executed from sibling examples.

> **Experimental only:** `N_LAYERS=2` requires
> `IMEX_TWO_LAYER_VALIDATION=1`. Functional regression tests passed for a
> previous solver revision; scientific TL13b approval is still pending.
> The default 200 s case has not been scientifically validated.

## 1. Open the example directory

From the repository root, enter this directory **once**:

```bash
cd EXAMPLES/EXAMPLE_ETNA_TWO_LAYER
```

**All subsequent commands below run here**, not from the repository root.
Requirements: Python 3, NumPy, pandas and a compiled
`../../src/IMEX_SfloW2D` executable.

## 2. Copy the DEM ZIP (once)

```bash
cp ../EXAMPLE_ETNA/DEM/Etna2014_crop.zip DEM/
```

This is the only dependency on `EXAMPLE_ETNA`. The ZIP is not tracked a
second time in Git. Do **not** manually copy the preprocessing scripts or
the original example input.

## 3. Prepare only the DEM (once)

```bash
python3 prepare.py
```

**This command never launches the Fortran solver.** It only extracts the
topography ZIP and runs the local Python ellipsoid preprocessor.

`prepare.py` uses only files in this example:
- `DEM/Etna2014_crop.zip`: original topography, extracted locally.
- `input_ellipsoid.py` and `create_input_ellipsoid.py`: the Etna ellipsoidal
  excavation, preserving the original `x1=500430`, `y1=4177500`, and
  `vol=100000` parameters. Outputs are `DEM/new_dem.asc`,
  `DEM/init_volume.asc` and `DEM/ellipsoid.asc`.
- `single_layer_seed.inp`: the original Etna **single-layer** configuration,
  used only by the explicit initialization step below.

## 4. Generate the conservative initial states (explicit, once)

```bash
python3 initialize_layers.py
```

**This is a separate user-invoked command.** It runs the compiled Fortran
solver once at `T_END=0` to calculate the canonical mass, momentum,
energy and composition of L1 from the generated DEM rasters. There is
**no time advancement** and **no two-layer simulation**, but the executable
is launched for initialization. L2 is constructed with exactly zero
conservative components.

This avoids a second Python implementation of the Fortran density and
energy formulas. You can skip this command if the two initial `.q_2d`
files already exist and are valid.

The initial-state generator writes:
- `INITIAL_CONDITIONS/etna_L1_0000.q_2d`: original dense avalanche,
  with an initial 100000 m³ excavation target, 50% solid volume fraction,
  373 K, Voellmy `mu=0.3`, `xi=300`.
- `INITIAL_CONDITIONS/etna_L2_0000.q_2d`: exactly zero conservative mass,
  particles, momentum and energy, with matching grid coordinates.

The initialization command will **not overwrite** existing snapshots.
If you change the DEM or excavation setup, remove the generated DEM rasters
and both snapshots explicitly before regenerating them.

## 5. Run the two-layer simulation

The tracked `IMEX_SfloW2D.inp` is the two-layer run configuration.
From the same directory:

```bash
IMEX_TWO_LAYER_VALIDATION=1 ../../src/IMEX_SfloW2D
```

The current input specifies `N_LAYERS=2`, `T_END=200 s`,
`DT_OUTPUT=5 s` and `CFL=0.24`. For a first exploratory run, consider
temporarily setting `T_END=10.0` and `DT_OUTPUT=1.0`. The 10 s horizon is
**not** an accepted scientific regression.

Outputs are written under `OUTPUT/exampleEtna_two_layer/`, including
L1/L2 conservative snapshots, ESRI fields, NetCDF, budget CSVs and binary
restart files.

**Physical setup:** L1 initially contains the same dense Etna avalanche
as the original example; L2 starts entirely dry. L1→L2 elutriation can form
the upper layer after shear develops; it carries both solids and carrier.
The preliminary run disables L2→L1 settling and Chezy drag to isolate
elutriation. L2 uses `rheology_model='none'`. The `LINEAR_USTAR` closure
uses lifted solid fraction 0.15, efficiency 0.01 and Rouse threshold 100.
These are exploratory coefficients, not calibrated predictions.

## 6. Optional local validation

The validation scripts are **already inside this directory**, and create
separate evidence directories without changing `IMEX_SfloW2D.inp` or the
prepared manual example.

Run lightweight Python unit tests (no solver or DEM required):

```bash
python3 -m unittest discover -s . -p 'test_*.py'
```

Run a controlled short smoke test (t=0.001 s):

```bash
python3 run_validation.py \
  --solver ../../src/IMEX_SfloW2D \
  --workdir ./validation-smoke
```

Run the more extensive conservation, restart and 1/2/4-thread checks:

```bash
python3 extended_validation.py \
  --solver ../../src/IMEX_SfloW2D \
  --workdir ./validation-extended
```

Both `--workdir` paths must be absent or empty. The scripts copy the
**local DEM ZIP only** into their isolated working directories. Their
preparation and solver runs use the local scripts and seed input.

The short smoke test checks both layer snapshots, ESRI/NetCDF output,
budget arithmetic and restart presence. The extended checker compares
continuous and resumed evolution at 0, 0.005 and 0.01 s, OpenMP 1/2/4
threads, component conservation, and a separate 2 s L2 formation
trajectory. The upper layer need not form during the first 0.01 s because
the avalanche initially rests. These are functional acceptance criteria,
not scientific TL13b validation.

## Files and ownership

| File | Purpose |
|---|---|
| `IMEX_SfloW2D.inp` | Editable two-layer simulation input |
| `single_layer_seed.inp` | Local original Etna single-layer t=0 seed template |
| `input_ellipsoid.py` | Local topography/excavation parameters |
| `create_input_ellipsoid.py` | Local terrain/initial-volume generator |
| `prepare.py` | DEM extraction and ellipsoidal raster preprocessing only |
| `initialize_layers.py` | Explicit, optional Fortran t=0 initialization of L1/L2 snapshots |
| `run_validation.py` | Short controlled solver regression |
| `extended_validation.py` | Restart, budgets and OpenMP regression |
| `test_*.py` | Script unit tests |
| `DEM/Etna2014_crop.zip` | Only file copied from `EXAMPLE_ETNA`, untracked |

The example's `DEM/*.asc`, `INITIAL_CONDITIONS/*.q_2d` and `OUTPUT/`
are generated artifacts and are not tracked. The original single-layer
example remains unchanged.

**Scientific limitation:** the code implements a coupled two-layer
experimental path, but coupled hyperbolicity/shear, longer evolution,
constitutive calibration and physically justified lift-off still require
TL13b analysis.
