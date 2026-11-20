# Etna two-layer: dense avalanche generating an upper dilute current

This is an **experimental validation example**, not an approved two-layer
production scenario. The solver rejects `N_LAYERS=2` by default; a controlled
run requires `IMEX_TWO_LAYER_VALIDATION=1`. Keep the production gate closed
until functional and scientific acceptance.

## Scientific initial condition

The lower layer **L1** is exactly the initial dense pyroclastic avalanche from
[`EXAMPLE_ETNA`](../EXAMPLE_ETNA/README.md): the same original Etna DEM,
ellipsoidal excavation, initial 100000 m³ target, 50% solid volume concentration,
373 K initial temperature, and Voellmy coefficients `mu=0.3`,
`xi=300`. The original preprocessing script and the single-layer solver
generate the conservative state: this example does **not** implement a second
interpolation, mixture-density, or energy-conversion algorithm.

The upper layer **L2 is entirely dry at t=0**. Every conserved component,
including carrier mass, solids, momentum and energy, is identically zero in
every cell. There is no manually prescribed dilute cloud, artificial initial
depth, atmospheric release or separate L2 initial volume. The only configured
mechanism that can initially generate upper-layer mass is **L1→L2
elutriation**, with physically required carrier co-lift. The elutriation
implementation supports a dry upper-layer cold start.

The first diagnostic configuration disables L2→L1 settling and Chezy
interlayer drag to isolate the development of the dilute layer. L2 is given
`rheology_model='none'` as a provisional no-basal-resistance choice.
The configured `LINEAR_USTAR` law uses 15% solids in uplifted parcels,
`elutriation_efficiency=0.01`, and `rouse_critical=100`. **These are
exploratory test coefficients, not a calibrated Etna eruption model**, and
the intended L2 closure and physically defensible parameter ranges still
require scientific validation. The original Etna L1 rheology, transport
properties and initial state are not changed.

## Prepare input data now

Requirements: Python 3, NumPy, pandas, and a compiled single-layer
`bin/IMEX_SfloW2D`. From the repository root:

```bash
python3 EXAMPLES/EXAMPLE_ETNA_TWO_LAYER/prepare.py
```

The script:

1. Extracts `Etna2014_crop.asc` from the **single** tracked
   `EXAMPLES/EXAMPLE_ETNA/DEM/Etna2014_crop.zip`. No DEM is committed
   a second time.
2. Invokes the **unchanged** `EXAMPLE_ETNA/create_input_ellipsoid.py`
   with its original `input_ellipsoid.py` parameters. The results
   `DEM/new_dem.asc` and `DEM/init_volume.asc` reproduce the same
   physical excavation.
3. Runs an isolated **single-layer** Etna initialization with `T_END=0`
   and conservative t=0 output enabled. No timestep of the dense avalanche
   is advanced; all original physical inputs remain unchanged.
4. Copies its initial conservative snapshot **byte for byte** to
   `INITIAL_CONDITIONS/etna_L1_0000.q_2d`. It produces
   `INITIAL_CONDITIONS/etna_L2_0000.q_2d` with identical cell coordinates
   and strictly zero conservative components. A format/cell-count/mass
   check prevents partial or accidentally empty initial conditions.

The temporary single-layer initialization run is deleted after extraction.
The generated DEM, snapshots, and eventual solver outputs are not tracked.

To prepare a separate directory or use a custom compiled binary:

```bash
python3 EXAMPLES/EXAMPLE_ETNA_TWO_LAYER/prepare.py \
  --directory /absolute/path/to/etna-two-layer \
  --solver /absolute/path/to/IMEX_SfloW2D
```

When using an isolated directory, copy the tracked
`EXAMPLE_ETNA_TWO_LAYER/IMEX_SfloW2D.inp` there before launching the
eventual two-layer solver. The preparator deliberately never overwrites
modified user input files or existing L1/L2 initial snapshots. If the
original Etna ellipsoid configuration changes, remove/recreate both
generated rasters and both initial snapshots before regenerating.

Run the lightweight preparation regression independently:

```bash
python3 -m unittest discover -s EXAMPLES/EXAMPLE_ETNA_TWO_LAYER \
  -p 'test_*.py'
```

This test requires only Python's standard library, uses synthetic small
snapshots and does not need the solver or DEM.

## Controlled two-layer execution — experimental TL13a

The production parser still rejects `N_LAYERS=2` by default. For a
controlled integration run only, set `IMEX_TWO_LAYER_VALIDATION=1`.
This opt-in is not an endorsement of the coupled scientific model.

To run the short prepared Etna end-to-end smoke test from the repository
root (requires compiled solver, Python 3, NumPy, pandas and the DEM archive):

```bash
python3 EXAMPLES/EXAMPLE_ETNA_TWO_LAYER/run_validation.py \
  --solver src/IMEX_SfloW2D \
  --workdir /tmp/imex-etna-tl13a
```

Use an empty directory for `--workdir` and inspect its solver outputs.
The test prepares the canonical single-layer Etna initial state and strictly
dry upper layer, runs to `t=0.001 s`, then checks snapshot time indices,
L1/L2 conservative files, upper-layer ESRI rasters, total mass budget
arithmetic, NetCDF presence and binary restart presence. It does **not**
yet validate restart continuation, physical conservation accuracy, L2
growth over longer times or OpenMP bitwise identity.

The `solver-validation.yml` GitHub Action exposes an optional
`two_layer_smoke` dispatch input. Set it to `true` to run this
experimental job; ordinary CI jobs continue without experimental activation.

## Two-layer reference input and future scientific validation

The checked-in `IMEX_SfloW2D.inp` requests `N_LAYERS=2`,
explicit `RHEOLOGY_L1/L2_PARAMETERS`,
`TEMPERATURE_L1/L2_PARAMETERS`,
`MASS_EXCHANGE_L1/L2_PARAMETERS`, each layer's four boundary NAMELISTs,
and one canonical pair-wide `INTERLAYER_PARAMETERS`.
It points to the two generated `.q_2d` files in `INITIAL_CONDITIONS`.
The model's common topography is the generated
`DEM/new_dem.asc`, as in the original dense avalanche.

**Do not remove the production safety gate.** In a disposable prepared
validation directory only, the controlled execution command is:

```bash
IMEX_TWO_LAYER_VALIDATION=1 ../../bin/IMEX_SfloW2D
```

Before functional or scientific acceptance, check: exact L2 dryness at t=0; onset of
elutriation only after flow motion; nonnegative L1/L2 mass and component
inventories; solid/carrier conservation across internal transfer; no
spurious energy/momentum production; separate L1/L2 outputs and shared
substrate budgets; restart continuation; and OpenMP reproducibility.
A passing input parser or a created L2 field alone is not scientific
validation.

## Extended TL13a functional regression

The optional `two_layer_extended=true` workflow input runs the separate
`extended_validation.py` checker. The script reuses the **same** canonical
Etna input, initial L1/L2 snapshots, production solver and snapshot/budget
reader. It runs a continuous trajectory with outputs at 0, 0.005 and 0.01 s,
restarts an otherwise identical trajectory from the midpoint binary checkpoint,
and compares the final conservative state and checkpoint. Independent runs
with 2 and 4 OpenMP threads are compared to the 1-thread reference.

Budget assertions cover L1/L2 mass and solids additivity, nonnegativity,
total mass and total solids conservation (relative tolerance 1e-6), and
strictly positive L2 mass growth. The cross-run conservative and budget
comparisons use 5e-10 relative and 1e-8 absolute tolerances. These are
**regression criteria**, not experimentally established error bounds.

Run locally (requires built solver):

```bash
python3 EXAMPLES/EXAMPLE_ETNA_TWO_LAYER/extended_validation.py \
  --solver src/IMEX_SfloW2D \
  --workdir /tmp/imex-tl13a-extended
```

The run preserves its evidence directory. Longer time and alternative
physics configurations, shear stability, convergence and quantitative
scientific validation remain separate work under TL13b.
