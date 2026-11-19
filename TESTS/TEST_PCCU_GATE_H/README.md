# Gate-H comparisons with the current KT-PYTHON operator

Run from the repository root with an executable built locally:

```sh
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE KEEP_TEST_WORKDIR=1 \
  TESTS/TEST_PCCU_GATE_H/run_test.sh \
  "$PWD/build-local/optimized/src/IMEX_SfloW2D"
```

The test uses three grids (`dx = dy = 0.4, 0.2, 0.1 m`) and both G=1 and
slope-plus-curvature modes. All runs use `N_RK=2` and
`KT_TRACE_RECONSTRUCTION='froude_blend'`.

The independent NumPy implementation lives in
`../TEST_SPATIAL_OPERATOR/reference/kt_python_2d_core.py`. It implements the
current Fortran operator: independently limited velocities, the four selectable
depth reconstructions, KT central-upwind transport, oriented hydrostatic paths,
dry-donor guards, copied boundary traces, outflow-only external fluxes and
curvature sources. It supports isothermal liquid with an optional solid species;
Gate-H uses pure liquid without rheology or sources. One-cell dimensions use
the same code as 2D with their inactive direction omitted.

The reference starts from the exact emitted ASCII initial state and independently
imports/interpolates the DEM. A test-only observational main records every
accepted timestep and its full-precision conservative state. Its canonical
output must remain byte-identical to the uninstrumented production run.

**Acceptance:** at *every* accepted timestep, KT-PYTHON advances once from
the **identical preceding Fortran conservative state**, with the actual accepted
`dt`. The maximum error over every grid cell and all three conservative
equations must remain below the fixed `32768 * float64_epsilon` normalized
per-step bound. The mass scale is the initial maximum conservative mass
(at least one), and both momentum scales use initial hydrostatic units.
This is a local operator check; its tolerance does not grow with timestep
count. Any missing or modified state stream, nonfinite field, or failed step
makes Gate-H fail.

**Independent trajectory:** KT-PYTHON also integrates from the common initial
state over the entire accepted timestep sequence. `reference_fields.npz` and
`field_evidence.json` retain its h, eta and all conservative momenta together
with the *original* stringent roundoff limits and their pass/fail indicators.
Independent trajectory drift is **diagnostic**, not an acceptance gate:
roundoff may amplify across many nonlinear steps even when the per-step
operator is correct. The mandatory `trajectory_evidence.json` records both
the cumulative drift and the authoritative `same_input_all_steps` result.

The separate aggregate checks retain nonnegativity, zero uphill leakage,
G=1 lateral symmetry and the original mass-budget tolerance `5e-10`,
including mass transported across open external boundaries. These remain
acceptance gates for all six grid/mode cases.

The original HP-PCCU cores, `reference_metrics.csv` and contracts v1-v3 remain
historical archives. Their different numerical scheme and SSPRK2 trajectories
are not acceptance oracles for the current KT operator. Other historical test
suites that explicitly exercise HP-PCCU are not migrated by this Gate-H test.

## Mandatory full-precision same-input audit

`diagnose_trajectory.py` is executed for **every** case, regardless of the
independent trajectory drift. It reads and verifies `accepted_states.bin` and
`accepted_times.bin`, computes one KT-PYTHON update for each accepted step
from the preceding Fortran state, checks all cells and conservative components,
and writes `trajectory_evidence.json`. A failed step yields a nonzero exit
status and fails the full Gate-H run. Additional independent trajectory
integration records the first drift exceeding the original roundoff threshold,
but does not decide acceptance. Optional short-prefix replay remains available
for investigating a failed mandatory operator check.

## Conservation of nearly dry cell momentum

The reference now preserves the accepted conservative state in positive-mass,
near-dry cells rather than zeroing their momenta at the reconstruction-depth
threshold. Fortran `update_spatial_cell` keeps conservative momenta unchanged
apart from the physical spatial update; regularization of dry velocities is
performed in primitive and interface reconstruction. The previous Python
post-step projections were not part of the numerical model, and the new
`test_kt_reference.py` case guards against their reintroduction. The Gate-H
comparison limits and all solver physics are unchanged.
