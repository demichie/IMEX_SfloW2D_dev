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
imports/interpolates the emitted DEM. It never uses evolved Fortran states as
reference input. A test-only observational main records the accepted timesteps
and must reproduce the uninstrumented canonical output byte for byte. The
Python time integration uses the same `n_RK=2` explicit tableau (forward Euler
for this source-free problem). Matching the accepted times isolates the spatial
comparison; this test does not assert autonomous adaptive-CFL equality.

`field_contract_v4.json` pins the active core and retains the original
32768-epsilon-per-step h/eta comparison bound. Every cell is compared, including
both momenta using independent initial hydrostatic mass/momentum scales.
`reference_fields.npz` and `field_evidence.json` retain reference fields, checksums,
errors, accepted times and the independently integrated boundary mass transport.

The aggregate checks retain nonnegativity, zero uphill leakage, G=1 lateral
symmetry and the existing mass-budget tolerance `5e-10`. The budget includes
mass crossing the open external boundaries: the new scheme can produce a small
outflow even at the coarse grid's short comparison time. Aggregate reference
values are calculated from the new Python fields, with the existing bounds.
A failed full-field comparison or invariant still fails the test, and all six
completed grid/mode cases retain their diagnostics.

The original HP-PCCU cores, `reference_metrics.csv` and contracts v1-v3 remain
historical archives. Their different numerical scheme and SSPRK2 trajectories
are not acceptance oracles for the current KT operator. Other historical test
suites that explicitly exercise HP-PCCU are not migrated by this Gate-H test.
