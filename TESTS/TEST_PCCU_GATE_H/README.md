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

## First-step failure localization

When a full-field comparison fails, `run_test.sh` automatically invokes
`diagnose_first_step.py`. It prepares a separate `first-step-diagnostic/`
run with exactly the original first accepted timestep as end time and
output event, and uses the already compiled observational main. It
**verifies that this replay accepted one timestep with matching `dt`**
before comparing the first conservative snapshot with a *single* independent
KT-PYTHON update. Thus no Fortran numerical kernel is modified.

The script retains `first_step_evidence.json` alongside the existing
`field_evidence.json`; it prints the worst mass cell, its three conservative
variables, the inferred first-step Fortran rate `(Q1-Q0)/dt`, the direct
Python spatial residual and a surrounding 3-by-3 error stencil. This is a
diagnostic, **not** a relaxation of `field_contract_v4.json`. If the
replay cannot be verified, the diagnostic fails independently and the
original Gate-H failure remains in force.

For a persistent local case, set `KEEP_TEST_WORKDIR=1` as above; otherwise
the temporary directory is removed by the runner.

The same observer also runs an exact **two-step replay** (`--steps 2`),
checking both accepted `dt` values against the full run. The
`prefix_2_step_evidence.json` report captures the earliest post-initial
discrepancies: if the first step matches but the second does not, inspect
state conversions and the reconstructed second-step traces. For multi-step
prefixes the reported inferred rate is a *mean rate*, not an instantaneous
spatial residual.

## Full-precision accepted-state trajectory audit

The Gate-H observer now additionally emits `accepted_states.bin`: one
Fortran stream-unformatted record of the three canonical L1 conservative
components over the complete grid, immediately after each accepted advance.
The original canonical snapshot is still checked byte-for-byte to establish
that recording was non-invasive. The SHA-256 of the state stream is kept
in `accepted_time_evidence.json`.

Upon a failed full-field comparison, `diagnose_trajectory.py` reads those
full-precision states and independently advances KT-PYTHON through the **same**
accepted `dt` sequence, comparing each step's conservative mass and momenta.
It records the first step exceeding the existing roundoff-scaled Gate-H
threshold, plus the maximum error reached and the worst cell, in
`trajectory_evidence.json`. This avoids drawing conclusions from the
12-decimal-digit ASCII snapshots. It does not modify any numerical object,
physical formula or acceptance limit.
