# Lie time-step transaction: migration roadmap

## Objective

The **time-step**, not an individual Lie source group, is the transaction.
A failure in any group invalidates the entire time-step. Preserve the
pre-spatial source evaluation for group 0 and maintain the Lie order.
Do not clip or limit the spatial conservative update.

## Rules

- The integrator owns the pre-step fluid snapshot, trial fluid state,
  substrate snapshot and geometry transaction.
- Source groups apply changes directly to the integrator-owned trial
  arrays. Distinct OpenMP workers update distinct cells. No per-group
  full-domain `q_group` or `qp_group` copies.
- The trial state is not visible in the published `state%q` or
  `state%qp` until the full step is accepted.
- A cell-local failure is collected outside the parallel region and
  returned to the integrator. The integrator restores the pre-step
  substrate and geometry, leaving published fluid untouched.
- Group 0 evaluates its source rates from `q_before` and
  `qp_before`, but applies increments to the post-spatial trial.
- Source-rate physical limitations remain in place. Spatial admissibility
  is checked without limiting the spatial rates.

## Milestones

- [x] Inspect the full time-step rollback and group-local candidate flow.
- [x] Eliminate full-domain per-group source-state copies; update trial
      cells directly using existing OpenMP partitioning.
- [x] Retain the time-step trial and pre-spatial snapshots at the integrator
      level and retain the existing rollback failpoint.
- [x] Remove obsolete per-group copy profiling counters; retain setup,
      evaluation and completion measurements.
- [x] Keep analytic and Newton results in cell-local candidates until validation.
- [x] Remove per-group inventory staging; update per-cell inventories in
      the time-step trial, retaining bed-rate buffers for projection.
- [ ] Run the Lie rollback/injected-failure regression Actions and PCCU tests.
- [ ] Re-run the OpenMP performance benchmark at 1, 2 and 4 threads;
      compare performance against commit 0e3a8a7.
- [ ] Investigate reducing remaining full-domain *time-step* copies only
      if profiling shows a benefit, preserving a single rollback boundary.

## Validation checklist

- Failure in an early group leaves published conservative and primitive
  state unchanged and restores substrate inventories and geometry.
- Failure after a topography-changing group also rolls back geometry.
- Successful updates match the previous implementation within the
  existing numerical regression tolerances.
- No OpenMP worker reads state from neighboring cells being updated
  by other workers. The active-cell index list must be unique.
- Rejected groups do not continue to subsequent groups.

## Notes

The time-step trial arrays are retained intentionally: this permits
in-place per-cell updates while keeping the published fluid state
unchanged until the entire time-step succeeds. Source-group-local
full-domain snapshots are unnecessary with this structure.

## Second implementation pass

Analytic and Newton operations now validate private cell candidates before
publishing each cell. Explicit mass-exchange updates the substrate per cell
inside OpenMP; a later failure rolls back the full time-step. Bed geometry
projection and its validation remain outside the parallel region. The
existing injected late-failure regression must still be run.

## Validation route identified

The existing `.github/workflows/lie-openmp-validation.yml` already contains
an end-to-end production-solver job, running
`TESTS/TEST_LIE_OPENMP_END_TO_END/compare_threads.py` against four fixtures:
standard, pore-pressure, Newton, and topography. It compares complete NetCDF
output at 1, 2, and 4 threads. The topography fixture also injects a late
failure after group 4 and checks `D20_FULL_STEP_ROLLBACK_PASS` for all three
thread counts. The strict-debug build and workspace-lifecycle regression
are separate jobs of the same workflow.

**Current status:** code committed, review performed, but compilation and
end-to-end Action execution have NOT yet been observed on this revision.
Run the validation workflow first, inspect all job results, and only then
run the performance benchmark. Do not label the migration verified on the
strength of a source-level review alone.

## Regression repair after first test run

The direct executor test `TEST_MASS_EXCHANGE/test_projection_stress.f90`
reported `Lie inventory rollback mismatch`. The old test invoked the
executor and rolled back geometry but still expected inventories to be
unchanged, which was the superseded *group-level* transaction contract.
The revised test now takes the same `substrate_snapshot_type` snapshot as
the production integrator and explicitly restores it after a deliberately
rejected group. It continues checking the fluid trial, geometry and all
inventory arrays. This corrects the test harness rather than restoring
unwanted per-group staging in the production executor.

Validation of this correction requires a new regression run; no pass has
been asserted yet.
