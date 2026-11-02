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
