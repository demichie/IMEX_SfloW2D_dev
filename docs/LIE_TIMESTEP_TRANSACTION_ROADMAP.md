# Lie timestep transaction — current contract

## Scope and ownership

The Lie integrator owns the **trial conservative and primitive fluid arrays**
and publishes them only after the spatial stage and every scheduled source
group have completed successfully. Source groups operate on cell-local
candidates and shared substrate inventories; they do not own independent
full-domain fluid states.

- `time_integration_2d` owns `q_trial`, `qp_trial` and the final fluid
  publication. An invalid spatial candidate or rejected Lie source causes
  a **fatal stop**. There is no timestep retry.
- `lie_source_executor_2d` orders explicit, Newton and analytic sources
  through the canonical source registry. It reports a failed group to
  the integrator; it does not promise whole-step substrate recovery.
- Group zero evaluates rates on the pre-spatial state `q^n` and applies
  accepted increments to the post-spatial candidate `q^*`.
- Shared substrate inventories may be updated during the cell sweep.
  There is **no recoverable global substrate snapshot** in production.
- On a rejected substrate candidate (D19, error 7), the executor prints a
  failure-only `D19_SUBSTRATE` witness for the first failed worklist cell:
  coordinates, Lie group, timestep, per-species before/candidate inventories,
  deposition/erosion rates and lava inventory. Its OpenMP capture is
  synchronized and does not alter the rejected source transaction.
- The geometry transaction defers bed-cache publication. Its internal
  rollback is a geometry-local mechanism, not a recoverable timestep.
- A two-layer trial uses the same Lie integration path, layer-specific
  physical configurations, a shared active workset and one fluid commit.
  The production model-layout gate still rejects `N_LAYERS=2`.

## Active invariants

1. Reject inadmissible spatial conservative states rather than clipping them.
2. Preserve the conservative limits, group order and dry-state behavior of
   every explicit, implicit and analytic physical source.
3. Do not make trial `q/qp` visible as authoritative fluid state before
   acceptance of the complete step.
4. Keep OpenMP cell-local writes disjoint and deterministic.
5. Do not describe a failed timestep as rolled back unless a future,
   explicitly reviewed design adds a complete substrate/geometry/OU
   transaction and verifies it end-to-end.
6. Preserve the single-layer numerical path and its evaluation order.

## Lie/OpenMP workload balance

The accepted source-group transaction is cell-major: each worklist cell
owns its local trial state, and the physical groups execute sequentially
inside that cell. Worklist indices are spatially ordered, but the physics
cost is heterogeneous (dry cells, settling, erosion and rheology).
The Lie cell loop therefore uses `SCHEDULE(DYNAMIC,16)` to distribute small,
contiguous chunks while preserving each cell's physical evaluation order.
No cell is processed more than once. Failure selection still uses the
smallest worklist index, and the existing synchronized substrate witness
remains deterministic. Per-worker scratch continues to be indexed by
the executing OpenMP worker.

The opt-in `IMEX_PROFILE_LIE=1` diagnostic reports
`LIE_SOURCE_DETAIL cell_balance_efficiency`, computed as the sum of all
instrumented per-cell worker times divided by the sum of
`team_size * max(worker_cell_time)` over timesteps. It is an approximate
cell-compute balance indicator, not a measurement of barrier time or
whole-simulation parallel efficiency: early-skipped cells, team creation,
setup and final synchronization are not included.

The change is motivated by the LAHAR `T_END=1800` baseline on macOS ARM64
(72.796 s / 48.395 s / 36.729 s wall time at 1/2/4 threads). These are
single-run times with the previous scheduling, not a verified prediction
for the new schedule. Repeat that exact benchmark and compare 1/2/4
thread reproducibility before claiming a speedup.

## Current verification

- `TEST_SPATIAL_OPERATOR` exercises the coupled Lie trial, wet/dry
  configurations, external erosion/deposition, interlayer transfers,
  drag, optional L1/L2 fields and OpenMP 1/4-thread comparisons.
- `TEST_MASS_EXCHANGE/test_projection_stress.f90` checks the direct Lie
  executor's rejected-source status without asserting recoverable
  substrate rollback.
- `TEST_WORKSPACE_LIFECYCLE` separately checks geometry-transaction
  restoration and canonical granular/lava substrate aggregates.
- Binary restart and NetCDF round-trip checks are in the spatial
  regression. These do not substitute for a production end-to-end run.
- Every new commit requires its own optimized/debug build and numerical
  regression results. Earlier successful workflows cannot certify
  later changes.

## TL13a/TL13b gate

Refer to [two-layer closure audit](TWO_LAYER_CODE_AUDIT.md) and
[two-layer architecture roadmap](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md).
The two-layer production gate remains closed pending an integrated
production-driver acceptance and subsequent scientific validation.

This document describes the **current** failure contract only. It does
not require retention of superseded transaction APIs or prior binary
restart formats.
