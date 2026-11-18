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
