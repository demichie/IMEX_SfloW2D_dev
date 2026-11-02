# Lie sources: cell-outer traversal and end-of-step geometry refresh

## Source-level audit (2026-11-02)

Inspected `src/numerics/lie_source_executor_2d.f90`, `src/numerics/time_integration_2d.f90`, `src/physics/source_terms_2d.f90`, and erosion, deposition and lava deposition kernels.

The exchange group computes `transfers%dB_dt` per cell, and updates per-cell `deposit`, `erosion`, `erodible`, and `lava_deposit`. Currently the executor projects `dB_dt` to vertices, calls `assign_topography_vertices`, and refreshes substrate totals after the mass-exchange group, before subsequent Lie groups. Other groups read `B_prime_x_geom`, `B_prime_y_geom`, curvature and gravity from shared geometry.

### Proposed numerical convention

Freeze the complete geometry (including slope, curvature, gravity, zone/arc data) throughout the spatial update and ordered Lie source sweep for each trial timestep. Evolve local fluid and substrate inventories through the ordered groups; accumulate cell-local `dB_dt` without publishing any geometry change. If all cells and groups succeed, project the final cell-rate once, refresh geometry once, and atomically publish accepted fluid, inventories and geometry. A rejected step must restore all authoritative fields.

**Numerical-semantic change:** groups after mass exchange currently see refreshed geometry; under the proposal they see the old geometry. This is not bitwise equivalent and must be covered by accuracy/regression tests.

### Implementation constraints

1. Preserve group-0 semantics: evaluate rates on immutable pre-spatial `qp_original`, apply to post-spatial conservative state.
2. Preserve ordered per-cell intermediate `q/qp`, admissibility checks and per-cell evolving erodible/deposit/erosion/lava inventories.
3. Do not publish geometry or substrate totals inside the cell loop. Keep the complete source group loop local per cell; merge `dB_dt` only after all groups succeed.
4. Geometry projection and `assign_topography_vertices` can fail, and therefore must execute **before** committing the authoritative timestep. Preserve rollback on failure.
5. Adapt `D20_INJECTED_LATE_LIE_FAILURE` test: its current `had_bed_change` assertion expects geometry changes *before* the late failure. With deferred projection, the fault-injection contract must check trial bed rates or inject failure after the deferred geometric update, without weakening rollback verification.
6. Preserve `LIE_GROUP_PROFILE` observability; distinguish per-group local evaluation from one global geometric commit. Avoid misleading setup/evaluate/commit accounting.
7. Test no-topo-change runs, topographic erosion/deposition, lava, multi-group geometric-source ordering, Newton and analytic paths, failing Newton, injected rollback and 1/2/4/6 thread consistency.
8. Preserve cell-local physics interfaces; never pass grid indices into `src/physics`.

### Status

The cell-major OpenMP Lie sweep and deferred geometry refresh are implemented in `src/numerics/lie_source_executor_2d.f90`. `LIE_CELL_SWEEP_PROFILE` and `LIE_GEOMETRY_COMMIT_PROFILE` now report aggregate wall times; `LIE_GROUP_PROFILE` retains group IDs, methods and call counts without falsely attributing wall time to individual groups. The injected late-step rollback fault now occurs after the deferred geometric refresh. **Compilation and D20 runtime validation are not yet verified; do not treat this as a validated numerical change.**
