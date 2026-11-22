# Active-domain OpenMP roadmap

**Status:** approved design, **not yet implemented in the production solver**.
**Branch:** `feature/two-layer`.
**Goal:** improve active-domain management with minimal changes and unchanged numerics.

This document covers only domain worksets. The
[two-layer roadmap](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md) remains authoritative
for physics; the [Lie transaction roadmap](LIE_TIMESTEP_TRANSACTION_ROADMAP.md)
governs accepted-state publication and fatal failure handling.

## 1. Agreed strategy

There is **one shared** `domain%solve_mask(j,k)` for all layers and
**one global ordered solve list**, stored in the existing
`domain%j_cent(:)`, `domain%k_cent(:)`, and `domain%solve_cells`.
A cell is wet when the accepted conservative total mass is positive
in **at least one** active layer. The final solve mask additionally
contains required source cells, boundaries and stencil halos.

The changes relative to the current implementation are limited to:

1. **Lie wet/dry update:** evaluate wetness in the existing final OpenMP
   Lie publication loop, only for old solve-list cells, and update
   `solve_mask` there. Today `check_solve` scans the full state.
2. **Parallel halo:** retain both existing snapshot-based, four-neighbor
   `.OR.` dilations, but replace the serial `!$OMP MASTER` work with
   an OpenMP loop that owns disjoint destination cells.
3. **Parallel global compaction:** replace the serial full-grid scan that
   builds `j_cent/k_cent` by OpenMP row counts, a small prefix sum and
   parallel output to the existing global list.

Keep the existing `solve_mask_temp(:,:)` as the read-only dilation
snapshot. Do **not** add per-layer masks, persistent linear-index copies,
thread-local cell lists, sparse merge, duplicate-removal passes or tiles.

The [standalone compaction tests](../TESTS/TEST_DOMAIN_COMPACTION/README.md)
compare serial scanning, OpenMP global scanning and sparse merging.
**Global OpenMP scanning is selected** for simplicity and generally
favorable measured performance, subject to end-to-end solver profiling.

## 2. Timestep algorithm

### A. Final accepted Lie loop: update the mask

Keep the existing `!$OMP PARALLEL DO` over the **old** global
`j_cent/k_cent` list. After all Lie source groups have succeeded:

1. Publish the accepted `state%q` and `state%qp` for the cell.
2. Check final conservative total mass (`idx_totMassEqn`), using
   an OR across the active layers (L1 only in single-layer runs).
3. Set the corresponding **single shared** `solve_mask(j,k)` to the
   resulting wet/dry flag. Do not expand the halo or build lists here.

Each worklist cell is visited exactly once, so updates are
race-free without atomics. Preserve source-group ordering, the accepted
transaction, geometry handling, and failure behavior.

**Important:** after Lie, `solve_mask` is provisional: it describes
wetness of the old worklist, not the new finalized solve domain.
The old list must remain valid until `check_solve` finishes.
Audit stochastic, diagnostics, output and restart accesses between
Lie and `check_solve`; no consumer may incorrectly rely on the
provisional mask as if it were the previous finalized workset.

### B. Next `check_solve`: complete the seed mask

Keep `check_solve` responsible for all non-Lie activations and their
current rules at the **new** simulation time:

- Enforce `solve_mask_time(j,k) <= t` as currently required.
- Restore the outer boundary membership maintained by the current code.
- Insert bottom radial/fissural seeds from `cell_source_fractions > 0`
  and radial seeds from `source_cell == 2`, when enabled.
- Recover wet cells newly eligible through delayed releases or changed
  external sources, even **outside the old solve list**.
  Use existing source/timing information and read their committed
  masses where needed; do not lose such cells by trusting the Lie
  worklist alone.
- Keep initialization, restart and `solve_all=.TRUE.` paths correct
  even when there was no preceding accepted Lie step.

Ordinary cells processed by Lie reuse its logical wetness result;
`check_solve` should not scan their real-valued state globally again.
A full-grid physical-state scan may still be required on exceptional
release/restart paths unless complete event enumeration is proven.
Do not duplicate source configuration or invent persistent activation
structures just to remove this scan.

The result is a **seed mask**. Do not interpret halo membership as
wetness after the next phase.

### C. Two OpenMP halo expansions, reusing `solve_mask_temp`

Keep `active_halo_layers=2` and the same four orthogonal neighbors,
interior-only destinations and boundary handling as today.

For each of the two expansion passes:

1. Copy `solve_mask` to the existing `solve_mask_temp`.
2. Parallelize the destination grid rows with OpenMP.
3. Read **only** `solve_mask_temp` and write to the thread-owned
   destination in `solve_mask`.
4. Synchronize before the next snapshot and expansion.

For interior destinations, the operation is:

```fortran
solve_mask(j,k) = solve_mask_temp(j,k) .OR. &
                  solve_mask_temp(j-1,k) .OR. &
                  solve_mask_temp(j+1,k) .OR. &
                  solve_mask_temp(j,k-1) .OR. &
                  solve_mask_temp(j,k+1)
```

The snapshot prevents propagation chains within one pass. Each
destination has a unique writer: **no atomic operations or thread
masks are required**. The result preserves the current
Manhattan-distance-two solution halo.

After both passes, apply the existing radial exclusion
(`source_cell == 1`) before building the solve worklist.
Do not reduce the halo width or alter boundary/source ordering.

### D. OpenMP global solve-list compaction

The final mask is now authoritative. Work only on this **logical**
array; do not read `q` or `qp` during compaction.

1. **OpenMP count:** `row_count(k)=COUNT(solve_mask(:,k))`
   independently for each `k=1..Ny`.
2. **Prefix sum:** `row_offset(1)=0`, followed by
   `row_offset(k+1)=row_offset(k)+row_count(k)`.
3. **OpenMP emit:** loop over rows, then increasing `j=1..Nx`,
   writing active cells directly to `j_cent(l)` and `k_cent(l)`,
   starting after `row_offset(k)`.
4. Set `solve_cells=row_offset(Ny+1)`.

Only temporary `O(Ny)` counts and offsets are needed. Each row
writes to a disjoint output segment: no atomics, per-thread lists,
duplicate removal, sorting or merging. Indices are ordered by
increasing `k`, then `j`, independent of scheduling/thread count.
Verify `solve_cells == COUNT(solve_mask)` and exact index order.

### E. Keep dependent worksets consistent

Preserve current face masks/lists, HP and reconstruction worksets,
including their membership and ordering. The reconstruction halo
extends **two additional read-only cells** around solve cells;
it does not enlarge the evolved solve domain.

Do not blindly move cell-driven updates of `solve_mask_x/y` into
the parallel compaction loop: adjacent cells can write the same
face location. Retain the current safe path initially, or
separately validate a face-owned parallel pass.

## 3. Implementation and verification

1. **Audit:** `check_solve` callers, wet/dry uses of `solve_mask`,
   time-dependent sources, deferred releases, restart, initialization,
   boundary rules, and stochastic/diagnostic consumers.
2. **Modify Lie:** add only the shared wet/dry mask assignment to
   accepted-state publication; check single-/two-layer behavior and
   safe timing of provisional membership.
3. **Parallelize the existing halo:** reuse `solve_mask_temp`
   without changing the two passes or stencil geometry.
4. **Parallelize list construction:** implement count / prefix sum /
   emit into the existing `j_cent/k_cent` arrays.
5. **Compare and profile:** check exact masks, ordered lists,
   face/HP/reconstruction worksets, boundaries, source exclusions,
   dry/wet transitions, delayed releases, empty/full domains and
   restarts against the old implementation. Run strict/debug and
   optimized builds, OpenMP 1/2/4/8, representative LAHAR cases,
   physical invariants and relevant experimental two-layer tests.

The standalone
[compaction benchmark](../TESTS/TEST_DOMAIN_COMPACTION/README.md)
validates only isolated list construction, **not** Lie mask updates,
activation completeness, halo equivalence or full solver behavior.

**Completion criterion:** one shared mask and the existing ordered
worklist remain consistent across all paths, with unchanged
numerical results and verified OpenMP behavior.

**Out of scope:** modifying physics, changing the two-ring solve
halo, per-layer masks, MPI/GPU/tiles, sparse production backends
and production activation of experimental two-layer physics.
