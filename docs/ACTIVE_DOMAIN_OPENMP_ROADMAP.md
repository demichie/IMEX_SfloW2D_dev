# Active-domain OpenMP roadmap

**Status:** agreed algorithm, **not yet implemented in the production solver**.
**Branch:** `feature/two-layer`.
**Scope:** active solve-cell mask, Lie publication, halo refresh and compact
solve worklist. Do not alter the physical model or numerics as part of this
performance refactor.

This document is the implementation specification for **active-domain
management only**. The [two-layer roadmap](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md)
remains authoritative for layer physics and Lie scheduling; the
[Lie transaction roadmap](LIE_TIMESTEP_TRANSACTION_ROADMAP.md) remains
authoritative for accepted-step publication and fatal failure semantics.

## 1. Decision and data ownership

**Selected approach:** rebuild the global solve worklist by an **OpenMP
full-grid scan of `solve_mask`**, with row counts, row offsets and parallel
output. The standalone test in
[`TESTS/TEST_DOMAIN_COMPACTION`](../TESTS/TEST_DOMAIN_COMPACTION/README.md)
compares serial full-grid, parallel full-grid and sparse merge compaction.
The global OpenMP method is selected for simplicity and generally favorable
measured performance; do not add a sparse/automatic switching backend.

- `domain%solve_mask(j,k)`: shared 2-D membership map for the **finalized**
  solve domain; one entry per cell, hence no worklist duplicates.
- `domain%j_cent(l)`, `domain%k_cent(l)`, `domain%solve_cells`: **existing**
  global, ordered worklist consumed by numerical kernels. Preserve these
  arrays; do not introduce a persistent, redundant `solve_index(:)` unless
  separately replacing all existing consumers.
- Temporary `row_count(:)` and `row_offset(:)` of size `O(Ny)` support
  worklist compaction. They are bookkeeping, not a second domain state.
- No thread-private cell lists, list merge, sorting, duplicate removal,
  tile framework or additional long-lived per-thread masks.
- The shared workset is the union of active cells across L1/L2. For a
  single-layer run the same algorithm uses L1 only.

The current code still reconstructs `solve_mask` from full-grid state,
expands two solution-halo layers inside `!$OMP MASTER`, and serially scans
the mask to populate `j_cent/k_cent`. The Lie publication currently updates
only `state%q/qp`. The sequence below describes the **target**, not the
current behavior.

## 2. Timestep sequence

### A. Accepted Lie publication: remove newly dry cells

After all Lie source groups have succeeded, retain the current parallel
publication loop over the **old** global solve worklist `S^n`:

1. Publish accepted conservative and primitive states for every cell in
   `S^n`, using the shared L1/L2 workset.
2. Evaluate final wetness from total conservative mass in **either** layer.
3. Set `solve_mask(j,k)=.FALSE.` for cells that became dry. Do not expand
   neighbors or construct a new worklist in this loop.

Each old worklist cell has exactly one loop iteration, so these deactivation
writes need no atomics. Publication is still performed **only after** the
whole accepted Lie transaction; never change partial/rejected-step behavior.
The old worklist remains readable until the new one is fully built.

The mask is *provisional* between this deactivation and the next domain
refresh: it must not be consumed as a finalized solve-domain membership
map during that interval. Audit any stochastic, diagnostic and other mask
users before moving the update into the Lie routine.

### B. Begin the next timestep: determine authoritative seeds

After time-dependent parameters and source geometry are refreshed, build the
seed set for the new domain at the **new** simulation time:

- Cells with positive committed conservative mass in either layer and
  `solve_mask_time(j,k) <= t`. Read final `state%q` to determine wetness,
  **not the mutable `solve_mask`**, because halo cells may have been marked
  active even if they are dry.
- The existing bottom radial/fissural and radial inlet seed rules:
  `cell_source_fractions > 0` (when enabled) and
  `source_cell == 2` (when enabled).
- Any boundary or delayed-release activations implied by the existing
  `check_solve` semantics, including cells outside the previous solve list.

Do **not** assume all seeds are in `S^n`: time-triggered releases and source
changes may activate previously inactive cells. Reuse existing source
geometry and scheduling information; do not duplicate source configuration.
Any new seed enumeration optimization must demonstrate complete coverage.
Initialization and the `solve_all=.TRUE.` path remain separately supported.

### C. Expand seeds into the solve halo, using OpenMP

Maintain the existing `active_halo_layers=2` and the existing
**four-neighbor (Manhattan)** geometry. Two layers around a seed cover
offsets `(dj,dk)` with `ABS(dj)+ABS(dk) <= 2`; do not accidentally create
a rectangular/diagonal halo or wrap across row boundaries.

Use an OpenMP loop over known wet/source seeds, with each iteration issuing
writes of `.TRUE.` to valid destination positions of `solve_mask`.
One straightforward option is to apply the complete radius-two offset set
directly to each seed: this avoids a second pass that could mistake a
newly activated *dry* neighbor for a wet seed. Alternatively, two
synchronized dilation passes from an immutable snapshot must yield exactly
the same mask.

**Concurrency contract:**

- A seed loop reads immutable committed mass/source information, **not**
  `solve_mask` for deciding whether to expand.
- Multiple threads may target the same neighbor. Use `!$OMP ATOMIC WRITE`
  (or a verified equivalent race-free method) for shared-mask writes.
  Repeating the value `.TRUE.` is harmless mathematically, but unsynchronized
  concurrent writes are **not** OpenMP-correct.
- Complete expansion and synchronize **before** reading the finalized mask
  or performing exclusions/compaction.
- Respect the exact boundary stencil of the existing implementation:
  outer boundary cells are initially kept active; expansion currently
  updates interior destinations only.
- Apply the existing `source_cell == 1` radial exclusion **after** halo
  creation. Do not make excluded cells part of the finalized solve worklist.

Before implementation, verify the required ownership and correctness of
intermediate mask values, particularly if the source geometry or its
activation time changes during a step. The target is equivalence to the
current two-layer halo logic, not a change in stencil width.

### D. Compact the finalized mask: global OpenMP row scan

For `Nx = comp_cells_x`, `Ny = comp_cells_y`, process each Fortran-contiguous
row `k` in increasing `j`. This is the **chosen production algorithm**:

1. **OpenMP count:** independently compute
   `row_count(k) = COUNT(solve_mask(1:Nx,k))` for `k=1..Ny`.
2. **Exclusive prefix sum:** set `row_offset(1)=0` and
   `row_offset(k+1)=row_offset(k)+row_count(k)`.
3. **OpenMP emit:** traverse each row again; start its private counter at
   `row_offset(k)` and write active coordinates into the existing global
   `j_cent/k_cent` arrays.
4. Set `solve_cells=row_offset(Ny+1)` after the parallel region and verify
   the worklist agrees exactly with the mask.

Each row owns a disjoint output segment; no atomics, private cell lists,
sort or duplicate removal are needed. The result is ordered by
`p=j+(k-1)*Nx` **independently of OpenMP scheduling or thread count**.
The compaction makes two full-grid mask passes, costing `O(Nx*Ny)`, plus
a small `O(Ny)` prefix sum. The row counters must be private while emitting.

### E. Build dependent worksets and make the new domain visible

Refresh the face worklists, HP workset and reconstruction worklist **from
the finalized new solve domain**, preserving their existing required
membership and ordering. In particular, the HP workset currently visits
solve cells before appending their neighbors; the reconstruction halo
extends two further cells for stencil reads and must not be conflated
with solve cells. A separate optimization may address remaining global
scans in face/reconstruction workset construction; the selected
solve-compaction algorithm does not itself remove those costs.

Ensure all workset updates are complete and synchronized before spatial,
CFL, stochastic or other consumers access the new domain. Keep any
necessary mask lookups `O(1)`; do not replace them with searches through
the list. Preserve initial/all-domain, restart, boundary and delayed-release
behavior.

## 3. Implementation and verification roadmap

| Phase | Change group | Required evidence |
| --- | --- | --- |
| 0 | Audit `check_solve` callers, mask consumers, release/source scheduling, initialization, face/HP/reconstruction dependencies and current OpenMP ownership. | Establish a baseline on the active HEAD, including optimized and strict builds. |
| 1 | Extend accepted Lie publication to deactivate dry cells without changing trial, failure or source semantics. | Single-/two-layer wet/dry equivalence; no shared writes to the same cell; no premature mask consumers. |
| 2 | Apply time-aware source/boundary seeds and the existing two-ring halo through race-free OpenMP writes. | Differential comparison of finalized `solve_mask` against current `check_solve`, including source exclusions, delayed releases, edges, empty and fully wet grids. |
| 3 | Replace the serial solve-list scan with count / prefix sum / parallel emit into `j_cent/k_cent`. | Exact ordered-list equality with serial compaction for 1, 2, 4 and 8 threads; unique coordinates; mask/count consistency. |
| 4 | Reconcile face, HP, reconstruction, stochastic and other worksets/consumers. | Unchanged memberships, ordering where required, and numerical outputs; REAL/COMPLEX and two-layer compatibility where applicable. |
| 5 | Profile representative LAHAR and other single-layer cases with output conditions matched, using 1, 4 and 8 threads. | End-to-end correctness, conservation/dry states, OpenMP reproducibility and measured total timestep cost; no unsupported performance claim. |

The [standalone compaction benchmark](../TESTS/TEST_DOMAIN_COMPACTION/README.md)
already tests the three candidate compaction algorithms independently.
It does **not** validate moving mask updates into Lie, time-dependent
source enumeration, the two-ring dilation, auxiliary worksets or
end-to-end solver equivalence. Add targeted integration/regression tests
for those behaviors before claiming completion.

## 4. Explicit exclusions and future work

- **No change to physics or numerical stencils.** Retain two solution-halo
  layers and the separate reconstruction stencil. Any proposal to reduce
  the halo width needs independent scientific authorization and regression
  evidence.
- **No current tile/MPI/GPU implementation.** A future tile design may
  replace row ownership with tile ownership while keeping the same rule:
  one owner writes each final worklist entry.
- **No sparse merge in production.** Keep its independent benchmark as
  comparative evidence, not as a second worklist-construction backend.
- **No claim of performance improvement until integration profiling.**
  Standalone benchmark results are hardware- and occupancy-dependent.

**Completion criterion:** the active mask and ordered solve worklist are
consistent across all supported execution paths, all dependent worksets
and numerical results are preserved, and relevant OpenMP/strict/optimized
regressions pass. Until then, this document is an approved design roadmap,
not a statement that the solver has already been migrated.
