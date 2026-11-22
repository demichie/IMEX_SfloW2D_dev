# Active-domain OpenMP roadmap

**Status:** implemented on `feature/two-layer`; **full solver CI and scientific acceptance pending**.
**Branch:** `feature/two-layer`.
**Goal:** remove redundant global work and memory transfers, then parallelize remaining loops.
Preserve the physical model, numerical stencil, source scheduling and accepted Lie transaction.

The [two-layer roadmap](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md) owns layer physics;
the [Lie transaction roadmap](LIE_TIMESTEP_TRANSACTION_ROADMAP.md) owns accepted-state
publication and fatal failure semantics. This document specifies only active-domain management.

## 1. Final strategy and data ownership

There is **one** shared `domain%solve_mask(j,k)` for all active layers and
**one** global ordered solve list in the existing `domain%j_cent(:)`,
`domain%k_cent(:)`, and `domain%solve_cells`. A cell is wet if total
accepted conservative mass is positive in at least one active layer.
The final solve domain also includes the current source, boundary and halo rules.

Reuse `domain%solve_mask_temp(:,:)` as the second ping-pong mask.
No per-layer solve masks, per-thread cell lists, sparse merging,
atomic writes, sorting or duplicate removal. Two small scalars in
`domain_type` track valid Lie wetness and the next deferred release:
`committed_mask_ready` and `next_release_time`; neither duplicates flow state.

The agreed sequence is:

1. **Accepted Lie loop:** update wet/dry membership for old worklist cells
   while publishing `q/qp`, without an extra global state scan.
2. **Start of next timestep:** `check_solve` completes the seed mask from
   existing source, time-gating, boundary, initialization and restart rules.
3. **First OpenMP dilation:** read `solve_mask` and write the first halo
   to `solve_mask_temp`, copying only the perimeter to preserve boundaries.
4. **Second OpenMP dilation:** read `solve_mask_temp`, write `solve_mask`,
   apply radial exclusions *after* dilation, and count final cells per row.
5. **Prefix sum and OpenMP compaction:** write the existing `j_cent/k_cent`
   list directly from the finalized logical mask.
6. **Dependent worksets:** preserve safe face, HP and reconstruction
   ownership and their current geometrical and ordering contracts.

Retain **two** four-neighbor/Manhattan halo layers
(`active_halo_layers=2`). The separate two-cell reconstruction-read halo
does not enlarge the evolved solve domain.

## 2. Cost comparison: production now versus the selected design

Let `N=Nx*Ny`, `M` be the old number of solve cells, and `Ny` the number
of grid rows. **Less work** means fewer operations or memory passes;
**OpenMP** means similar work in parallel, potentially reducing elapsed
time. Neither alone proves faster end-to-end simulations.

| Operation | Current production code | Target | Benefit type |
| --- | --- | --- | --- |
| Wet/dry mass checks | OpenMP on `N` cells | Lie OpenMP on `M` cells | **Less work** if `M << N` |
| Reset interior mask | Full-grid reset | Lie clears old dry cells | **Potential saving**; audit needed |
| Source/release activation | Global checks | Preserve coverage | **No guaranteed saving** |
| Halo mask copies | Two full `O(N)` copies | Border-only `O(Nx+Ny)` copy | **Less work / traffic** |
| Two halo OR passes | Two serial stencil scans | Two OpenMP stencil scans | **OpenMP wall time** |
| Radial exclusions | Separate global scan | Fuse into second halo | **Less work / traffic** |
| Count cells by row | No separate row counts | Count inside second halo | Avoid **new** global count scan |
| Global list scan | Serial `O(N)` | OpenMP `O(N)` | **OpenMP wall time** |
| Prefix sum | Absent | Serial `O(Ny)` | **New small cost** |
| Mark face masks | In serial list scan | Separate safe `O(M)` pass | **Possible new cost** |
| Worklist `COUNT` assertion | Full-grid count | Retain initially | **No saving** |
| Face/HP/reconstruction lists | Existing passes | Unchanged initially | **No claimed gain** |

**Do not conflate the baseline with an intermediate design.**
The pre-refactor production list scan counted and emitted cells
together; the fused `row_count(k)` saves a full-grid count pass
**relative to the earlier proposed OpenMP count/emit scheme**, not
relative to the current serial scan. The retained
`COUNT(solve_mask)` consistency assertion is another full-grid read.

The two eliminated full-mask copies are an actual memory-traffic
reduction. Parallelizing the two halo OR passes does **not** remove
their neighbor tests. Thread startup/barriers, memory bandwidth,
additional face marking and deferred releases can offset speedups.
The decisive metric is the complete timestep, not kernel timing alone.

## 2a. Dependent workset optimization

The solve halo and ordered solve list remain the only authoritative active
domain representation. After row-wise solve compaction, `check_solve` now builds
the dependent worksets within one additional OpenMP team:

- Face masks are computed directly from the OR of their adjacent solve cells.
  Every face is written by exactly one thread, eliminating serial overlapping
  writes from the solve list and full-mask clearing.
- The reconstruction read halo is expanded twice with the existing
  `reconstruction_mask` and `solve_mask_temp`, without two full-grid copies.
- Each HP cell belongs to the first adjacent solve cell in the original
  `(k,j)` traversal. Per-row counts and prefixes preserve the historical HP
  ordering (solve cells first, then unique W/E/S/N neighbours) while emitting
  in parallel. No per-thread lists, atomics or sorting are introduced.
- Reconstruction and face index arrays are emitted in parallel into disjoint
  row-owned segments, preserving the exact lexicographic order.
- Existing full-grid validation scans are replaced by per-row compacted-count
  checks. Small local row counts/offsets are scratch, not persistent state.

The production spatial test independently reconstructs HP neighbours, both
reconstruction dilation rings and face lists using the former serial
algorithms, then compares the exact masks, counts and ordered indices for
initial, deferred-release, incremental, full-domain and shared-layer cases.
OpenMP numerical regression and end-to-end performance must be measured on
the actual solver before claiming a speedup.

## 3. Algorithm and concurrency details

### A. Accepted Lie publication

Use the existing OpenMP loop over the **old** unique solve worklist.
After every Lie source group has succeeded, publish `state%q/qp` and
assign the wet/dry value of the one shared `solve_mask(j,k)` from
positive conservative mass of either active layer. For single-layer
runs only L1 participates. Do not change the wetness threshold to
primitive depth or alter trial/publication ordering.

Each old worklist coordinate has one owner, requiring no atomic write.
The resulting mask is **provisional** until the next `check_solve`.
The old list remains available until the new list is complete.
Audit consumers between these points: stochastic updates, diagnostics,
output, checkpoint/restart and any other mask readers must not
misinterpret this provisional mask as the finalized old workset.

### B. Seed completion at the new timestep

Preserve the current `solve_mask_time(j,k)<=t` activation rule;
boundary values; source rules for `cell_source_fractions>0` and
`source_cell==2`; initialization, restart, and `solve_all=.TRUE.`.
Recover delayed releases or newly activated source cells **outside**
the old list. Do not treat previously active dry halo cells as wet.

The normal path should reuse Lie's logical wetness updates rather than
rereading real mass across the full grid. However, **this saving depends
on proving complete activation of cells outside the previous list**.
Keep the necessary state checks on exceptional paths if required.
Do not duplicate persistent wetness or source-state representations.

**Release scheduling:** `next_release_time` is initialized from the
read-only `solve_mask_time` raster and refreshed when an event is crossed.
Normal committed Lie steps reuse the previous worklist's wetness flags and
check release eligibility there without a global mass scan. When an event
is due, `check_solve` performs a complete positive-mass check to recover
newly eligible cells anywhere on the grid. Without an accepted Lie
publication (including direct tests and restart initialization), the
current conservative state is scanned to reconstruct wetness.
The outer perimeter keeps its original membership rules; Lie updates
interior wetness only. These are correctness conditions, not parallel
legacy implementations.

### C. First dilation: `solve_mask` -> `solve_mask_temp`

Copy **only** the perimeter values from `solve_mask` into
`solve_mask_temp`. The current numerical halo modifies interior
destinations only. In an OpenMP loop with exclusive destination
ownership, calculate the interior first halo:

```fortran
solve_mask_temp(j,k) = solve_mask(j,k) .OR. &
                       solve_mask(j-1,k) .OR. &
                       solve_mask(j+1,k) .OR. &
                       solve_mask(j,k-1) .OR. &
                       solve_mask(j,k+1)
```

Synchronize all threads before starting the second dilation.
The preserved perimeter is required because cells next to the
boundary read its values in the second pass. Degenerate grids
with no interior must not read nonexistent neighbors.

### D. Second dilation: exclusions plus per-row counts

For interior cells read **only** `solve_mask_temp`:

```fortran
active = solve_mask_temp(j,k) .OR. &
         solve_mask_temp(j-1,k) .OR. &
         solve_mask_temp(j+1,k) .OR. &
         solve_mask_temp(j,k-1) .OR. &
         solve_mask_temp(j,k+1)
```

Perimeter destinations retain the previous `solve_mask(j,k)`.
**After** obtaining the two-ring halo, apply the current exclusion
`source_cell(j,k)==1` when `radial_source_flag` is active;
then assign the final `solve_mask(j,k)`. The exclusion must not
remove a seed **before** it contributes to the second dilation.

Assign each complete `k` row to one OpenMP iteration and increment
a row-private counter for each final `.TRUE.` cell, including border
cells and accounting for their exclusions. Write `row_count(k)`
once per row. Each thread owns different destinations and counts;
all stencil reads come from the immutable temporary mask.
No atomics, no intermediate activation chains and no separate
global exclusion or counting pass.

### E. Prefix sum and parallel ordered worklist

Compute `row_offset(1)=0` and, for increasing `k`,
`row_offset(k+1)=row_offset(k)+row_count(k)`.
Set `solve_cells=row_offset(Ny+1)`.

In an OpenMP loop over rows, visit `j` in increasing order and write
the existing `j_cent/k_cent` arrays only in that row's interval
`(row_offset(k),row_offset(k+1)]`. This is a single global logical-mask
scan, with no `q/qp` accesses, sorting, duplicates or atomic writes.
Its ordering exactly matches the current serial scan: increasing `k`,
then increasing `j`, independent of thread scheduling.

### F. Faces, reconstruction and verification scan

The pre-refactor serial scan also marked `solve_mask_x` and
`solve_mask_y` for adjacent faces. Placing those writes directly
in the parallel cell emitter would be unsafe: neighboring cells
can write the same face. Initially mark faces in a safe separate
pass over the compact global solve list, and measure that `O(M)` cost.
Preserve the existing HP, face-list and reconstruction-list order
and membership. Do not change the separate reconstruction halo.

The existing unconditional full-grid `COUNT(solve_mask)` assertion
remains in place. Any removal or movement of this runtime
check is a separate review decision; the integrated `row_count`
does not silently replace its independent consistency test.

## 4. Implementation and acceptance

1. **Dependency audit:** all consumers of `solve_mask`, source/release
   schedule, restart/initialization, stochastic and diagnostics behavior,
   boundary rules and shared L1/L2 workset semantics.
2. **Differential correctness:** compare old two-copy dilation with
   ping-pong dilation and fused count, checking masks, row counts and
   the **exact** ordered worklist for small/large, empty/full,
   isolated, corridor and border cases, including source exclusions.
3. **OpenMP correctness:** test 1/2/4/8 threads, disjoint writes,
   phase barriers and identical deterministic output. Verify that
   provisional masks cannot be read as finalized worksets.
4. **Coherent production changes:** Lie wet/dry update, seed completion,
   ping-pong halo, fused exclusion/count, prefix sum, ordered emission
   and safe face-mask construction; no second production algorithm.
5. **End-to-end verification:** strict/debug and optimized builds,
   single-layer LAHAR examples, experimental two-layer configurations,
   dry/wet transitions, delayed releases, restart, mass conservation,
   stochastic behavior and face/HP/reconstruction workset equivalence.
   Profile complete `check_solve` and the total timestep for
   1/4/8 threads with matched inputs and output settings.

The [compaction and halo suite](../TESTS/TEST_DOMAIN_COMPACTION/README.md)
now includes a strict OpenMP reference-versus-ping-pong halo test and
exact ordered compaction checks. The production spatial integration suite
also compares incremental and full domain refresh, including deferred
releases, shared-layer wetness, sources, exclusions and dependent worksets.
These checks are included in the manually dispatched
`physics-regressions.yml` spatial job. They do **not** replace end-to-end
optimized/debug simulation and performance verification.
Do not claim production speedups until those measurements exist.

**Completion criterion:** exactly one shared solve mask and one existing
ordered global list, unchanged numerical behavior, correct edge/source
semantics and verified OpenMP performance.

**Out of scope:** new physics, halo-width changes, per-layer solve masks,
sparse production backends, tiles/MPI/GPU, or two-layer production approval.
