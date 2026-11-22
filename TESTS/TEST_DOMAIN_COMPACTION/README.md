# Active-domain compaction experiment

Standalone Fortran/OpenMP differential test and microbenchmark; this directory does **not** change
`domain_2d`, `time_integration_2d`, physical formulas, timestep policy, or the production build.

## Algorithms

- **Reference:** scan `solve_mask(nx,ny)` in `(k,j)` order, producing sorted linear indices
  `p = j + (k-1)*nx`.
- **Sparse parallel candidate:** take sorted, unique `old_indices` and `external_seeds`.
  Merge virtual offset streams `(0,-1,+1,-nx,+nx)` from both inputs, restricted to disjoint
  destination-index blocks. Each block counts its selected mask entries; an exclusive prefix
  sum provides disjoint output positions; a second parallel pass emits ordered, unique indices.
  Row-wrapping offsets are rejected. The ten streams are **not** materialized as arrays.

The required **completeness contract** is that every true entry of the updated mask lies in
`halo_1(old_indices union external_seeds)`. External seeds stand for new, known source/release/
boundary activations; they must be provided by the caller. This test deliberately does not
validate the physical correctness of how the solver generates these seeds, the width of the
numerical stencil, or whether a cell becomes wet. Exclusions can remove candidates.

No thread-local candidate/output lists, global-grid mask scan inside the sparse implementation,
or atomic writes to the output list are needed. The implementation uses small per-block count
and offset arrays, and its output order is deterministic across thread counts.

## Run

From the repository root:

```sh
TESTS/TEST_DOMAIN_COMPACTION/run_test.sh
TESTS/TEST_DOMAIN_COMPACTION/run_test.sh --benchmark
TESTS/TEST_DOMAIN_COMPACTION/run_test.sh --benchmark-large
```

Set `FC` to change the Fortran compiler. The correctness-only command uses `-O0`, `-fcheck=all`,
`-Wall -Wextra -Werror`. The benchmark commands compile with `-O3` and repeat correctness
checks before timing. All commands need a compiler supporting OpenMP and `omp_lib`.

The correctness suite compares exact index arrays and reconstructed masks across 1, 2, 4 and
8 threads. It covers empty and full domains, 1-cell dimensions, edges, thin regions, islands,
external activation, source-like exclusions and candidate block boundaries. The benchmark
uses compact-square and corridor geometries on 500x500 or 2000x3000 grids with nominal active
fractions of 0.1%, 1%, 10% and 100%. It reports median elapsed times of three repetitions;
mask preparation, reconstruction worksets and solver execution are excluded. The percentages
are geometry-generation targets, not necessarily exact old-list occupancies.

The reference is a serial scan; benchmark results **must not** be interpreted as a comparison
against an optimized parallel full-grid scan. Results are hardware-specific and are intended
for algorithm selection, not production-speedup claims. Test kernels are intentionally kept
independent of the production modules and are not yet wired into `make check` or GitHub Actions.
