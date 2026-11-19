# Solver integration CI: runtime and diagnostic policy

The `Solver integration and OpenMP validation` workflow runs `make check`
with optimized (four OpenMP threads) and strict-debug (one thread) profiles.
The job limit is 90 minutes. Each individual solver fixture has a shorter
wall-clock limit: **20 minutes optimized**, **30 minutes strict-debug**.
These bounds are for diagnosis; they do not constitute solver convergence
criteria or silently shorten simulated physical time.

## Long-running geographical fixtures

Only the test copies of Campi Flegrei and Taal were coarsened. The
corresponding `EXAMPLES/EXAMPLE_CF` and `EXAMPLES/EXAMPLE_TAAL` inputs
and the compressed original DEMs are unchanged.

| Fixture | Old grid / cell size | New grid / cell size | Domain extent |
| --- | --- | --- | --- |
| `TEST_CF` | 700 × 520 / 25 m | 350 × 260 / 50 m | 17,500 × 13,000 m |
| `TEST_TAAL` | 500 × 700 / 50 m | 250 × 350 / 100 m | 25,000 × 35,000 m |

Both grids have **one quarter the number of cells** while retaining their
initial domain origins, extents, physical durations, source parameters,
rheology, gas/solid properties and DEM input. Unlike a pure performance
optimization, coarsening **changes the spatial discretization and numerical
solution**: it must not be presented as equivalent to the original fine-grid
scientific baselines. Existing higher-resolution example configurations are
available independently for full scientific validation.

`VERBOSE_LEVEL=0` in the coarsened CI fixtures allows the solver to emit
time, timestep and budget progress into the per-case log. The runner prints
START/END records with timestamps, elapsed time, exit status and the
last logged `t`/`dt`. Under the per-case timeout it emits a heartbeat about
once per minute, even if the Fortran log is buffered. A missing timestep in
the heartbeat is not by itself evidence of a deadlock.

On failure or timeout, the runner shows the final 160 log lines, preserving
the **complete** solver log and case timing in the job artifact:

- `integration-optimized-<run_id>`
- `integration-strict-debug-<run_id>`

These artifacts are uploaded on failure using `if: always()`, and retained
for 30 days. The temporary simulation directory can still be cleaned safely.

## Focused reproduction

Build the solver first. For example, from the repository root:

```sh
mkdir -p /tmp/imex-integration-evidence
IMEX_TEST_LOG_DIR=/tmp/imex-integration-evidence \
IMEX_TEST_CASE_TIMEOUT=20m \
OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE \
sh TESTS/run_tests.sh src/IMEX_SfloW2D TEST_CF
```

The supported focused cases are `TEST_CF`, `TEST_TAAL`, and
`TEST_DYNAMIC_FINH`, `TEST_OFF_FINH`, `TEST_STATIC_FINH`,
`TEST_STEP_FINH`. The final four use the existing DAMBREAK input files.
Without a focused selector, the original complete integration sequence is
run and individual failures still make `make check` fail.

Use the first timed-out case's full log to distinguish a slow but advancing
simulation, a collapsing `dt`, an iteration/convergence failure, or a
genuine stalled solver. Do not weaken numerical tolerances or change the
physical model just to clear CI.
