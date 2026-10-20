# CFL policy and accepted timestep tests

This production test verifies `directional_min` (including omission) and
`combined` on the same Q1 excavation. It covers `N_RK=2`, the prescribed
`KT_TRACE_RECONSTRUCTION='froude_blend'`, the `CFL=-1` sentinel, unchanged 2D CFL limits and invalid policy names. It does
not approve an N9 scientific baseline or establish universal IMEX stability.

The temporal fixture adds eight dry cells per side to avoid treating open-boundary
outflow as a closed-domain mass defect. The frozen Gate-H fixture remains unchanged
(`--pad-cells` defaults to zero). The initial compact temporal fixture failed its
unchanged budget bound: mass reached the open north/east/south boundary, so that
result is retained as a fixture diagnosis, not waived or used to relax the bound.

Run with a freshly compiled solver:

```sh
KEEP_TEST_WORKDIR=1 OMP_NUM_THREADS=4 \
  sh TESTS/TEST_CFL_POLICY/run_test.sh /absolute/path/to/IMEX_SfloW2D
```

The test generates 18 valid trajectories and three invalid-input controls.
A test-only main records both the spatial proposal and the accepted step;
all canonical conservative outputs and the final binary checkpoint must be
byte-identical to the uninstrumented solver. An independent checker verifies
the existing growth limiter, all three output events and the final event.
Limits fixed before execution are `64*eps64` for time algebra and
`8192*eps64*max(1,sum(abs(initial_mass)))` for the closed-domain mass budget.
Default/omitted and all sentinel-policy comparisons require exact byte identity.

Set `IMEX_CFL_BASELINE` to the executable built at `691a3e21`, using the same
compiler profile, to additionally compare omitted-policy trajectories against
the preceding implementation for all three stage counts. The evidence records
whether this optional historical-executable comparison was actually performed.
After that run, `compare_baseline.py CASE_ROOT BASELINE_EXECUTABLE REPOSITORY`
also checks spatial proposals and accepted times byte-for-byte against a
numerically neutral observer of the historical executable.
The full acceptance runner executes this test with real one/four-thread teams
in strict and optimized profiles and compares canonical outputs/checkpoints.

Supplementary original-duration restart runs can explicitly request
`IMEX_AUDIT_CFL_POLICY=combined` through the acceptance capture adapter. The
override is written to each effective input and invocation manifest; omission
does not change any original audit fixture.

`TEST_SPATIAL_OPERATOR` separately checks both formulas against the exact same
face-speed payload, the clean Python combined bound, zero propagation, caps,
sentinel bypass, unequal x/y spacings and rotated 1D x/y/point grids. Historical
Python cores and the existing scientific comparison limits remain unchanged.
