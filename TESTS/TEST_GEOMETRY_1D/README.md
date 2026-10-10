# One dimensional geometry regression

This checks refreshed Q1 topography, filtered slopes/curvatures and center/face
gravity factors along both Cartesian axes, including short lines and a single
cell. It also retains a two-dimensional control with a mixed derivative.

```sh
KEEP_TEST_WORKDIR=1 sh TESTS/TEST_GEOMETRY_1D/run_test.sh
```

The [frozen contract](contract.json) defines 112 cases per compiler profile:
line lengths 1, 2, 3, 4, 5, 9 and 17 along x/y, counting the common single-cell
grid once, plus a 7-by-9 control; four flat/linear/quadratic/rough beds and
slope correction OFF/ON. All caches are poisoned before refresh. Strict and
optimized builds must agree with an independent polynomial/Q1 oracle within
`2048*epsilon(float64)` times each field scale. Every raw snapshot must match
bitwise between actual one/four-thread teams.

Lines of at least five cells retain the original five-point filter and nearest
interior extrapolation. Short lines fit degree `min(2,n-1)` over their available
centers; a single cell has no resolved derivative. Inactive and mixed 1D
derivatives are zero. Nodal elevations and Q1 face/center values are not filtered.

`--baseline-only` selects the 32 previously valid x/2D controls for a comparison
with the preceding module. It is not a complete geometry acceptance run. The
full regression does not claim support for arbitrary two-dimensional grids
with fewer than five cells along a resolved axis.
