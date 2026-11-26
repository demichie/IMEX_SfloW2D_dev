# Geometry and shared active-domain worksets

`src/geometry/` manages Cartesian terrain geometry and the common
active cell/interface worksets used by finite-volume kernels.
Cell-local physical laws do not traverse the grid.

## Modules

| File | Responsibility |
| --- | --- |
| [`geometry_2d.f90`](geometry_2d.f90) | Grid, elevations, slopes, curvature and refresh after terrain changes |
| [`domain_2d.f90`](domain_2d.f90) | Shared wet/active masks, source/release activation and compact worklists |

## Active-domain implementation

The one- and two-layer solver uses **one** `domain%solve_mask(j,k)` and
one globally ordered cell list `j_cent/k_cent`. A cell is wet when its
accepted conservative mass is positive in at least one active layer.

- Accepted Lie publication updates provisional wetness on the old
  worklist. `check_solve` rebuilds the finalized workset, including
  boundary, source, initialization, restart and deferred release rules.
- Two four-neighbour/Manhattan halo layers
  (`active_halo_layers=2`) are formed with OpenMP row-owned passes over
  `solve_mask` and `solve_mask_temp`. Perimeter handling and the order
  of radial source exclusions are correctness conditions.
- Per-row counts and offsets compact the final mask in parallel into
  the **same globally ordered** `j_cent/k_cent` worklist. Face,
  hydrostatic-path and reconstruction worksets preserve exact
  membership and ordering. The separate reconstruction read halo
  does not expand the evolved solve domain.
- The production implementation does not require per-layer solve lists,
  per-thread merge/sort stages, atomics or duplicate persistent wetness.

## Verification and performance

`TEST_DOMAIN_COMPACTION` and `TEST_SPATIAL_OPERATOR` check halo and
ordered workset equivalence, including deferred releases, sources,
shared-layer wetness and boundary cases. Changes require optimized
and strict/debug tests plus representative end-to-end simulations,
OpenMP reproducibility and restart checks.

Reducing mask copies and parallelizing stencil passes do **not**
automatically establish a whole-simulation speedup. Compare complete
timings with matched cases, outputs and 1/2/4/8-thread settings.
See [the numerical integration guide](../numerics/README.md).
