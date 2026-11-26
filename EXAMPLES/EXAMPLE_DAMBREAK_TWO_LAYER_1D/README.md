# Two-layer 1D dam-break on flat topography

This case isolates the two-layer spatial operator on flat topography. The mesh has
**200 cells in x and exactly one cell in y** (`COMP_CELLS_Y=1`).
The topography is constant zero; no transverse evolution is intended.

Initially L1 is 2 m deep for `0 <= x < 60 m`, with zero velocity, and
dry elsewhere. L2 starts completely dry. The initial L1 volume is 120 m3.
Both layers use the existing primitive-input initializer. Solid volume
fraction in L1 is 50%, temperature 300 K. Rheology, erosion, settling,
elutriation and interlayer drag are disabled to isolate spatial transport.

## Generate inputs

From this directory:

```bash
python3 prepare.py
```

This generates only `DEM/topography.asc` and `DEM/l1_thickness.asc`.
There is no need to run the solver to prepare initial conditions.

## Run

With the executable built at `../../src/IMEX_SfloW2D`:

```bash
../../src/IMEX_SfloW2D
```

The simulation ends at 3 s and writes output every 0.25 s in
`OUTPUT/dambreak_two_layer_1d/`. Because L2 remains dry, this test
does not exercise two wet interacting layers. For the coupled
characteristic limits, see
[two-layer hyperbolicity](../../docs/TWO_LAYER_HYPERBOLICITY.md).

## Interpretation

Verify that L1 momentum develops from rest and that its wet front propagates
in the positive x direction. Since all L2 transfer processes are disabled,
L2 should remain dry. Check total L1 mass conservation, zero transverse
momentum, and the timestep history. For a separate L1 single-layer baseline, use `N_LAYERS=1` and remove
`RESTART_L2_PARAMETERS` and the other L2 namelists. Single-layer input
accepts either `RESTART_PARAMETERS` or `RESTART_L1_PARAMETERS` and either
`INITIAL_PARAMETERS` or `INITIAL_L1_PARAMETERS`, with one canonical set
of L1 parameters;
do not change the physics or numerical parameters between comparisons.

The single y cell has no y-direction interior flux and is intentional.

## Wet/dry-front mass diagnostic

Enable the same optional diagnostic in both a single-layer reference
and the two-layer run:

```bash
IMEX_FRONT_MASS_DIAGNOSTICS=1 ../../src/IMEX_SfloW2D > front.log 2>&1
grep FRONT_MASS front.log | head -n 30
```

At each accepted timestep the solver finds the first x-face whose L1
left cell is wet and right cell dry. It prints the face index, time and
timestep, the oriented numerical mass fluxes delivered to each adjacent
cell, and the right cell's L1 mass before and after the spatial update,
after source groups, and after publication. A positive `flux_to_right`
means rightward mass transport; the right-cell spatial mass residual must
then be negative when this is its only nonzero mass flux. If the spatial
candidate gains mass but the post-source value becomes zero, investigate
the source executor. If the post-source value is positive but the
published value is zero, investigate the final state commit. The diagnostic
reads existing flux arrays and does not alter the physical state.

The companion `FRONT_TRACE` records show whether the selected x-face is
active, the conservative masses in its adjacent cells, the reconstructed
left/right depths and velocities, free-surface values, Froude-blend depths,
and characteristic speed bounds. These distinguish a missing active face,
dry reconstructed states, and vanishing coupled characteristic speeds.
If the front mass flux remains identically zero, inspect `FRONT_TRACE`
rather than the later Lie source or commit stages.

### Inspecting the flux kernel

With `IMEX_FRONT_MASS_DIAGNOSTICS=1`, the same L1 wet-to-dry
x-face is also logged **inside** the numerical flux kernel. Use
`grep FRONT_KERNEL front.log | head -n 40` to inspect one or two
timesteps. The records are identical for `N_LAYERS=1` and
`N_LAYERS=2`; run both with the same initial L1 state to identify
the first divergent stage.

- `FRONT_KERNEL enter`: the face was actually visited by `eval_flux_PCCU`.
  If the later `FRONT_TRACE` reports the face as active but this entry is
  absent, inspect the compact interface worklist.
- `input`, `speeds`, `after_pccu`: conservative face masses, normal
  velocities, characteristic bounds, endpoint and raw PCCU mass fluxes.
- `after_component_limiter`, `after_dry_guard`,
  `after_rest_guard`: successive flux values, with a
  `suppressed_by_...` reason when a zeroing branch executes.
- `oriented`, `final`, `exit`: oriented face contributions after
  hydrostatic path splitting, dry-side treatment and return to the caller.

Diagnostics are read-only, are disabled by default and do not change the
physical model. For the 1D example, there is only one right-facing wet/dry
front; in a general 2D domain, diagnostic lines from multiple x-faces can
interleave under OpenMP.
