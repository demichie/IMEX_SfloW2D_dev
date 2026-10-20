# TEST_SOURCE_MASS

Checks source mass injection and rejection of a source setup with invalid
buoyancy.

`make_cases.py` writes nine mass-balance cases on a flat bed with settling,
erosion and entrainment off, so the mass in the domain is the mass the source
put in:

| Case | Source |
|---|---|
| `fissures` | two fissures at 800 K and 700 K, one oblique to the grid; linear ramps, a phase and a second window |
| `fissure_map_coordinates` | one fissure always on, at map coordinates of order 1e6 m |
| `lateral_window` | lateral radial source, periodic window with cosine ramps, cut off by `TIME_PARAM(4)` inside a window |
| `lateral_boundary_{west,east,south,north}` | lateral boundary sources on all four domain sides, with prescribed positive depth and constant mass-flow rate |
| `bottom_window` | elliptical bottom source, periodic window with linear ramps and a phase |
| `bottom_intervals` | circular bottom source with `N_INTERVALS` velocity steps |

An additional expected-failure case uses a source mixture lighter than the
ambient gas with a prescribed Richardson number. The solver must reject it
before deriving the source thickness or creating output. Thus the suite
contains ten cases total: nine mass checks and one input-validation check.

`check_mass.py` sums the solid mass (the mixture mass for lateral sources,
whose input is a mixture mass flow rate) in each `.q_2d` output and compares
it with the time integral of the prescribed source, computed independently of
the solver. It reads snapshots from `OUTPUT/<case>/snapshots/`. Every output
must agree to `1e-9` relative.

`lateral_window` also checks the impermeable inner boundary of the radial
source. Injection is applied as a volume source in the surrounding cells;
the excluded source interior must not receive numerical mass flux. The KT
free-surface and selected transport-depth traces must therefore preserve the
reflecting ghost states after reconstruction. Recomputing those scalar traces
without restoring the wall previously caused a mass deficit of about 0.437%
with `froude_blend`, including mass loss after injection stopped.

The four `lateral_boundary_*` cases check that `froude_blend` preserves the
prescribed source depth at each domain side when estimating wave speeds and
computing the advective and hydrostatic interface contributions.

```sh
TESTS/TEST_SOURCE_MASS/run_test.sh src/IMEX_SfloW2D
```

`KEEP_TEST_WORKDIR=1` keeps the cases and their output. Only `python3` and
`numpy` are needed.
