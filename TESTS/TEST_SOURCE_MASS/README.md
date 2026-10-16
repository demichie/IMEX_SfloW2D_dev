# TEST_SOURCE_MASS

Checks that every source type injects the mass it prescribes.

`make_cases.py` writes five small cases on a flat bed with settling, erosion
and entrainment off, so the mass in the domain is the mass the source put in:

| Case | Source |
|---|---|
| `fissures` | two fissures at 800 K and 700 K, one oblique to the grid; linear ramps, a phase and a second window |
| `fissure_map_coordinates` | one fissure always on, at map coordinates of order 1e6 m |
| `lateral_window` | lateral radial source, periodic window with cosine ramps, cut off by `TIME_PARAM(4)` inside a window |
| `bottom_window` | elliptical bottom source, periodic window with linear ramps and a phase |
| `bottom_intervals` | circular bottom source with `N_INTERVALS` velocity steps |

`check_mass.py` sums the solid mass (the mixture mass for the lateral source,
whose input is a mixture mass flow rate) in each `.q_2d` output and compares
it with the time integral of the prescribed source, computed independently of
the solver. It reads snapshots from `OUTPUT/<case>/snapshots/`. Every output
must agree to `1e-9` relative.

```sh
TESTS/TEST_SOURCE_MASS/run_test.sh src/IMEX_SfloW2D
```

`KEEP_TEST_WORKDIR=1` keeps the cases and their output. Only `python3` and
`numpy` are needed.
