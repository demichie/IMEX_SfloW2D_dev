# Fortran 1-D inclined excavation counterpart

This case is the Fortran setup corresponding to `python_diagnostics/kt_v2_1d.py`
with `--integrator fortran2`. The Python solver and its spatial discretization
are left untouched.

The matching setup is 120 x 1 cells, `dx=0.25 m`, slope `10°` down toward +x,
excavation `10 <= x <= 20 m` with depth `2 m`, initially zero velocity, pure
liquid (`rho=1000 kg/m^3`), `g=9.81 m/s^2`, `CFL=0.14`, `dt0=1e-4 s`,
`max_dt=0.02 s`, output interval `0.1 s`, and `n_RK=2`. Both use the KT
spatial scheme and disabled slope/curvature corrections. The lateral extent is
one cell only, so the Fortran solver takes its 1-D execution path.

`create_case.py` generates the ASCII DEM, initial `.q_2d` restart and input
file. The DEM is sampled from the exact face-elevation profile used in the
Python case; the restart thickness uses the exact cell-average face weights
used by Python. Run the generator from this directory, then run the Fortran
executable with this directory as the working directory:

```sh
python3 create_case.py
../../IMEX_SfloW2D
```

The expected Fortran behavior is not assumed to match the Python result yet:
the purpose is to make any difference reproducible with identical initial
geometry, state, and principal time-integration settings. Compare the console
failure time and the saved conservative snapshots against:

```sh
python3 ../kt_v2_1d.py --integrator fortran2 --t-end 1.5
```

The Python diagnostic includes its own explicit dry-donor check and raw
negative-thickness stop; these are part of that reference implementation and
are not altered by this Fortran case setup.
