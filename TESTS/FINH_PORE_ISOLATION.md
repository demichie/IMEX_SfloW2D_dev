# FINH isolation experiment: disable pore pressure

The original `TESTS/TEST_DAMBREAK/TEST_OFF_FINH` fixture is unchanged.
`compare_finh_pore_toggle.py` runs from the local fixture, without requiring Git metadata.
If the fixture is original, it copies it twice into temporary directories:

- `original`: original input, restart and numerical configuration.
- `no_pore_pressure`: changes `PORE_PRESSURE_FLAG=T` to `F` and removes
  the `PORE_PRESSURE_PARAMETERS` namelist **only in the temporary copy**.

If the local fixture has already been manually simplified (`PORE_PRESSURE_FLAG=F`
and no pore-pressure namelist), the script runs only the `no_pore_pressure`
variant and records its output; it cannot reconstruct the original input
from that modified file. To compare both in one invocation, restore the
original fixture first.

The namelist must be removed because the single-layer input parser rejects
`PORE_PRESSURE_PARAMETERS` when pore pressure is disabled. Disabling pore
pressure also disables gas loss and dynamic permeability. No numerical
algorithm, CFL limiter, grid, restart file or rheology parameter is changed.

Build the current branch's optimized solver, then run from repository root:

```sh
python3 -m unittest discover -s TESTS -p 'test_compare_finh*p*.py'
python3 TESTS/compare_finh_pore_toggle.py \
  --solver src/IMEX_SfloW2D \
  --evidence /tmp/imex-finh-pore-isolation \
  --timeout-seconds 120
```

Outputs: `finh_pore_comparison.csv`, and for each variant,
`IMEX_SfloW2D.inp`, `solver.log`, and `timesteps.csv`.

Compare the first 20 `TIMESTEP_DIAG` and `CFL_X_FACE` records, especially
step 6, characteristic speed maxima, and time-step reduction. A timeout
does not establish success or failure of the physical model. A run that
exits successfully must still be checked against the requested 0.30 s
physical horizon and expected solver output.

Interpretation: if the no-pore-pressure variant avoids the early CFL
collapse, subsequently isolate pore-pressure coupling, source terms, and
state conversion; if both variants collapse, next disable rheology in
another temporary copy before changing reconstruction. Do not introduce
a timestep floor or modify the FINH reference fixtures as a workaround.
