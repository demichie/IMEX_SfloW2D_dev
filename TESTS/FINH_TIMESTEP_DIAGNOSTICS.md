# FINH timestep-collapse investigation

The four single-layer granular dam-break fixtures in
`TESTS/TEST_DAMBREAK` share the same 300 x 15 grid, 0.01 m cells, 0.30 s
physical horizon, regularized mu(I) rheology, gas loss, and dynamic
permeability. The difference is
`PORE_PRESSURE%F_INHIBIT_MODE`: `OFF`, `STEP`, `STATIC`, or `DYNAMIC`.

A prior optimized CI run of `TEST_DYNAMIC_FINH` timed out after 1200 s
while the rounded progress log still showed `t=0.0006` and
`dt=0.00000`. The rounded `F9.5` output cannot tell whether the actual
step is nonzero, collapsing, or no longer advancing floating-point time.
The origin of the small step is not established. No timestep floor,
limiter alteration, pore-pressure formula change, or change to the
fixtures is authorized by this observation.

## Reproduce with numerical evidence

Build the existing solver. Run from the repository root:

```sh
python3 TESTS/compare_finh_timestep.py \
  --solver src/IMEX_SfloW2D \
  --evidence /tmp/imex-finh-diagnostics \
  --timeout-seconds 90
```

This executes each of the four **unchanged** fixture inputs through
`TESTS/run_tests.sh`, with one OpenMP thread, the original physical end
time of 0.30 s, and a wall-time cap of 90 s **per case**. The cap is
diagnostic only and is not an acceptance criterion. Successful tests
still require the full physical horizon and all ordinary output checks.

To reproduce only one mode:

```sh
python3 TESTS/compare_finh_timestep.py \
  --solver src/IMEX_SfloW2D \
  --evidence /tmp/imex-finh-dynamic \
  --timeout-seconds 90 --case TEST_DYNAMIC_FINH
```

The report contains `finh_comparison.csv` and
`finh_comparison.json`, plus per-case `timesteps.csv`,
`TEST_*_FINH.solver.log`, `TEST_*_FINH.timing.txt` and
`runner.log`. Compare `last_t`, `min_cfl_dt`,
`min_accepted_dt`, `max_a_x`, `max_a_y` and the per-step trends
before attributing a failure to dynamic inhibition.

The solver writes `TIMESTEP_DIAG` records when
`IMEX_TIMESTEP_DIAGNOSTICS=1`. The records are printed for the first
20 steps and every 100 steps thereafter. They include the exact
scientific-notation `t`, unmodified CFL proposal, applied step after
the existing output/previous-step caps, and maxima of the x/y face wave
speeds. The diagnostic uses the CFL work already computed during the
normal evolution. When diagnostics are disabled, neither this output
nor a second CFL evaluation occurs.

The driver also rejects any selected nonfinite, nonpositive or
non-advancing timestep (`t + dt <= t`) with `INVALID_TIMESTEP`.
This is an explicit failure, not a silent lower-bound substitution.
The check is performed before updating the two previous-step values
or advancing the state.

Run the diagnostic-parser unit test independently:

```sh
python3 -m unittest discover -s TESTS -p 'test_compare_finh_timestep.py'
```

GitHub Actions workflows remain **manual-only** (`workflow_dispatch`);
the comparison script does not dispatch any workflow.

## Interpretation

If `cfl_dt` becomes tiny together with rapidly growing face-speed
maxima, inspect the reconstructed states, dry/wet transition, and pore
pressure/mixture conversion that feed the HP-PCCU characteristic speeds.
If `cfl_dt` stays reasonable while `accepted_dt` shrinks, inspect
the existing previous-step growth cap and output-time clipping.
If `DYNAMIC` differs from `OFF`, `STEP`, and `STATIC`, isolate
the first divergent cell/state before changing the constitutive law.
A timeout alone is not evidence of a deadlock or proof of a faulty
inhibition formula.
