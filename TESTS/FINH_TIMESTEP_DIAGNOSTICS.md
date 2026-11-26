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

For sampled steps where `max_a_x > 100 m/s`, the solver additionally
prints `CFL_X_FACE` records. The first active x-face attaining the exact
global maximum is selected deterministically, after the OpenMP reduction.
The records include its `j,k` face indices, signed characteristic bounds,
the two full physical face states actually passed to the speed-bound kernel
(including selected depths and the dry-state velocity zeroing), cached
hydrostatic coefficients, and the neighboring conservative/primitive cell
states. The left cell is `j-1,k` and the right cell is `j,k`, when
within grid bounds. Printing is restricted to the first 20 steps and every
100th step thereafter, and only when `IMEX_TIMESTEP_DIAGNOSTICS=1`.
This diagnostic does not recompute waves or alter reconstruction, CFL,
or source integration.

For the observed OFF/STEP failure, inspect steps 5 through 9: the first
abnormal x characteristic bound occurs at step 6 (`t≈6.22e-5 s`).
Compare `CFL_X_FACE left_speed_input` and `right_speed_input` with the
neighboring cell states to determine whether the huge bound originates
from a reconstructed velocity, selected depth, or another local wave-speed
input. Do not interpret the records as proof of a specific cause until the
actual limiting face has been inspected.

The driver also rejects any selected nonfinite, nonpositive or
non-advancing timestep (`t + dt <= t`) with `INVALID_TIMESTEP`.
This is an explicit failure, not a silent lower-bound substitution.
The check is performed before updating the two previous-step values
or advancing the state.

Run the diagnostic-parser unit test independently:

```sh
python3 -m unittest discover -s TESTS -p 'test_compare_finh_timestep.py'
```

This diagnostic script does not dispatch any GitHub Action. Development
workflows remain manually selectable; the separate stable-release pipeline
executes source/numerical integration checks automatically.

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
