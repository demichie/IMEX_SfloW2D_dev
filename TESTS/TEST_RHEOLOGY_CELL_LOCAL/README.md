# D02 cell-local rheology characterization

This suite measures production sources before extraction into rheology modules.
It deliberately reports discrepancies as failures. A completed run is not
automatically a scientific acceptance: consult the JSON `acceptance` field and
the exit status. There are no expected-failure allowances or tolerance increases
for known production defects.

Requirements: GNU Fortran, LAPACK, Python 3 and NumPy. The integration/parser
checks additionally require a configured GNU out-of-source solver build and
NetCDF-Fortran (`nf-config` and `nc-config`). No GitHub Actions are launched.

## T1 values, derivatives and isolated hazards

From the repository root:

```sh
D02_WORK_DIR="$PWD/build-local/d02" sh TESTS/TEST_RHEOLOGY_CELL_LOCAL/run_test.sh
```

`D02_WORK_DIR` is optional; the default is a retained `/tmp/imex-d02.*`
directory. All objects, modules, executables, fixture namelists and logs remain
in the work directory. Each profile writes `results.json` and `cases/*/probe.log`.
The wrapper runs both profiles even when characterization fails, then exits
nonzero if either has failures. It does not modify production inputs or source.

The profiles are strict (`-O0 -g -fcheck=all -finit-real=snan`) and optimized
(`-Ofast -g -funroll-all-loops`). Both enable OpenMP and
`-ffpe-trap=invalid,zero,overflow`; driver-side finite differences stay at `-O0`.
Fixtures run serially in isolated processes, with frozen stochastic Z and no
RNG draws. Four optimized SI rest probes are repeated without enabling traps
in the main program, using the identical production module objects. They are
diagnostics, not replacements for the trapping checks.

The independent NumPy reference builds conservative states from prescribed
phase volume inventories, then decodes and evaluates the equations. It never
uses production primitive states as oracle inputs. Models 1–6 and 9–12 are
covered, with disabled rheology and a separate legacy model-8 ODE check.
Pore-pressure diffusion, all inhibition modes, mono/multi-class permeability,
Sutherland viscosity and conservative Darcy gas loss have dedicated fixtures.

For regular states, complex-step uses `1e-25*max(abs(q_j),1)` in physical
conservative coordinates. Central finite differences use relative steps
`1e-5,1e-6,1e-7,1e-8` with the same scale. At least two successive steps must
agree with CS. The column infinity-norm relative tolerance is `3e-5`, with a
floor `max(1e-12,1e-5*max(abs(S))/max(abs(q_j),1))` to handle near-zero
columns. Value checks use `2e-10` relative to `max(1,abs(reference))`; gas volume
loss instead uses a `1e-12` denominator floor. Every measured error is retained.

A second CS probe uses the actual Newton normalization:
`delta=n_vars*epsilon(wp)`, `q_org=max(abs(q),1e-3)`. The scaled column must
agree with the physical column times q_org to `2e-10`. The Newton Jacobian
assembly remains private and unchanged; full stage solves provide a separate
integration check. Derivatives of explicit mass/composition columns are
identified separately from columns selected by the current Newton mask.

At zero pressure the pressure derivative is checked with a positive one-sided
FD, matching the implemented `P>=0` guard. At rest, exact dry states and named
depth/temperature/velocity switching points, value/limit probes replace a
claim of a unique derivative. The exact STEP packing switch accepts either
one-sided value because independent mass-to-volume rounding can change the
selected side; it never accepts an intermediate value. Points on either side
are checked normally.

Newton residuals are checked against the production criterion
`1e-5 + 1e-5*abs(initial_residual)` for every component, independently of the
returned `converged` flag. The quadratic stage also has an independent scalar
root oracle; its state tolerance is `1e-5`, matching production solve accuracy.

## T2 actual integration and parser checks

```sh
python3 -B TESTS/TEST_RHEOLOGY_CELL_LOCAL/check_stages.py \
  build-local/optimized build-local/d02/stages-optimized
python3 -B TESTS/TEST_RHEOLOGY_CELL_LOCAL/check_stages.py \
  build-local/strict-debug build-local/d02/stages-strict
python3 -B TESTS/TEST_RHEOLOGY_CELL_LOCAL/check_inputs.py \
  build-local/optimized/src/IMEX_SfloW2D build-local/d02/inputs
```

`check_stages.py` updates the supplied configured build with `make -j2`, then
links the existing IMEX observer against its actual production objects.
It exercises one step with `n_RK=2` at 1 and 4 threads on uniform flat fixtures,
checking the central cell against quadratic DIRK, Coulomb stopping and coupled
lava cooling/drag oracles. It checks known stage assembly, final assembly,
effective increment after corrections and unchanged final projection. The
coupled thermal root is obtained by independent bisection. This does not import
or certify the historical Python spatial solvers.

The observer's optional eighth argument selects model 1 or 3 for D02 only;
omission preserves its existing behavior. Lava selects `h_x`; other fixtures
select `froude_blend`. No production tableau or rheology selector is changed.

`check_inputs.py` uses a fresh source-mass fixture and the real solver parser.
It distinguishes input rejection text from the process exit code: plain Fortran
STOP can return zero. Valid control runs must reach the end of a short run.
Five invalid parameter configurations are currently accepted and are reported
as failures, not interpreted as successful negative tests.

## Acceptance boundary

The initial characterization finds unresolved derivative, domain, configuration
and convergence issues; **D02 is not accepted**. See
`docs/acceptance/D02_CELL_LOCAL_CHARACTERIZATION_2026_10_21.md` for evidence,
scope and the exact work required to close the gate. This suite is distributed
as a standalone check, not added to passing `make check`/CI gates while those
production findings remain unresolved. Scientific formula changes, ownership
migration and model-8 removal belong to separately validated follow-ups.
