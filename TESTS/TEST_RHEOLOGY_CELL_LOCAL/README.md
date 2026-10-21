# D02 cell-local source characterization

This suite measures production sources against independent references while
the rheology laws, pore-pressure terms and general thermal sources are being
extracted into dedicated modules. It also checks canonical input-name
resolution and reverse mapping.
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
Fixtures run in isolated processes, with frozen stochastic Z and no RNG draws.
The dedicated `parallel-voellmy-workspace` fixture repeats the same four frozen
forcing values in a verified four-thread OpenMP team and compares each result
with its serial reference. Four optimized SI rest probes are repeated without enabling
traps in the main program, using the identical production module objects. They
are diagnostics, not replacements for the trapping checks.

The independent NumPy reference builds conservative states from prescribed
phase volume inventories, then decodes and evaluates the equations. It never
uses production primitive states as oracle inputs. The parser uses the names
`voellmy`, `plastic`, `lava`, `lahar`, `empirical`, `quadratic`, `froude`,
`lucas`, `mui`, and `mui_voellmy`; `none` denotes disabled rheology. The
historical model 8 is no longer supported.
Pore-pressure diffusion, all inhibition modes, mono/multi-class permeability,
Sutherland viscosity and conservative Darcy gas loss have dedicated fixtures.
The source-module checks independently exercise inlet heat, mass-exchange heat,
vaporization, conservative pore-pressure inlet/mass-exchange terms, and the
D09 friction-only momentum projection (alignment, stopping, and no magnitude
increase). The exact lava-cooling checks exercise convection-only closed-form
decay, radiative-convective flow composition, sub-ambient approach to ambient,
and the initial ODE derivative.

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
assembly remains private; full stage solves provide a separate
integration check. Derivatives of explicit mass/composition columns are
identified separately from columns selected by the current Newton mask.

At zero pressure the pressure derivative is checked with a positive one-sided
FD, matching the implemented `P>=0` guard. At exact rest, smooth flat-bed viscous laws 3/4/5 also receive derivative
checks. At nonsmooth rest, exact dry states and named depth/temperature/velocity
switching points, value/limit probes replace a claim of a unique derivative. The exact STEP packing switch accepts either
one-sided value because independent mass-to-volume rounding can change the
selected side; it never accepts an intermediate value. Points on either side
are checked normally.

Newton residuals are checked against the production criterion
`1e-5 + 1e-5*abs(initial_residual)` for every component, independently of the
returned `converged` flag. The quadratic stage also has an independent scalar
root oracle; its state tolerance is `1e-5`, matching production solve accuracy.
A large-state stagnation fixture requires a failed convergence flag when the
independently recomputed residual exceeds its bound. It tests accurate failure
reporting, rather than accepting an inaccurate stage solution.

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
The 32 cases cover accepted Voellmy and Bursik--Woods inputs, case-insensitive
name handling, rejection of unknown names and legacy numeric IDs, missing/
zero/negative denominators, valid zero Coulomb coefficients, nonfinite
coefficients and ambient cooling configuration. Generated fixtures allow
overwrite only inside their supplied output directory. Lava uses `h_x`; other
cases use `froude_blend`.

The parser requires finite active rheology coefficients, positive `xi`, `Fr_0`
and `U_w`, and nonnegative Coulomb coefficients. Integer exponent inspection
validates IEEE binary32/64 inputs under `-Ofast` without NaN comparisons.
Cooling uses the configured positive finite `T_env`, defaulting to 300 K when
omitted. The test reads its startup diagnostic because the historical backup
writer omits the temperature namelist.

## Acceptance boundary

The current suite contains 255 fixtures per profile; the three legacy model-8
ODE probes were removed with D04. The original 249-fixture
characterization and its findings remain historical evidence in
`docs/acceptance/D02_CELL_LOCAL_CHARACTERIZATION_2026_10_21.md`; follow-up fixes,
contracts and validation are recorded in
`docs/acceptance/D02_NUMERICAL_FIXES_2026_10_21.md`.

Pressure diffusion contributes zero for dry states, negative excess pressure or
absent carrier pore volume. Dynamic permeability additionally requires a solid
skeleton. Guards suppress only pressure relaxation and retain other sources.
Smooth complex algebra preserves composition and height sensitivities, including
inhibition transitions. The existing diffusion sphericity weights remain
unchanged; D06 preserves and tests them separately from the Darcy gas-loss
weights. Any harmonization is a separate scientific decision, not an implicit
consequence of extraction.

D06 places general inlet and mass-exchange heat/vaporization in
`thermal_sources_2d`, including lava cooling, while `rheology_lava_2d` owns
the lava rheological terms. `pore_pressure_2d` owns pressure diffusion, Darcy
gas loss and the conservative pressure-source couplings. The diffusion (`alpha*d^2`) and
gas-loss (`alpha/d`) sphericity providers remain distinct and are tested
against their separate independent formulas; this extraction does not silently
change either scientific convention.

This suite remains a standalone check, distributed with the project. Scientific
formula changes, ownership migration, model-8 removal and broader compiler or
platform acceptance belong to separately validated follow-ups. No test tolerance
has been increased to close the initial findings.
