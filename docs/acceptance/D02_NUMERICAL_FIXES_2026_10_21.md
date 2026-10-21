# D02 numerical fixes and acceptance follow-up — 21 October 2026

Baseline: `3e07c6aa57273a3564bc25b53cdc3e437e38d04d` (initial D02 characterization).
The [initial report](D02_CELL_LOCAL_CHARACTERIZATION_2026_10_21.md) and its
failure counts remain historical evidence. This follow-up records production
corrections and the additional regression checks needed before D03.

Status: **LOCAL ACCEPTANCE COMPLETE**. Focused checks and the full local regression matrix passed.
Machine-readable evidence: [results](results_d02_numerical_fixes_2026_10_21.json).
Reproduction commands: [suite README](../../TESTS/TEST_RHEOLOGY_CELL_LOCAL/README.md).

## Corrections and explicit contracts

1. **Complex source sensitivities.** Dynamic permeability, Sauter diameter,
   diffusion sphericity and smooth inhibition now retain complex composition
   and height perturbations. Only branch/clamp selection uses real base states.
   Complex `SUM` products avoid `DOT_PRODUCT` conjugation. Integer polynomial
   exponents avoid an unnecessary logarithm at zero. Flat-bed viscous laws 3/4/5
   retain their smooth velocity derivatives at exact rest. These checks do not
   assert a unique derivative of Coulomb friction or a switching surface.
2. **Pressure diffusion domain.** Decode the mixture when pore pressure is
   enabled independently of rheology. Dry states, negative excess pressure and
   absent carrier pore volume contribute zero pressure relaxation. Dynamic
   permeability additionally requires a solid skeleton. A named local block
   exits before singular coefficient evaluation and preserves all other
   assembled source terms. This is an explicit no-pore-network policy, without
   a fitted permeability floor. Diffusion's existing `alpha*d**2` sphericity
   weights remain unchanged; reconciling them with Darcy gas-loss weights
   `alpha/d` remains D06 work.
3. **Input validation.** Reject nonfinite active rheology coefficients before
   floating-point comparisons. Require positive denominator parameters `xi`,
   `Fr_0` and `U_w`, while allowing zero nonnegative Coulomb coefficients. A
   private IEEE binary32/binary64 exponent-bit check survives GNU `-Ofast`
   finite-math assumptions, which can fold `IEEE_IS_FINITE` to true. Unsupported
   representations fail closed. Invalid inputs are tested through the real
   parser, including Fortran `STOP` implementations that return exit status 0.
4. **Cooling configuration.** Remove the source kernel's write of `T_env=300`.
   Input initialization supplies 300 K when omitted, while a supplied positive
   finite value is honored. The startup log exposes the effective value; the
   historical backup namelist writer is unchanged. This intentionally changes
   the thermal source when a user had configured an ambient other than 300 K.
5. **Real SI rest.** Return zero at exact zero horizontal velocity before norm
   and division evaluation. This prevents optimized trapped builds from
   evaluating a reciprocal at rest; moving-state formulas are unchanged.
6. **Newton convergence.** A small correction or roundoff stagnation cannot
   substitute for the residual criterion. Form the conservative difference in
   a local `VOLATILE` array before adding forcing, so `-Ofast` reassociation
   cannot erase a small residual between equal large states. The added extreme
   state regression requires failure reporting for an unrepresentable update;
   valid stage solves must still satisfy the original residual bounds. The
   legacy model-8 ODE bypass is preserved pending D04.

No tracked input file, tableau, solver tolerance or expected-failure allowance
was changed. Builds and generated fixtures stay outside `src`.

## Focused evidence

GNU Fortran (Homebrew GCC 15.2.0) on macOS ARM64; production sources are compiled
with strict `-O0 -g -fcheck=all -finit-real=snan` and optimized
`-Ofast -g -funroll-all-loops` profiles. Both cell-local profiles enable
`-ffpe-trap=invalid,zero,overflow`. Finite-difference driver arithmetic stays
strict. The test suite contains **257 fixtures per profile**, up from 249:
pressure-only/empty/all-solid domains, three smooth viscous rest derivatives,
and one Newton stagnation regression add eight fixtures. Independent NumPy
oracles and tolerances are retained. Each profile checks 793 derivative columns;
the maximum CS/FD error is `1.306e-5`, below the unchanged `3e-5` threshold.
The extreme stagnation case has an independently nonzero residual (1092) and
correctly returns a failed convergence flag in both profiles. See the JSON
evidence for exact hashes.

The real parser has **28 cases per build**, including valid zero coefficients,
missing/zero denominators, NaN/Infinity and default/explicit/invalid ambient
cooling temperature. Actual `n_RK=2` integration checks cover quadratic DIRK,
Coulomb stopping and coupled lava cooling/drag at 1 and 4 threads, **six cases
per build**. Lava fixtures use `h_x`, other fixtures `froude_blend`.

Before/after comparison of common finite `REAL` and `SI` source records finds
435 unchanged records in strict and 427 in optimized, within scaled differences
`6.53e-12` and `8.16e-12`, respectively. The sole changed finite record is the
intentional ambient-temperature correction (`ambient-overwrite`, about 0.16743
scaled difference). Nonfinite/trapped historical records cannot establish
value equivalence; those cases instead use the explicit domain contract and
independent oracle. Raw comparison: `build-local/d02-fixes/value-equivalence.json`.

## Local regression matrix and acceptance decision

Completed: **26/26 commands passed** (13 per profile). For each profile, the
build, parser and actual-stage checks passed; `make check`, binary restart,
stochastic restart, Gate H and source-mass passed at both 1 and 4 threads.
Detailed run records: `build-local/d02-fixes/validation-results.json`; raw logs:
`build-local/d02-fixes/validation-logs/`. No GitHub Actions were run.

D02 acceptance is limited to this local source characterization and regression
matrix. It does not certify ifx/Linux/Docker, global thread safety, the old
Python spatial cores or a broad physical calibration campaign. Shared
`xi_temp` ownership remains D03 work; model-8 removal remains D04; differing
permeability providers remain D06; final cross-model scientific acceptance
remains D15. D03 is the next roadmap step after this gate closes.
