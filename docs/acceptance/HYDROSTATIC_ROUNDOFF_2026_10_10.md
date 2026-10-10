# Stationary hydrostatic roundoff correction

This records the earlier stationary-only milestone. The subsequent
[geometry and transport package](N7_B_GEOMETRY_AND_TRANSPORT_2026_10_10.md)
fixes the y-directed geometry defect and closes N7-B with the full isolated
contact matrix. The original results and scope below retain their provenance.

The variable-composition hydrostatic equilibrium subset now passes: all 90
original equilibria in both strict and optimized builds retain exactly zero
state changes, component/carrier budget errors and force-scaled momentum
residuals. N7-B remains open for its moving-contact test. CFL policy, physical
laws, reconstruction and the clean reference kernels are unchanged.

The [machine-readable evidence](results_hydrostatic_roundoff_2026_10_10.json)
pins the source and numerical fixtures, contracts, focused results and full
audit. The base revision is `e4333a5cf0575d23287a32a082fe18fbd939bcbe`.
Only `src/hyperbolic_2d.f90` changes in production.

## Correction and exclusions

In the preceding equilibrium investigation, cancellation roundoff produced a
velocity of about `1e-17 m/s` and a nonzero discharge. At the next stage that discharge
disabled the existing exact-rest scalar-flux guard, allowing finite CU diffusion
of composition and thermal energy. The [historical diagnosis](N7_B_DIAGNOSIS_2026_10_10.md)
retains the failing inputs and observations.

After normal spatial assembly, the correction requires both momenta of the
cell and every adjacent face endpoint to equal zero exactly. Only then may a
normal-momentum residual be zeroed, using the fixed arithmetic allowance

```text
abs(residual) <= 64 * epsilon(wp) * scale / directional_spacing
scale = max over adjacent directional endpoints of
        abs(G * Gamma) * h * max(h, abs(eta))
```

The pressure/path scale accounts for the subtraction scale of the absolute
free surface; a dry endpoint contributes zero. This is a tested arithmetic
allowance, not a universal certified error bound or an equilibrium theorem for
every thermodynamic/topographic configuration. Its multiplier was frozen before
execution and was not tuned after failures.

There is no velocity epsilon, absolute force floor, scalar-flux modification or
timestep modification. Any nonzero participating momentum bypasses the
correction. A genuinely nonzero pressure force above the allowance remains
active even when initially at rest. The helpers use automatic local scalars,
read-only reconstruction/thermodynamic data and the existing cell-owned output;
there are no new workspace allocations or shared scratch writes.

## Verified tests

Both compiler profiles use actual one/four-thread teams. The focused inventory
has 648 solver runs, counted separately from the general 38-script audit.

| Inventory per profile | Result |
| --- | --- |
| 90 original composition equilibria | Exactly zero equilibrium, force-residual, boundary-mass-flux and component/carrier budget errors; all original raw-stage, thermal, positivity and fraction bounds pass. |
| 36 moving/force controls | Complete `result.bin` and `composition.bin` match the old hyperbolic operator bitwise, and match between actual one/four-thread teams. |
| Existing 102 N7-A cases | Original equilibrium, transport, Ritter and excavation assertions pass. |
| Existing IMEX-stage suite | Supported stage counts and explicit/implicit/mixed assertions pass. |

The equilibria cover liquid, gas and gas-liquid closures, each with two solids;
the gas closures include an added gas. They include constant-density Q1-bed
states at constant free surface and flat-bed states with varying `Gamma` and
depth at constant `Gamma*h**2`. The original 1D/2D, slope/curvature and RK2/3/4
matrix and its `32768*epsilon(float64)` criteria were not changed.

The controls resolve both directions on flat 20-by-5 or 5-by-20 grids. Normal
velocities are signed `1e-20` and `1e-12 m/s`; tangential velocity is `1e-20 m/s`.
Initially resting pressure-imbalance controls accelerate to at least
`7.82687e-9 m/s`, above the frozen `1e-10 m/s` minimum. The old hyperbolic module
is a checksum-pinned test-only snapshot linked with the same other current
objects and observer. The main build never compiles it; this isolates the
modified operator rather than introducing a production facade.

The complete acceptance run passes **38/38** scripts. All six full
Gate-H/deterministic-restart/stochastic-restart cross-thread comparisons pass,
with original restart durations and existing numerical limits. It does not
approve N9. Source and numerical fixture hashes refer to that executed snapshot;
documentation and closure metadata were reconciled afterward.
The final metadata/checker checks also pass: 13 acceptance-tool unit tests and
five composition-checker unit tests. They do not add solver runs to the archived
38-script total.

## Comparison with the preceding full audit

The comparison uses the full captured manifests, not their compact test-summary
records. All 16 paired full-solver entries have identical effective input bytes.
Ten have identical canonical `.q_2d` and `restart.bin` fingerprints: both
profiles of the lakes and deterministic/stochastic restarts. Six differ: both
thread configurations of Gate H and the embedded-thread inclined excavation,
in both profiles. Initial fields agree; some final fields/checkpoints differ.

The current Gate-H metrics, symmetry, positivity and zero-uphill criteria still
pass. Fingerprint differences alone do not quantify a field norm; no new
scientific baseline or blanket output-neutrality claim follows from them.
Long ETNA and lava runs were not repeated here and should be repeated before
adopting the corrected scientific baseline.

## Failed test constructions retained

The archived control contracts `v1/v2/v3` retain why their original fixtures
could not isolate the intended predicate. Inactive 1D momentum handling erased
tangential velocity; the old one-column y geometry fitter then failed bounds
checks, even with a flat bed. Contract v4 resolves both directions with at least
five transverse cells for the retained fit stencil. Speeds, pressure criteria,
the production predicate and the 64-epsilon multiplier were not relaxed.

The first general audit failed only because the new test runner required `.git`
inside an exported source tree. The checksum-pinned operator reference removed
that infrastructure dependency; the final exported-source audit passes. The
evidence index retains the earlier audit path and explanation.

## Remaining work

H04 is satisfied for the defined stationary composition matrix. N7-B still
requires an independently valid moving-contact fixture/reference: the original
velocity taper influences its purported uniform-velocity comparison interval.
The failing draft remains archived and is not counted as a passing gate.

The pre-existing out-of-bounds x-slope extrapolation in `geometry_2d.f90` for
`nx=1, ny>1` is not fixed here. Address it as a separate bounded change before
using rotated one-dimensional y fixtures. N7-C, N7-D, D-N7-CFL and N9 remain
open; this correction does not authorize layer or implicit-group development.

## Reproduction

From the repository root:

```sh
KEEP_TEST_WORKDIR=1 sh TESTS/TEST_HYDROSTATIC_ROUNDOFF/run_test.sh
PYTHONDONTWRITEBYTECODE=1 python3 TESTS/TEST_N7_COMPOSITION/test_checker.py
```

To run only the original equilibrium subset:

```sh
KEEP_TEST_WORKDIR=1 sh TESTS/TEST_N7_COMPOSITION/run_test.sh --equilibrium-only
```

Without that option, the composition script attempts the unresolved contact
draft and must not be used as a passing acceptance gate. See the respective
[roundoff regression](../../TESTS/TEST_HYDROSTATIC_ROUNDOFF/README.md) and
[composition diagnostic](../../TESTS/TEST_N7_COMPOSITION/README.md) documentation.
