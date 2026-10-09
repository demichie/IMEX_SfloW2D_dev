# Froude friction correction and stochastic restart validation

The strict stochastic restart failure is resolved by restricting rheology model
9 to positive reduced gravity. Both clean production builds and all 33 audit
test runs pass. All 14 previously completed production cases retain identical
canonical conservative-output and checkpoint SHA-256 fingerprints. N9 remains
open for the other scientific coverage requirements.

## Cause of the strict failure

The unmodified executable stops in `eval_nh_semi_impl_terms` at
`Fr = mod_vel / SQRT(r_red_grav * r_h)`. On this macOS arm64 toolchain the
floating-point trap is reported as `SIGILL`. The failing operation is in the
semi-implicit REAL rheology evaluation, not checkpoint decoding or the OU random
generator.

A diagnostic build reports a transported gas tail with thickness
`3.3491779338939385e-12 m`, temperature `300 K`, no solids, and mixture density
equal to ambient density (`1.1763298740177413 kg/m3`). Its horizontal speed is
`3.7378768273725189e-6 m/s`. Reduced gravity is exactly zero, so a moving wet
state enters the old branch and divides by zero. Dry-cell detection alone does
not prevent this: the thickness exceeds machine epsilon.

The [dedicated regression](../../TESTS/TEST_RHEOLOGY_FROUDE/run_test.sh) also
fails with `SIGILL` against the unmodified numerical sources, without requiring
OU evolution or restart. Successful optimized execution did not make this
intermediate arithmetic valid under the strict floating-point contract.

## Correction and physical domain

The production change is confined to
[`equation_terms_2d.f90`](../../src/equation_terms_2d.f90): model 9 evaluates its
Froude-dependent basal friction only when both speed and reduced gravity are
positive. Neutral mixtures have no buoyancy load in this law; buoyant mixtures
are outside its compressive basal-contact domain. Both retain the initialized
zero friction source instead of evaluating a zero or negative Froude
denominator.

The positive-gravity formula, stochastic Froude shift and clipping, slope
factor, state conversion and other rheologies are unchanged. No gravity floor,
dry threshold, noise amplitude, timestep, reference tolerance or floating-point
trap has been relaxed. The REAL/COMPLEX implicit API is unchanged; model 9 has
no implicit friction contribution there.

## Verification results

| Check | Result |
| --- | --- |
| New regression on unmodified numerical sources | Fails with `SIGILL`, exit 132. |
| New regression on corrected sources | Passes with strict and optimized numerical-module builds, trapping enabled in both. |
| Clean strict and historical optimized production builds | Both pass. |
| Existing unit/documentation scripts plus the new regression | 17 scripts pass; each retains its own compilation policy. |
| Full production-executable audit | All 16 runs pass, eight per profile. |
| Deterministic and stochastic continuation | Exact checkpoints and final conservative snapshots in both profiles, with 1 and 4 actual threads. |
| Cross-thread Gate-H and restart comparisons | All six completed pairs match exactly. |
| Comparison to the previous numerical baseline | All 14 previously completed cases have identical canonical q/checkpoint hashes. |
| Audit-tool self-tests | Seven pass. |

The new regression covers moving neutral gas, the measured tiny transported
tail, warm buoyant gas, dry cells, dense cells at rest, and moving dense cells
with/without slope correction. It checks deterministic friction, positive
stochastic shifts and negative shifts clipped to zero, and confirms finite zero
implicit sources through both REAL and complex-step APIs. Its assertion driver
is compiled strictly even when numerical modules use `-Ofast`, so NaN/Inf
checks cannot be optimized away by the driver.

The unchanged-output comparison includes both optimized stochastic runs. The
old strict stochastic trajectories failed and therefore cannot supply an
old/new full-trajectory comparison; their replacement evidence is completion,
exact restart and exact one/four-thread agreement.

## Candidate provenance and reproduction

The tested candidate is an uncommitted working-tree snapshot based on
`5051d017ab5ecf9010e29ce80e12cdf70b294d39`, not an executable built from that
commit alone. Its numerical sources differ from the previous `752b47fe`
baseline only in `src/equation_terms_2d.f90`. The candidate file SHA-256 is
`f5d4cd0512fe00db79c8774b785c493a474d951da4766a2627dfc9fcf1ff7625`.

[`results_froude_fix_2026_10_09.json`](results_froude_fix_2026_10_09.json)
retains candidate source/test/tool hashes, compiler flags, input/output
fingerprints, actual thread counts and the old/new comparison. Full build and
invocation logs remain in
`/tmp/imex-stochastic-fix.zOOJ2c/acceptance-after-final`. Observed output hashes
are not approved golden references, and not all raw field payloads are archived.

To test an uncommitted candidate, choose a new output directory:

```sh
python3 TESTS/ACCEPTANCE/run_acceptance.py \
  --working-tree --revision HEAD \
  --output /tmp/imex-froude-repeat \
  --netcdf /opt/homebrew/opt/netcdf-fortran --jobs 3
```

After committing, use the resulting commit SHA and omit `--working-tree` for a
tracked-source reproduction. The audit tooling now records a candidate
snapshot explicitly and retains an archive file before extraction, avoiding the
observed BSD-tar pipe termination (`tar` succeeds but Git receives `SIGPIPE`).
Neither tooling change alters a solver input's physical parameters.

## Remaining acceptance work

The strict stochastic failure gate is closed for this candidate. The
[original N7/N8/N9 matrix](N7_N8_N9_AUDIT_2026_10_09.md) remains the pinned record
of the earlier baseline, not a current failure claim. Next come isolated
production inventory-limit and evolving-bed geometry tests, followed by
production residual/CFL and mixed-composition comparisons. Long-run provenance,
source/output contracts and explicit scientific sign-off remain separate N9
conditions.
