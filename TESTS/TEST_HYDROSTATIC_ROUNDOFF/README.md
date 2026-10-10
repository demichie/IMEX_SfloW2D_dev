# Stationary hydrostatic roundoff regression

This verifies the bounded stationary correction in `hyperbolic_2d`, not the
entire N7-B contact/transport gate. Run from the repository root:

```sh
KEEP_TEST_WORKDIR=1 sh TESTS/TEST_HYDROSTATIC_ROUNDOFF/run_test.sh
```

Both strict and historical optimized profiles run 90 original composition
equilibria with actual one/four-thread teams. They retain the original N7-B
contract and all raw budgets, positivity, thermal closure, TVD and force criteria.
The separate 36-case control inventory includes signed normal speeds of
`1e-20` and `1e-12 m/s`, tangential speed `1e-20 m/s`, and pressure imbalances
that must produce nonzero acceleration. Its full raw payloads must be
bitwise identical to the pre-correction hyperbolic operator, as well as across
thread counts.

`reference/hyperbolic_2d.f90` is the exact module from commit `e4333a5c`, with
SHA-256 `1d5af9f9d5e1936acd1f37003b448f2cd94005de82ad573e505d7b7ff9342ee3`.
It is a test-only operator snapshot, not a production facade. Both operators
link against the same other current module objects and observer, isolating the
changed hyperbolic module. The main build never compiles this reference. The
runner verifies its checksum and also works in a source export without `.git`.

The [guard contract](../TEST_N7_COMPOSITION/roundoff_guard_contract.json) fixes
the predicate, 64-epsilon pressure allowance, controls and byte-equality criteria.
Earlier control constructions are retained as `roundoff_guard_contract_v1/v2/v3.json`:
1D inactive-direction handling zeroed tangential motion, and the old one-column
y geometry fitter went out of bounds even for a flat bed. Corrected controls
resolve both axes on flat 20-by-5 or 5-by-20 grids; velocities, force criteria
and production conditions were not weakened. The original 90-case Q1/constant
pressure equilibrium matrix was never changed.

The retained build contains all inputs, complete binary payloads and
`roundoff_guard_evidence.json` / `equilibrium_evidence.json` for each profile.
The [correction report](../../docs/acceptance/HYDROSTATIC_ROUNDOFF_2026_10_10.md)
records results and remaining work.
