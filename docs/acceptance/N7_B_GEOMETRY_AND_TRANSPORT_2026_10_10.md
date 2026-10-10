# One dimensional geometry and variable composition transport

N7-B is satisfied in its defined equilibrium and material-contact fixtures.
Both compiler profiles pass the full 162-case inventory with actual one/four
OpenMP teams. The production change is confined to `geometry_2d.f90`: correct
the missing y-directed 1D fit and safely handle short lines and a single cell.
Fluxes, thermodynamic laws, IMEX tableaux and CFL policy are unchanged.

The [evidence index](results_n7_b_2026_10_10.json) records frozen contracts,
source/test hashes, independent comparisons, retained failures and the complete
acceptance audit. The base revision is `a8148e4e`, after the stationary-roundoff
correction. N7 remains open for N7-C, N7-D and D-N7-CFL; N8 stays satisfied and
N9 approval remains separate.

## Geometry correction

The preceding module reads the nonexistent third x column on a one-column
y grid; its strict build reproduces an index-out-of-bounds error. The corrected
module filters along the resolved y direction using the same five-point
polynomial-fit coefficients and volatile scalar accumulation as the x kernel.
It extrapolates only along an axis with a valid interior stencil. All inactive
and mixed 1D derivatives are explicitly refreshed to zero.

Lines of at least five cells retain the original filter and nearest-interior
boundary values. For one through four cells, the available centers determine a
least-squares polynomial of degree `min(2,n-1)`: a single cell has zero resolved
derivatives, two resolve a line, and three/four resolve a quadratic. This fallback
does not filter or alter the authoritative vertices or their Q1 faces/centers.
The helper is private and has no shared scratch or persistent workspace.

The dedicated regression runs 112 cases per profile, or 448 solver invocations:
lengths 1, 2, 3, 4, 5, 9 and 17 along both axes, the common single-cell grid once,
and a 7-by-9 control, each with four beds and slope correction OFF/ON. It poisons
every cache before refresh and independently checks Q1 geometry, derivatives
and gravity coefficients. Maximum relative field errors are below `2.6e-14`,
within the frozen `2048*epsilon(float64)` allowance. All 224 profile/case
one/four-thread pairs have identical complete snapshots.

The 32 preceding valid x/2D controls per profile were captured before the source
change. Every candidate snapshot matches its preceding counterpart bitwise.
This test does not extend the production contract to arbitrary short 2D grids
with both axes resolved.

## Material contact reference and isolation

The scalar reference is valid in a uniform-velocity, constant-depth region
because fractions vary in the specific-volume nullspace: density and the
hydrostatic coefficient remain constant at uniform temperature. Independent
specific-volume and heat-capacity sums lift the scalar profile into every
canonical conservative component, including thermal energy. The scalar CU/TVD
operator uses independently specified current explicit IMEX weights, not a
production flux call or a one-solid gas surrogate.

The original contact's taper reached its purported uniform-velocity comparison
interval. That failed construction remains in `contract.json` and its historical
diagnosis. The v2 extension isolated the central comparison but left insufficient
quiet buffers at physical boundaries; its nonzero boundary flux failed the
unchanged criterion and remains archived.

The frozen v3 fixture translates the same pulse and complete comparison interval
into a 480 m domain. Its grids have 480/960/1920 longitudinal cells along x and
960 along y. The original spacings, timestep rule, final time, limiter, pulse
width, 16 m comparison width and all error limits are unchanged. No comparison
cells were removed and no tolerance was enlarged.

Before execution, both the central comparison and physical boundaries must lie
outside a conservative local dependency cone of `4*n_RK*steps` cells. The
prescribed timestep avoids global timestep coupling. Minimum remaining guard
widths are 37 cells centrally and 29 at the boundaries. Full-domain raw-stage
budgets, positivity, thermodynamic/fraction closure and boundary flux checks
remain mandatory; only the independent uniform-velocity field comparison uses
the isolated interior. This is a finite-time isolated transport test, not a
general certification of every boundary-condition policy.

## Full N7-B results

Each profile covers three closures: liquid/two-solid, gas/two-solid/added-gas
and gas/liquid/two-solid/added-gas. All supported stage counts are tested:
`N_RK=2,3,4`; one stage remains inadmissible.

| Inventory per profile | Verified result |
| --- | --- |
| 90 original variable-composition equilibria | Original states, flags, times and limits unchanged; exactly zero state changes and force-scaled residuals. |
| 72 moving contacts | Both directions and axes pass independent conservative fields, analytic transport/refinement and whole-domain component/carrier budgets. |
| Raw stages and reconstructed faces | Positive mass/components/carrier, thermal closure and TVD fractions pass before final repair. |
| Actual one/four-thread teams | Complete base and composition payloads agree bitwise in all 324 profile/case pairs. |

Maximum component-scaled contact reference errors are `1.072e-15` in strict and
`8.939e-16` in optimized builds, within the unchanged step-scaled roundoff limit.
Maximum analytic profile errors are below `0.0078` in L1 and `0.046` in Linf;
the largest refinement ratio is below `0.304`, against the frozen `0.9` limit.
Each conservative component and residual carrier retains its own raw-stage and
final budget; the largest inventory-scaled error is below `6.55e-15`. Boundary
mass flux is exactly zero. Signed centroid motion agrees with the prescribed
direction.

The focused N7-B matrix has 648 solver invocations. Together with the new
geometry regression it has 1096 invocations and 548 bitwise thread pairs;
these counts do not include repeats inside the general acceptance audit.

## Closure and next batch

The complete general acceptance audit passes **40/40** scripts, including
original-duration deterministic/stochastic continuation and all six full
one/four-thread comparisons. All 16 preceding full-case canonical output and
checkpoint fingerprint groups are identical, with unchanged effective input
bytes. The numerical source/test snapshot is pinned in the evidence index;
documentation and closure metadata were reconciled after execution. Final
metadata checks pass 14 acceptance-tool unit tests and six composition-checker
unit tests; these do not add solver runs to the 40-script total.

H04 and H05 are satisfied. N7-05 closes; N7-03 and N7-07 retain only their
N7-C requirements. Existing restart and thread contracts remain active. The
historical stationary correction report retains its original partial status
and links to this subsequent completion rather than changing archived evidence.

The next user-approved batch combines N7-C and N7-D: finite-time rough/front,
oblique-ramp, circular symmetry and dynamic-curvature checks together with the
reproducible Gate-H thickness/free-surface comparator. D-N7-CFL remains a
separate explicit decision; neither this batch nor a passing script approves
N9 or new golden scientific outputs. Long ETNA/lava acceptance remains in N9.

## Reproduction

From the repository root:

```sh
KEEP_TEST_WORKDIR=1 sh TESTS/TEST_GEOMETRY_1D/run_test.sh
KEEP_TEST_WORKDIR=1 sh TESTS/TEST_N7_COMPOSITION/run_test.sh
PYTHONDONTWRITEBYTECODE=1 python3 TESTS/TEST_N7_COMPOSITION/test_checker.py
```

`--equilibrium-only` deliberately omits moving contacts and cannot close N7-B.
`--diagnose` retains the original historical contact construction and is not an
acceptance gate. The [geometry](../../TESTS/TEST_GEOMETRY_1D/README.md) and
[composition](../../TESTS/TEST_N7_COMPOSITION/README.md) test READMEs specify the
contracts and retained binary observation formats.
