# Nodal projection and geometry acceptance completion

N8-B completes the original N8 refinement and repeated-update requirements.
N8 is now satisfied within the supported uniform Cartesian production contract.
The user accepted the separate unequal-area algebra test on October 10, 2026;
this does not enable nonuniform production meshes. N7 and N9 remain open.

The [projection evidence](results_n8_b_projection_2026_10_10.json) records fixed
fixtures, refinement metrics, stress envelopes, source hashes and all-team
payload fingerprints. The [full regression](results_n8_b_acceptance_2026_10_10.json)
passes 36/36 runs and six cross-thread Gate-H/restart pairs. All 16 full-solver
cases retain the preceding canonical output/checkpoint fingerprints, and all
16 N8-A payload fingerprints remain unchanged.

## Fixtures and acceptance limits

The new driver runs alongside the unchanged N8-A driver in `TEST_MASS_EXCHANGE`.
Both strict and optimized production builds use actual one/four-thread teams;
assertion arithmetic is compiled strict. The grids are 16/32/64 cells in 1D and
16x12/32x24/64x48 in 2D, on a fixed unit domain with rectangular cell spacing.

| Fixture | Check |
| --- | --- |
| Smooth prescribed erosion | Area-L1 mismatch decreases on both refinements; measured ratios are approximately 0.25 in both dimensions. |
| Wet/dry step or quadrant | L1 ratios are approximately 0.50; mismatch is confined to a one-cell band on each side of the front. Linf remains finite, not necessarily convergent. |
| Signed cell-scale rough rates | 32 consecutive public-projector updates per grid with seed 104729; global bed volume and finite refreshed geometry are required, not rough-field convergence. |
| Retained-state mass exchange | 64 consecutive real erosion transactions per grid with seed 130363; flow and bed persist, only the prescribed substrate inventory is replenished. |

The frozen refinement requirement is a positive L1 ratio below 0.8. The
roundoff multiplier remains 2048 times machine epsilon. Derivative checks use
the predeclared operation scale: maximum Q1 bed elevation times the LS kernel
coefficient L1 norm divided by its dx/dy normalization. This accounts for
cancellation and refinement without scaling by a nearly zero derivative.
No bound was increased after observing the thread-comparison failure.

An independent serial AC/4 scatter supplies nodal proposals and quadrature
weights, including dry contributors and physical boundaries. Every update
checks the nodal bed proposal, global volume, signed mismatch and independently
recomputed Q1 centers/faces, five-point LS slopes/Hessian and center/face G.
Actual transactions additionally check every conservative equation and inventory,
and require dry flow states to remain exactly zero, even where their nodal bed
changes through neighboring proposals. Wet-plus-halo worksets leave other dry
cells inactive. These are controlled source/projection tests, not advection
convergence tests or a new erosion law.

Across four configurations, 2352 bed updates pass: 1584 production exchange
transactions and 768 signed projection updates. Full intermediate conservative,
bed, geometry and inventory fields are serialized and compared bitwise between
teams, rather than comparing only final states. Four checker self-tests cover
valid evidence and reject incomplete teams/cases, nonfinite or inconsistent volume/localization evidence
and absent L1 refinement.

## Correction exposed by the tests

The initial optimized 1D comparison differed in 6060 cached geometric values,
with a maximum absolute difference of approximately 7.8e-13. Conservative state,
authoritative bed, inventories and all 2D fields already matched. Strict runs
matched throughout. Optimized vector/scalar stencil accumulation depended on
the OpenMP chunk layout, violating the retained bitwise comparison contract.

`topography_reconstruction` now uses two explicitly thread-private `VOLATILE`
scalar accumulators only in its 1D five-point stencil. OpenMP, coefficient
values, filter definition and the 2D kernels are retained. The final strict and
optimized intermediate fields agree bitwise with one/four-thread teams. This
bounded reproducibility correction is the only production source change;
the complete regression was rerun after it.

## Unequal areas and completion scope

The separately approved synthetic algebra uses six unequal cell areas and
constant, signed and isolated-source fields. AC/4 weighted gathering agrees
with both nodal quadrature and independent cell-wise Q1 integration. A negative
control using unweighted arithmetic gathering fails the unequal-area identity.
The uniform-area specialization is verified against the actual production
projection throughout the numerical fixtures.

Run from the repository root:

```sh
KEEP_TEST_WORKDIR=1 sh TESTS/TEST_MASS_EXCHANGE/run_test.sh
```

N8-A, N8-B and D-N8-AREA are retained as completed actions in the
[closure checklist](N7_N8_CLOSURE_2026_10_10.md). Next is N7-A: the remaining
1D dynamic and equilibrium cases. Universal positivity, other physical-law
hierarchies, nonuniform production grids and N9 scientific sign-off are not
claimed by this package.
