# N7 dynamic geometry and excavation completion

N7-C and N7-D pass their defined contracts after the explicitly authorized
dry-state corrections. The corrected candidate passes all 41 audit groups
and all six full-solver cross-thread comparisons. N7 remains open only for
the separate CFL interpretation, and N9 scientific baseline approval remains
open. Neither milestone is silently closed by passing this package.

The [machine-readable evidence](results_n7_cd_2026_10_10.json) pins all
production sources, versioned contracts, fixture hashes, case inventories,
reference errors and observer fingerprints. The base revision is `aed31540`;
the candidate is identified by its source hashes, not by that historical commit.

## Authorized numerical changes

Auxiliary centre velocities used in direct and tangential reconstruction are
zero at or below the existing shared dry-depth threshold. Original volumetric
reconstruction candidates are retained. The explicit local curvature source
is evaluated only above that threshold, selecting resolved cells before any
quadratic velocity contraction. Other local sources remain unchanged.

After observing and checking the raw final IMEX assembly, the accepted-state
projection sets only the two conservative momenta to zero in unresolved dry
cells. Positive mass, thermal energy and every transported component remain
intact. The corresponding primitive momenta and auxiliary velocities are
updated without repeating the thermodynamic closure. This pure, allocation-free
mapping uses the already recovered mixture depth and the shared cutoff.
Raw IMEX stage states and implicit Newton solves are not projected.

These are bounded numerical policy corrections, not output-neutral
optimizations. CFL policy, IMEX tableaux, dry threshold, flux/path kernels,
historical clean-core checksums and all numerical error limits are unchanged.
Accepted time sequences may change because the computed states change; that
is not a change to the timestep selection formula.

## Dynamic geometry results

Each strict/optimized profile passes 69 finite-time cases with actual one-
and four-thread teams: 24 seeded rough/drainage/compact/inclined fronts, 24
mixed-Hessian dynamic/equilibrium cases, 12 axial/oblique ramp and rotation
cases, and nine circular dry-front refinement cases. All supported stage
counts, `N_RK=2,3,4`, are included.

All 138 trajectory payload pairs are bitwise equal across thread teams.
Every raw solved stage of the last step also has an independent normal and
tangential face audit: 828 observer runs, or 414 bitwise team pairs. These
check face thickness, free surface, momenta, thermal/density closure and
characteristic bounds. The maximum final thickness difference from the
independently assembled current-policy reference is `1.333e-15 m`; the
largest relative mass drift is `6.053e-15` in the defined fixtures.

The six curvature ON/OFF differentials per profile have nonzero mixed-Hessian
contributions and scaled differences approximately `6.54e-6` to `6.65e-6`,
above the frozen `1e-8` non-vacuity bound. Resting counterparts remain balanced.
Circular refinement uses the user-approved change of the signed fourth
spatial harmonic from initial to final state, retaining the absolute final
bound. Finest/coarsest evolution ratios are approximately `0.0493`, `0.0477`
and `0.0459`. The failed former final-only measurement is retained explicitly;
this spatial harmonic is not a physical angular momentum.

## Gate-H field results

All 24 grid/mode/profile/team combinations pass complete thickness and
free-surface comparisons from common initial fields at `t=0.1 s`, using
the current `N_RK=2` IMEX assembly. Accepted times are recorded by a test-only
main linked against the tested numerical objects. Its canonical final output
must be byte-identical to the uninstrumented executable; the reference uses
those times, never intermediate actual states.

For slope plus curvature, the largest thickness errors across both profiles
and thread teams are:

| Grid spacing [m] | Maximum thickness error [m] | Unchanged bound [m] |
| ---: | ---: | ---: |
| 0.4 | 4.993e-12 | 2.066e-9 |
| 0.2 | 6.180e-12 | 2.183e-9 |
| 0.1 | 1.518e-11 | 1.208e-9 |

The bound remains
`32768 * float64_epsilon * accepted_steps * max(1, initial_h_max)`;
accepted step counts are 142, 150 and 83 respectively. All free-surface and
G=1 comparisons also pass. Existing aggregate, positivity, mass and
zero-uphill requirements pass unchanged.

Separately regenerated historical SSPRK2 fields pass the user-approved
normalized L1 limits `2e-4` for thickness and `1e-4` for free surface and
the `0.003 m` Linf limit. Historical slope-plus-curvature Linf differences
are approximately `0.503`, `0.398` and `0.734 mm`. The historical algorithm
retains its own temporal and dry-cell policies; it is not relabeled as the
current IMEX oracle.

## Dry momentum diagnosis and remaining CFL decision

Before the accepted-state projection, the intermediate common-time grid
failed by approximately `3.256e-5 m`. A cell of depth `2.196e-15 m` retained
and amplified momenta while auxiliary velocities and the local curvature
source were already zero. Rewetting exposed that retained momentum to the
wet-neighbour reconstruction. The [diagnostic history](N7_CD_DIAGNOSTICS_2026_10_10.md)
preserves the failing audits, local replay evidence and earlier face-switch
mechanism; none of these failures was waived by changing a tolerance.

The corrected observational trace has byte-identical canonical output and
zero accepted dry momentum over all 150 steps, while recording 1770 instances
of positive retained dry mass. The formerly problematic cell still retains
`2.196349e-12 kg/m2` without its old momenta. Unit tests also protect exact
mass/thermal/component preservation, cutoff equality, wet momenta,
idempotence and rewetting.

The supplementary autonomous-timestep comparison now passes on the strict
intermediate slope-plus-curvature case, with a thickness error of
`4.052e-10 m`. This single-case result does not approve the hierarchy-wide
adaptive-control contract. Earlier autonomous failures remain archived;
`D-N7-CFL` remains open under the user-approved separation. No retry mechanism
or CFL change has been introduced.

The complete corrected-candidate audit retains the original-duration
deterministic and stochastic restarts, geometry, variable-composition,
thermodynamic, mass-exchange and source regressions. All 41 groups pass in
the recorded compilation policies. Source hashes match the live numerical
files; a subsequent reference-adapter module-docstring edit is separately
verified to leave its numerical AST unchanged. Closure metadata is reconciled
after execution and checked independently. Long ETNA/lava scientific baseline
acceptance and broader model characterization remain in N9.
