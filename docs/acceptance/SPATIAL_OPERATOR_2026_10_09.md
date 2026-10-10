# Direct HP-PCCU spatial-operator validation

The production spatial operator agrees with the unchanged reference cores for
the covered residuals, geometry, reconstructed traces and directional CFL
bounds. All 128 direct comparisons and all 35 acceptance runs pass. No production
source or numerical tolerance was changed. The combined 2-D timestep remains
different from the explicit reference: this is documented, not accepted as
agreement. N7 is still partial and N9 remains open.

The subsequent [IMEX stage and order audit](IMEX_STAGES_2026_10_10.md)
validates the covered autonomous two/three/four-stage cases, reproduces a
single-stage assembly defect and measures tightening stage CFL bounds. Its
characterization results do not close the combined-CFL or N9 decision.

## Scope and independent reference

[`TEST_SPATIAL_OPERATOR`](../../TESTS/TEST_SPATIAL_OPERATOR/run_test.sh) exercises
the public production `spatial_operator_type` API, actual domain worksets and
the production curvature-source evaluator. The Python comparator imports the
unmodified 1-D/2-D cores from the September 26 drySafe reference archive. Their
SHA-256 fingerprints are checked before import and retained in the
[evidence artifact](results_spatial_operator_2026_10_09.json). Neither reference
formulas nor production kernels are patched to obtain agreement.

Eight fixtures comprise four lakes at rest (flat, smooth, seeded rough and
strong ramp beds) and four dynamic states (smooth, rough, drainage and compact
wet support). Each runs with all four slope/curvature flag combinations, two
compiler profiles and actual teams of one and four threads: 128 comparisons,
or 64 thread pairs. The rectangular grid has 18 by 16 cells, with
`dx=0.75 m`, `dy=1.25 m`; random beds use NumPy PCG64 seed `20261009`.
The smallest observed slope gravity factor is about `0.2688`.

The closure is pure liquid, with one zero solid component, constant temperature
`300 K`, liquid density `1000 kg/m3` and ambient density `1.2 kg/m3`. This avoids
conflating the prototype composition layout with production mass fractions.
Variable temperature/composition, rheology, stochastic evolution and actual
IMEX advancement are outside this focused comparison. A constant boundary
collar supplies common boundary semantics; general boundary-condition
equivalence is not claimed.

The compact fixture uses the real solve-list builder: 208 active cells include
196 dry halo cells, while 80 cells are inactive. The independent full-grid
reference must have zero residual outside this workset. Assertions do not
silently crop a nonzero reference residual out of the comparison. Undefined
inactive output storage is normalized only by the test driver before export.

## Assertions and results

The comparator checks complete active-cell arrays for `R = -spatial + sources`,
including both momentum components and the separate curvature source. Thermal
transport is checked against mass transport at the fixed heat capacity and
temperature. A closed-collar global mass budget and lake-at-rest zero residual
are checked independently.

Additional assertions cover Q1 centers/faces, retained LS slopes and Hessian,
center/face gravity factors, both conservative and primitive face states,
normal and tangential velocities/momenta, thermal and zero-solid traces,
surface elevation, nonnegative final depths, compatible mass/normal-momentum
means, one-sided characteristic bounds and both directional timestep limits.
The driver also checks facade/backend equality, input-state immutability,
refresh of poisoned reconstruction caches, repeated CFL calls, the timestep
cap and the existing CFL sentinel contract.

The fixed limit is `8192*epsilon(float64)*scale`. Normally
`scale=max(1,max(abs(reference)))`; lake residuals use
`rho_c*g*max(1,hmax^2)/min(dx,dy)`. No limits were relaxed.

| Check | Observed result |
| --- | --- |
| Largest normalized array/residual discrepancy | `2.79e-15` |
| Largest scaled lake-at-rest residual | `1.13e-16` |
| Largest scaled global mass-budget error | `3.88e-17` |
| Complete test payloads, one versus four threads | Bitwise identical within each profile, all 64 pairs |
| Curvature toggle at fixed slope flag | All spatial fields and timestep bitwise unchanged |
| Comparator negative controls | Seven deliberate corruptions detected in each profile |
| Complete acceptance suite | 19 unit/documentation scripts plus 16 full-solver runs: 35/35 |
| Audit-tool self-tests | 7/7 |
| Existing full-solver output/checkpoint comparisons | All 16 unchanged against the N8 baseline |
| Cross-thread Gate-H/restart comparisons | All six pairs pass |

At rest, curvature sources are exactly zero in value; signed `+0/-0` is not
claimed to be byte-identical. Strict and optimized production modules both
enable OpenMP and invalid/divide-by-zero/overflow traps; the assertion driver
is compiled strictly in both profiles. Different profiles are not required to
produce identical bytes.

## Combined timestep: remaining scientific decision

With identical directional characteristic bounds, the retained production
policy is `min(max_dt, dt_x, dt_y)`, while the explicit SSPRK2 reference uses
`1/(1/dt_x + 1/dt_y)`. At the same configured CFL, the measured production/reference
ratio is `1.5253` to `1.6540`. For the flat lake it is `1.6`:
`0.04066153193134508 s` versus `0.025413457457090675 s`.

The optional `--require-reference-timestep` gate exits with an assertion for
this mismatch; its command, failure and log fingerprint are archived. The
default suite validates the retained policy explicitly and reports the gap;
its success does not certify combined timestep agreement.

Production input already limits 2-D CFL to `0.25`. Consequently the retained
directional-minimum policy bounds the combined Courant number by `2*CFL <= 0.5`.
The mismatch alone is not evidence of instability. Conversely, an SSPRK2
positivity argument does not establish admissibility of the production IMEX
stages. Before changing this policy, characterize stage admissibility,
temporal order and convergence for the actual tableaux. The raw Euler-mass
diagnostic retained in the artifact is not an IMEX-stage proof.

## Provenance and reproduction

This is an uncommitted test/documentation snapshot based on
`8af163fe11c0ab3b4b5f18837cef4f54db0783c7`. The production source hashes are
unchanged against the [N8 baseline](N8_MASS_EXCHANGE_2026_10_09.md).
The actual repository is `/Users/demichie/Codes/GIT/IMEX_SfloW2D_dev`.

From that repository root:

```sh
sh TESTS/TEST_SPATIAL_OPERATOR/run_test.sh
```

Set `KEEP_TEST_WORKDIR=1` to retain generated fixtures, binaries and payloads.
Dependencies are gfortran, LAPACK, Python 3 with NumPy, `nf-config` and
`nc-config`. To isolate the combined-CFL mismatch, run the comparator with a
retained strict executable, a fresh working directory and
`strict --case lake_flat --require-reference-timestep`; a nonzero exit is
expected for this unresolved gate.

For the full uncommitted candidate, choose a new output directory:

```sh
python3 TESTS/ACCEPTANCE/run_acceptance.py \
  --working-tree --revision HEAD \
  --output /tmp/imex-spatial-repeat \
  --netcdf /opt/homebrew/opt/netcdf-fortran --jobs 3
```

After committing, use its SHA without `--working-tree`. The JSON artifact
retains source/test/tool hashes, input/output fingerprints, 128 direct records,
per-field maximum errors, complete audit results and the baseline comparison.
Full logs remain in `/tmp/imex-spatial-acceptance-20261009-final`; temporary
directories are not permanent archives or approved golden outputs.

Mixed-composition dynamics, the remaining 1-D reference hierarchy, IMEX stage
and order acceptance, source-overlap/I/O contracts and explicit baseline
approval remain tracked in the
[acceptance plan](../../TESTS/ACCEPTANCE/acceptance_plan.json). This step does not
authorize the next structural/layer milestone.
