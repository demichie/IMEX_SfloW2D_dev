# N8 mass exchange and evolving bed validation

The [N7/N8 closure checklist](N7_N8_CLOSURE_2026_10_10.md) reconciles these
passing transaction assertions with the original required refinement and
repeated-update tests. Broader physical-model coverage boundaries below do not
automatically add new mandatory N8 tests.

Production erosion/deposition tests now verify limited cell inventories,
conservative sources, global bed volume, refreshed geometry and exact restart.
They exposed five defects corrected in this candidate. All 34 acceptance runs
pass, and all 16 existing full-solver cases retain identical conservative-output
and checkpoint fingerprints. This establishes the controlled N8 test scope;
the broader N7/N9 acceptance requirements remain open.

## Production corrections

The changes are confined to
[`equation_terms_2d.f90`](../../src/equation_terms_2d.f90) and
[`mass_exchange_2d.f90`](../../src/mass_exchange_2d.f90).

| Defect | Correction |
| --- | --- |
| Below `alphastot_min`, erosion consumed substrate while conservative and bed sources remained zero. | Apply the existing cutoff before proposing any exchange rates, keeping the transaction consistently zero. |
| Positive liquid loss alone did not activate the mass-exchange updater. | Include allocated, positive liquid `loss_rate` in the activation condition, with nested allocation guards. |
| A negative available carrier reserve turned a nonnegative requested loss into carrier gain. | Clamp the availability cap to zero before applying the existing loss limit. |
| Re-erodible deposits changed shared `T_erodible` inside the parallel cell loop. | Select the substrate temperature locally for each source evaluation, preserving the flow-temperature policy without mutating configuration. |
| The bottom-source mask applied to radial sources but not fissural sources. | Apply the shared `cell_source_fractions` mask for either geometry, before updating inventories and projecting bed proposals. |

The cutoff's physical policy is unchanged: it still suppresses the entire
transaction, including carrier terms, below the configured solid fraction.
Signed pore-pressure gas exchange is not clipped to a new nonnegative law.
Source coverage geometry, settling laws, porosity conventions, HP-PCCU fluxes,
CFL calculation and reference tolerances are unchanged.

## Test cases and independent checks

[`TEST_MASS_EXCHANGE`](../../TESTS/TEST_MASS_EXCHANGE/run_test.sh) compiles the
production modules, including the real checkpoint writer and reader. Its
rectangular 8 by 7 fixture contains two solid classes and a liquid carrier,
spatially varying temperature, two dry cells on the solve list and one wet
inactive corner.

Nine independent cases cover saturated deposition, saturated erosion, combined
exchange, full/partial radial masking, full/partial fissural masking, the solid
cutoff, liquid loss alone, an exhausted carrier reserve and local thermal
sources. A restart case checks three further transactions: the initial update,
continuous second update and resumed second update. The second update is
nonzero, so restart equivalence cannot pass through an inactive continuation.

Expected conservative masses, momenta, thermal energy and each
`deposit`/`erosion`/`erodible` inventory are assembled independently from known
phase volumes and limited increments. The reference reuses only terminal
settling velocity; it does not call the production exchange evaluator or nodal
projection. Explicit assertions ensure that the intended erosion and
deposition caps are saturated.

Independent nodal gathering includes dry and inactive zero proposals in the
geometric denominator. Q1 quadrature uses quarter-area corner weights,
half-area edge weights and full-area interior weights. Its global bed-volume
change must equal `dx*dy*SUM(cell_delta)`. Inventory increments have units of
solid volume per area [m], while conserved solid masses have units [kg/m2].

Each transaction checks Q1 centers/faces, the retained five-point LS slopes and
Hessian, and center/face gravity factors. Checkpoint reading is tested after
poisoning all these derived geometry caches. The read must restore the exact
nodal bed and inventories and regenerate the caches before continuation.

The diagnostic checker compares production cell/geometric volume rates and
signed/L1/Linf mismatch against independent reductions. The first four fields
have units [m3/s]; Linf has units [m/s]. Nonzero local mismatch is required in a
dry cell and the wet inactive corner when the bed evolves. A zero global signed
mismatch therefore does not imply zero local geometric mismatch.

## Verification results

| Check | Result |
| --- | --- |
| Finalized test driver against original numerical sources | Two cases pass; eight fail at assertions reproducing the five defects. |
| New tests with strict numerical modules | All 12 transactions pass with actual teams of 1 and 4 threads. |
| New tests with optimized numerical modules | All 12 transactions pass with actual teams of 1 and 4 threads. |
| Production bed diagnostic comparisons | 48 comparisons pass across the four configurations. |
| New test snapshots and checkpoints | Bitwise identical between 1 and 4 threads within each profile; continuous/resumed conservative state, bed and inventories match bitwise. |
| Acceptance suite | 18 unit/documentation scripts and 16 full-solver runs pass, 34 in total. |
| Existing full-solver output comparison | All 16 cases retain identical canonical q/checkpoint SHA-256 fingerprints against the Froude-fix baseline. |
| Cross-thread Gate-H and full restart comparisons | All six pairs pass. |
| Audit-tool self-tests | Seven pass. |

The new test uses strict `-O0 -g -fcheck=all` and optimized
`-Ofast -funroll-all-loops` numerical modules. Both enable OpenMP and
invalid/divide-by-zero/overflow traps. The assertion driver remains strict in
both builds. State/geometry assertions use
`2048*epsilon(wp)*max(1, reference_scale)`; parsed diagnostic assertions use
`2e-12*max(1, abs(reference))`. Exact comparisons use bytes or integer bit
patterns, not a relaxed floating-point tolerance. Diagnostic reduction order
may change low bits between thread counts without changing the state.

## Candidate provenance and reproduction

The tested sources are an uncommitted working-tree snapshot based on
`5c47656f19b35d300a7ce1054266e12d8a06b98e`; that commit alone does not include
these corrections. The actual source files are in
`/Users/demichie/Codes/GIT/IMEX_SfloW2D_dev`.

[`results_n8_mass_exchange_2026_10_09.json`](results_n8_mass_exchange_2026_10_09.json)
retains source/test/tool hashes, build flags, original/effective input hashes,
output fingerprints, actual teams, the 48 diagnostic records, pre-fix assertion
failures and comparison with the
[preceding baseline](FROUDE_FRICTION_FIX_2026_10_09.md). Raw build/invocation logs
remain in `/tmp/imex-n8-acceptance-20261009-final`; pre-fix logs remain in
`/tmp/imex-n8-before.9YofFU`. Observed output hashes are not approved golden
references, and these temporary directories are not permanent archives.

Run the focused test from the repository root:

```sh
sh TESTS/TEST_MASS_EXCHANGE/run_test.sh
```

Set `KEEP_TEST_WORKDIR=1` to retain its temporary binaries, snapshots and
checkpoints. The test requires gfortran, LAPACK, Python 3 and both `nf-config`
and `nc-config`, using their library paths to link the production I/O module.

For the complete uncommitted candidate, choose a new output directory:

```sh
python3 TESTS/ACCEPTANCE/run_acceptance.py \
  --working-tree --revision HEAD \
  --output /tmp/imex-n8-repeat \
  --netcdf /opt/homebrew/opt/netcdf-fortran --jobs 3
```

After committing, use the resulting SHA without `--working-tree` for tracked
source reproduction. Test outputs are generated outside the working repository.

## Remaining acceptance work

This fixture validates a liquid/two-solid exchange regime, not all physical
models. Gas exchange, pore pressure, vaporization, entrainment, source overlap
geometry and end-to-end injection budgets require additional characterization.
The existing cutoff policy for initially solid-free flow remains a separate
physical-model question, not a silently changed behavior in this correction.

The [current criterion definitions](../../TESTS/ACCEPTANCE/acceptance_plan.json)
retain these boundaries. Production residual/CFL and dynamic mixed-composition
comparisons are the next numerical acceptance work. Long ETNA/lava provenance,
NetCDF/source contracts, approved references and explicit scientific sign-off
remain N9 conditions. The
[original acceptance matrix](N7_N8_N9_AUDIT_2026_10_09.md) remains the historical
record for its earlier baseline.
