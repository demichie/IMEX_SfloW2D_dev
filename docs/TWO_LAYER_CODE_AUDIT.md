# Two-layer closure audit — TL13a

**Branch:** `feature/two-layer`. **Audit baseline:** `7f8189a76e`, 18 November 2026; output-schema wording verified against `src/io/inpout_2d.f90` on 19 November 2026.
**Scope:** architectural closure and functional gate; scientific model validation remains TL13b.
This document is the current audit, not a historical compatibility matrix.

## Architectural ownership

| Domain | Authoritative owner | Current invariant |
| --- | --- | --- |
| Model configuration | `model_layout_2d` and per-layer physical descriptors | State arrays have shape `(component,layer,x,y)`; layout admits one or two layers, while the production parser requires explicit TL13a opt-in for two. |
| Fluid and substrate inventories | `state_fields::state_type` | `q/qp` are per-layer; bed, deposit, erosion and erodible inventories are shared. |
| Contact exposure | `layer_contacts_2d::evaluate_layer_contacts` | Four two-layer wet/dry occupations define atmosphere/substrate contacts, not layer-dependent copies of physics. |
| Source scheduling | `source_registry_2d` | L1, L2 and INTERLAYER scopes are validated in one registry. |
| Spatial and Lie integration | `spatial_operator_2d`, `time_integration_2d`, `lie_source_executor_2d` | Shared OR workset, per-layer operators, one cell-local Lie transaction and one accepted fluid publication. |
| Coupled physics | `interlayer_mass_exchange_2d` and `interlayer_momentum_exchange_2d` | Conserved gross transfers are limited jointly; drag remains a separate analytic process. |
| Binary checkpoint and NetCDF | `inpout_2d` | One canonical binary checkpoint with a layer-count header stores full conservative/primitive state; visualization NetCDF writes named scalar physical fields per layer, without `q_layers` or `qp_layers` tensors. |

The single-layer NAMELIST names without `_L1` are **intentional input design**, not obsolete interfaces. No compatibility with previously generated binary checkpoint or NetCDF schemas is required. `INTERLAYER_PARAMETERS` is the only pair-wide NAMELIST, covering mass exchange and momentum drag; the redundant input spelling and its parser branch have been removed.

## TL13a cleanup in this change

The prior `substrate_snapshot_type` plus `state_type%snapshot_substrate` and
`state_type%restore_substrate` were not used in the current production timestep.
They existed only to demonstrate a caller-managed recovery in two test programs.
This API has been removed. The workspace regression retains the granular/lava
aggregate-inventory checks and independent geometry-transaction checks.
The direct Lie rejection regression now verifies the executor's failure
status and diagnostic, without claiming that a failed timestep restores
shared substrate state.

**Failure contract:** A source or spatial rejection terminates the timestep.
No retry and no recoverable whole-step substrate rollback are supported.
`q/qp` are published only after successful Lie execution. The geometry
transaction remains a separate mechanism for deferred bed-cache updates;
its internal rollback fixture is still valid.

## Canonical I/O continuation and end-to-end verification

The existing structural `TEST_SPATIAL_OPERATOR` fixture now appends a second
NetCDF frame using `restart=.TRUE.`, then checks that the canonical named L1/L2
physical fields are readable, unchanged and that the time dimension has exactly
two records. This extends the round-trip from initial serialization to
restart continuation without introducing a second I/O implementation.

The existing `solver-validation.yml` is the production-driver acceptance
route: optimized and strict/debug builds, `make check`, binary and stochastic
restarts, source-mass and HP-PCCU tests, canonical two-layer structural
regressions, and complete single-layer NetCDF/OpenMP comparison across
1/2/4 threads. The OpenMP job is enabled by default for manual dispatch.
This is intentionally a *single-layer production driver* test plus a
*structural two-layer* test; it is not a claim of two-layer production execution.

## Evidence and remaining functional gates

- On baseline `4bc2dce0`, optimized/debug build
  [#44345533346](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/44345533346)
  and numerical regressions
  [#44345568177](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/44345568177)
  both concluded **success**.
- No workflow outcome exists yet for the subsequent TL13a integration
  extensions. The new SHA must pass the full solver-validation workflow
  before functional acceptance claims can be strengthened.
- The production `N_LAYERS=2` gate remains closed. Structural two-layer
  fixtures are not equivalent to a complete production run.
- A **two-layer** production end-to-end scenario still requires an approved,
  controlled layout activation and must exercise parsing, initial conditions,
  boundaries, spatial/Lie evolution, substrate budgets, output, continuation
  and OpenMP determinism through the actual driver.
- On the final cleanup SHA, run optimized and strict/debug builds, all
  numerical/source/state regressions, single-layer examples and the
  representative lava/dry variants. Compare scientific budgets, not only exit codes.
- TL13b must separately validate coupled hyperbolicity/shear, density
  contrasts, pressure and closure limits before any production release.

**TL13a gate remains open:** structural two-layer integration and real
single-layer end-to-end validation are distinct. A two-layer integration run is now available only through explicit TL13a
validation opt-in;
scientific approval and TL13b remain separate.


## Definitive interlayer input interface

The canonical `source_input_2d::read_interlayer_parameters` and
`write_interlayer_parameters` handle one pair-wide
`INTERLAYER_PARAMETERS` NAMELIST that refers directly to the
existing authoritative mass- and momentum-exchange configuration
objects. The additional exchange-only NAMELIST and its conditional
reader were removed; no shadow configuration or compatibility
parser remains. Single-layer short NAMELISTs are intentionally retained.


## Interlayer configuration regression

`TEST_SOURCE_REGISTRY` now reads one canonical pair-wide NAMELIST with
distinct mass and drag groups and checks the Chezy coefficient. It verifies
that the backup writes the canonical group, rejects duplicate occurrences,
and rejects interlayer input when only one layer is configured. The test
uses temporary Fortran scratch units and restores both configuration
objects, so it does not create a parallel configuration or modify other
test inputs.


## Canonical scientific naming audit (18 November 2026)

The nomenclature migration is deliberately **semantic**, not a general
case/style conversion. It updates the production Fortran owners, their
REAL/COMPLEX consumers, shared `.inc` kernels and directly dependent tests.

| Quantity | Sole current name | Ownership and distinction |
| --- | --- | --- |
| Reciprocal bulk mixture density | `inv_rho_m` | REAL/COMPLEX `physical_state` and mixture kernels; distinct from `inv_rho_c`. |
| Pair-wide mass-exchange configuration | `interlayer_mass_exchange` | Owned by `interlayer_mass_exchange_parameters_2d`; directly referenced by `INTERLAYER_PARAMETERS`. No pointer alias is retained. |
| Rheology configuration pointer | `layer_physics_type%rheology` | Non-owning reference to per-layer canonical rheology, parallel to the thermal/pore/mass references. |
| Equation and transported-variable counts | `n_eqns`, `n_vars` | Shared by global index conventions, model layout and source registry. One `n_vars` sizes both `q` and `qp`; the redundant `n_physical_variables` was removed. |
| Optional component counts | `n_add_gas`, `n_stoch_vars`, `n_pore_vars` | Layout metadata; no second spelling for those counts. |
| Primitive temperature index | `model_layout%idx_T` | Explicitly initialized in the layout, independent of global-index initialization. |
| Conservative energy equation index | `model_layout%idx_engyEqn` | Semantically distinct from primitive temperature even when both index values are four. |
| Kernel flow depth and total solid fraction | `flow_depth`, `total_solid_fraction` | Erosion and deposition interfaces use names shared with their downstream mathematical kernels. |
| Settling drag inputs | `particle_diameter`, `solid_density`, `carrier_density` | Identical physical meanings throughout the terminal-velocity kernel. |

Further intentional distinctions:

- `q` (conservative), `qp` (primitive) and their `trial`/`original`
  stages are not interchangeable.
- `idx_totMassEqn`, `idx_uEqn`, `idx_vEqn`, `idx_engyEqn`
  describe conservative equations; `idx_h`, `idx_u`, `idx_v`,
  `idx_T` describe primitive components. Equal integer values are not
  synonymous physical quantities.
- Horizontal speed, total speed and squared speed remain distinct.
  The local letter `T` in cell-level physics represents temperature.
- `r_h`, `rho_c`, `rhoc` and `config` outside the audited kernel
  interfaces are not all proven synonyms. Further broad renaming would
  require checking each procedure's units, scope, formulation and
  keyword-argument callers; a superficial global replacement is unsafe.

Verification for this nomenclature migration must be performed against
its final branch SHA. Previously green CI and currently running workflows
on preceding commits do **not** validate these edits. The scientific
formulations, source ordering, input NAMELIST group names and production
`N_LAYERS=2` gate are unchanged.

## TL13a — ESRI two-layer output integration (20 November 2026)

The production `output_esri` routine previously returned before its
`output_two_layer_esri` call, making all L2/combined ESRI raster output
unreachable. The unconditional `RETURN` has been removed. The existing
`TEST_SPATIAL_OPERATOR` two-layer I/O round-trip now enables ESRI output
through `output_solution` and checks creation of L2 depth, L2 temperature,
combined depth and upper-surface rasters. This is a structural regression of
the real output entry point; it is not a production two-layer run. Run strict
and optimized CI on this exact commit to verify the newly exercised path.
The `N_LAYERS=2` production gate remains closed pending controlled end-to-end
functional validation and TL13b scientific approval.

## TL13a — controlled production-driver activation (20 November 2026)

The structural layout now accepts one or two layer slots. The **only**
activation gate in the production parser remains closed by default:
`N_LAYERS=2` is rejected unless `IMEX_TWO_LAYER_VALIDATION=1` is explicitly
set in the environment. This is a diagnostic opt-in, **not scientific approval**.
The single-layer path retains its original layer count and numerical formulae.
The driver now forms one OR workset at each step, prepares L1/L2 CFL traces,
and obtains the preliminary coupled bound using the existing coupled spatial
operator before choosing the timestep. The coupled Lie integrator then uses
those traces through its existing production call. No alternate numerical
integrator, state representation or physical source law is introduced.

This commit has not yet completed an end-to-end two-layer run or scientific
TL13b gate. Run both CI profiles, the single-layer regression set and a short
controlled two-layer run with output/restart/OpenMP before claiming TL13a
functional acceptance. The preliminary total-depth wave bound is not a
validated coupled eigenvalue/hyperbolicity analysis.

## TL13a — new controlled smoke-test evidence route

The optional `two_layer_smoke` input in the manual solver-validation
workflow starts an isolated strict/debug Etna run. The standalone script
`EXAMPLES/EXAMPLE_ETNA_TWO_LAYER/run_validation.py` first verifies that
`N_LAYERS=2` is rejected by default, then prepares the canonical L1 initial
state and exactly dry L2, enables the diagnostic opt-in and checks the first
advanced snapshots, L2 ESRI fields, manifest, mass-budget arithmetic,
NetCDF file and checkpoint file. It preserves the solver's `runner.log`
under the evidence directory. This test cannot replace long-time
conservation, physical calibration, coupled wave-speed validation, restart
continuation or OpenMP comparisons. Its first CI result is pending.

## TL13a — Etna initial snapshot separators (20 November 2026)

The production Fortran conservative snapshot writer emits one blank record
per grid row **and two additional blank records after the final row**.
The original Etna preparator counted every blank record as a grid-row
separator and therefore rejected valid 250 x 150 snapshots with 152
blank records. The preparator now validates the 150 row separators and
allows up to the two canonical trailing records, rejecting additional
records or internal misplaced separators. L1 remains a byte-for-byte
copy of the single-layer snapshot; all L2 conservative components remain
exactly zero at initialization. The synthetic preparation regression now
includes the two trailing records and rejects excessive trailing blanks.
This correction does not alter any physical initial condition or numerical
kernel. The controlled two-layer execution remains pending on the new SHA.

## TL13a — ESRI fixture initialization (20 November 2026)

The newly enabled `output_solution` ESRI regression reached the canonical
`output_esri` writer but failed at its first `grid_output=-9999` assignment:
`TEST_SPATIAL_OPERATOR` initializes its own geometry arrays without calling
`geometry_2d:init_grid`, which normally allocates the shared real/integer
ESRI grids. The test's canonical `initialize_fixture` now allocates
`grid_output` and `grid_output_int` with the same `nx,ny` dimensions as
`init_grid`. This fixes the missing test setup, without altering the
production writer or its physical results. Strict/debug and optimized
spatial regressions on this commit still require execution.

## TL13a — L2 rheology and pore-pressure input completion (20 November 2026)

The first experimental Etna run passed the explicit validation gate but
failed on the parser's unconditional rejection of `RHEOLOGY_L2_PARAMETERS`.
The layer count was correctly two; the parser contained unfinished L2
NAMELIST handling. The L2 rheology group is now required and read into the
canonical `rheology_l2` configuration for two-layer runs, while single-layer
input still rejects L2 groups. The existing early and later pore-pressure
scans now reject L2 only with `N_LAYERS=1` and load the optional typed L2
pore-pressure configuration when present, before computing component counts.
The L2 method/group checks are included in the existing source scheduling
validation. No physical formulas, state schema or input-file aliases changed.
Strict/debug CI and the experimental Etna run remain unverified on this SHA.

## TL13a — extended functional regression (20 November 2026)

An explicit `two_layer_extended` dispatch option now tests the real Etna
two-layer driver beyond the initial smoke stage: 0/0.005/0.01-second
outputs, L2 emergence, mass/solid budget additivity and approximate
conservation, initial L2 dryness, numerical restart equivalence and
full-checkpoint equality, and independent 1/2/4-thread conservative
snapshot comparisons. All runs use the canonical Etna input and explicit
experimental activation; existing production tests are unaffected.
This is a new unverified **functional** gate, not acceptance of the
coupled characteristic spectrum or hyperbolicity with shear. The extended
job must be executed and its numerical outcomes inspected before closure.

The extended-checker Python logic has solver-independent tests for the
binary-restart NAMELIST edit, three-record budgets, numerical snapshot
comparison and rejection of negative inventories. These run before the
optional extended Etna simulation in its workflow job.

## TL13a — separate Rouse-driven L2 activation from initial short-time regression

The first extended run failed the assertion that L2 must grow by t=0.01 s,
although the initial Etna avalanche is at rest. In the canonical
`LINEAR_USTAR` law, positive flux depends on shear-generated `ustar`, the
settling/Rouse limiter and available donor solids; `ELUTRIATION_FLAG=T`
alone does not produce positive flux. The short 0.01-second trajectories
continue to check conservation, restart and OpenMP without assuming a
nonzero physical source. A separate unchanged-physics Etna trajectory to
2 seconds (configurable `--growth-horizon`) now requires L2 growth and
reports time/L2-mass budgets on failure. This is a pending diagnostic,
not evidence that positive mass has already been observed. No source
law, coefficient or canonical example input is changed.

## TL13a — binary in-place restart test input (20 November 2026)

The extended Etna validation failed at the second invocation with
`run output directory already exists`. The Python fixture supplied the
explicit path `OUTPUT/exampleEtna_two_layer/restart/restart.bin`, which is
correctly treated as an external binary restart by the production parser.
In-place continuation instead requires the canonical token
`RESTART_FILES='restart.bin'`; the solver resolves this token inside the
existing output directory and applies its normal output manifest and
snapshot-index safety checks. The fixture and its solver-independent
regression now exercise that contract. No `OUTPUT_OVERWRITE` escape hatch
was added and the output-path implementation was left unchanged.
The end-to-end continuation is not verified on this revision yet.

## TL13a final closure audit — shared-domain ownership (20 November 2026)

Source and Fortran test-tree audit found exactly one active-cell/face workset
in use: `simulation%domains(1)`. The allocated second `domain_type` was
initialized/finalized but never read or populated by the driver or Lie
integrator. `simulation_context_type` now owns **one** `domain` for the
union of L1/L2 wet and boundary-injection cells. The production driver and
time integrator use that authoritative object; separate per-layer
`spatial_operators` and stochastic workspaces remain as required.
No physical quantities, numerical formulas, source scheduling, time bounds,
or output/restart schemas were modified. The single-layer path still uses
the same domain implementation.

### Verified evidence and pending CI

- `847f7906`, workflow run `44710265332`: the experimental Etna job
  finished **success**, including short smoke, independent Python checker
  tests, and extended conservation, L2 emergence, segmented binary restart,
  and OpenMP 1/2/4 comparison. The overall run was later **cancelled**.
- Optimized structural integration on that workflow finished successfully,
  but strict/debug and remaining integration checks were interrupted.
- The earlier `2e64c86` workflow run `44705008879` completed optimized
  and strict/debug general integration successfully; this earlier SHA does
  **not** verify the subsequent TL13a changes or domain cleanup.
- **TL13a stays open until both integration profiles and the two-layer
  functional job are green on the final commit.** No automated execution
  of the new domain ownership code has yet been observed. TL13b scientific
  validation and production acceptance remain separate and open.

## TL13a functional closure — verified CI (20 November 2026)

**Verified baseline:** `89a9f2b46e52f718417e28b7c4db7ee11bba65b9` on
`feature/two-layer`, GitHub Actions workflow run
[44712747883](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/44712747883):
**completed / success**. Every job finished successfully:

- `gfortran / optimized / integration` — success.
- `gfortran / strict-debug / integration` — success.
- `End-to-end OpenMP comparison (1, 2, 4 threads)` — success.
- `Experimental TL13a two-layer Etna smoke` (short smoke and optional
  extended conservation, continuation and two-layer OpenMP checks) — success.

The full workflow verifies the definitive single shared-domain ownership
change as well as the earlier Etna, budget and checkpoint integrations.
**TL13a is functionally complete.** This statement refers strictly to the
verified baseline SHA, not to future documentation-only commits, and does
not authorize a production two-layer run. The experimental environment gate
remains intact. **TL13b remains open:** principal-symbol hyperbolicity with
shear, density contrasts, coupled pressure, wet/dry, refinement and model
closure validation require separate scientific evidence.

## TL13 follow-up — static deposition restricted to L1 (20 November 2026)

The contact evaluator continues to report the actual geometric ground contact
for an isolated wet L2. This does **not** make L2 eligible for static
sediment deposition. The existing Lie group scope identifies L1, and that
boolean is passed through `source_terms_2d` to the canonical
`mass_exchange_2d` coordinator. Static substrate deposition and its
coupled steam-vaporization process now require **both** terrain contact and
L1 scope. Erosion, liquid carrier loss, lava immobilization and atmospheric
air entrainment retain their original geometric eligibility, pending their
separate physical review. No numerical/constitutive formulas were altered.

The paired Lie transfer still supports L2-to-dry-L1 settling and
L1-to-dry-L2 elutriation; neither transfer is conditioned on receiver
presence. The existing atomic/conservation/wet-dry tests remain in
`test_interlayer_lie_transaction.f90`, extended with explicit external
substrate-deposition eligibility checks. Buoyancy deletion is **not**
redefined here as a physical lift-off model; its density threshold and
budget accounting still require a separate scientific decision.

**Verification status:** source review and new regression assertions are
recorded, but strict/optimized builds and integrated Actions on this change
must succeed before declaring the implementation fully verified.

**Dependency follow-up:** `TEST_MASS_EXCHANGE/test_projection_stress.f90`
exercises real terrain/atmosphere contact masks directly and now supplies
`static_deposition_allowed=(layer==1)` to the common mass coordinator;
its L2 terrain-carrier-loss and atmospheric-entrainment assertions remain
unchanged. The paired transaction regression also verifies that a sole
wet L2 cannot trigger vaporization coupled to substrate deposition.
