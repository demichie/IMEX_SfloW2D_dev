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
