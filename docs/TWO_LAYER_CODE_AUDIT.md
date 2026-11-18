# Two-layer closure audit — TL13a

**Branch:** `feature/two-layer`. **Audit baseline:** `7f8189a76e`, 18 November 2026.
**Scope:** architectural closure and functional gate; scientific model validation remains TL13b.
This document is the current audit, not a historical compatibility matrix.

## Architectural ownership

| Domain | Authoritative owner | Current invariant |
| --- | --- | --- |
| Model configuration | `model_layout_2d` and per-layer physical descriptors | State arrays have shape `(component,layer,x,y)`; production `supported_layer_count` still permits only one layer. |
| Fluid and substrate inventories | `state_fields::state_type` | `q/qp` are per-layer; bed, deposit, erosion and erodible inventories are shared. |
| Contact exposure | `layer_contacts_2d::evaluate_layer_contacts` | Four two-layer wet/dry occupations define atmosphere/substrate contacts, not layer-dependent copies of physics. |
| Source scheduling | `source_registry_2d` | L1, L2 and INTERLAYER scopes are validated in one registry. |
| Spatial and Lie integration | `spatial_operator_2d`, `time_integration_2d`, `lie_source_executor_2d` | Shared OR workset, per-layer operators, one cell-local Lie transaction and one accepted fluid publication. |
| Coupled physics | `interlayer_mass_exchange_2d` and `interlayer_momentum_exchange_2d` | Conserved gross transfers are limited jointly; drag remains a separate analytic process. |
| Binary checkpoint and NetCDF | `inpout_2d` | One canonical binary format with a layer-count header; `q_layers` and `qp_layers` are the only NetCDF layer-state tensors. |

The single-layer NAMELIST names without `_L1` are **intentional input design**, not obsolete interfaces. No compatibility with previously generated binary checkpoint or NetCDF schemas is required.

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
NetCDF frame using `restart=.TRUE.`, then checks that both canonical L1/L2
state tensors are readable, unchanged and that the time dimension has exactly
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
single-layer end-to-end validation are distinct. A production two-layer
integration run is unavailable while the intentional layout gate is closed;
scientific approval and TL13b remain separate.
