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

## Evidence and remaining functional gates

- The optimized/debug build on baseline `5f00eec55` succeeded.
- On the later `7f8189a76` regression run
  [#44343739896](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/44343739896),
  the spatial and state jobs were observed successful. Other jobs must be
  checked for their final status before declaring that run fully green.
- No build or regression evidence exists yet for the cleanup commit that
  updates this document. Previous CI runs do **not** validate new code.
- The production `N_LAYERS=2` gate remains closed. Structural two-layer
  fixtures are not equivalent to a complete production run.
- A production end-to-end scenario must exercise parsing, initial conditions,
  all four boundary sides, spatial/Lie evolution, substrate budgets, output,
  continuation and OpenMP determinism through the actual driver.
- On the final cleanup SHA, run optimized and strict/debug builds, all
  numerical/source/state regressions, single-layer examples and the
  representative lava/dry variants. Compare scientific budgets, not only exit codes.
- TL13b must separately validate coupled hyperbolicity/shear, density
  contrasts, pressure and closure limits before any production release.

**TL13a gate is open:** this cleanup removes an obsolete interface, but
does not by itself authorize two-layer production or establish scientific validity.
