# D19 — Lie splitting migration: implementation guide and live tracker

**Status:** IN PROGRESS — architecture rewritten in production modules; no final compilation or tests yet.
**Branch:** `master`
**Companion architecture:** [D18](D18_IMEX_RK_SOURCE_ARCHITECTURE_PRINCIPLES.md). D18 describes the numerical contract; **D19 tracks the implementation**, file-by-file changes, dependencies, evidence, and remaining work.
**Last reviewed:** 2026-11-01
**Rule for maintenance:** commits may deliberately be non-compilable during the large migration. Do not require intermediate builds; compile and run the complete acceptance suite once the end-to-end architecture and references are finished. Record actual commits and checks; do not mark a task DONE without evidence.

## 1. Frozen target contract

The **only global temporal integration method** is first-order Lie splitting. Retain configurable physical source groups, their order (from GROUP 1 onward), method selection, per-model numerical interfaces and physical ownership of limiters. Do not hard-code a source family to EXPLICIT, NEWTON or ANALYTIC. Preserve semi-implicit interfaces as supported model/NNewton auxiliaries, rather than RK channels.

The spatial finite-volume operator always belongs to **GROUP 0**, which is always EXPLICIT and executes first. Additional GROUP-0 sources must also be EXPLICIT. Crucially, distinguish **the state used to evaluate rates** from **the state to which their increments are applied**:

1. Preserve the original conservative and required physical state `q_n` for GROUP 0.
2. Evaluate spatial flux divergence on `q_n` and **advance only the spatial operator** to `q_sp = q_n + dt*H(q_n)`.
3. Check admissibility of `q_sp`. **No physical source-rate limiter is applied to the spatial operator or its rate.** Its FV stability, conservation, boundary and CFL requirements remain independent.
4. If GROUP-0 sources exist, evaluate **every source rate using the same original `q_n`** (and corresponding original primitive/inventory inputs), not `q_sp`.
5. Apply any *physical joint rate limiter* to **only the individual GROUP-0 source rates**. The candidate source increments will be applied to `q_sp`, so check physical availability and consistency against the **actual recipient state/inventories** where appropriate. Define and test this interaction explicitly; do not presume that using `q_n` for rate evaluation ensures admissible increments on `q_sp`.
6. Apply the sum of limited GROUP-0 rates to `q_sp` **once**: `q_0 = q_sp + dt*sum(S_p_limited(q_n))`. Commit any associated inventory transfer coherently. Check admissibility of `q_0`.
7. Execute groups `1,2,...` in ascending configured order. Each group takes as input the validated complete state coming out of its predecessor and advances for the **full same `dt`** using its configured supported method. Apply appropriate physical limiters within the group and check modified state/inventories after each completed group.
8. Commit accepted state and advance the global clock **once** to `t_n+dt`. Never increase the clock once per group. On failure, fail consistently or reject/retry the whole step without partial fluid/bed/topography commitment.

There is no `STAGE`/`STEP` source cadence, global IMEX-RK tableau, DIRK stage diagonal, interstage RK rates or Strang option.

## 2. Scope: KEEP / CHANGE / DELETE inventory

Based on static inspection of the branch at the planning stage. **All code operations below are pending.** File-specific names must be audited for remaining references before deleting.

| Path | DELETE (RK/cadence-only) | CHANGE / KEEP |
| --- | --- | --- |
| `src/numerics/time_integration_2d.f90` | `n_RK` handling, `validate_runge_kutta_stages`, Butcher tableaux (`a_tilde_ij`, `a_dirk_ij`, weights), `q_rk`, `qp_rk`, `EX_terms`, `IM_terms`, stage-indexed `spatial_terms`, stage loops, `imex_RK_solver`, RK observer events | Implement Lie driver, preserve `q_n` through GROUP 0, one spatial update, independent spatial admissibility, group-0 source-only limiter/update, then group loop and post-group checks. Keep only needed spatial/source workspaces. |
| `src/numerics/nonlinear_solver_2d.f90` | `solve_imex_stage`, RK-only `solve_rk_step` adapters, `effective_im_term`, `a_diag`, `known_stage_residual`/`Rj_not_impl`, `n_RK==4` tolerance special-case | Retain reduced conservative-variable Newton, complex-step Jacobian, line search, diagnostics and model interfaces. Backward-Euler group residual `F(q_out)=q_out-q_in-dt*S(q_out)`. Preserve correct triangular treatment or expand unknown set when excluded components have nonzero coupled sources. |
| `src/numerics/source_catalog_2d.f90` | `SOURCE_CADENCE_STAGE/STEP`, `default_cadence` | Retain physical source catalog, capabilities and `SOURCE_METHOD_*` metadata. Model/method compatibility remains validated. |
| `src/numerics/source_registry_2d.f90` | Group `cadence` field; separate `stage_entries`/`step_entries`; `stage_groups()`/`step_groups()`; stage/step-only validation and placeholders | Single ordered group registry with group ID, method, source names, numerical unknown equations as needed; GROUP-0 EXPLICIT validation and no repeated source application. |
| `src/numerics/source_input_2d.f90` | Every `*_cadence` variable, namelist entry, parsing and validation | Retain `METHOD`, `GROUP`, physical coefficients and active-model selection; revise input documentation. |
| `src/physics/source_configuration_2d.f90` | Two-list `stages`/`steps` builder, stage ordering, post-RK ordering, forced empty explicit-stage descriptor, cadence validation | Build one normalized group list; keep independent model method validation and physical grouping restrictions. |
| `src/numerics/step_source_executor_2d.f90` | Post-RK `execute_step_source_groups` scheduling and special-case RK adapter calls | Refactor into generic group executor, possibly rename. **Preserve** atomic mass-exchange transactions, `apply_mass_exchange_increment`, deposits/erosion/erodible inventory and topography updates. |
| `src/physics/source_terms_2d.f90` | STAGE/STEP-specific wrapper names and coupling where redundant | Keep `explicit`, `semi_implicit`, `implicit`, `analytic` model interfaces; preserve geometry/injection and substrate context; unify only dispatch layers without discarding physical kernels. |
| `src/core/parameters_2d.f90` | `n_RK` | Keep all non-RK physical and numerics settings. |
| `src/io/inpout_2d.f90` | `N_RK` and `CADENCE` I/O, RK stage validation and obsolete messages | Preserve input/restart consistency. Define backward compatibility policy for obsolete namelist keys and review restart metadata independently. |
| `src/core/simulation_2d.f90` | RK-then-post-RK two-phase flow; `TIMESTEP_EVENT_RK_COMPLETE` and `TIMESTEP_EVENT_STEP_SOURCES_COMPLETE` names | One Lie call; retain source-time-window preparation, stochastic preparation, diagnostics, and clock advance; replace events with meaningful spatial/group events. |
| `src/IMEX_SfloW2D.f90` | Stage registry lookup for boundary injection and other RK-only main-driver references | Preserve physical initialization, CFL computation, domain refresh, restart and output logic. |
| `src/numerics/equation_metadata_2d.f90` | **TO_VERIFY:** partition machinery used only to distinguish global RK explicit/implicit channels | Keep any group-local reduced-Newton equation selection and Jacobian indexing. |
| `src/spatial/spatial_operator_2d.f90` and `src/spatial/*` | **None planned** except obsolete RK-facing adapters discovered by audit | Preserve FV reconstruction, face fluxes, boundary handling and CFL; its output must not go through physical source limiter. |
| Physics families: `src/physics/rheology_laws/*`, `mass_exchange/*`, `pore_pressure/*`, `mass_injection_terms/*`, `thermal_laws/*` | **No broad deletion planned** | Keep model physics and all supported numerical interfaces, rheology anti-reversal ownership, analytic Darcy with full conservative changes, and physical inventory limiters. |

## 3. Phased implementation checklist

Use statuses `TODO`, `IN_PROGRESS`, `BLOCKED`, `DONE`. Add date, commit and test evidence when changing status.

| ID | Phase / deliverable | Status | Evidence / dependencies |
| --- | --- | --- | --- |
| L00 | Baseline: record build/CI result and search all callers/namelists/restarts/tests for RK, STAGE, STEP and CADENCE | TODO | Static inventory above is not a complete call-graph/test audit. |
| L01 | Single method-capable source-group registry and validation, including always-explicit GROUP 0 | IN_PROGRESS | Single Lie registry implemented; final validation pending. |
| L02 | Remove CADENCE from catalog, source configuration, source input and I/O | IN_PROGRESS | Production configuration and input moved away from STAGE/STEP; examples and tests pending. |
| L03 | Implement RK-independent group-local backward-Euler Newton API and analytic dispatcher | IN_PROGRESS | Renamed internal group/input-state symbols and removed misleading stage/DIRK documentation; numerical coupling review pending. |
| L04 | Generic explicit/NEWTON/ANALYTIC source-group executor | IN_PROGRESS | Lie source executor in place; shared rate-limiter and coupled inventory verification pending. |
| L05 | Lie time driver with GROUP-0 spatial update, source rates evaluated on `q_n`, source-rate limits applied to `q_sp`, then GROUP 1+ | IN_PROGRESS | Lie spatial-first driver in place, GROUP 0 rates evaluated from original state; spatial and group state reject non-finite values; final acceptance pending. |
| L06 | Simplify production timestep, observers, diagnostics and stochastic/source-window plumbing | IN_PROGRESS | Removed RK_COMPLETE observer and unused global equation partition from Lie driver; tests and events still require final review. |
| L07 | Delete IMEX-RK tableau, stage arrays, RK stage loops, adapters and `N_RK` across build/input/tests | IN_PROGRESS | RK tableau and workspaces removed; obsolete references/tests pending. |
| L08 | Exhaustive reference cleanup, formatting and documentation updates, D18 alignment | TODO | No obsolete symbols or contradictory architecture descriptions in current-source docs. |
| L09 | Comprehensive acceptance: unit, integration, conservation, stability, restart, examples and performance | TODO | Explicit executed suite names, commits, failures and resolutions. |

Recommended implementation order is L00–L07, followed by L08–L09; isolated cleanup can be moved forward only if builds remain bisectable and tested. Avoid simultaneous destructive changes in multiple modules without an independently compiling checkpoint.

## 4. Acceptance tests required

- **GROUP 0 common evaluation state:** spatial advance uses `q_n`; every GROUP-0 source rate uses `q_n`, *even though source increments are applied to `q_sp`*. Construct a test where `q_sp != q_n` and prove this distinction.
- **Spatial-only validation:** compute `q_sp`, inspect admissibility, and assert no source-rate limiter receives the FV spatial rate or changes a face-flux divergence. Spatial CFL/positivity and conservation tested independently.
- **GROUP-0 limited source update:** multiple processes competing for the same inventory jointly limited using separate contributions; apply once, verify acceptance and no double counting; handle source demands whose full `dt` increment cannot fit in `q_sp`.
- **Sequential groups:** each GROUP 1+ sees its immediate predecessor's accepted state; full `dt` per group; clock increment once. Changing group order changes only the expected splitting result.
- **Method freedom:** test EXPLICIT, NEWTON and ANALYTIC for the same physical family when each is supported; reject unsupported combinations, including non-explicit GROUP 0.
- **Newton:** residual `q_out-q_in-dt*S(q_out)`, complex-step Jacobian, reduced unknowns, failure reporting, no unphysical post-solve rate projection.
- **Analytic:** physical output applied once with no RK weighting, including complete Darcy conservative gas/mass/momentum/energy/pore-pressure updates when the model prescribes them.
- **Mass exchange:** same physically limited transfers applied to fluid, deposit, erosion, erodible substrate and optional changing topography, with coherent failure/commit.
- **Physical regressions:** dry/wet transitions, hydrostatic balance, gravity/friction equilibrium, injection timing, stochastic forcing, OpenMP and boundary conditions.
- **I/O and repository:** builds, tests, examples, input generators, restart, output diagnostics, memory/performance comparison with old RK; no runtime path requires an RK stage or STAGE/STEP cadence.

Do not claim that first-order Lie is unconditionally stable or well-balanced for arbitrarily stiff source/spatial couplings.

## 5. Work log (append newest entries here)

| Date | Change | Commit | Checks executed | Notes |
| --- | --- | --- | --- | --- |
| 2026-11-01 | Decoupled Lie driver from unused global equation partition and replaced misleading RK-completion observer with Lie-completion event. | Pending commit | Static inspection only; compilation intentionally deferred | Full step-level atomicity still open, including geometry state and failure rollback. |
| 2026-11-01 | Strengthened source transaction preconditions: finite timestep, valid source coverage, complete mass-exchange transfer finiteness, valid group limiter and post-transfer erodible inventory checks. | Pending commit | Static review, not compiled/tested | Final coherent rollback and model-specific physical tests remain pending. |
| 2026-11-01 | Removed remaining stage-group terminology from the cell-local Newton solver, clarified Lie-group input-state naming and checked explicit rates/transfers for non-finite values before committing inventory updates. | Pending commit | Static inspection only; no compilation/tests | Additional coupled-source and rollback review required. |
| 2026-11-01 | Strengthened spatial and source-group admissibility: reject infinities as well as NaNs; enforce solid-mass <= mixture-mass after every completed source-group cell update. | Pending commit | Static review only; no compilation/tests (migration policy) | No post-update clipping; all failures are explicit. |
| 2026-11-01 | Large in-place migration: replaced RK tableau with Lie driver, removed STAGE/STEP catalog and registry, removed N_RK from primary input, routed sources through group executor, removed RK Newton wrappers. | 5dfa59b3, a47087c1, b7ca9c28, bad22a6d, 8f757be7, 67edbf78, 9ce49244 | Deliberately not run | Major integration and physics-consistency work remains. |
| 2026-11-01 | Replaced transitional Lie configuration bridge with a direct METHOD/GROUP builder; spatial GROUP 0 remains mandatory and explicit. | 46ba9be3 | Not run | Production remains on the legacy integrator; further validation of defaults and collision rules required. |
| 2026-11-01 | Added opt-in source-configuration bridge to unified Lie registry; conflict detection and mandatory spatial group. | This commit | Not run | Transitional: depends on legacy metadata; replace with direct physical METHOD/GROUP builder in L02. |
| 2026-11-01 | Introduced opt-in Lie source group registry, mandatory explicit GROUP 0, ordered IDs, method capability and duplicate-source checks. | This commit | Not run | Legacy STAGE/STEP callers remain; L01 IN_PROGRESS. |
| 2026-11-01 | Added first RK-independent public backward-Euler group entry point and migrated existing STEP Newton caller; full Lie driver not yet present. | This change | Not executed (remote edit) | Compile/test verification required before marking L03 DONE. |
| 2026-11-01 | Created migration guide and pending checklist from static code inventory. | This document's initial commit | None (documentation only) | All code phases remain TODO. |

## 6. Open engineering details — resolve during implementation

1. **GROUP-0 limiter inventory:** rate evaluation is fixed at `q_n`; candidate source update is applied to `q_sp`. Establish which available conservative quantities and substrate inventories cap each source and how to reject/update consistently, without applying a limiter to spatial rates.
2. **Post-group admissibility:** identify which checks are diagnostic, which corrections belong inside a physical source model, and when to reject/retry a full timestep. Never silently alter accepted Newton solutions or conservative FV fluxes.
3. **Group metadata vs physical process:** GROUP 0 spatial operator is mandatory but need not be represented as a fake physical source catalog item. Decide simple representation and avoid a second grouping layer.
4. **Implicit grouping:** retain group-specific reduced unknown equations; evaluate auxiliary semi-implicit rates consistently with the Newton group without reintroducing global IMEX stages.
5. **Analytic groups:** preserve the current one-physical-process analytic-group constraint unless a mathematically justified combined exact flow is supplied.
6. **Input compatibility and restart:** decide handling of legacy `N_RK` and `*_CADENCE` keys and existing restart files before updating examples.
7. **Non-FV physical source time:** evaluate all GROUP-0 sources from the same `q_n` but pass source time consistent with the timestep, time-window averaging and any stochastic sample.

**Closure rule:** D19 is done only when all L00–L09 are DONE, acceptance evidence is linked, and D18 matches the implemented architecture. D19 must never itself substitute for executable evidence.
