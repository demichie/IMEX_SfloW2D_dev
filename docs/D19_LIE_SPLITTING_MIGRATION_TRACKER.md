# D19 — Lie splitting migration: implementation guide and live tracker

**Status:** IN PROGRESS — architecture rewritten in production modules; no final compilation or tests yet.
**Branch:** `master`
**Companion architecture:** [D18](D18_IMEX_RK_SOURCE_ARCHITECTURE_PRINCIPLES.md). D18 describes the numerical contract; **D19 tracks the implementation**, file-by-file changes, dependencies, evidence, and remaining work.
**Last reviewed:** 2026-11-01
**Rule for maintenance:** commits may deliberately be non-compilable during the large migration. Do not require intermediate builds; compile and run the complete acceptance suite once the end-to-end architecture and references are finished. Record actual commits and checks; do not mark a task DONE without evidence.

## 0. Frozen ownership and group-transaction rules (2026-11-01)

The following architectural decisions are **binding for subsequent production refactoring**. They supersede any transitional code organization documented below. Implement them completely before compiling or migrating tests.

### Module ownership

- `spatial_operator_2d` and its spatial/flux collaborators assemble **only spatial terms**. No spatial rate, flux or FV increment passes through a physical source limiter.
- `physics/source_terms_2d.f90` assembles **only physical source rates** (including group-combined contributions, mass-exchange transfers and source-related physical geometry inputs) and **owns their rate limiters**. Evaluating a source against a bed slope does not make that source part of the spatial operator.
- `state_2d` owns conservative/primitive **state validity and admissibility**, without silently limiting physical source rates.
- `geometry_2d` owns substrate/deposit/erodible and topography update operations, geometric cache refresh (slopes, curvature, etc.) and the rollback-capable geometry transaction.
- `lie_source_executor_2d` sequences groups and commits already assembled, physically consistent results. It must not independently compute physical rate limiters or silently apply partial transfer scalings.
- `nonlinear_solver_2d` solves fixed-geometry, fixed-inventory backward-Euler residuals on the complete structurally affected conservative equation set. Newton residual/Jacobian evaluations and line-search trials **never commit** substrate, deposit, erodible or topography changes.

### One coherent limiter for one coupled explicit group

After assembling **all** explicit physical contributions of group `G`, `source_terms_2d` determines a single admissible factor `lambda_G` between zero and one. The availability constraints are evaluated against the actual **recipient** fluid state and substrate inventory. For GROUP 0 specifically, evaluate rates on `q_n`, but limit their application against the post-spatial state `q_sp` and the current inventory. Never apply this limiter to the spatial rate itself.

The **same** factor must scale every coupled component of the group transaction: conservative source vector, solid and carrier deposition/erosion, change in erodible material, gas/mass and energy transfers, lava inventory and topography/bed rate. Symbolically: `S_G_limited=lambda_G*S_G`, `Ddot_limited=lambda_G*Ddot`, `Edot_limited=lambda_G*Edot`, `Bdot_limited=lambda_G*Bdot`. No component may be multiplied twice or left unscaled. The source assembler should return **an already limited, complete transaction**, rather than requiring the Lie executor to apply a separately returned factor. The physical coupling must be preserved even if different processes share a group. This applies to explicit-source rate limiting; do not invent a post-hoc Newton or analytic increment limiter that would invalidate the chosen integrator.

### Group boundaries are the only substrate/geometry commit points

For each ordered Lie group:

1. Freeze its input fluid state, deposit/erodible inventories, and current geometry while assembling sources and, when selected, evaluating Newton iterations or analytic operators. The group uses the geometry left by **the preceding completed group**.
2. Build the group's candidate fluid and substrate/topography changes consistently. Check state admissibility through `state` and physical transfer validity through the responsible physics/geometry modules.
3. **Only after the whole group succeeds**, commit its deposit, erosion, erodible and bed changes consistently. If topography changes, apply the projected update and refresh all derived geometry before starting the next group; later groups must see the updated slopes and curvatures.
4. If the group makes **no** substrate/topographic changes, skip unnecessary geometry/inventory updates. Do not refresh topography per cell or during a Newton iteration.
5. Retain a transaction spanning the **entire Lie timestep**. A failure in any later group restores the pre-timestep fluid state, substrate inventories, authoritative bed and geometric caches, even if earlier groups had updated geometry for subsequent groups.

The group boundary is the visibility point for these coupled changes. Per-cell candidate preparation may be temporary, but must not expose partially committed substrate/topography to another cell or to the solver while the same group is being evaluated. Group-0 spatial admissibility remains an independent check without clipping or rate limiting.

### Required production changes, before tests

- Refactor `source_terms_2d` to return a complete, jointly limited explicit-source transaction (not raw physical and exchange rate fragments plus a limiter factor for the executor to combine).
- Remove source-rate limiter calculations and mass-exchange-specific coverage/transfer scaling from the Lie executor, giving them to physical source assembly.
- Stage deposit/erodible and topography changes across a group, then publish/refresh geometry at group completion. Keep Newton's geometry and inventories immutable across its evaluations.
- Preserve whole-step rollback, and retire obsolete RK-era source APIs/comments as part of the production cleanup. **Only once this production structure is complete** should build and tests be resumed.

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

## 1.1. Frozen Newton unknown-set contract (decision 2026-11-01)

For every Lie source group whose configured method is **NEWTON**, define the conservative equation index set

`I_G = { i : the group's source rate S_G,i(q) is structurally nonzero }`.

**Structural** means that the source process can produce a nonzero contribution to that conservative equation for an admissible physical state; a rate that happens to evaluate to zero in the current cell does **not** remove the equation from `I_G`. For a group containing multiple sources, `I_G` is the union of all equations potentially modified by any constituent source, including induced total mass, momentum, thermal energy, gas, solid and transported-variable changes. Do not infer the set by sampling rate values or by looking only at the process's primary physical quantity.

1. Solve **all and only** equations in `I_G` simultaneously using backward Euler, with complete conservative candidate `q` supplied to the existing physical source-rate evaluator:
   `F_i(q) = q_i - q_in,i - dt*S_G,i(q) = 0` for every `i in I_G`.
2. All other conservative equations remain exactly at their incoming group values: `q_j=q_in,j` for `j not in I_G`, even if the rates in `I_G` depend on these fixed values.
3. Assemble the complex-step Jacobian on **rows and columns indexed by `I_G`**: `J_ij = delta_ij - dt*dS_G,i/dq_j`, `i,j in I_G`. Extract these entries with the equation-index maps already used in the solver; do **not** physically permute the global state/Jacobian.
4. **Do not eliminate additional nonzero-rate equations** merely to reduce the Newton dimension. In particular, no physics-specific reconstructed-variable equations, assumed zero off-diagonal blocks, general Schur-complement or block-factorization framework are required. A coupled process that changes solid mass, total mass, momentum and energy solves all affected conservative equations, not only the solid-mass unknown.
5. The numerical Newton solver remains **physics-agnostic**: it accepts a validated group, its complete structurally affected equation indices, the generic real/complex source-rate evaluators and numerical parameters. The catalog/registry/configuration layer supplies and validates the affected-equation set. It must include coupled secondary conservative changes and must not be defined by whichever variables appear convenient as Newton unknowns.
6. Check convergence and full backward-Euler residual on every equation in `I_G`; verify that equations outside `I_G` remain unchanged. Preserve existing small dense linear solvers, pivoting, line search, complex-step differentiation, failure reporting and physical admissibility checks.

**Example:** if sedimentation can remove a solid constituent and therefore affect total mixture mass, both momenta and thermal energy, all those conservative equations belong to `I_G`, even though the process is driven by solid loss. This contract avoids deriving process-specific analytic reconstruction rules and keeps the solver independent of physical models.

**Implementation audit:** check that each source catalog's declared affected-equation indices cover **all** structurally nonzero conservative rates, including indirect coupled contributions. Confirm that `source_group%equation_ids` represents this full union before simplifying `nonlinear_solver_2d.f90`.

**Metadata implementation (2026-11-01):** The catalog exposes `source_newton_affected_equations(name)`; group configuration forms a unique union and the registry rejects missing/extra indices. Currently Newton-capable rheology changes both momenta, and cooling changes energy; pore-pressure dissipation is analytic-only and therefore never entered in a Newton group. Any new Newton-capable process with mass or multi-equation exchange must declare its full conservative rate footprint here before enabling its capability.

**Cleanup 2026-11-01:** Removed duplicated simulation-global and executor-local equation partitions from the production Lie/Newton call path. `source_group%equation_ids` now provides the index set directly; `equation_metadata_2d.f90` remains for the final legacy/test reference audit. No intermediate compilation or tests.

**Reference audit following partition cleanup:** `src/physics/equation_terms_2d.f90` contained a now-unused `USE equation_metadata_2d`, removed. The registry's one-use `in_list` helper was inlined in `same_set`. The legacy `equation_metadata_2d.f90` module is still referenced by old test targets and `src/Makefile.am`; do not delete it until those tests/build targets are migrated or removed as part of L07–L08. No intermediate compilation or tests were run.

## 2. Scope: KEEP / CHANGE / DELETE inventory

Based on static inspection of the branch at the planning stage. **All code operations below are pending.** File-specific names must be audited for remaining references before deleting.

| Path | DELETE (RK/cadence-only) | CHANGE / KEEP |
| --- | --- | --- |
| `src/numerics/time_integration_2d.f90` | `n_RK` handling, `validate_runge_kutta_stages`, Butcher tableaux (`a_tilde_ij`, `a_dirk_ij`, weights), `q_rk`, `qp_rk`, `EX_terms`, `IM_terms`, stage-indexed `spatial_terms`, stage loops, `imex_RK_solver`, RK observer events | Implement Lie driver, preserve `q_n` through GROUP 0, one spatial update, independent spatial admissibility, group-0 source-only limiter/update, then group loop and post-group checks. Keep only needed spatial/source workspaces. |
| `src/numerics/nonlinear_solver_2d.f90` | `solve_imex_stage`, RK-only `solve_rk_step` adapters, `effective_im_term`, `a_diag`, `known_stage_residual`/`Rj_not_impl`, `n_RK==4` tolerance special-case | Retain reduced conservative-variable Newton, complex-step Jacobian, line search, diagnostics and model interfaces. Backward-Euler group residual `F(q_out)=q_out-q_in-dt*S(q_out)`. The Newton unknown set is exactly the union of structurally nonzero conservative source-rate equations; do not eliminate coupled nonzero-rate components or introduce model-specific reconstruction. |
| `src/numerics/source_catalog_2d.f90` | `SOURCE_CADENCE_STAGE/STEP`, `default_cadence` | Retain physical source catalog, capabilities and `SOURCE_METHOD_*` metadata. Model/method compatibility remains validated. |
| `src/numerics/source_registry_2d.f90` | Group `cadence` field; separate `stage_entries`/`step_entries`; `stage_groups()`/`step_groups()`; stage/step-only validation and placeholders | Single ordered group registry with group ID, method, source names, full set of structurally affected conservative equations for Newton; GROUP-0 EXPLICIT validation and no repeated source application. |
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
| L03 | Implement RK-independent group-local backward-Euler Newton API and analytic dispatcher | IN_PROGRESS | Catalog defines structurally affected indices and validates exact unions; removed redundant global and local equation partitions from production Lie Newton. Final reference and coupled-rate audit remain. |
| L04 | Generic explicit/NEWTON/ANALYTIC source-group executor | IN_PROGRESS | Newton and analytic solver failures now propagate via status and trigger rollback of geometry and inventories; other fatal validation paths and retry policy remain. |
| L05 | Lie time driver with GROUP-0 spatial update, source rates evaluated on `q_n`, source-rate limits applied to `q_sp`, then GROUP 1+ | IN_PROGRESS | Fluid q/qp now advance in trial workspaces and publish after all Lie groups return; substrate/topography transaction and final acceptance remain pending. |
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
- **Newton:** residual `q_out-q_in-dt*S(q_out)` on the **complete structurally affected conservative equation set**; verify coupled mass, constituent, momentum and energy rates are not omitted; inactive equations remain fixed; complex-step Jacobian on the same index rows/columns; failure reporting and no unphysical post-solve rate projection.
- **Analytic:** physical output applied once with no RK weighting, including complete Darcy conservative gas/mass/momentum/energy/pore-pressure updates when the model prescribes them.
- **Mass exchange:** same physically limited transfers applied to fluid, deposit, erosion, erodible substrate and optional changing topography, with coherent failure/commit.
- **Physical regressions:** dry/wet transitions, hydrostatic balance, gravity/friction equilibrium, injection timing, stochastic forcing, OpenMP and boundary conditions.
- **I/O and repository:** builds, tests, examples, input generators, restart, output diagnostics, memory/performance comparison with old RK; no runtime path requires an RK stage or STAGE/STEP cadence.

Do not claim that first-order Lie is unconditionally stable or well-balanced for arbitrarily stiff source/spatial couplings.

## 5. Work log (append newest entries here)

| Date | Change | Commit | Checks executed | Notes |
| --- | --- | --- | --- | --- |
| 2026-11-01 | Froze the physics-agnostic Newton unknown set: every conservative equation with structurally nonzero source rate, indexed directly without permutation or Schur complement; untouched equations fixed. | Pending commit | Design decision documented; implementation and final tests pending | Audit source catalog's affected-equation maps for coupled changes before Newton simplification. |
| 2026-11-01 | Added source-group entry validation for primitive states, stochastic values, terrain slopes/curvatures/gravity and erodible inventories before physical kernel calls; kept unsupported model/method errors fatal configuration errors. | Pending commit | Static review only; no build/tests | Kernel-internal arithmetic errors remain possible and require physical-specific status interfaces or final validation. |
| 2026-11-01 | Propagated detailed source-group rejection messages (including Newton/analytic diagnostics and bed projection) through Lie rollback, with finite analytic-rate validation. | Pending commit | Static review only; compilation deferred | Underlying physical kernels may still terminate directly; exhaustive failure-path audit remains. |
| 2026-11-01 | Snapshot and restore exact derived topography caches alongside B_vertex, avoiding reconstruction during rollback when reconstruction itself may be invalid. | Pending commit | Static review only; tests/build deferred | Fatal errors within physical kernels and full retry remain outside this change. |
| 2026-11-01 | Checked all derived bed fields after Lie bed reconstruction, so a non-finite geometry result requests rollback; also reject inadmissible conservative composition before entering real constitutive decoding. | Pending commit | Static review only, no build/tests | Remaining constitutive routines may produce physically invalid but finite values; full error-path coverage needs final test suite. |
| 2026-11-01 | Added optional status-bearing conservative-to-primitive conversion and routed Lie group conversion failures to timestep rejection; spatial conversion reports explicit rejection before geometry transaction starts. | Pending commit | Static review; no build/tests | Constitutive helper failures and derived geometry still need review; this does not prove universal rollback. |
| 2026-11-01 | Reject negative transported additional-gas masses before primitive conversion in spatial and source candidate states; source-side failures use rollback. | Pending commit | Static only; no build/tests | Other physical conversion failures can still terminate internally; no blanket claim of full recoverability. |
| 2026-11-01 | Made candidate nodal-bed update return recoverable rejection to Lie executor and step rollback, preserving existing fatal API for callers without status. | Pending commit | Static review only; no build/tests | Geometry reconstruction and primitive conversion can still fatal-stop internally; full recovery requires status-aware lower-level APIs. |
| 2026-11-01 | Route explicit Lie cell validation failures through the same step_failed return path as Newton/analytic failures; the caller rolls back the authoritative bed and substrate snapshot. | Pending commit | Static only; no build/tests | Configuration errors and internal conversion ERROR STOP remain fatal. Automatic retry is not implemented. |
| 2026-11-01 | Propagate analytic/NEWTON group failure to Lie driver; rollback geometry and substrate snapshot before fatal rejection, without publishing q/qp trial fields. | Pending commit | Static review only; no tests or build | Not yet universal: fatal validators still stop directly and automatic retry is not implemented. |
| 2026-11-01 | Added geometry/inventory pre-step snapshot with explicit commit and rollback operations, called around the Lie group sweep. Intermediate bed updates remain visible to successor groups. | Pending commit | Static review only, no compilation/tests | Fatal errors still terminate instead of invoking rollback; full recoverable all-state rejection is not complete. |
| 2026-11-01 | Validate complete candidate B_vertex and projected bed rates before changing the authoritative nodal topography; reject non-finite cell proposals. | Pending commit | Static review only, no build/tests | This is geometry-update staging, not full-step isolation; derived geometry refresh and substrate globals remain nontransactional. |
| 2026-11-01 | Moved all spatial and Lie-source fluid updates into dedicated timestep q/qp trial arrays; publish both only after the full Lie executor returns. | Pending commit | Static review; no build/tests | Substrate and geometry remain global during execution, so full-timestep atomicity is NOT yet guaranteed. |
| 2026-11-01 | Stage and validate complete cell-local explicit Lie fluid, primitive and substrate candidates before committing; remove partial candidate writes on failed inventory checks. | Pending commit | Static review only; compilation intentionally deferred | This is cell-local staging, not full-timestep rollback; bed geometry remains globally committed per mass-exchange group. |
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

### Legacy equation-partition retirement (2026-11-01)

The legacy `src/numerics/equation_metadata_2d.f90` module and its dedicated `TESTS/TEST_EQUATION_METADATA` unit test were removed, with object/source entries deleted from `src/Makefile.am`. The production Lie Newton solver and source registry already use `source_group%equation_ids` directly.

**Outstanding obsolete test callsites (must migrate/remove before final compile):** `TESTS/TEST_SOURCE_REGISTRY/test_source_registry.f90`, `TESTS/ACCEPTANCE/test_d17_production_newton.f90`, `TESTS/TEST_RHEOLOGY_CELL_LOCAL/test_grouped_implicit_stage.f90`, `TESTS/TEST_RHEOLOGY_CELL_LOCAL/probe_sources.f90`, `TESTS/TEST_MASS_EXCHANGE/test_mass_exchange.f90`, and `TESTS/TEST_IMEX_STAGES/*`. They reference obsolete equation partitions, source configuration signatures, or stage/step entrypoints. Do not reintroduce the removed module to keep them working; replace RK-specific tests with Lie-specific validation in L07–L09. **No build/tests run**.

### IMEX-RK test retirement (2026-11-01)

Removed all `TESTS/TEST_IMEX_STAGES` files (stage-count/RK tableau/reconstruction-stage checks); removed its physics-regressions CI matrix entry. Retired obsolete `test_grouped_implicit_stage.f90` and `probe_sources.f90` plus their invocation/characterization from the rheology runner. Kept independent rheology model, analytical velocity, cooling, and pore-pressure tests. The remaining `TEST_SOURCE_REGISTRY/test_source_registry.f90`, `TEST_MASS_EXCHANGE/test_mass_exchange.f90` and `ACCEPTANCE/test_d17_production_newton.f90` still contain legacy partition and IMEX calls and **must be replaced with Lie tests**, not left in final CI. No intermediate compilation/tests.

### Lie test conversion (2026-11-01)

Replaced the RK-era source-registry regression with Lie group-order, method and Newton footprint checks, including co-grouped rheology/cooling. Replaced the production Newton acceptance fixture with a group-local backward-Euler solve and residual validation using `solve_backward_euler_group`. Both require final compilation and acceptance. The large N8 mass-exchange harness still contains obsolete RK/STEP APIs and requires separate test restructuring.

### Mass-exchange legacy harness retirement (2026-11-01)

Removed `TESTS/TEST_MASS_EXCHANGE/test_mass_exchange.f90` because its RK-stage solver, cadence, legacy equation partition and step-source executor are deleted. The mass-exchange runner now executes the retained independent projection/limiter/process-momentum checks. **Important coverage gap:** the old Fortran harness also exercised inventory budgets, checkpoint restart and 1-/4-thread reproducibility. Those checks are **not replaced yet**; L09 must add direct Lie mass-exchange, topography transaction, restart and multithreading tests. Retained `test_projection_stress.f90` for subsequent migration to `execute_lie_source_groups` (it still refers to the removed step-source executor and is not invoked). No compilation/tests performed.

### Restored independent projection regression (2026-11-01)

Refactored `test_projection_stress.f90` to test only the independent production geometry projector: smooth, wet/dry and signed random fields, with conservative-flow invariance, bed/cached-geometry consistency and nodal integral checks. Removed calls to deleted STEP source executor and historical RK count. The mass-exchange runner now compiles and runs that projection test in strict/optimized profiles on 1 and 4 OpenMP threads and compares their binary snapshots, in addition to retained physics source-check scripts. **This does not restore mass-exchange inventory/restart transaction coverage**, which remains L09. No intermediate compilation or execution performed.

### Lie mass-exchange transactional regression (2026-11-01)

Added a cell-grid Lie GROUP-0 erosion regression to `test_projection_stress.f90`: initializes a real registry, evaluates original conservative-to-primitive states, calls the actual `execute_lie_source_groups` with `q_original` and `qp_original`, checks conservative solid-mass / erosion inventory transfer, erodible depletion, nonnegative inventories, finite results and refreshed geometry. The runner already executes this fixture in strict/optimized and 1-/4-thread configurations. Final compilation, runtime acceptance, and restart/rollback tests are still pending. No intermediate compilation/tests.

### Rejected Lie mass-exchange regression (2026-11-01)

Added a late-cell nonfinite-input rejection fixture in `test_projection_stress.f90`, with a real group-0 erosion registry and the production Lie source executor. The test opens a geometry transaction, forces the second cell's original primitive input invalid, expects recoverable `step_failed` and a diagnostic, invokes `rollback_lie_geometry_transaction`, and checks bed, cached geometry and all inventory snapshots exactly. Fluid trial state is discarded by the test as in `advance_lie`. **Actual restart serialization/continuation remains outstanding**; do not count this as restart coverage. Final compile and runtime verification remain deferred.

### Production Newton separation (2026-11-01)

Removed the obsolete multi-group implicit/analytic sweep and semi-implicit momentum seed from `nonlinear_solver_2d.f90`. The exported `solve_backward_euler_group` now validates the single NEWTON group's unique equation index set and directly calls the existing complex-step Newton kernel over precisely that set; all non-group conservative components remain fixed and failed solves restore the group input. The Lie executor retains the EXPLICIT / NEWTON / ANALYTIC dispatch; no model-dependent reconstruction, RK stage, or separate equation partition is introduced. **No tests or compilation until architecture completion.**

### Lie executor cleanup: coherent coverage and admissibility (2026-11-01)

In `lie_source_executor_2d.f90`, consolidated source-only conservative-mass admissibility checks in a side-effect-free predicate, used before committing explicit candidates and after group advances. Centralized source-coverage scaling of the *complete* `mass_exchange_rates_type` transfer record (including coupled mixture, solid/gas mass, energy, carrier erosion/deposition, withdrawal volume, lava and bed rates), while scaling the conservative exchange-rate vector exactly once. The group-0 spatial operator remains independently advanced and **never rate-limited**. This is production-only refactoring; tests/build deferred.

### State-validation ownership (2026-11-01)

Moved `conservative_mass_admissible` from `numerics/lie_source_executor_2d.f90` to the public state API in `state/state_2d.f90`. Both the Lie spatial-only admissibility check in `time_integration_2d.f90` and source-group candidate/post-group checks invoke the same side-effect-free state predicate. Numerical orchestration no longer embeds the conservative mass-budget validity policy. No rate limiting or clipping was introduced for the spatial update. Build/tests intentionally deferred.

### Lie state checks: primitive-state ownership and duplicate removal (2026-11-01)

Added `physical_cell_state_finite` to `state_2d` and routed the Lie source executor's original/predecessor primitive-state check through this state-owned predicate. Removed a duplicated conservative-candidate finite-value check already covered by `conservative_mass_admissible`. Rate finiteness, geometry finiteness, and physical substrate inventory checks remain distinct from state validation. No intermediate build or tests.

### Physical rate limiter ownership (2026-11-01)

Moved the coupled explicit-source rate limiter out of `numerics/lie_source_executor_2d.f90` into `physics/source_terms_2d.f90` as `limit_explicit_source_rates`. The physical assembler determines the single shared scale factor from the post-spatial reference conservative inventory, mass components, remaining carrier, and erodible substrate inventory. The Lie executor only applies that factor to the already assembled source rates and synchronized substrate transaction. No spatial-rate limiter is introduced. State admissibility remains in `state_2d`. No tests/build.

### Production explicit Lie transaction assembly (2026-11-01)

The physical `source_terms_2d` module now defines `assemble_lie_explicit_transaction`: it combines all explicit contributors, applies the source-coverage correction to **all** coupled mass-exchange rates, determines the common physical limiter against the *recipient* conservative state/inventory, and scales both conservative and substrate rates exactly once. It does not modify fluid, inventory or geometry fields. `lie_source_executor_2d` requests that fully limited transaction and applies it without re-implementing physical rate assembly or applying a second factor. Group-end geometry projection remains separate. **Still pending:** staging all candidate inventory changes until group completion, geometry-only ownership of substrate updates, and the final build/test phase.

### Group-boundary staging of substrate inventories (2026-11-01)

Production commit `2885c92e` stages every mass-exchange cell's candidate deposit, accumulated erosion, erodible inventory, lava inventory and already-limited bed rate in group-local arrays. Only after **every** cell succeeds are these substrate values published together and bed rates recorded for the group-end geometry projection. No per-cell inventory change is visible during evaluation of the remaining cells; later groups see completed group geometry including refreshed slopes/curvatures. Non-explicit mass-exchange groups are rejected as invalid configuration. Existing whole-timestep geometry snapshot/rollback covers deposits, erosion, erodible, lava, bed vertices and geometric caches. This is a **production-only intermediate step**: subsequent refactoring should move the final inventory commit interface into `geometry_2d`, without moving source-rate limitation out of `source_terms_2d`. No build/tests yet.
