# D18 — IMEX-RK source integration: design principles and target contract

**Status:** agreed design requirements; implementation and validation are **not** complete.  
**Scope:** source-term integration in IMEX-RK stages (`STAGE`) and post-RK source splitting (`STEP`).  
**Date:** 2026-10-31.

This document is the **design contract** for the next source-integration refactoring. It describes the desired behavior independently of the existing Fortran types, registry, and solver. Historical D17/D18 implementation notes are not reasons to preserve obsolete abstractions.

## 1. Design process: architecture first

1. Specify the required mathematics, data flow, interface contracts, execution order, and failure behavior **before** altering existing modules.
2. Design the **simplest sufficient target architecture from scratch**, without assuming that existing wrappers, partitions, registries, or adapters must survive.
3. Compare the target to the current implementation. Identify **keep / change / delete**, including obsolete tests and compatibility layers.
4. Prefer removing redundant layers, duplicated state, transformations, and indirection over introducing bridge wrappers to fit old assumptions.
5. Preserve only deliberately established public contracts (e.g. four numerical interfaces per physical model), not incidental implementation details.
6. Implement the target consistently across all callers; add focused tests for the target behavior and delete tests whose only purpose was to preserve removed internal interfaces.
7. Require build, numerical regression, conservation/limiter, sequencing and temporal-accuracy checks before declaring the refactor complete. Do not treat a committed change as validated without evidence.

## 2. Physical sources, groups, and conservative rates

- **Groups partition the *application of physical source processes*, never the conservative equations.**
- Every physical source interface returns a rate vector of length `n_eqns`; components not affected by that process are zero.
- A group identifies its source process(es), numerical `METHOD`, `CADENCE` (`STAGE` or `STEP`), and execution `GROUP` order. The group does **not** grant or prohibit writes to a subset of conservative equations.
- Different groups may act on the **same** conservative components. This is necessary for sequential source splitting.
- **Within one EXPLICIT group**, evaluate every member source on the **same group-entry state** `q^(g-1)`, apply each process's own limits as defined for that common state, and sum the resulting full conservative rate vectors. **Do not update the state between members of the same explicit group.** Update the conservative state **once** using the group's summed rate before executing the next group.
- A Newton solver may choose a **smaller set of unknowns** to compute a group's full conservative rate. Those solver unknowns are a numerical detail, **not** group membership or a restriction on source-rate components.
- Avoid representing a general source group through an exclusive equation membership map. If a Newton unknown map is useful, keep it internal to the relevant solve.

### Within-group versus between-group execution

For an EXPLICIT group `g` containing sources `k=1,...,m`:

```text
q_in = q^(g-1)
S_g = sum_k S_(g,k)(q_in)         ! All processes see identical input q_in
q^(g) = q_in + tau * S_g         ! One update, after summing
```

**Explicit GROUP 0 is reserved for spatial discretization:** the spatial flux/divergence operator (including its associated spatial contributions) belongs to explicit group `0` and is evaluated on the **initial state of the explicit stage evaluation**, before any modifications by subsequent explicit source groups. Its rate is not recomputed on the intermediate states of groups `1,2,...`. Other processes, if assigned to explicit group `0`, must also be evaluated on this same initial state; no within-group state update is permitted. After completing group `0`, subsequent explicit source groups follow sequential updated-state semantics.

**Distinct EXPLICIT groups are sequential:** evaluate each explicit group's rates using the conservative state already modified by all preceding groups, including preceding EXPLICIT groups. This applies to both `STAGE` and `STEP`; explicit groups are not evaluated independently against a common original state. The correct intermediate update coefficient at `STAGE` must be derived from the explicit IMEX-RK tableau rather than assumed equal to the implicit diagonal coefficient.

The sequential updated-state rule applies **between** groups. An explicit group's member processes are **not** a sequential sub-splitting. If their combined increment needs additional admissibility control, its behavior must be defined explicitly rather than silently changing the common-input evaluation contract.

## 3. Sequential splitting contract (both CADENCE values)

For each cell and each ordered source group `g`, let `q^(g-1)` be the **current conservative state** after the preceding group. All source evaluations and limiters of group `g` use this state (reconstructing primitive/physical variables as needed).

The group computes its limited conservative rate `S_g` and updates:

```text
q^(0) = starting conservative state
for g in configured order:
    S_g = compute_limited_group_rate(q^(g-1), tau)
    q^(g) = q^(g-1) + tau * S_g
S_total = sum_g S_g
```

For a common nonzero interval `tau`:

```text
S_total = (q^(last) - q^(0)) / tau
```

This identity is an **invariant to test**. Do not evaluate each group against the same original state. Do not reset previously updated components when moving to a subsequent group. Do not apply the same increments again during assembly.

- `STAGE`: execute the source-group sweep as part of the IMEX-RK stage, using the appropriate diagonal-stage interval (`tau = a_ii * dt` in the currently discussed diagonal stage design).
- `STEP`: execute an ordered source-group sweep after the RK step using `tau = dt`.
- Both modes have the **same sequential splitting semantics**, differing in placement and interval.
- If `tau = 0`, implementation must avoid division by zero and define a zero-increment result.
- Numerical order of the overall IMEX/RK plus sequential splitting is **not assumed**: verify using convergence tests.

### Shared inventories and limiter ownership (agreed rule)

- **Within the same group:** all member processes evaluate on the same group-entry conservative state. If two or more processes consume the **same finite inventory** (for example, carrier mass), their requested withdrawals must be **limited jointly** against that shared available inventory before the group's effective full-vector rate is assembled and applied. Independent process limits alone are insufficient.
- **Between distinct groups:** each group first evaluates its own rates, applies all necessary process and shared-inventory limiters, and updates the conservative state. The next group **must reconstruct its physical quantities, available inventories, rates and limits from that already limited, updated state**, not the previous group's input state. No additional cross-group shared-inventory limiter is needed **provided that every group guarantees an admissible update**.
- These rules apply to both `STAGE` and `STEP`, and to every permitted numerical treatment. A group may contain multiple explicit sources; an analytic group remains single-process. Shared limits are defined where the physical inventory is owned, not via a generic rate-clipping wrapper.
- This rule concerns **sequential group updates**. It does not prove admissibility of the final combination of different Runge–Kutta stages; stage-combination admissibility requires separate verification.

## 4. Methods and group exclusivity

- **ANALYTIC** groups contain **exactly one physical process** for either `STAGE` or `STEP`. They cannot contain or be merged with NEWTON or EXPLICIT processes under the same group identifier/cadence.
- The analytic physical model evolves the received state over `tau` (possibly with frozen coefficients), applies its own process-specific limiters and returns the **full conservative mean rate** `(q_after-q_before)/tau`.
- NEWTON groups solve their configured physical source contributions with the appropriate coupled unknowns, using the current sequential input state, and expose the complete limited rate.
- EXPLICIT and semi-implicit behavior must be specified consistently with the configured group/order and existing IMEX tableau; their detailed placement relative to the RK explicit stage is an **implementation/design verification point**, not permission to violate sequential state propagation.
- No artificial disjointness requirement may be imposed between different groups' affected equations or between a group's affected equations and another group's Newton unknowns.

## 5. Physical model interface contract

Keep the four standard model entry points even when a method is unavailable:

- explicit;
- semi-implicit;
- implicit/Newton;
- analytic.

Unsupported entries **fail explicitly** (e.g. `fatal_error`), not a silent zero source. Registry/preflight validation must reject configurations using unavailable methods before integration.

### Darcy pore-pressure process

- Darcy pore-pressure diffusion and optional interstitial carrier-gas withdrawal belong to the pore-pressure physical model, not to mass exchange.
- Darcy supports **ANALYTIC/STAGE** and **ANALYTIC/STEP**, with `ANALYTIC/STAGE` as the intended default. The three other numerical methods remain present but unsupported.
- Its analytic operator can change total mass, momenta, energy, gas components and conservative pore pressure. It always returns the **complete conservative rate**; no `equation_ids` filtering of its update.
- Limiting (including gas inventory and any physical admissibility bounds) is owned by the Darcy process, not retrofitted by the group registry.

## 6. Minimal target responsibilities (conceptual, not mandatory module names)

| Responsibility | Required behavior |
| --- | --- |
| Source descriptor/configuration | Declares available methods, cadence, group order and source membership; validates incompatible configurations. |
| Physical source model | Computes full conservative rates and process-specific admissibility/limiting. |
| Source-group executor | Applies ordered groups to the **updated** conservative state for STAGE/STEP; sums effective rates and prevents double application. |
| Newton numerical solver | Resolves only necessary unknowns for its source group and returns an effective **full** rate; cannot define group membership by equation ownership. |
| IMEX-RK integrator | Provides stage state and interval, consumes the stage's effective total source rate exactly once with correct RK bookkeeping. |

**No extra wrapper or metadata layer is justified merely to retain the old equation-partition architecture.** Decide actual modules/types only after tracing the minimal required data through callers.

## 7. Refactoring review and acceptance checklist

Before further structural code changes, write a concise inventory of current source executors, group metadata, equation partitions, stage rate assembly, and call chains. Draw a minimal target call/data flow and a **keep / change / delete** map.

Acceptance requires tests of:

1. Two sequential groups changing the same conservative component: second sees the first's updated state.
2. Both NEWTON → ANALYTIC and ANALYTIC → NEWTON at STAGE, and corresponding ordered STEP combinations where supported.
3. Full-vector rate and limiter conservation/admissibility, notably Darcy with and without gas loss.
4. Per-group and total-rate identities; **no duplicate increments** in RK stage or post-step assembly.
5. Rejection of mixed/duplicate ANALYTIC groups and unsupported numerical methods; retained unsupported entry points fail.
6. Pure explicit, pure Newton, pure analytic and mixed runs; relevant legacy physical regressions.
7. Temporal-convergence behavior of chosen IMEX tableau **with splitting**.
8. Build and automated tests, documenting precisely what was run and what remains unverified.

## 8. Decisions still to make by reading the full RK implementation

These are deliberately **not** frozen as architecture until independently checked:

- Exact semantics of EXPLICIT and semi-implicit stage groups within the Butcher tableau.
- How the RK integrator carries and combines effective per-stage source rates without altering order or applying updates twice.
- Whether any current `equation_partition_type` or intermediate wrapper serves a necessary independent role; prefer removing it if not.
- Internal representation of Newton-group unknowns when two sequential groups solve overlapping conservative variables.
- Zero-diagonal stages (`a_ii=0`) and consistent zero-interval behavior.
- Interaction between sequential splitting, stage-order accuracy and non-source RK flux updates.

**Rule of interpretation:** the numbered requirements above take precedence over old implementation constraints; unresolved points must be analyzed from the RK mathematics and current call graph before deciding a minimal code structure.

## 9. Architectural decisions to resolve before implementation (review 2026-10-31)

These are **open design decisions**, not implementation approvals. Resolve them from the IMEX-RK equations, then record each choice and its tests **before** restructuring code. In particular, do not retain an old type merely because an existing caller expects it.

| ID | Decision | Minimal option to assess | Required verification |
| --- | --- | --- | --- |
| A1 | **One STAGE execution path:** how does spatial discretization execute first in EXPLICIT GROUP 0, before sequential explicit source groups and the other methods? | One group executor with method-specific evaluations, rather than separate explicit accumulation and implicit loop. | Butcher-tableau consistency, sequential observation of updated state, complete rate bookkeeping. |
| A2 | **EXPLICIT at STAGE:** what stage-consistent intermediate state increment allows sequential explicit groups with a distinct explicit RK tableau? | Distinguish RK explicit evaluation from an actual state-advancing split operator where mathematically necessary; do not silently assign the same `a_ii dt` to both. | Derive stage and end-of-step formulas; check RK order on noncommuting and coupled source examples. |
| A3 | **Newton group unknowns:** what is the smallest representation of active unknowns when successive groups touch the same equations? | A group-local unknown-index array passed only to Newton; eliminate or radically shrink `equation_partition_type` if no independent user remains. | Overlapping-group solves in both orders; convergence and residual correctness. |
| A4 | **Single rate contract:** should group evaluation return only complete conservative `rate(n_eqns)`, or also a candidate state/status? | Prefer `rate(n_eqns)` plus admissibility/status; a single executor owns the `q += tau*rate` update. | Conservation, limiter semantics, no duplicate updates, error paths. |
| A5 | **STAGE effective rate:** how are all sequential increments represented in `IM_terms` and `EX_terms`? | Compute one effective increment from final minus starting stage state for the relevant tableau contribution, not a second independently applied source. | Stage reconstruction, no double counting, correct RK recurrence. |
| A6 | **STAGE vs STEP executor:** one shared loop or two? | Share the minimal `apply_group` operation and sequential state update; keep distinct orchestration only where RK/stage context or bed updates require it. | Equivalent operator semantics; avoid wrappers and duplicate dispatch. |
| A7 | **Group identity and order:** are `METHOD`, `CADENCE`, `GROUP`, and process names sufficient? | Remove source-group `equation_ids` for non-Newton methods; group identity comes only from source membership and configuration. | Exclusive analytic groups; shared conservative components; deterministic order. |
| A8 | **Semi-implicit method:** is it a standalone group method, a modifier to a Newton solve, or an independent numerical option? | Decide its actual mathematical role before representing it as an additional layer. | Existing friction/limiter behavior and correctness under sequential splitting. |
| A9 | **Zero-diagonal RK stages:** how do groups behave when `a_ii=0`? | No division by `a_ii dt`; distinguish zero-duration analytic/implicit increments from explicit RK source evaluations. | Schemes with explicit first stage, dry states, no invalid divisions. |
| A10 | **STEP mass exchange and bed state:** how does a conservative full-vector source group coexist with source transactions that also update bed/inventories? | Explicitly define a transactional group result when non-conservative-state side effects are necessary; do not hide bed mutation behind a rate-only interface. | Bed + fluid conservation, sequential visibility and rollback/failure behavior. |
| A11 | **Primitive-state reconstruction and limits:** who owns re-evaluation of physical state and limiters between groups? | Executor passes current `q`; the physical process reconstructs only necessary state and owns process-specific limiting. | Group B uses post-A temperature/mass/velocity; limiting under overlapping processes. |
| A12 | **Lifetime of compatibility types and wrappers:** which existing modules are truly needed? | Build a call graph and explicit **KEEP / CHANGE / DELETE** inventory before coding, including `equation_partition_type`, dispatch coordinators and old tests. | No dead paths, no duplicate interface conversion or stale API compatibility. |

### Evidence motivating these questions

- In `src/numerics/time_integration_2d.f90`, explicit STAGE source groups are currently evaluated against the same `qp_rk` and their rates added into `EX_terms`; the spatial operator is evaluated separately. **Target contract:** spatial discretization must be first, in explicit GROUP 0, using the original explicit-stage state; only subsequent explicit groups see sequentially updated states. The implicit stage uses `solve_imex_stage` and returns an effective `IM_terms` separately.
- `src/numerics/nonlinear_solver_2d.f90` currently reconstructs the effective implicit rate as `(q_stage-q_known)/(dt*a_diag)` and uses Newton unknown lists in the same stage coordinator as analytic groups.
- `src/numerics/step_source_executor_2d.f90` already loops sequentially through STEP groups, but contains method-specific paths and a distinct mass-exchange transaction updating the bed.
- `src/physics/source_configuration_2d.f90` and `src/numerics/equation_metadata_2d.f90` still carry historical equation-partition metadata. Its continued existence must be justified from actual numerical requirements, not historical usage.

### Proposed decision order

First resolve **A2 and A5** from the RK mathematics: these control whether a fully unified STAGE sweep is meaningful without changing the integration scheme. Then resolve **A1, A3, A4, A6, A7, A8 and A9** as the smallest implementation consistent with that mathematics. Finally resolve the non-fluid-state transaction **A10**, physical-state ownership **A11**, and keep/change/delete inventory **A12**.

**Naming convention agreed:** numerical model interfaces use `eval_<model>_explicit`, `eval_<model>_semi_implicit`, `eval_<model>_implicit`, and `eval_<model>_analytic`. This is a naming contract; it does not by itself certify that older analytic implementations already return complete conservative rates.

**Acceptance gate:** do not start a broad source architecture rewrite until the selected options above have been recorded, the RK formula has been checked, and the minimum set of needed types and executors has been drawn.

## 10. Confirmed RK integration decisions (2026-10-31)

The following decisions **supersede the corresponding open questions A1, A2, A5 and A8 in section 9**. They are design requirements, not claims that the current code already implements them.

1. **Explicit integration interval:** between sequential EXPLICIT source groups at `STAGE`, advance the intermediate conservative state with the interval associated with that **explicit RK stage**. For `STEP`, use the **complete time step** `dt`. The actual explicit-stage interval is to be taken from the implemented explicit RK discretization; do not substitute the diagonal implicit coefficient automatically.
2. **Distinct RK rate channels:** EXPLICIT source groups (including spatial discretization in EXPLICIT GROUP 0) contribute **only to the explicit RK rates** (`EX_terms`, with spatial terms handled consistently). NEWTON and ANALYTIC groups contribute **only to the implicit RK rates** (`IM_terms`). Do not mix explicit groups into an implicit-group solve or vice versa. Splitting is sequential *within each channel* in configured group order.
3. **NEWTON versus ANALYTIC groups:** both belong to the implicit channel but **must remain different groups**; neither method may be combined with the other in the same group. Execute the ordered implicit-group sequence with each group receiving the conservative state already updated by its predecessors.
4. **Effective rate accumulation:** initialize the complete rate accumulator for each RK channel to zero, then **add each group's effective, limited full-vector rate** to it after evaluating and applying that group. Thus `S_exp_total = sum_g S_exp_g` and `S_imp_total = sum_g S_imp_g`, independently. Never reapply the accumulated increments when assembling stage rates.
5. **Semi-implicit treatment:** semi-implicit contributions are **always evaluated together with their implicit/Newton group**, as in the existing implementation. They are not independent source groups, and do not introduce a third RK rate channel.

**Consequence for architecture:** do not force one mixed EXPLICIT/NEWTON/ANALYTIC STAGE sweep. An implementation may share a minimal group-evaluation primitive, but the explicit and implicit tableau contributions remain separate. A single complete conservative state and full-vector rates are used within each sequence; source groups do not partition equations.

**Remaining numerical checks:** establish the precise value of the explicit RK stage interval for each supported tableau; verify how the separately accumulated `EX_terms` and `IM_terms` enter the RK recurrence; and test the temporal order and admissibility. These checks do not reopen the decisions above.
