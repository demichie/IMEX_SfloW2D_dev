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
- A Newton solver may choose a **smaller set of unknowns** to compute a group's full conservative rate. Those solver unknowns are a numerical detail, **not** group membership or a restriction on source-rate components.
- Avoid representing a general source group through an exclusive equation membership map. If a Newton unknown map is useful, keep it internal to the relevant solve.

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
