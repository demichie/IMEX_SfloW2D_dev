# D18 — IMEX-RK source architecture: target specification

**Status:** target architecture agreed; full Fortran implementation and regression validation pending.  
**Scope:** source contributions during IMEX-RK stages (`STAGE`) and source operators applied after each RK step (`STEP`).

This document specifies **only the current target design**. It is not a chronology of discussions or superseded approaches. Where a physical limiter's detailed algorithm has not been selected, that is stated explicitly; no unagreed algorithm is implied.

## 1. Design and ownership

- Refactor to the **simplest architecture that satisfies the mathematical contract**. Avoid compatibility wrappers, redundant partition types, duplicate state conversion, and extra orchestration layers without independent purpose. Determine KEEP / CHANGE / DELETE from the actual call graph.
- Physical processes are organized by their physics, not by exclusive ownership of conservative equations. Distinct processes and groups may contribute to the **same** components.
- A physical source provides its complete conservative rate vector, length `n_eqns`, with zeros only on components it does not affect. A group similarly returns a complete vector; group metadata does not restrict which equations may receive a rate.
- A source configuration identifies physical process name(s), available numerical METHOD, CADENCE (`STAGE` or `STEP`), and GROUP identity/order. Unsupported methods or physically incompatible configurations are rejected at configuration/preflight, and unsupported model entry points fail explicitly rather than silently returning zero.
- Keep four model-specific interfaces, consistently named `eval_<model>_explicit`, `eval_<model>_semi_implicit`, `eval_<model>_implicit`, and `eval_<model>_analytic`. Numerical responsibilities and real/complex state types remain appropriate to the interface. Do not invent a universal wrapper state merely to preserve legacy calls.
- The RK integrator owns the Butcher tableau, stage states, and final weighted assembly. Source physics owns constitutive calculations and physical limiting; generic Newton owns only nonlinear/linear numerical solution and line search.

## 2. Methods, grouping, and time intervals

| METHOD / CADENCE | Execution | State propagation | Relevant interval |
| --- | --- | --- | --- |
| EXPLICIT / STAGE | All explicit processes, including the spatial operator, evaluated at one common RK stage state | **No advancement between groups** | No splitting interval; explicit rheology limiter horizon is `dt_global` |
| NEWTON / STAGE | Ordered implicit groups, each acting on current implicit working state | Sequential between implicit groups | `tau_imp=a_ii*dt_global` |
| ANALYTIC / STAGE | Ordered exact/analytic implicit group | Sequential with implicit groups | `tau_imp=a_ii*dt_global` |
| EXPLICIT / STEP | Ordered post-RK groups | Advance after each group | `dt_global` |
| Other supported STEP methods | Apply their configured method as post-RK operators | Sequential between groups | `dt_global` |

- The explicit stage channel contributes **only** to the explicit RK weights. Newton and analytic STAGE groups contribute **only** to the implicit RK weights. Semi-implicit rates are internal to their associated Newton group; there is no third RK rate channel and no standalone semi-implicit group.
- An ANALYTIC group contains **exactly one physical process**. It is never merged with EXPLICIT or NEWTON. NEWTON and ANALYTIC groups remain distinct even where they have the same cadence or affect overlapping equations.
- For each implicit RK stage, skip NEWTON/STAGE and ANALYTIC/STAGE when `a_ii=0`, set the full implicit rate to zero, and avoid division by zero. Explicit rates must nevertheless be evaluated whenever later explicit tableau coefficients or final weights require them, even if the implicit diagonal coefficient or an explicit tableau entry is zero.
- The argument name `dt` in source-family interfaces denotes the **interval relevant to that call**; the caller supplies `dt_global` at STEP, `a_ii*dt_global` for applicable implicit/semi-implicit STAGE operations, and `dt_global` as the explicit rheology limiter **horizon**. The latter does not authorize a stage-local state update.
- The currently supported tableau stage counts are 2, 3, and 4. The existing RK tableaux and weights are **not changed** by this architecture.

## 3. EXPLICIT/STAGE: simultaneous evaluation

For each RK stage `i`, let `q_i` be the conservative stage state constructed by the IMEX integrator.

1. Evaluate the spatial operator `F(q_i)` and **every** EXPLICIT/STAGE process `S_p(q_i)` from exactly this same, **unmodified** stage state.
2. Explicit GROUP 0 includes the spatial operator and may include source processes. Other explicit GROUP values remain useful for source organization and configuration, **not** for propagation of intermediate states at STAGE.
3. Preserve the **separate process-rate vectors**, across physical families and group boundaries, until simultaneous physical joint limiting is complete (section 6). Apply any justified process/family-specific constitutive constraints as specified there.
4. Aggregate the physically adjusted contributions exactly once:
   ```text
   R_exp,i = F(q_i) + sum_p S_p,limited(q_i)
   ```
   Assemble this rate with the existing **explicit** RK tableau. No `q_i += dt*S` update is performed between explicit STAGE groups, and no local split update is added a second time by the integrator.
5. The spatial operator is evaluated on the common stage state, not reevaluated in each group. Its discretized intercell fluxes retain conservation; do **not** locally rescale a flux-divergence vector as though it were an independently disposable source contribution.
6. Different explicit STAGE groups are **simultaneous for evaluation and joint inventory limiting**, even if their configuration GROUP IDs differ.

**Numerical rationale:** a sequential stage-local Euler update (using `dt_global`, `a_ii*dt_global`, or another fixed nonzero multiple) modifies the stage RHS by an order-`dt` term and can reduce a nominal second-/third-order RK scheme to first order. The target therefore uses a **common stage state** and an unchanged explicit RK rate.

## 4. STEP: genuinely sequential source operators

After completing the RK step, start from `q^(0)=q_RK,out` and apply configured STEP groups in their specified order.

```text
q_work = q_RK,out
for each STEP group g in order:
    evaluate each member process on q_work (same state for all members)
    jointly limit competing contributions within this group
    rate_g = sum of adjusted full conservative rates
    q_work = q_work + dt_global * rate_g
    commit associated bed/substrate/deposit inventory changes atomically
q_out = q_work
rate_total = (q_out - q_RK,out) / dt_global
```

- All sources **inside one group** see the same group-entry state and do not update it individually. A subsequent group reconstructs its physical state from the **already updated** conserved state and inventories.
- Competing processes within one STEP group must be jointly limited against the available inventory. Distinct STEP groups are **not** one simultaneous limiter domain: each sees the remaining inventory after the preceding committed update.
- STEP processes that affect bed elevation `B`, erodible thickness, or deposit inventories must update those quantities consistently with their conservative fluid increments; maintain a single coherent transaction per group, including failure handling. These substrate quantities are **not additional RK-integrated conservative variables**.
- STEP rates are not inserted into `EX_terms` or `IM_terms` and are not reapplied by RK final assembly.
- Ordered first-order STEP composition is **not claimed** to preserve the order of the full RK+STEP method for dynamic sources. Temporal convergence of combined problems must be measured.

## 5. Implicit STAGE sequence and reduced Newton

- At an implicit RK stage with `tau_imp=a_ii*dt_global>0`, build the proper stage-known state from the IMEX tableau. Apply configured NEWTON and ANALYTIC groups **in order**; each receives the full conservative state resulting from its predecessor. Rates may overlap in conservative components; no global equation-exclusive partition is imposed.
- NEWTON iterates on **conservative `q`**, never primitive `qp`. Each group has one local list `newton_equations(:)`, used for **both** residual rows and Jacobian columns; do not maintain separate unknown and equation lists merely to describe the same reduced system.
- For a rheology-only Newton group, the reduced list comprises conservative momenta `[idx_uEqn,idx_vEqn]`, **not** the primitive velocity indices. Complex-step differentiation builds the required reduced Jacobian and the solver handles convergence/line search. Other conservative quantities can affect rate calculations while remaining fixed because the present group's sources do not evolve them.
- A group's *physical source* can return nonzero rates outside its Newton list. Where outside components E depend only on solved implicit components I, and the I residual is independent of the **updated** E, solve I first and then reconstruct E **exactly** from the complete source evaluated at the solved state. This is the triangular block case
  ```text
  J_source = [ 0  B ]
             [ 0  D ]
  J_resid  = [ I -tau*B ]
             [ 0 I-tau*D ]
  ```
  with ordering `(E,I)`. No extra Newton unknowns or generic Schur solver are needed for this case. The actual dependency conditions must be checked, not inferred from which columns a partial Jacobian routine happens to evaluate.
- If the source on E depends on unknown updated E, or the I residual depends on updated E, the simple reconstruction is **not justified**. Use a physically justified analytic elimination or enlarge the Newton system; never silently omit nonzero coupled residuals.
- **Current source dependency audit:** the existing implicit dispatcher writes rheology rates only in momenta `(2:3)` and thermal cooling only in energy `(4)`, with corresponding reduced Newton lists. Consequently, rates outside the selected components are currently zero for these groups. Rheology rates can nevertheless depend on other conserved quantities (mass, composition, energy) held fixed through that group. Darcy pore pressure is currently ANALYTIC, not NEWTON; mass-exchange models are not currently dispatched in the Newton implicit-source path. The generalized reconstruction is a **future-required behavior when a real nonzero outside-list rate appears**, not a reason to introduce unused block machinery now.
- ANALYTIC solves the specified physical evolution over `tau_imp` (including relevant physical admissibility) and returns the complete effective conservative rate. Do not apply an extra Newton projection to an analytic result.
- The implicit rate for the RK tableau represents the full effective change from the stage-known state over the implicit diagonal interval, **exactly once**. Validate effective-rate bookkeeping when multiple implicit groups are applied sequentially.

## 6. Physical limiters: ownership and simultaneous contributions

**Rule:** the limiter of a simultaneously applied set of processes receives their **separate full-vector rate contributions**, the common input state, and relevant inventory/physical context. Preserve process identity through limiting; **sum only the adjusted rates** to form the final conservative contribution. Rates from different physics families and distinct explicit STAGE GROUP IDs can compete for the same quantity.

- **Constitutive constraints** specific to one law belong in its physical model/family and may be enforced during model evaluation; do not duplicate them at the generic RK/Newton level.
- **Shared-inventory limits** (mass, phase inventory, erodible substrate, etc.) must consider **all simultaneous demands**, not only each process in isolation. Distinguish production from consumption and retain process provenance; limiting a net sum alone can hide an inadmissible withdrawal.
- For EXPLICIT/STAGE, cross-family and cross-group joint limiting is coordinated **after collecting all common-state source contributions** and before returning the rate to RK. A limiter cannot mutate the common input stage state. The spatial contribution may need to be considered by a **separate conservation-preserving admissibility strategy**; do not arbitrarily clip local spatial flux-divergence rates.
- For EXPLICIT/STEP, joint limiting occurs among members of the **current** group, followed by the group's state/inventory update. The next group performs a fresh physical-state and inventory evaluation.
- The assembler coordinates the limiter call and keeps individual rates available; the **physical domain defining the constraint** owns the rule for how individual contributions are adjusted. The generic RK or Newton solver must not impose a universal scalar clipping of all conservation laws.
- Generic dry-state handling, roundoff repair, conversion constraints, and final RK admissibility checks remain distinct from physical source limiters. **The weighted RK combination may still need its own validation:** limiting instantaneous stage rates does not prove final-state admissibility.

**Detailed physical limiter algorithms not yet fully specified:** which joint constraints require cross-family adjustment, how competing withdrawals are apportioned, and exactly how source admissibility interacts with the conservative spatial flux discretization. These are the remaining limiter-design tasks; no generic algorithm is assumed.

### 6.1 Rheology momentum reversal

- The common rheology anti-reversal limiter is owned by `src/physics/rheology_laws/rheology_2d.f90`, at the **family aggregation** level, after summing rheology-model contributions on the same state, not redundantly inside every rheological law.
- It applies to **EXPLICIT** and **SEMI_IMPLICIT** rheology only. The explicit limiter uses the **full global `dt` as its limiter horizon**, even when RK requires a rate at a stage with zero explicit tableau coefficient. This does not advance the STAGE state. The semi-implicit limiter uses `tau_imp=a_ii*dt_global` for its predictor **before Newton**.
- For the incoming momentum vector `m0` and candidate isolated rheology-updated momentum `m1`, prevent **actual reversal**, `dot(m1,m0)<0`, by removing only the reversing component along `m0`. Preserve transverse momentum and impose **no upper bound on momentum magnitude** solely for this rule. One possible expression is `m_after=m1-min(0,dot(m1,m0)/dot(m0,m0))*m0` when `dot(m0,m0)>0`; the effective rheology rate reflects that limited increment.
- When both incoming conservative momentum components are **exactly zero**, return zero explicit or semi-implicit rheology momentum rate immediately; this is separate from dry-state tolerances. Do not generalize the shortcut to unrelated sources.
- Do **not** apply this extra reversal projection after a fully implicit NEWTON solve, or to ANALYTIC rheology. A post-Newton projection could invalidate the converged residual; analytic flows must satisfy their own model physics.
- A rheology-specific limiter is **not** a blanket limiter on the assembled sum of rheology plus gravity, injection, spatial fluxes, etc. Any shared total-update constraint would need its own justified physical rule.

### 6.2 Other physical limiter domains

- Mass exchange owns shared carrier-withdrawal/packing and erosion/deposition inventory logic, including its STEP substrate transaction.
- Darcy pore-pressure physics owns its gas withdrawal and pore-pressure admissibility constraints, independently of mass-exchange source ownership.
- State conversion, generic admissibility verification and Newton line search are numerical safeguards, not substitutes for physically identified shared-inventory limiters.

## 7. Single source-executor contract

The design uses **one shared source-group execution mechanism**, with the orchestration appropriate to cadence/method. Conceptually:

```fortran
CALL execute_source_groups(groups, q_in, dt, q_out, rate, ...)
```

- `q_in`, `INTENT(IN)`: immutable conservative cell state.
- `q_out`, `INTENT(OUT)`: resulting conservative cell state where an operator actually advances state (STEP or an ordered implicit group operation); it must **not** replace the common STAGE state between explicit groups.
- `rate`, `INTENT(OUT)`: full conservative effective rate, already limited and counted once. For EXPLICIT/STAGE it is the sum of simultaneously evaluated rates; for STEP with positive `dt` it matches `(q_out-q_in)/dt`.
- Maintain local cell workspaces; reconstruct primitive/physical state when the input conservative state actually changes. No domain-wide copy or bespoke wrapper type is needed just to satisfy the interface.
- The same evaluation/dispatch primitives can be reused between cadences, but the **explicit STAGE common-state path and sequential STEP/implicit paths have different advancement semantics**. Do not force one mixed explicit/implicit stage sweep that violates distinct tableau channels.
- Bed and inventory context are provided only where physically required, notably mass exchange at STEP. Errors must not leave a partially committed coupled fluid/bed transaction.
- For EXPLICIT/STAGE, `q_out` is **not** an intermediate state handed to another group. Where the interface demands a candidate state, it may be derived for diagnostics, but the RK update remains owned exclusively by the tableau.

## 8. Verification and implementation requirements

**Required tests:**

1. EXPLICIT/STAGE: spatial operator and all source groups evaluate at the identical stage state; no stage-local inter-group update; RK weights and order (pure explicit 2-/3-/4-stage tests) remain unchanged.
2. EXPLICIT/STAGE: two physical processes in different groups removing the same inventory; verify a **joint** limiter receives both separate contributions; verify no false admissibility from individually limited rates.
3. EXPLICIT/STEP: sequential groups modifying overlapping equations; later groups see prior limited fluid and substrate state; within-group processes see the same group-entry state.
4. Reduced rheology/thermal Newton: correct residuals, frozen non-evolving conserved components, correct complex-step Jacobian, and overlapping sequential implicit groups.
5. NEWTON/ANALYTIC combinations in both permissible orders, full-vector effective rates, skipped implicit stages with `a_ii=0`, explicit rate evaluation as RK requires, no double application.
6. Rheology explicit/semi-implicit anti-reversal including transverse motion, zero momentum, absence of post-Newton projection, and interactions with other momentum sources.
7. Conservation and admissibility for fluid, gas, solids, bed erosion/deposition, Darcy with and without gas loss; stage and final RK admissibility.
8. Temporal convergence of the **combined** IMEX-RK + STEP process, especially when STEP contains dynamic operators; do not claim RK order transfers automatically.
9. Build/CI regressions and performance comparisons against relevant v2 cases; report what was actually compiled/run.

**Implementation state:** this is a target specification, not a certification of completed refactoring. The code already contains parts of the source registry, reduced Newton and STEP transaction, but the unified executor, joint limiter wiring, and explicit rheology limiter interface still require implementation/verification. The remaining semi-implicit limiter in the generic nonlinear solver must be moved to the rheology-family owner; the post-Newton reversal projection has already been removed. No build or regression test is claimed by this documentation rewrite.
