# D18 — Lie-splitting source architecture: target specification

**Status:** agreed target architecture; Fortran refactoring and verification remain pending.  
**Scope:** first-order Lie splitting of the finite-volume spatial operator and physical source groups. This document specifies the **current design only**, not the history of previous designs.

## 1. Time integrator: Lie only

The global time integrator is **first-order Lie operator splitting**. The IMEX Runge–Kutta method, its tableau, its RK-stage arrays, and the configuration distinction `STAGE`/`STEP` are **not part of the target architecture**. There is **no Strang option** in this target.

For a step of duration `dt_global>0`, start from the complete state at `t_n` and execute each configured operator in its specified order, passing the **updated** state to its successor:

```text
state_0 = state_n                       # immutable pre-spatial/pre-source state
spatial_rate = evaluate_spatial_operator(state_0, t_n)
for each additional configured EXPLICIT source p in GROUP 0:
    source_rate_p = evaluate_explicit_source(p, state_0, t_n)
jointly_limit_group0_source_rates_and_check_combined_admissibility()
state_work = apply_group0_once(state_0, spatial_rate, source_rates, dt_global)
validate_group0_updated_state()

for group g = 1, 2, ... in configured order:
    state_before = state_work
    state_work = advance_configured_group(g, state_before, t_n, dt_global)
    validate_updated_components_and_coupled_inventories(g, state_work)
state_(n+1) = state_work
```

**GROUP 0 is mandatory and EXPLICIT**, and always contains the finite-volume spatial operator. It is executed first. Any additional physical source configured in GROUP 0 must also be EXPLICIT and is **evaluated on exactly the same pre-spatial state** as the spatial operator; no evaluation of a GROUP-0 source may use a partly advanced spatial state. Spatial and source contributions are assembled into **one** forward-Euler update, subject to physically justified joint limits and conservative spatial-flux treatment. After completing GROUP 0, the integrator loops over groups **1 onward** in configured order, with each seeing the validated output state from its predecessor. **After every completed group**, verify admissibility of the conservative components and physical inventories it has changed; do not project a converged implicit solution or alter intercell flux conservation with a generic post-update clip.

The complete state includes conservative fluid unknowns and, for processes that require them, topography, deposit, erodible material, and related physical inventories. A group call may be implemented with local input/output states, preserving an immutable input for evaluation. Do not copy the entire simulation domain unnecessarily.

**Illustrative configuration (not prescribed):**

| Order | Physical operator/group | Local integrator |
| --- | --- | --- |
| 0 | Finite-volume spatial operator + injection processes | Forward Euler explicit |
| 1 | Rheology/friction | Backward Euler implicit (reduced Newton) |
| 2 | Mass exchange, including deposition/erosion and substrate transaction | Forward Euler explicit with physical joint limits |
| 3 | Darcy pore-pressure process | Analytic evolution |

The spatial operator must reside **only in GROUP 0**; GROUP-0 source processes must be EXPLICIT. All subsequent groups select any model-supported configured numerical method (EXPLICIT, NEWTON, ANALYTIC), without a fixed method per physical family. Group 0 is never repeated in the later group loop.

The example assigns methods only to illustrate one valid configuration. **No source family is permanently assigned EXPLICIT, NEWTON, or ANALYTIC by this architecture:** rheology, for example, can be explicit, implicit/Newton, or analytic where the selected model supports it. The spatial finite-volume operator itself remains explicit. Group order and method remain configured and validated; placing friction immediately after spatial acceleration is an optional configuration choice. Lie order affects the result.

Each group advances its own physical subsystem over the **same full interval** `dt_global`, even though the operators are called successively. Their calls are not consecutive independent physical-time intervals; one overall Lie step advances time from `t_n` to `t_n+dt_global`. **Time-dependent source laws**, including injection schedules, must receive the appropriate physical-time argument and be treated consistently with the chosen first-order discretization (e.g., left endpoint for forward Euler). Never increment the physical clock by `dt_global` after each group.

## 2. Group contract and minimal responsibilities

- A group consists of one or more physical processes configured with a compatible numerical treatment and application transaction. Each source family retains its supported EXPLICIT, NEWTON/implicit, and ANALYTIC choices, rather than having one hard-coded treatment. The semi-implicit interface is retained as an auxiliary part of its associated NEWTON solve, not a standalone global Lie method. Different groups may modify **overlapping conservative equations**. There is no global exclusive equation partition.
- Every physical process returns a complete conservative rate vector `rate(1:n_eqns)`; unaffected components are zero. **Keep individual process contributions distinct until their joint limiter runs**.
- All processes **within one explicit group**, including spatial and source contributions in GROUP 0, evaluate the same immutable group-entry fluid and inventory state. Assemble and jointly limit physical source demands, then apply one limited group update. A subsequent group reconstructs its primitive and physical variables from the updated state. GROUP 0 is special only because it contains the domain-coupled spatial operator and is executed before the loop over groups 1 onward.
- **Two distinct limiter responsibilities:** before an update, physical constitutive/joint-inventory limiters act on separate proposed process rates or transfers; **after each group**, a state-admissibility check covers the components and coupled inventories updated by that group. If the candidate is inadmissible, use a physically/conservatively justified adjustment or reject/retry the step with a smaller `dt`; do not silently clip mass, momentum, or Newton results after the fact.
- Conceptually, `advance_group(group, state_in, time, dt, state_out)` advances the group in one cell or an appropriate domain/spatial transaction. The actual Fortran signature should be minimal. `state_in` is never mutated, `state_out` includes all affected fluid/bed inventories. Separate `q_in`/`q_out` are preferred; additional bed/inventory arguments are introduced only where physically required.
- The **spatial operator is domain-coupled**, unlike cell-local source groups. Evaluate face fluxes consistently across the domain, update conservative cells with the FV discretization and handle boundaries before subsequent cell-local groups. Group 0 may contain injection terms evaluated on the same pre-update cell states as the spatial fluxes. Do not arbitrarily clip a local flux-divergence vector if doing so would destroy intercell conservation.
- A step must not be partially committed after an error in a coupled fluid/bed transaction; implement validation and coherent commit/rollback or fail-before-commit handling as needed.
- Preserve the existing four physics-model entry points named `eval_<model>_explicit`, `eval_<model>_semi_implicit`, `eval_<model>_implicit`, `eval_<model>_analytic` where they are an established interface contract. Unsupported methods fail explicitly and source configuration rejects invalid assignments. A semi-implicit contribution, where used, belongs inside the associated NEWTON group; it is not a global time channel or an independent RK stage.
- Avoid extra wrapper state types, duplicate conversions, equation-partition compatibility layers, or unused solvers. Derive a KEEP / CHANGE / DELETE list from the actual source call graph before implementation.

## 3. EXPLICIT group update (any supported physical family)

For an explicit group with process rates `S_p`:

```text
for every p in group:
    rate_p = physical_rate_p(state_in, t_n)
adjust all competing rate_p with physical joint limiters
rate_g = sum_p rate_p_adjusted
q_out = q_in + dt_global * rate_g
```

Any associated inventory or topographic update must use the **same limited physical transfers**. Source input states must not be updated between member processes. GROUP 0 additionally includes conservative spatial face-flux divergence, assembled consistently with its simultaneous injection rates evaluated on the **pre-spatial** state. Its spatial evaluation can be performed first operationally, but must not commit a spatial state change before the GROUP-0 source rates have been evaluated. A single combined update follows. The explicit state update must respect the spatial CFL and physical admissibility constraints; time-step selection is not solely controlled by the implicit friction method.

A group configured as EXPLICIT uses the first-order forward Euler method regardless of whether the source is stiff; **stiffness imposes stability/positivity restrictions** on any source left explicit. A source requiring unconditional or large-step stability must be assigned a suitable implicit/analytic treatment or reduce the step.

## 4. NEWTON group: generic implicit treatment

- **Any group configured as NEWTON** is advanced over the full `dt_global` by backward Euler solved with Newton. This is a method-specific contract, **not a requirement that rheology always use Newton**. Rheology may instead use EXPLICIT or ANALYTIC when the selected law supports them; no global IMEX tableau remains.
- Newton works on **conservative** variables. For a rheology-only group **configured as NEWTON**, the reduced unknown/residual index list is `[idx_uEqn,idx_vEqn]`, not primitive velocity indices. One list `newton_equations(:)` is used for both residual rows and Jacobian columns; groups do not claim exclusive equation ownership.
- The group-entry state is fixed as the backward-Euler reference. The full physical source is evaluated at the implicit outgoing state, and the solved residual must be respected; the Newton solver manages numerical convergence and line search, not rheology-specific physical clipping.
- Current supported NEWTON sources affect only the selected momenta for rheology or energy for cooling, so all other components stay at group-entry values. For a future source that modifies conservative components **outside** the reduced Newton list, reconstruct them after solving only when the exact triangular dependence conditions hold: outside-component rates must not depend on updated outside unknowns, and Newton-equation residuals must not depend on their updated values. Otherwise enlarge the Newton system or use an exact justified elimination. Do not discard nonzero residual components.
- The existing semi-implicit reological predictor, if still used by a model, is an **internal part of the implicit group** over `dt_global`, and its physical limiter is applied before the nonlinear solve. Do not apply an extra rheology momentum projection after convergence because it invalidates the solved residual.

## 5. ANALYTIC group: generic analytic treatment

- An ANALYTIC group contains **exactly one physical process**, integrated by its analytic physical evolution operator over `dt_global`. It updates the **complete relevant conservative state**, without an equation-ID filter or a second Euler application of its increment.
- The currently supported Darcy pore-pressure model is **ANALYTIC** in the illustrative configuration; its model-specific method capabilities do not constrain those of other source families. Its model owns pore-pressure admissibility and any interstitial gas withdrawal. The process can affect total mass, momentum, energy, gas and pore-pressure conservative variables as dictated by the physical equations.
- Analytic evolution may require input-state assumptions or frozen coefficients; document and test their validity. The exactness claim applies to the stated subproblem, not automatically to the full split coupled evolution.

## 6. Physical limiter ownership and joint rates

- Physical constitutive constraints belong in the corresponding model or family. Limiting shared finite inventories belongs to the physical domain that owns the inventory, **coordinated across all simultaneous processes in one group**.
- A joint limiter receives (a) the immutable incoming state and applicable inventories, (b) **each individual physical process's complete rate**, and (c) the time interval and any additional physical context. It can adjust the individual process contributions consistently; aggregate them **only after limiting**. Net rate alone is insufficient when independent withdrawals compete.
- Between groups, there is **no simultaneous cross-group limiter**: a group first computes joint physical limitations, performs one admissible update and any coupled inventory transaction, and then **checks all changed conservative components and associated inventories**. Only after acceptance does the next group use the new state/inventories. Inadmissibility requires a justified conservative correction at the proper ownership level or a failed/retried step, never arbitrary post-update clipping.
- The FV flux contribution is subject to its own conservation-preserving discretization/admissibility policy. A generic cellwise source limiter must not rescale face-flux divergence independently of neighboring cells.
- State-conversion roundoff repair, dry-state consistency, final-state admissibility checks and Newton line search are numerical safeguards, not replacements for physical inventory limits.
- **Detailed algorithms to finalize before code integration:** inventory ownership across processes in one group, allocation/priorities for concurrent withdrawals, and consistency of source transfers with bed/deposit changes. Do not introduce an unphysical universal per-component clipping factor.

### 6.1 Rheology anti-reversal

- The rheology-family coordinator `src/physics/rheology_laws/rheology_2d.f90` owns anti-reversal logic, **not** the general Newton solver, and individual constitutive routines need not duplicate it.
- For an explicit rheology update (if configured), limit the combined rheological momentum rate over **the actual full group interval `dt_global`**. Likewise limit an isolated semi-implicit predictor before Newton if it is used. Do not apply the rate-level reversal limiter to an implicit converged solution or to an analytic rheology flow.
- The anti-reversal rule prevents only actual reversal: for momentum `m0` and candidate `m1`, forbid `dot(m1,m0)<0` by removing the reversing component parallel to `m0`, preserving transverse momentum and not imposing a general cap on magnitude. For `dot(m0,m0)>0`, one expression is `m_after=m1-min(0,dot(m1,m0)/dot(m0,m0))*m0`.
- If both incoming conservative momentum components are **exactly zero**, the isolated reological momentum rate is zero; do not invoke the constitutive friction model unnecessarily. This is separate from near-dry-state handling.

### 6.2 Mass exchange and Darcy

- Mass exchange physically coordinates erosion, deposition, entrainment, carrier withdrawal, packing limits and the corresponding fluid–substrate transfers; keep competing processes affecting one substrate/inventory in a **joint group transaction** unless a genuinely conservative alternative split is demonstrated.
- Darcy owns interstitial-gas depletion limitations and any conservative change implied by its pore-pressure evolution; it is not silently folded into mass exchange.

## 7. Configurable group dispatch and implementation cleanup

All existing source-group membership and ordering rules remain, but without the obsolete STAGE/STEP cadence. The group METHOD is independently configured and validated against the selected model's capabilities. The source-family dispatcher must invoke the selected numerical interface, not infer it from names such as rheology or mass exchange. No group is required to own exclusive conservative equations.

## 8. Data flow and implementation cleanup

Target global loop:

```text
while t < t_final:
    dt = choose_admissible_dt(current_state, geometry, CFL, source_limits)
    dt = min(dt, t_final - t)

    state_before = current_state
    spatial_rate = spatial_operator(state_before, t)
    source_rates_0 = evaluate_all_group0_explicit_sources(state_before, t)
    limited_rates_0 = jointly_limit_group0_sources(source_rates_0, state_before, dt)
    current_state = apply_spatial_and_sources_once(
        state_before, spatial_rate, limited_rates_0, dt)
    validate_group_updated_components(group=0, current_state)

    for group in configured_groups_with_ID_at_least_1:
        state_before = current_state
        current_state = advance_group_by_configured_method(
            group, state_before, t, dt)
        validate_group_updated_components(group, current_state)

    commit_validated_fluid_and_bed_state()
    t = t + dt
```

The group-0 spatial rate may be computed first, but the **input** to every group-0 source evaluation is the same original pre-spatial state. No separate spatial update is committed before group-0 source assembly; do not count any contribution twice.

Remove RK-stage bookkeeping used only by IMEX integration (`q_rk`, `qp_rk`, `EX_terms`, `IM_terms`, per-stage `spatial_terms`, tableau matrices and weights, RK stage loops, and stage-specific implicit diagonal handling). Preserve any **workspace genuinely required by the FV spatial discretization, reconstruction, diagnostics, or source evaluation** even if its old name included an RK term, until its remaining uses are audited.

Eliminate source `CADENCE=STAGE/STEP` from the target configuration and merge the currently separate STAGE/STEP orchestrators into an ordered Lie group executor. Preserve a distinct domain-wide spatial advance primitive and atomic STEP-like substrate transactions as required by the physical methods. Use a **single full-step time variable**, with local interval parameter `dt`; there are no RK stage weights or rates to accumulate into final RK assembly.

## 9. Verification and acceptance

1. Pure spatial FV+injection forward Euler: conservation, hydrostatic well-balancing where applicable, CFL and positivity, time-dependent source scheduling, dry/wet transitions and boundary-source handling.
2. **Mandatory GROUP 0**: spatial operator evaluated first and every configured group-0 explicit source evaluated from exactly the same original pre-spatial state; one combined update only. Groups 1 onward see their predecessor's validated updated fluid/bed/inventory state; the integrator checks updated components **after every group**, with no double application and deterministic order.
3. Shared finite-inventory limiter: retain distinct individual process rates until joint limiting, including competing consumption and production, conservation across fluid/substrate, and proper mass-exchange bed transaction.
4. Configurable family methods: test rheology EXPLICIT, NEWTON, and ANALYTIC when supported, and reject invalid method/model combinations. For NEWTON rheology check reduced residual/Jacobian, convergence, zero-momentum shortcut, and no post-Newton momentum projection; test stiff linear friction against `u' = a-ku`.
5. Analytic Darcy: correct complete conservative update, pore-pressure-only versus interstitial gas-loss behavior, and physical admissibility.
6. Multiple operators on the same conservative components: Lie-order sensitivity, temporal convergence **to first order** on coupled noncommuting test problems, and relevant legacy v2 physical regressions. Do not mistake first-order splitting for a well-balanced or uniformly stiff-accurate method.
7. Whole-step state/geometry consistency and failure behavior, including updated topography before subsequent spatial flux evaluation; reject invalid method configurations and unsupported entries.
8. Compare wall-clock performance and memory with the prior IMEX-RK version; run relevant CI, and report exact tests actually executed.

**Known limitations:** Lie splitting generally has first-order global temporal accuracy even when component operators are exact. Configuring friction implicitly avoids an explicit friction stability restriction but does **not** guarantee accurate gravity–friction equilibrium for arbitrarily large `dt`; test on stiff balanced regimes. No numerical convergence or production performance claim is made by this architecture document.

**Implementation status:** D18 is the target specification, not evidence that the Fortran refactoring is complete. No source edits or tests are claimed as part of this documentation revision.
