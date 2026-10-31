# D18 — Lie-splitting source architecture: target specification

**Status:** agreed target architecture; Fortran refactoring and verification remain pending.  
**Scope:** first-order Lie splitting of the finite-volume spatial operator and physical source groups. This document specifies the **current design only**, not the history of previous designs.

## 1. Time integrator: Lie only

The global time integrator is **first-order Lie operator splitting**. The IMEX Runge–Kutta method, its tableau, its RK-stage arrays, and the configuration distinction `STAGE`/`STEP` are **not part of the target architecture**. There is **no Strang option** in this target.

For a step of duration `dt_global>0`, start from the complete state at `t_n` and execute each configured operator in its specified order, passing the **updated** state to its successor:

```text
state_work = state_n
for group g in configured Lie order:
    state_work = advance_group(g, state_work, t_n, dt_global)
state_(n+1) = state_work
```

The complete state includes conservative fluid unknowns and, for processes that require them, topography, deposit, erodible material, and related physical inventories. A group call may be implemented with local input/output states, preserving an immutable input for evaluation. Do not copy the entire simulation domain unnecessarily.

**Reference configuration:**

| Order | Physical operator/group | Local integrator |
| --- | --- | --- |
| 0 | Finite-volume spatial operator + injection processes | Forward Euler explicit |
| 1 | Rheology/friction | Backward Euler implicit (reduced Newton) |
| 2 | Mass exchange, including deposition/erosion and substrate transaction | Forward Euler explicit with physical joint limits |
| 3 | Darcy pore-pressure process | Analytic evolution |

This is an **example ordering**, not a hard-coded universal source order. The key recommended adjacency is spatial acceleration followed immediately by implicit rheology. Source-group order remains configurable where physics and method capabilities allow. With Lie, swapping groups may change results by an order-`dt` global splitting error.

Each group advances its own physical subsystem over the **same full interval** `dt_global`, even though the operators are called successively. Their calls are not consecutive independent physical-time intervals; one overall Lie step advances time from `t_n` to `t_n+dt_global`. **Time-dependent source laws**, including injection schedules, must receive the appropriate physical-time argument and be treated consistently with the chosen first-order discretization (e.g., left endpoint for forward Euler). Never increment the physical clock by `dt_global` after each group.

## 2. Group contract and minimal responsibilities

- A group consists of one or more physical processes sharing a numerical integration method and an application transaction. Different groups may modify **overlapping conservative equations**. There is no global exclusive equation partition.
- Every physical process returns a complete conservative rate vector `rate(1:n_eqns)`; unaffected components are zero. **Keep individual process contributions distinct until their joint limiter runs**.
- All processes **within one explicit group** evaluate the same immutable group-entry fluid and inventory state. Assemble and jointly limit them, then apply one limited group update. A subsequent group reconstructs its primitive and physical variables from the updated state.
- Conceptually, `advance_group(group, state_in, time, dt, state_out)` advances the group in one cell or an appropriate domain/spatial transaction. The actual Fortran signature should be minimal. `state_in` is never mutated, `state_out` includes all affected fluid/bed inventories. Separate `q_in`/`q_out` are preferred; additional bed/inventory arguments are introduced only where physically required.
- The **spatial operator is domain-coupled**, unlike cell-local source groups. Evaluate face fluxes consistently across the domain, update conservative cells with the FV discretization and handle boundaries before subsequent cell-local groups. Group 0 may contain injection terms evaluated on the same pre-update cell states as the spatial fluxes. Do not arbitrarily clip a local flux-divergence vector if doing so would destroy intercell conservation.
- A step must not be partially committed after an error in a coupled fluid/bed transaction; implement validation and coherent commit/rollback or fail-before-commit handling as needed.
- Preserve the existing four physics-model entry points named `eval_<model>_explicit`, `eval_<model>_semi_implicit`, `eval_<model>_implicit`, `eval_<model>_analytic` where they are an established interface contract. Unsupported methods fail explicitly and source configuration rejects invalid assignments. A semi-implicit contribution, where used, belongs inside the associated implicit group; it is not a global time channel or an independent RK stage.
- Avoid extra wrapper state types, duplicate conversions, equation-partition compatibility layers, or unused solvers. Derive a KEEP / CHANGE / DELETE list from the actual source call graph before implementation.

## 3. Explicit group update

For an explicit group with process rates `S_p`:

```text
for every p in group:
    rate_p = physical_rate_p(state_in, t_n)
adjust all competing rate_p with physical joint limiters
rate_g = sum_p rate_p_adjusted
q_out = q_in + dt_global * rate_g
```

Any associated inventory or topographic update must use the **same limited physical transfers**. Source input states must not be updated between member processes. Group 0 additionally includes conservative spatial face-flux divergence, assembled consistently with these simultaneous injection rates. The explicit state update must respect the spatial CFL and physical admissibility constraints; time-step selection is not solely controlled by the implicit friction method.

An explicit group can use the first-order forward Euler method regardless of whether the source is stiff; **stiffness imposes stability/positivity restrictions** on any source left explicit. A source requiring unconditional or large-step stability must be assigned a suitable implicit/analytic treatment or reduce the step.

## 4. Implicit rheology group

- Rheology/friction is advanced for the **full `dt_global`** using backward Euler solved by Newton; no global IMEX tableau remains.
- Newton works on **conservative** variables. For a rheology-only group, the reduced unknown/residual index list is `[idx_uEqn,idx_vEqn]`, not primitive velocity indices. One list `newton_equations(:)` is used for both residual rows and Jacobian columns; groups do not claim exclusive equation ownership.
- The group-entry state is fixed as the backward-Euler reference. The full physical source is evaluated at the implicit outgoing state, and the solved residual must be respected; the Newton solver manages numerical convergence and line search, not rheology-specific physical clipping.
- Current supported NEWTON sources affect only the selected momenta for rheology or energy for cooling, so all other components stay at group-entry values. For a future source that modifies conservative components **outside** the reduced Newton list, reconstruct them after solving only when the exact triangular dependence conditions hold: outside-component rates must not depend on updated outside unknowns, and Newton-equation residuals must not depend on their updated values. Otherwise enlarge the Newton system or use an exact justified elimination. Do not discard nonzero residual components.
- The existing semi-implicit reological predictor, if still used by a model, is an **internal part of the implicit group** over `dt_global`, and its physical limiter is applied before the nonlinear solve. Do not apply an extra rheology momentum projection after convergence because it invalidates the solved residual.

## 5. Analytic group

- An ANALYTIC group contains **exactly one physical process**, integrated by its analytic physical evolution operator over `dt_global`. It updates the **complete relevant conservative state**, without an equation-ID filter or a second Euler application of its increment.
- The Darcy pore-pressure process remains an **analytic** group, not an implicit-Newton group. Its model owns pore-pressure admissibility and any interstitial gas withdrawal. The process can affect total mass, momentum, energy, gas and pore-pressure conservative variables as dictated by the physical equations.
- Analytic evolution may require input-state assumptions or frozen coefficients; document and test their validity. The exactness claim applies to the stated subproblem, not automatically to the full split coupled evolution.

## 6. Physical limiter ownership and joint rates

- Physical constitutive constraints belong in the corresponding model or family. Limiting shared finite inventories belongs to the physical domain that owns the inventory, **coordinated across all simultaneous processes in one group**.
- A joint limiter receives (a) the immutable incoming state and applicable inventories, (b) **each individual physical process's complete rate**, and (c) the time interval and any additional physical context. It can adjust the individual process contributions consistently; aggregate them **only after limiting**. Net rate alone is insufficient when independent withdrawals compete.
- Between groups, there is **no simultaneous cross-group limiter**: a group first performs an admissible update and any coupled inventory transaction, then the next group uses the new state/inventories to evaluate fresh physical constraints.
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

## 7. Data flow and implementation cleanup

Target global loop:

```text
while t < t_final:
    dt_global = choose_admissible_dt(current_state, geometry, CFL, source_limits)
    dt_global = min(dt_global, t_final - t)

    advance_spatial_and_injection_explicit(dt_global)
    advance_rheology_implicit(dt_global)
    advance_mass_exchange_explicit(dt_global)
    advance_pore_pressure_analytic(dt_global)
    # actual ordered groups follow validated configuration, not hardcoded calls

    validate_and_commit_complete_state()
    t = t + dt_global
```

Remove RK-stage bookkeeping used only by IMEX integration (`q_rk`, `qp_rk`, `EX_terms`, `IM_terms`, per-stage `spatial_terms`, tableau matrices and weights, RK stage loops, and stage-specific implicit diagonal handling). Preserve any **workspace genuinely required by the FV spatial discretization, reconstruction, diagnostics, or source evaluation** even if its old name included an RK term, until its remaining uses are audited.

Eliminate source `CADENCE=STAGE/STEP` from the target configuration and merge the currently separate STAGE/STEP orchestrators into an ordered Lie group executor. Preserve a distinct domain-wide spatial advance primitive and atomic STEP-like substrate transactions as required by the physical methods. Use a **single full-step time variable**, with local interval parameter `dt`; there are no RK stage weights or rates to accumulate into final RK assembly.

## 8. Verification and acceptance

1. Pure spatial FV+injection forward Euler: conservation, hydrostatic well-balancing where applicable, CFL and positivity, time-dependent source scheduling, dry/wet transitions and boundary-source handling.
2. Lie ordering: later groups use the actually updated fluid/bed/inventory state; all members of one simultaneous explicit group use the same initial state; no double application; configured group order is deterministic.
3. Shared finite-inventory limiter: retain distinct individual process rates until joint limiting, including competing consumption and production, conservation across fluid/substrate, and proper mass-exchange bed transaction.
4. Backward-Euler rheology: reduced Newton residual/Jacobian, convergence and failed-solve handling, zero-momentum shortcut, no post-Newton momentum projection; test stiff linear friction against `u' = a-ku` and more representative nonlinear rheology.
5. Analytic Darcy: correct complete conservative update, pore-pressure-only versus interstitial gas-loss behavior, and physical admissibility.
6. Multiple operators on the same conservative components: Lie-order sensitivity, temporal convergence **to first order** on coupled noncommuting test problems, and relevant legacy v2 physical regressions. Do not mistake first-order splitting for a well-balanced or uniformly stiff-accurate method.
7. Whole-step state/geometry consistency and failure behavior, including updated topography before subsequent spatial flux evaluation; reject invalid method configurations and unsupported entries.
8. Compare wall-clock performance and memory with the prior IMEX-RK version; run relevant CI, and report exact tests actually executed.

**Known limitations:** Lie splitting generally has first-order global temporal accuracy even when component operators are exact. Treating friction implicitly avoids an explicit friction stability restriction but does **not** guarantee accurate gravity–friction equilibrium for arbitrarily large `dt`; test on stiff balanced regimes. No numerical convergence or production performance claim is made by this architecture document.

**Implementation status:** D18 is the target specification, not evidence that the Fortran refactoring is complete. No source edits or tests are claimed as part of this documentation revision.
