# D17.15 — Time-step ordering

## Audit result

The current driver and numerical kernels follow the ordering in the D17 source
contract without requiring a change to the external orchestration:

1. `IMEX_SfloW2D.f90` computes and limits `dt`, prepares source-window averages,
   then prepares the stochastic timestep.
2. `time_integration_2d.f90` constructs each known stage state and converts it
   before the ordered implicit-group sweep.
3. The implicit coordinator returns the effective IM increment; requested EX
   rates are evaluated on the solved stage, followed by the spatial operator
   and its boundary hook.
4. The final RK assembly applies the existing positivity/dry-state handling
   and primitive reconstruction before returning to the driver.
5. The driver applies the registered atomic mass-exchange batch after RK and
   before advancing runtime time. That batch applies the source-coverage mask,
   updates state/inventories and reconstructs the primitive state before the
   accumulated bed update.
6. Runtime time advances before existing diagnostics and output scheduling.

No source evaluation was moved between stages, and no blanket EX evaluation
was added. The existing `t` argument and tableau-driven `need_explicit_stage`
behavior are unchanged.

## Local evidence

`TESTS/TEST_IMEX_STAGES/run_test.sh` passed in strict and optimized profiles,
covering observed stage trajectories and final assembly for N_RK=2/3/4, mixed
ANALYTIC/NEWTON group trajectories, and boundary/source reconstruction.
`TESTS/TEST_MASS_EXCHANGE/run_test.sh` passed in strict and optimized profiles,
covering aggregated post-RK exchange, source masks, state and inventory
budgets, bed updates, and checkpoint continuation.

## Remaining acceptance gap

These tests exercise the stage integrator and post-RK exchange separately.
There is no automated single-timestep test of the production driver that
captures the full sequence, including source-window averaging and OU
preparation. Therefore D17.15 is **partially verified**, not complete. No
cross-platform or CI acceptance is claimed.
