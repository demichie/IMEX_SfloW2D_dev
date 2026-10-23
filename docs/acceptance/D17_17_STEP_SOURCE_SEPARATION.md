# D17.17 — Post-RK step-source separation

## Change

The post-RK mass-exchange transaction has been moved from
`time_integration_2d` to `mass_exchange_step_2d`. The IMEX integrator now
contains the RK tableau, stage progression and final assembly, while the
dedicated step-source module owns the registered mass-exchange transaction:
registry validation, cell-local evaluation through `equation_terms_2d`,
source-coverage masking, conservative and inventory updates, and coordination
with the existing geometry-owned bed update.

`simulation%advance_timestep` continues to call the executor after RK and
before advancing the runtime clock. The step remains atomic and post-RK; no
rate formula, limiter, mask, geometry operation or physical ordering was
changed. `apply_mass_exchange_increment` moved with the transaction, since it
has no other production caller.

The `time_integration_2d` module no longer exports or owns mass-exchange
transaction routines. Production and regression-test callers import the
executor directly. Its compilation dependencies are declared in `src/Makefile.am`.

## Verification

`TEST_MASS_EXCHANGE/run_test.sh` passed in strict and optimized profiles,
including transaction, inventory, bed-projection and checkpoint-continuation
cases with 1 and 4 OpenMP threads. `TEST_IMEX_STAGES/run_test.sh` passed in
strict and optimized profiles, verifying the production sequence remains RK,
post-RK step sources, then clock advancement. The full strict-debug solver
build also passed.

This change does not complete D17.17–D17.18: the broader audit of numerical and
physical module dependencies, term-selective mass-exchange evaluation, and the
full D17.20 acceptance matrix remains open.
