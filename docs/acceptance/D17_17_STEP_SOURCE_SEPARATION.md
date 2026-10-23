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

## Dependency-boundary audit

`TESTS/ACCEPTANCE/test_acceptance_tools.py` checks the complete `numerics` and
`spatial` source directories for direct imports of physical law modules. It
also checks that rheology, thermal, pore-pressure and mass-exchange law modules
do not import numerical or spatial orchestration. Spatial modules may use the
small approved `equation_terms_2d` flux, hydrostatic-coefficient, speed-bound,
boundary-source and component-flux APIs; physical source-law imports remain
behind the equation bridge. The boundary test passes against the current tree.

This is a static module-boundary audit, not a complete call-graph or
term-by-term runtime audit. The full D17.20 term/physics matrix and the
term-selective mass-exchange evaluation remain open.

## Verification

`TEST_MASS_EXCHANGE/run_test.sh` passed in strict and optimized profiles,
including transaction, inventory, bed-projection and checkpoint-continuation
cases with 1 and 4 OpenMP threads. `TEST_IMEX_STAGES/run_test.sh` passed in
strict and optimized profiles, verifying the production sequence remains RK,
post-RK step sources, then clock advancement. The full strict-debug solver
build also passed. The focused static boundary test passes. The complete
Python acceptance-tool suite has one unrelated failure:
`test_n7_b_and_geometry_completion_evidence` expects a different pinned
digest for `TESTS/TEST_GEOMETRY_1D/run_test.sh`; the observed digest is
`b650775807bc8b3afc9378577d03a57744dc832da33184ac608396f3c85aa56e`, while
the test expects `b7471ae4c748750ec99cc96a41bf078ad5264beb056d2ae213fdaf6f5bb3a1df`.
The baseline digest was not changed.

This change does not complete D17.17–D17.18: call-graph/runtime verification
and term-selective mass-exchange evaluation remain open, as does the full
D17.20 acceptance matrix.
