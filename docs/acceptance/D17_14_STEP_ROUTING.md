# D17.14 — Post-RK `step` routing

## Implemented scope

The production driver now passes the validated source registry to the
post-RK mass-exchange fractional step. That routine checks the ordered `step`
schedule, dispatches the existing aggregated mass-exchange transaction once
per step (one aggregate evaluation per solve cell), and returns without
changing state when no mass-exchange terms are registered. The batch remains
after RK completion and before the simulation clock advances; it is not
evaluated in RK stages.

The current source catalog exposes mass exchange as the sole post-RK batch.
Non-mass `step` terms are rejected explicitly until their evaluation kernels
are implemented. The legacy direct-call interface remains available for
standalone callers that do not pass a registry.

## Regression evidence

`TESTS/TEST_MASS_EXCHANGE/run_test.sh` passes in strict and optimized
profiles, with actual OpenMP teams of one and four threads. The production
mass-exchange cases now build and pass the configured source registry through
the post-RK dispatcher. Existing budgets, inventories, topography, projection
and checkpoint-continuation comparisons remain unchanged. An additional
regression verifies bitwise that an initialized registry with an empty `step`
group leaves the conservative state, inventories and bed unchanged even when
legacy erosion flags are active.

No cross-platform or CI acceptance is claimed.
