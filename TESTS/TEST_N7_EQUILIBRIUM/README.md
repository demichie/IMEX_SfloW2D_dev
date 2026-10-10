# N7 equilibrium and one dimensional regression tests

This package implements the bounded N7-A inventory: finite-time smooth
transport, a wet/dry Ritter dam break, fully wet parabolic lakes, inclined
one-cell excavation and the complete 1D/2D lake flag matrix. It uses the
unchanged production raw-stage driver from `TEST_IMEX_STAGES`, not a replacement
time integrator. The script does not edit production sources or frozen reference
kernels. Required tests exposed two production defects corrected in this package:
near-rest Richardson overflow and a fused-roundoff HP discharge at a clipped
zero velocity bound. The [acceptance report](../../docs/acceptance/N7_A_EQUILIBRIUM_2026_10_10.md)
records their regressions and the changed optimized dynamic outputs.

Run from the repository root:

```sh
KEEP_TEST_WORKDIR=1 sh TESTS/TEST_N7_EQUILIBRIUM/run_test.sh
```

The script compiles strict and historical optimized production modules,
with a strict observer in both executables. Every case runs with actual one
and four-thread teams, and the complete binary payloads must match bitwise.
The retained temporary directory contains inputs, logs, raw final/stage arrays,
all-step diagnostics and `evidence.json` for each profile.

## Frozen cases and criteria

[`contract.json`](contract.json) fixes grids, states, times, supported stages,
flags and error limits before execution. Each profile has 102 cases: 72 lakes
(three beds, two dimensions, four flag pairs and three stage counts), nine
transport refinements, nine Ritter refinements and 12 excavation cases.

Lake tests compare every conservative component against its unchanged initial
value, test all raw last-step states and verify the residual inferred from the
second known stage. The independent frozen spatial core must also return zero
residual. Momentum errors never use the much larger thermal-energy scale.

Smooth transport uses a Gaussian temperature pulse in a uniform moving liquid
on a flat bed. With limiter zero, its independent scalar central-upwind
semidiscrete operator is linear; the same explicit IMEX weights give a direct
full-field numerical oracle. An analytic translating Gaussian additionally
checks spatial refinement. This is a passive thermal-transport case, not a
variable-composition test.

The Ritter gate compares thickness and momentum against exact cell integrals
of the shallow-water rarefaction, using the configured ambient-density correction
to gravity. A separate full-field reference trajectory applies the current
IMEX tableaux to the checksum-pinned 1D kernel. Both the analytic error limits
and refinement-ratio limits must pass.

The excavation is cut into a plane descending in positive x, with shared nodal
walls spanning exactly one cell. It must accelerate downslope while preserving
total mass and resolving no uphill mass outside the western crest. Full-field
reference trajectories cover all four slope/curvature flag pairs.

## Raw states and dry canonicalization

Every step checks known, solved and raw-final minimum mass, raw temperature,
component/carrier fractions, implicit status and used versus computed CFL bound.
All-step raw and final mass budgets are checked separately.

The existing rest-transport suppression can initially leave a nonzero pressure
impulse in a cell with exactly zero mass. Production resets this dry momentum
at final assembly; the reference adapter applies the same explicit mapping only
after capturing its raw state. The test compares each raw field and each step's
repair magnitude before accepting that dry cleanup. It does not allow clipping
negative mass to establish positivity. Reference stages are checked before the
kernel's own sanitization is called.

The comparison uses the current IMEX weights, not the reference's SSPRK2 stepper
or retry mechanism. This package does not resolve D-N7-CFL or approve N9.
