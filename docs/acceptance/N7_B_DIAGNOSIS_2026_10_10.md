# N7 variable composition equilibrium diagnosis

This is the historical pre-correction diagnosis. N7-B is not complete. A defined variable-composition equilibrium failed in both
strict and optimized builds of `e4333a5c`, with identical results between actual
one/four-thread teams within each build. Production sources were unchanged
during this diagnosis. The subsequently approved correction and its passing
equilibrium regressions are recorded in the
[correction report](HYDROSTATIC_ROUNDOFF_2026_10_10.md); the contact reference
remains open. The CFL and clean reference kernels are unchanged.

The [evidence index](results_n7_b_diagnosis_2026_10_10.json) records the unchanged
production revision, frozen contract, executed test hashes, case findings,
binary fingerprints and retained raw directories. The
[test README](../../TESTS/TEST_N7_COMPOSITION/README.md) gives reproduction
commands and distinguishes diagnostics from acceptance.

## Equilibrium failure and mechanism

The first witness is a fully wet 20-cell liquid/two-solid mixture over a
continuous Q1 bed. Temperature is 300 K, free surface is constant and varying
mass fractions are constructed in the specific-volume nullspace so mixture
density and Gamma are constant mathematically. Both solid densities and heat
capacities are distinct. Hydrostatic balance is therefore defined independently
of the production conversions.

With RK2 and `dt=0.001 s`, the first step leaves mass, thermal energy and component
masses unchanged. It generates a velocity of about `1.27e-17 m/s` from the
roundoff-level hydrostatic residual. At the second step, the exact-zero
normal-discharge condition in `hyperbolic_2d` no longer suppresses scalar flux.
The central-upwind dissipative term depends on endpoint jumps and acoustic
speeds, rather than on that tiny physical velocity. The component-flux limiter
also activates on some faces. Independently reconstructing these scalar fluxes,
including the limiter and exact-rest guard, reproduces the observed spatial
term within `1.5e-17` of the corresponding component scales in both profiles.

After 20 steps (`t=0.02 s`), the strict witness has maximum fraction changes
`3.082e-5` and `4.161e-5`, and a temperature change of `0.006459 K`. Optimized
changes are `3.081e-5`, `4.159e-5` and `0.004996 K`. These are small absolute
changes but many orders above the frozen equilibrium tolerance. The maximum
component-scaled conservative-state discrepancy is approximately `4.805e-4`
versus the fixed limit `7.276e-12`.

This witness preserves total mass exactly and global thermal/component budgets
within roundoff, with no final repair. Global conservation therefore does not
establish local well balancing. Other gas-containing and constant-pressure
diagnostic trajectories also fail, and several additionally violate the closed
boundary budget/flux criteria. Those failures are retained; they are not hidden
by the first witness's successful global budget.

## Contact reference limitation

The moving-contact fixture prescribes a central uniform velocity with smooth
outer transitions to zero. Its scalar reference assumes uniform velocity. At
80 cells, RK3 and `t=0.1 s`, the comparison's outer cells contain coupled effects
absent from that reference. The largest scaled discrepancy is about `4.08e-8`
for the liquid case. As an explanatory probe, the central pulse interval
`16 <= x <= 24` agrees to about `6.4e-16` across the three strict closures.
That probe does not replace the frozen `12 <= x <= 28` criterion or turn the
failed comparison into a pass.

A new transport fixture/reference contract needs demonstrable isolation from
the velocity transition or an independent oracle that includes the transition.
The current contact discrepancy must not be classified as a solver defect based
on this scalar comparison alone. No tolerance or frozen state has been tuned.

## Executed coverage and existing regressions

The partial diagnostic command completes 144 solver runs: 36 cases per profile,
each with actual one/four-thread teams. This is 18 full-duration 1D equilibrium
probes, 12 shorter prefixes and six moving contacts per profile. Both profiles
also replay the first frozen acceptance case and reject it. The planned
162-case/profile inventory has not completed.

The existing IMEX-stage suite passes in both profiles. Rebuilding its old HEAD
driver against the same production objects and replaying the suite gives
bitwise-identical payloads in all 396 old/new comparisons. The optional
composition observer therefore preserves the existing default test output.
N7-A also passes all 102 cases per profile, with actual one/four-thread payload
agreement, for 408 solver runs. Four new checker unit tests and shell/diff
checks pass. These results do not satisfy the new failing N7-B gates.

## Correction and unresolved decisions

The numerical issue was amplification of a roundoff-level departure from rest
into finite scalar diffusion. The approved correction now applies a fixed
pressure-scaled arithmetic allowance only to exactly stationary stencils;
it does not introduce a velocity cutoff. All original 90 equilibria, including
the boundary flux and budget checks after HP finalization, pass in both profiles.

After fixing the contact reference, complete the full transport matrix before
closing N7-B. The bounded numerical correction has its own regression package
and does not depend on declaring the flawed contact oracle valid.
N7-C, N7-D, D-N7-CFL and N9 remain pending; the
architecture/implicit-groups refactor has not started in this package.
