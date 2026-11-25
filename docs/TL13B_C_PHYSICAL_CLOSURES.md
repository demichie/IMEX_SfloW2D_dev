# TL13b-C — Interlayer physical closures and budgets

**Status (25 November 2026):** source-level audit and targeted regression
extensions started on `feature/two-layer`. **Scientific validation remains
open.** The reference code at the start of the audit is `2ccc70f`.
No governing equation, coefficient, numerical method or physical closure
has been changed by this work.

This document is the continuing scientific acceptance checklist for the
existing physical-source implementation. It supplements, rather than
duplicates, [the two-layer roadmap](TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md).

## 1. Source ownership and current implementation

| Process | Authoritative implementation | Current behavior |
| --- | --- | --- |
| L1/L2 settling and elutriation | `src/physics/interlayer_mass_exchange/interlayer_settling_2d.f90`, `interlayer_elutriation_2d.f90` | Independent proposed solid-volume rates, per solid class. |
| Paired solid/carrier limits and sources | `src/physics/interlayer_mass_exchange/interlayer_mass_exchange_2d.f90` | Common donor-inventory and packing limiter; antisymmetric mass, momentum and sensible-thermal contributions. |
| L1 substrate deposition and L2 air entrainment | `src/physics/mass_exchange/` | External sources, not interlayer transfers. Deposition is restricted to L1. |
| Chézy interlayer momentum drag | `src/physics/interlayer_momentum_exchange/` | Analytic drag-only momentum exchange at fixed layer masses. |
| Source order and atomic commit | `src/numerics/source_registry_2d.f90`, `lie_source_executor_2d.f90` | Ordered cell-local Lie transactions, rejecting an inadmissible candidate. |

The current Etna input enables L1 deposition, L2 atmospheric entrainment,
L1-to-L2 elutriation with the dynamic local L1 solid fraction, and Chézy
interlayer drag. It disables L2-to-L1 settling, both erosion sources,
L1 entrainment and direct L2 substrate deposition. Thus **the Etna run alone
cannot validate simultaneous bidirectional interlayer mass exchange**.

## 2. Implemented interlayer mass budgets

Let `S_i >= 0` denote accepted L2-to-L1 solid-volume settling rates and
`E_i >= 0` accepted L1-to-L2 elutriation rates [m/s]. Let `D_c >= 0`
be the accepted downward L2 carrier-volume rate and `U_c >= 0` the
upward L1 carrier-volume rate. The paired source uses:

```text
solid i, L1:  +rho_s(i) * (S_i - E_i)
solid i, L2:  -rho_s(i) * (S_i - E_i)

carrier, L1: +rho_c,2 * D_c - rho_c,1 * U_c
carrier, L2: -rho_c,2 * D_c + rho_c,1 * U_c
```

The densities multiplying carrier rates are the **source donor**
densities; recipient carrier-volume changes are obtained from transferred
mass using the recipient density. Consequently, carrier **mass** is
antisymmetric, but carrier **volume** need not be the same on the two
sides when densities differ.

For each interlayer source evaluation, L1 and L2 receive opposite
conservative rates for total mass, each solid species, any additional
gas species, x/y momentum and sensible thermal energy. Settling solids
carry an interpolated interfacial velocity and L2 temperature;
elutriated solids carry L1 velocity and temperature. Each gas stream
uses the properties of its own donor.

These equalities are **source-level invariants**, not statements about
the full-system energy including kinetic, gravitational, terrain and
boundary contributions.

## 3. Elutriation and gas co-lift

The proposed solid rates come from `LINEAR_USTAR` or `SHEAR_POWER`,
with availability and shear/Rouse factors. They are not prescribed by
the carrier-flow model.

For a positive, explicitly supplied
`ELUTRIATION_SOLID_VOLUME_FRACTION=alpha_lift`, the requested
carrier rate is:

```text
U_c = [(1 - alpha_lift) / alpha_lift] * sum_i E_i
```

With `ELUTRIATION_SOLID_VOLUME_FRACTION=-1` (also the default when
omitted), `alpha_lift` is taken from the **current L1 donor inventory**:

```text
alpha_lift = sum_i V_s,i,1 / (sum_i V_s,i,1 + V_c,1)
```

This is the **mean L1 mixture** assumption; it is not evidence that the
actual interfacial entrained parcel has that composition. The prescribed
alternative remains available. `LIMIT_SOLIDS` preserves the target
solid/carrier ratio by restricting the solid flux when necessary;
`ALLOW_CONCENTRATED` permits less co-transported carrier.

**Scientific question C1:** compare these assumptions in otherwise
identical tests and quantify L2 solid concentration, mass flux, gas
inventory, packing and runout. Do not select a physically preferred law
on the basis of one run merely completing without failure.

## 4. Settling, packing and simultaneous gross fluxes

The L2-to-L1 solid settling rate reuses the basal deposition-volume
kernel with L2 fluid properties. Unlike basal deposition, it transfers
solid material *into L1*, not into the substrate. The paired limiter
computes the **minimum downward carrier stream** needed to keep the
updated L1 composition below the configured packing limit. A downward
solid flux can therefore carry zero gas if L1 already has capacity.

Both gross solid streams and both gross carrier streams are limited
against their **original donor inventories for the Lie substep**:
material received simultaneously cannot immediately be spent. The
limiter also checks both final layer packings.

**Scientific question C2:** validate the rule assigning additional
carrier to settling solids, especially for a concentrated/dilute pair,
and study simultaneous sedimentation and elutriation rather than
assuming that their net solid flux is sufficient.

## 5. Compressible carrier and thermal exchange

The paired limiter converts accepted carrier mass into recipient volumes
through the donor-to-recipient density ratios. With unequal
temperatures and carrier densities, there is no general assumption
that the **sum of carrier volumes** is conserved.

The physical source uses donor-state sensible heat:
`rho_c * U_c * cp_c * T` for carrier, and
`rho_s(i) * E_i * cp_s(i) * T` for each solid. Source contributions
to the L1 and L2 thermal equations remain exactly opposite. Extra
carrier gas species, when active, follow the donor gas composition.

**Scientific question C3:** check thermodynamic admissibility after
source integration, consistency between ideal-gas density and mixture
conversion, and the physical interpretation of transferred volumes
at fixed prescribed pressure. Source antisymmetry alone does not
establish correct compressible-flow thermodynamics.

## 6. Drag: momentum conservation and mechanical dissipation

The analytic Chézy source holds both layer masses fixed, conserves their
sum of x/y momenta and relaxes their relative velocity. Its tests
verify a non-increasing kinetic energy and the exact quadratic
relative-velocity decay. The drag source leaves the existing
**sensible thermal energy** of each layer unchanged.

Thus the decrease in resolved mechanical kinetic energy is **not**
converted into a thermal source in the present implementation. This
is an explicit model property, not a test failure. It must not be
silently changed to enforce a purported total-energy balance.

**Scientific question C4:** specify whether drag dissipation is
unresolved (lost from the resolved mechanical system), is assigned
to one or both thermal budgets, or requires a separate physical
closure. Any change requires scientific approval and regression
evidence.

## 7. Lie ordering and external source accounting

For equal group IDs the registry orders sources by scope:
`L1`, then `L2`, then `INTERLAYER`. The current Etna example
assigns group 4 to L1 deposition, L2 ambient entrainment and
interlayer elutriation. Therefore **L1 deposition precedes
elutriation** in the accepted cell-local state and the dynamic lifted
fraction uses the updated L1 inventory. L2 air entrainment in that
group is evaluated **before** the same-group interlayer source.

In a cell with initially dry L2, the air-entrainment source cannot
operate on the newly elutriated material until a subsequent timestep
under this ordering. This is a Lie-splitting effect, not an additional
physical entrainment law.

The full mobile-layer mass balance includes external L2 ambient-air
input and external L1 deposition output, besides domain fluxes and any
other enabled sources. Interlayer rates cancel in the sum of L1+L2.
For **solid species** (with no other active solid sources and absent
domain fluxes), L1+L2 mobile solid mass plus mass deposited into the
substrate is the appropriate closure check. Treat any programmed
buoyancy-removal mechanism separately in integrated budgets.

**Scientific question C5:** test source-order sensitivity and
conservation over complete timesteps with deposition and entrainment.
Do not change group ordering without an explicit numerical and
scientific decision.

## 8. Verification matrix and acceptance criteria

| Case | Existing evidence or new targeted coverage | Still required |
| --- | --- | --- |
| Prescribed composition, dry/wet pair | Cell-local and paired Lie regressions | Fresh strict/optimized execution on current revision. |
| Dynamic composition, L1 donor inventory | Cell-local 10%/20% cases and post-spatial donor fixture | Compare physical output with prescribed mode. |
| Opposite gross solid/carrier flows, unequal densities | `test_interlayer_mass_exchange` now also covers dynamic mode for both carrier policies | Check physical closure validity over representative parameter ranges. |
| Dynamic, simultaneous settling/elutriation | `test_interlayer_lie_transaction` now includes full two-layer conservative/primitive update | Run both build profiles and a sustained representative scenario. |
| Source-level mass/species/momentum/thermal antisymmetry | Cell-local tests with distinct donor densities, cp and temperatures | Integrated budget over multiple Lie groups and timesteps. |
| Chézy momentum and kinetic dissipation | Exact analytic unit test | Scientific decision about unmodelled dissipation heating. |
| L1 deposition + L2 atmosphere + interlayer co-lift | Etna input and example documentation | Integrated layer/bed/atmosphere budget and refinement study. |

**Technical exit:** strict/debug and optimized source regressions,
paired Lie admissibility, finite conserved-source budgets, cold-start
and wet/dry tests, and representative two-layer OpenMP/end-to-end
checks on the **same commit**. Run at least
`sh TESTS/TEST_SOURCE_REGISTRY/run_test.sh` and the existing
`solver-validation.yml` workflows.

**Scientific exit:** resolve C1-C5, assess the empirical parameters
and compare with suitable analytical/benchmark/laboratory evidence.
This is separate from TL13b-A hyperbolicity and TL13b-B spatial
positivity and does **not** authorize removing
`IMEX_TWO_LAYER_VALIDATION=1`.

The last verified solver-integration Action currently predates the
dynamic co-lift change; the new test extensions have not yet been
compiled or executed. No claim of scientific acceptance is made.
