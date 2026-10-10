# Cell exchange acceptance completion

N8-A and criterion N8-01 are satisfied by the extended production exchange
tests. No production source, exchange law or acceptance tolerance changes.
N8 remains open for nodal refinement and repeated-update tests in N8-B,
including the explicit D-N8-AREA interpretation. N7 and N9 remain open.

The [execution evidence](results_n8_a_2026_10_10.json) records the candidate
base revision, compiler, flags, source/fixture hashes, diagnostic values and
payload fingerprints. The [closure checklist](N7_N8_CLOSURE_2026_10_10.md)
retains the original requirements and moves N8-A to completed actions.

## Flat cell and vertex evaluations

Before each of the 18 production transactions, the same physical cell state
and erodible inventory are passed once to `eval_mass_exchange_terms` with zero
slope, then to four identical vertex evaluations in an OpenMP loop. All outputs
are compared bitwise: solid erosion/deposition, carrier erosion/loss,
conservative equation sources and the bed proposal. Their quarter-weighted sum
also agrees within the unchanged roundoff tolerance.

Each comparison is repeated with a dry state, requiring every output to be
zero. Input state/inventory and the shared temperature, viscosity, permeability,
substrate density, erosion coefficient, porosity and packing configuration must
remain bitwise unchanged. This checks constitutive evaluation identity; it does
not replace the cell-local production update with a vertex-based exchange law.

## Gas carrier transactions

The added fixture has two solid classes and an ambient-gas carrier, with one
conserved pore-pressure variable. It uses the existing OFF-inhibition,
fixed-permeability law; Sutherland viscosity, vaporization and entrainment are
disabled. The optional liquid-loss parameter is deliberately unallocated.
Initial conservative states are constructed from known phase volumes and
temperature. Expected gas density, signed Darcy requests and volume limits
are recovered independently, without calling the production exchange evaluator.

| Case | Required assertion |
| --- | --- |
| gas_packing | Positive pressure-driven loss reaches the packing cap. |
| gas_reserve | Positive loss reaches the available-carrier reserve before the packing cap. |
| gas_exhausted | An exhausted reserve produces zero loss, not a spurious carrier gain. |
| gas_inflow | Negative excess pressure preserves signed gas inflow. |
| gas_combined | Saturated solid deposition/erosion and capped carrier loss form one conservative transaction and bed proposal. |
| gas_masked | Full/partial masks scale that transaction, inventories and bed proposal before nodal gathering. |

All seven conserved equations are checked individually, so the thermal-energy
scale cannot hide a momentum or carrier error. Assertions also cover
nonnegative carrier/solid inventories, unchanged dry/inactive flow states,
deposit/erosion inventories, the mass-lumped nodal proposal, weighted bed volume
and refreshed Q1/LS/Hessian/G geometry. Negative gas loss remains an inflow term;
no nonnegative clamp or new gas law is introduced.

## Verification

Run from the repository root:

```sh
KEEP_TEST_WORKDIR=1 sh TESTS/TEST_MASS_EXCHANGE/run_test.sh
```

Both strict and optimized production modules pass with actual one/four-thread
teams: 72 production transactions and diagnostic comparisons, 144 wet/dry
single-versus-four evaluation comparisons, and 576 individual vertex bitwise
comparisons. Within each profile, conservative/geometry/inventory snapshots,
checkpoint files, gas snapshots and flat-evaluation payloads agree bitwise
between teams.

All eight historical snapshot/checkpoint fingerprints from the original 12
transactions still match [the previous N8 evidence](results_n8_mass_exchange_2026_10_09.json).
The diagnostic checker rejects five deliberately invalid logs: a missing flat
assertion, a missing gas case, a wrong actual team, a corrupted volume
diagnostic and a nonfinite reference. Existing Fortran and diagnostic
tolerances remain unchanged.

All archived production source hashes still match the preceding 36-run audit.
This bounded test-only change reruns the extended exchange suite, not that full
audit. The results do not approve N9 or establish other inhibition laws,
additional gases, arbitrary rheologies or a nonuniform production mesh.
