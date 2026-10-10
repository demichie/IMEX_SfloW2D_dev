# Gate-H full-field excavation comparisons

The compact four-sided one-cell excavation uses `dx = dy = 0.4, 0.2, 0.1 m`.
The nominal 8 m by 6 m depression is translated by half a cell in each
direction and built from one globally continuous Q1 vertex field. Initial
thickness is its exact cell average, giving the unexcavated inclined plane
as the initial free surface.

The test retains the existing three grids, G=1 and slope-plus-curvature
modes, aggregate limits and zero-uphill requirement. It additionally checks
every final thickness and free-surface cell against two independent initial
value trajectories: current-policy IMEX at common accepted times and the
unchanged historical SSPRK2 prototype under its separately approved limits.
`reference_metrics.csv` remains the aggregate regression: positivity,
conservation, uphill transfer, G=1 lateral symmetry, mass partitions, moments,
total variation and extrema. Passing these aggregates alone is not proof of
full-field agreement. Both clean-core checksums are frozen in the contracts.

The test-only observational main records accepted times immediately before
the production IMEX call. It reuses the tested numerical objects and compiler
flags when available, changes no production source, and must produce
byte-identical canonical conservative output. The accepted-time stream and
canonical output hashes are checked before comparison. The reference uses
those times and common initial fields, never actual intermediate states.

`field_contract.json` preserves the original autonomous adaptive comparison.
`field_contract_v2.json` records the user-approved separation from the
common-time spatial/stage gate; all error limits are unchanged. The autonomous
CFL comparison remains non-satisfied under `D-N7-CFL`, not a passing N7-D
assertion. Add `--autonomous` to `compare_fields.py` to reproduce and record
that diagnostic separately; its failure is retained in field evidence.

The current `field_contract_v3.json` also records the user-authorized dry-only
momentum projection of the accepted final IMEX state. Only the two momenta
are reset at the unchanged shared depth threshold; conservative mass, thermal
energy and components remain intact. The current-policy reference explicitly
matches this mapping, while the historical SSPRK2 core remains unchanged.

Run from the repository root to retain fields, accepted times and evidence:

```sh
KEEP_TEST_WORKDIR=1 OMP_NUM_THREADS=4 sh TESTS/TEST_PCCU_GATE_H/run_test.sh /path/to/IMEX_SfloW2D
```

The harness respects the requested thread team;
the full acceptance runner checks actual one/four-thread runs in both compiler
profiles. Compilation of the observation-only main requires gfortran and the
configured NetCDF development libraries, as do the other solver tests.

For provenance, the earlier port audit reported the following historical
slope-plus-curvature differences; these are not current acceptance limits:

| dx [m] | L1(h) | L1(eta) | Linf(h) [m] |
| ---: | ---: | ---: | ---: |
| 0.4 | 7.1697e-5 | 2.7927e-5 | 1.0116e-3 |
| 0.2 | 6.5239e-5 | 2.5412e-5 | 8.6013e-4 |
| 0.1 | 6.4436e-5 | 2.5099e-5 | 1.2830e-3 |
