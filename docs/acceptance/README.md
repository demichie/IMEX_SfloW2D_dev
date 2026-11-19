# Historical scientific acceptance records

This directory retains dated reports and machine-readable evidence
for earlier scientific acceptance decisions. **Some `results_*.json`
files are read directly by active regression programs in
`TESTS/ACCEPTANCE/` and must remain in the checkout.** Do not remove
or relocate them without migrating and verifying all those tests,
their plans and scientific checks.

Dated reports and their JSON files are evidence of the explicitly
recorded source revisions and configurations, not current
compatibility promises or results for the present branch. Retained
report links to raw data also point to the immutable
[pre-cleanup acceptance snapshot](https://github.com/demichie/IMEX_SfloW2D_dev/tree/d91843a0faf8c8491cf50c240cdd50a7cd006407/docs/acceptance),
which fixes the historical revision and prevents accidental
reinterpretation of later test artifacts.

The former very large D03 source-state manifest, which has no active
test consumption in the audited paths, is recoverable from that
snapshot. The current executable acceptance definitions live in
`TESTS/` and the associated workflows. Always check the exact SHA.
