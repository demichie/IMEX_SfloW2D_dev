# D17 — Audit of examples, input generators and launch infrastructure (2026-10-25)

**Repository:** `demichie/IMEX_SfloW2D_dev`; **branch:** `feature/d17-source-refactor`; **audited HEAD:** `d5d4ad25083b3b2320e151fb9b736406bcad6bfd`.

## Post-audit update — 25 October 2026

On branch `feature/d17-source-refactor`, commit `1d60cc8f` hardens `create_inputfiles.py`, `launch_jobs.py` and `templatedir/run.sh`: selective text placeholder substitution, preservation of binary/symlink/unchanged files, atomic replacement, preflight executable checks and failure propagation. **User-reported functional verification: working** (25 October 2026). No accompanying log, member count, input/output comparison, or negative-case evidence was supplied. Therefore the fix is **implemented and user-verified**, not an automated or complete end-to-end acceptance. The open work is targeted failure-case tests, generated-input/parser evidence, and smoke-test execution of other D17 example families. The original findings below describe the state at the prior audited HEAD and should not be read as unresolved source defects without rechecking.

## Status and method

**STATIC/PARTIAL — NOT AN END-TO-END ACCEPTANCE.** Queried the recursive Git tree (618 entries, not truncated), inventoried the 122 tracked files in `EXAMPLES/`, compared the D17 roadmap and earlier `D17_EXAMPLES_TESTS_INPUT_GENERATORS_AUDIT_2026_10_25.md`, checked the documented example release verification and inspected the principal generator/runner source code. No D17-HEAD solver execution or generated-input parser run was performed. The previous C06 example acceptance refers to earlier pinned revisions and must not be transferred to the D17 HEAD. Do not interpret a clean source scan as an executable PASS.

## Inventory

20 example families are documented by the existing second-pass guide audit (including the four BUMP regimes as separately standardized cases); 19 legacy `README.txt` files are retained. The current tree contains generated-input templates alongside generators for 2D, 2D_PP, BOX, four BUMP variants, BW, BW2D, DAMBREAK, ENGWELL2016, INCLINED_EXCAVATION_2D, LATERAL_SOURCE and SLOPE. CF, ETNA, ETNA_LAVA, TAAL and TAAL_BOTTOM contain tracked `.inp` files; the ensemble Etna template contains its own `.inp`, launch script and preprocessing script. Template absence alongside a tracked input is not by itself a defect.

Separate facilities: `tools/docs/build_examples.py`, `tools/release/build_archives.py`, visualization utilities and multiple acceptance generators under `TESTS/`. Their existence was inventoried but no full dynamic regression was run in this audit.

## Confirmed source-level risks requiring action

1. **P0, ensemble generator text/binary safety:** `EXAMPLES/EXAMPLE_ENSEMBLE_ETNA/create_inputfiles.py::replace_strings` loops through every regular file in each copied member directory, opens each with `open(fname,'r')` and then unconditionally rewrites it with `open(fname,'w')`, even if no `ENSEMBLE_...` placeholder was found. This is neither restricted to known text templates nor an atomic write. Non-UTF-8/binary regular files would trigger decoding failures; unchanged text can be rewritten unnecessarily; partial failure can leave a member in an incomplete state. Restrict replacement to an explicit allowlist or detectable text templates, preserve unchanged files and fail on unresolved placeholders; maintain original physical parameters and substitution semantics. `shutil.copytree(..., symlinks=True)` makes symlink handling part of the required safety review. This is a **source-level confirmed hazard**, not evidence that any particular ensemble run is currently corrupted.
2. **P1, ensemble launch contract:** `EXAMPLES/EXAMPLE_ENSEMBLE_ETNA/launch_jobs.py` launches `./run.sh` with `cwd=directory`, passing an extra member-name argument; the current `templatedir/run.sh` does not use that positional argument and expects `./IMEX_SfloW2D` in the member's current directory. The tracked template does not include that executable. Verify how provisioning/links are created and ensure a clear preflight failure when missing. The launcher does check per-member return codes and limits active processes to five; preserve both guarantees. Do not silently change the supported deployment model.
3. **P1, ensemble workflow error propagation:** `templatedir/run.sh` performs preprocessing followed by solver execution without `set -e` / an explicit preprocessing-status check. A preprocessing failure can be masked by a successful subsequent command. Preserve log destinations while ensuring failure propagation.
4. **P1, dam-break CLI documentation drift:** `EXAMPLE_DAMBREAK/create_dambreak.py` is the committed generator, yet its docstring/usage examples still mention `create_dambreak_static.py`. Update text-only instructions without altering its ten-argument interface or initial-condition construction.
5. **P1, template substitution validation:** multiple generators use raw `str.replace` on `IMEX_SfloW2D.template` placeholders. Before marking D17 compatibility, generate the actual inputs, check unresolved tokens and native/legacy namelists, and run the production parser's negative/positive cases. The expected fallback for absent `SOURCE_GROUPS` and `THERMAL_COOLING_FLAG` is intentional; do not force-convert existing legacy examples.
6. **P1, model-specific acceptance:** exercise `EXAMPLE_ETNA_LAVA` with cooling absent/present configuration preserving the D17 legacy default; check the lateral source's `ALPHAS_SOURCE` volume fraction semantics, the deposition/carrier loss coupled STEP transaction and restart/inventory balances. These are validation obligations, not confirmed failures.

## Historical evidence, not D17 approval

`EXAMPLES/README_SECOND_PASS_AUDIT.md` documents complete runs of the 20 standardized examples on an older synthetic `9.9.25` release candidate, including the 20-member Etna ensemble, and a narrower `9.9.27` candidate verification. These are valuable baselines but not runs of the currently audited D17 branch. The separate six-job historical physics CI run 40175887639 passed on `cb129420`; those tests do not substitute for a D17-HEAD execution of all example families.

## Required acceptance matrix

| Scope | Current finding | Dynamic D17 status |
|---|---|---|
| Legacy-generated inputs: 2D, 2D_PP, BOX, BUMP (4), BW, BW2D, DAMBREAK, ENGWELL2016, INCLINED_EXCAVATION_2D, LATERAL_SOURCE, SLOPE | Templates and generators present; inspect outputs, dependencies and parser results | NOT RUN |
| Tracked input examples: CF, ETNA, ETNA_LAVA, TAAL, TAAL_BOTTOM | Namelists present; previous C06 runs belong to older revisions | NOT RUN |
| Etna ensemble sample/input generation, per-member launch, probability aggregation | Generator/launcher hardened in `1d60cc8f`; user reports workflow works; negative-case and full-member evidence not supplied | USER-VERIFIED (scope unspecified); FULL GATE OPEN |
| TESTS-generated inputs, scripts, acceptance harness | Inventoried; no blanket compatibility claim | NOT RUN |
| Docs/release catalogue builders and UTILS postprocessing | Inventoried; prior packaging evidence version-scoped | NOT RUN |

## Next coherent change set (without violating D17)

(1) Fix ensemble generation and runner reliability in one coordinated reviewable change, with unit tests for placeholder replacement, non-text preservation, failure propagation and missing binary; (2) add a manual-only example-smoke CI or reproducible local runner that checks generated inputs, solver exit, finite diagnostics, NetCDF outputs, mass/restart where relevant and records each case/revision; (3) use existing 20-example C06 results strictly as reference baselines; (4) update this report and roadmap only after actual D17 executions. Never reintroduce numeric physical source IDs, removed sub-solvers, legacy adapters, or alter physics/tolerances to obtain green tests.
