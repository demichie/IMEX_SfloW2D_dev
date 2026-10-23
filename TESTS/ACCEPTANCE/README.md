# HP-PCCU acceptance evidence audit

This audit records N7, N8 and N9 evidence for a pinned solver revision. It builds
tracked sources from `git archive`, runs existing tests, records actual OpenMP
teams and retains failures. A successful command does **not** approve N9: the
coverage gaps in `acceptance_plan.json` require separate tests and scientific
sign-off. No production source, numerical formula or reference tolerance is
changed by the audit.

## Run the audit

From the repository root, choose a new output directory outside the repository:

```sh
python3 TESTS/ACCEPTANCE/run_acceptance.py \
  --revision 752b47fe853f68798d2ec84ca4cb905a53e24300 \
  --output /tmp/imex-acceptance-new-run \
  --netcdf /opt/homebrew/opt/netcdf-fortran \
  --jobs 3
```

Replace the revision, output path and NetCDF prefix as appropriate. The output
directory must not already exist. The optional `--roadmap FILE.docx` records the
roadmap hash; `--reference-zip FILE.zip` verifies the two dry-safe core hashes
listed in the pinned Gate-H README. That hash check does not execute the external
Python scientific-reference suite.

For a candidate not yet committed, add `--working-tree` and use the current HEAD
as `--revision`. This explicitly overlays tracked and nonignored untracked files
in the temporary export; ignored build artifacts are excluded. The manifest
labels the Git SHA as a **base revision only**, retains the candidate patch and
hashes the exported source/test files. It does not create a Git commit. Without
this option, uncommitted edits are deliberately not tested.

The audit requires Git, tar, Autotools, make, a configured Fortran compiler,
NetCDF Fortran, LAPACK, Python 3 and NumPy. M12c made the build compiler
configurable and added CI integration; this driver tests the selected local
toolchain but does not by itself establish portability to other compilers.

It returns nonzero for failed builds, failed test runs, failed available
one/four-thread comparisons or a supplied reference-package mismatch. The JSON
manifest retains failures and sets `n9_accepted` to false regardless of the
command's return code. The per-command timeout defaults to 900 seconds and can
be adjusted with `--timeout`; a timeout is a failure, not a skipped test.

## Evidence and interpretation

The [N7/N8 closure checklist](../../docs/acceptance/N7_N8_CLOSURE_2026_10_10.md)
and [machine-readable plan](n7_n8_closure_plan.json) distinguish satisfied
original contracts, mandatory missing tests and broader coverage limitations.
`closure_status` and `remaining_actions` in the execution criteria mirror that
plan. `satisfied` does not assert universal model coverage or N9 sign-off;
`execution_status=pass_for_tested_scope` alone never closes an open requirement.
The checklist retains historical audits and separately records each bounded
completion run. N8-A is completed by the
[cell-exchange report](../../docs/acceptance/N8_A_CELL_EXCHANGE_2026_10_10.md).
The [projection report](../../docs/acceptance/N8_B_PROJECTION_2026_10_10.md)
completes N8-B and the user-approved area interpretation. N8 is satisfied;
N8-A/N8-B/D-N8-AREA are in `completed_actions`. The
[equilibrium and 1D report](../../docs/acceptance/N7_A_EQUILIBRIUM_2026_10_10.md)
completes N7-A with 408 production runs and a new 37/37 full audit. It explicitly
records two bounded production corrections and changed optimized dynamic
fingerprints; output neutrality is not claimed for that historical package.
The subsequent N7-B/N7-C/N7-D packages are completed. The
[CFL completion report](../../docs/acceptance/CFL_POLICY_COMPLETION_2026_10_11.md)
records the authorized dual-policy contract, 45/45 full-audit groups, bitwise
retained-default compatibility and additional original-duration `combined`
restart checks. N7 and N8 are satisfied. N9 scientific-baseline approval is a
separate owner decision, now recorded for the pinned v5 scope in the
[v5 report](../../docs/acceptance/N9_RESTART_V5_2026_10_11.md) and
[machine-readable evidence](../../docs/acceptance/results_n9_restart_v5_2026_10_11.json).

`IMEX_AUDIT_CFL_POLICY=combined` explicitly selects the optional policy in
supplementary restart fixtures through the capture adapter. Every override is
recorded in effective inputs and invocation manifests; leaving it unset
preserves the normal audit fixtures. It is not silently applied by the runner.

`manifest.json` contains the solver Git SHA, tool probes, compiler flags,
executable and source checksums, test results and coverage criteria. Each test
has a log; each solver invocation has its own input snapshots, full solver log,
actual thread count, output SHA-256 fingerprints and conservative-field
finiteness checks. These files remain in the output directory after completion.
The existing shell tests may remove their temporary solver fixtures; the audit
does not retain every raw field or restart payload. Reproducing the fingerprints
requires rerunning the pinned revision. The scientific owner approved checksum-
and metadata-based evidence retention for the recorded v5 baseline; bulky raw
payloads are not required to remain archived.

The two full-executable profiles are:

- `strict_debug`: `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`;
- `historical_optimized`: `-Ofast -funroll-all-loops`, preserving the existing
  optimization policy for this audit, not endorsing it as the final strict policy.

Unit-test scripts compile their own executables with their existing flags. They
run once and are labeled `unit_script`; they must not be counted as testing both
production profiles.

The adapter changes only `SERIAL_FLAG` to `F` in temporary test inputs and
disables dynamic OpenMP teams. It stores the original/effective input hashes
and checks the solver's reported team and conservative fields on successful
invocations. Expected preflight rejections are retained as evidence without
being mistaken for completed solver runs. Lake-at-rest and short-excavation
scripts retain their internal one/four-thread comparisons. Gate-H and both
restart scripts run separately with one and four actual threads; canonical q
snapshots and checkpoints from successful invocations are compared across
completed pairs. Floating-point flags and physical input values are not
altered by the adapter.

Conservative output parsing accepts Fortran `D` exponents and three-digit
exponents with an omitted `E`. It does not rewrite or relax the field data.

## Characterize the N9 long cases

The selected cases are ETNA (200 s) and ETNA_LAVA (3600 s). The lava mass check
uses saved outputs only through 3180 s, the last output before boundary contact.
This does not shorten its full-field or restart checks. The
[historical v4 report](../../docs/acceptance/N9_LONG_CASES_2026_10_11.md)
preserves the originally observed CSV continuation failures. The
[v5 report](../../docs/acceptance/N9_RESTART_V5_2026_10_11.md) records the
subsequent correction and passing 8/8 full/split comparisons, while the
[versioned v5 evidence](../../docs/acceptance/results_n9_restart_v5_2026_10_11.json)
retains their source, executable, input and output hashes and records the
separate scientific owner approval for its pinned revision and scope. A fresh
audit run does not by itself extend that approval to a different revision.

Using executables and source hashes from a retained clean acceptance build:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 TESTS/ACCEPTANCE/n9_long_cases.py \
  --repo /absolute/path/to/IMEX_SfloW2D_dev \
  --build-audit /absolute/path/to/clean-acceptance-build \
  --lava-evidence /absolute/path/to/retained-lava-evidence.json \
  --reference-zip /absolute/path/to/HP_PCCU_clean_reference_cores_2026-09-26_drySafe_rebase.zip \
  --output /tmp/imex-n9-new-run
```

The lava evidence must still point to retained full-run payloads, with matching
NetCDF/checkpoint hashes. The driver independently runs ETNA and both split
trajectories in each build/team configuration. Numeric ASCII parsing uses NumPy
in place of the preprocessing script's pandas table reader; its geometry and
formatters are unchanged. Generated-data checksums remain in the evidence.
Native NetCDF is decoded sequentially through the installed C library.

When source hashes differ from the signed baseline and no current lava
revalidation is available, a fresh current-source diagnostic can be run
without reusing that approval:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 TESTS/ACCEPTANCE/n9_long_cases.py \
  --repo /absolute/path/to/IMEX_SfloW2D_dev \
  --build-audit /absolute/path/to/current-acceptance-build \
  --diagnostic-current-source \
  --output /tmp/imex-n9-current-diagnostic
```

This mode requires every current `src/` Fortran, include and C source hash to
match the supplied build audit. It records fresh full/split, thread,
checkpoint, NetCDF, raster and diagnostic observations, but does not compare
against the historical scientific baseline or approve N9; `n9_decision`
remains `open`. A scientific owner must review and sign off any extension of
the prior baseline approval separately.

To independently recharacterize retained completed full/split payloads without
rerunning or overwriting any of them, add
`--characterize-completed /absolute/path/to/retained-case-root` and choose a new
`--output` directory. `--reuse-first-halves` is a separate, explicitly mutating
recovery mode for continuing verified completed first legs in that temporary
case root; do not combine the two options. The completed characterization
checks final effective inputs, actual teams, schedules, intermediate checkpoint
time/schema, full decoded fields, rasters and final checkpoints. It reports
lost CSV history and wrong common-time values as failures and returns nonzero.
That nonzero result must not be relabeled as scientific baseline acceptance.

## Test the audit tools

The separate M12c build-system work and its current verification status are
tracked in the [M12c progress report](../../docs/acceptance/M12C_BUILD_SYSTEM_PROGRESS_2026_10_11.md).
M12c is complete and accepted for the tested Linux configurations; see its
[reproducibility manifest](../../docs/acceptance/M12C_REPRODUCIBILITY_MANIFEST_2026_10_11.json)
for the validated compiler, Docker and CI scope.

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s TESTS/ACCEPTANCE -p 'test_*.py' -v
```

These tests cover parsing, log/exit-code recording, actual-team evidence,
missing/different output fingerprints, plan links, and selected D17 module
dependency boundaries. They do not test the physical model. The dated
acceptance matrix is maintained separately in
[`docs/acceptance`](../../docs/acceptance/N7_N8_N9_AUDIT_2026_10_09.md).
