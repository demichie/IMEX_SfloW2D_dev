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

The audit requires Git, tar, Autotools, make, gfortran, NetCDF Fortran, LAPACK,
Python 3 and NumPy. Compiler portability and CI integration belong to M12c;
this driver does not remove the existing hard-coded compiler settings.

It returns nonzero for failed builds, failed test runs, failed available
one/four-thread comparisons or a supplied reference-package mismatch. The JSON
manifest retains failures and sets `n9_accepted` to false regardless of the
command's return code. The per-command timeout defaults to 900 seconds and can
be adjusted with `--timeout`; a timeout is a failure, not a skipped test.

## Evidence and interpretation

`manifest.json` contains the solver Git SHA, tool probes, compiler flags,
executable and source checksums, test results and coverage criteria. Each test
has a log; each solver invocation has its own input snapshots, full solver log,
actual thread count, output SHA-256 fingerprints and conservative-field
finiteness checks. These files remain in the output directory after completion.
The existing shell tests may remove their temporary solver fixtures; the audit
does not retain every raw field or restart payload. Reproducing the fingerprints
requires rerunning the pinned revision. Archived approved golden outputs remain
a separate N9 requirement.

The two full-executable profiles are:

- `strict_debug`: `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`;
- `historical_optimized`: `-Ofast -funroll-all-loops`, preserving the existing
  optimization policy for this audit, not endorsing it as the final strict policy.

Unit-test scripts compile their own executables with their existing flags. They
run once and are labeled `unit_script`; they must not be counted as testing both
production profiles.

The adapter changes only `SERIAL_FLAG` to `F` in temporary test inputs and
disables dynamic OpenMP teams. It stores the original/effective input hashes and
checks the solver's reported team. Lake-at-rest and short-excavation scripts
retain their internal one/four-thread comparisons. Gate-H and both restart
scripts run separately with one and four actual threads; canonical q snapshots
and checkpoints are compared across completed pairs. Floating-point flags and
physical input values are not altered by the adapter.

Conservative output parsing accepts Fortran `D` exponents and three-digit
exponents with an omitted `E`. It does not rewrite or relax the field data.

## Test the audit tools

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s TESTS/ACCEPTANCE -p 'test_*.py' -v
```

These tests cover parsing, log/exit-code recording, actual-team evidence,
missing/different output fingerprints and plan links. They do not test the
physical model. The dated acceptance matrix is maintained separately in
[`docs/acceptance`](../../docs/acceptance/N7_N8_N9_AUDIT_2026_10_09.md).
