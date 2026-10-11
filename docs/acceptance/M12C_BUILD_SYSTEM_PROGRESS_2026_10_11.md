# M12c build-system progress

Status: **complete; accepted**.

This report records M12c build/test infrastructure work after the N9-approved
solver baseline `3f7329c85d2e30e33cf68348721329f7726b01c1`. The work includes
test-harness corrections and robust rejection of overflowing diagnostic CSV
values in strict-debug builds; it does not change the numerical contract.

## Changes in this increment

- Removed the `FC`, `MPFC` and `LINKER` assignments that forced GNU Fortran
  from `src/Makefile.am`. The configured `FC` now controls compilation and
  linking.
- Kept user-provided `FCFLAGS`; the default remains the historical optimized
  profile, and `--enable-strict-debug` selects runtime checks and floating-point
  traps when `FCFLAGS` is not explicitly set.
- Split OpenMP flags, NetCDF include/link flags and NetCDF/LAPACK libraries
  into overridable configure variables.
- Kept the explicit Fortran module/include dependency graph in
  `src/Makefile.am`; clean out-of-tree parallel builds pass in both profiles.
- Wired `make check` to the legacy example runner. It copies its fixtures into a
  temporary directory and checks expected outputs, so cleanup scripts cannot
  delete or rewrite files in the source tree.
- Added a GitHub Actions matrix for optimized/strict-debug builds and actual
  one/four-thread restart and HP-PCCU Gate-H validation.

## Local verification

Using GNU Fortran 15.2.0, Autoconf 2.72, Automake 1.18 and the Homebrew
NetCDF-Fortran installation:

- optimized and strict-debug out-of-tree `configure` and `make -j2`: passed;
- `make check` in both profiles: passed all nine legacy integration cases;
- optimized binary restart/runout/CSV regression: passed, including all 18
  diagnostic continuation checks and evolving-topography restart, with
  requested OpenMP settings of one and four threads;
- stochastic restart: passed in the optimized profile with actual one- and
  four-thread teams, and in strict-debug with an actual four-thread team;
- optimized HP-PCCU Gate-H: full six-case numerical matrix passed with
  requested thread counts one and four;
- strict-debug HP-PCCU Gate-H: full six-case numerical matrix passed with
  requested thread count four.

## Remote CI verification

The GitHub Actions run
[`38154444911`](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/38154444911)
completed successfully for all four GNU Fortran `solver` matrix jobs:

- optimized, one thread;
- optimized, four threads;
- strict-debug, one thread;
- strict-debug, four threads.

Each job completed the clean Linux build, `make check`, binary/runout restart
validation, stochastic restart validation, and HP-PCCU Gate-H validation.
This closes the remote-CI gate and confirms both restart profiles/thread
counts, including strict-debug CSV overflow rejection.

## Docker and compiler portability

The Docker validation run
[`38158816015`](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/38158816015)
passed the optimized and strict-debug container build/run matrix with one and
four threads on revision `fdef1356c00ee0273d3649ebdbabf07201bdcbde`. The user
also confirmed manual Docker verification for the current repository state
`9ccb349e2dd2b82ebe6b218fc655f3881b966bff`.

Intel ifx portability was verified by the successful
[`38168602274`](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/38168602274)
run on revision `9ccb349e2dd2b82ebe6b218fc655f3881b966bff`. Its Ubuntu 24.04
matrix built with Intel oneAPI ifx 2025.1 and a matching NetCDF-Fortran 4.6.2,
and passed integration, binary/runout restart, stochastic restart, and Gate-H
validation in optimized and diagnostic profiles with one and four threads.

The tested support scope is GNU Fortran and Intel ifx on Linux; it does not
claim compatibility with every Fortran compiler. The Docker base uses the
`ubuntu:24.04` tag and system package repositories rather than immutable image
and package digests. The exact tested source revision, profiles, toolchains,
and evidence references are frozen in the
[M12c reproducibility manifest](M12C_REPRODUCIBILITY_MANIFEST_2026_10_11.json).

All defined M12c gates are complete. This milestone is accepted for the tested
Linux toolchain/configuration scope above.
