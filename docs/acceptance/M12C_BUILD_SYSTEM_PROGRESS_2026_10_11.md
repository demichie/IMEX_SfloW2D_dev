# M12c build-system progress

Status: **in progress; not accepted**.

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
[`38151328516`](https://github.com/demichie/IMEX_SfloW2D_dev/actions/runs/38151328516)
completed successfully for all four `solver` matrix jobs:

- optimized, one thread;
- optimized, four threads;
- strict-debug, one thread;
- strict-debug, four threads.

Each job completed the clean Linux build, `make check`, binary/runout restart
validation, stochastic restart validation, and HP-PCCU Gate-H validation.
This closes the remote-CI gate and confirms both restart profiles/thread
counts, including strict-debug CSV overflow rejection.

Docker is not installed in the local environment, and a non-GNU Fortran
compiler has not been tested. Therefore M12c remains **in progress; not
accepted**, and these CI results do not establish compiler portability.

## Remaining M12c gates

1. Verify Docker build/run, or record and explicitly accept a concrete
   environment limitation.
2. Confirm compiler portability and document the supported toolchains and
   dependency versions.
3. Freeze a reproducibility manifest for the accepted build configuration,
   then update the roadmap status only after all acceptance criteria pass.
