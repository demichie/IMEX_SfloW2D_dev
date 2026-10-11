# M12c build-system progress

Status: **in progress; not accepted**.

This report records the initial M12c build/test infrastructure work after the
N9-approved solver baseline `3f7329c85d2e30e33cf68348721329f7726b01c1`.
It changes build and test orchestration only; no solver source or numerical
contract is changed.

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

The new GitHub Actions workflow has not yet run remotely. Docker is not
installed in the local environment, and a non-GNU Fortran compiler and
strict-debug/one-thread Gate-H run have not been tested locally. Consequently
this increment is **not** evidence of complete M12c acceptance or compiler
portability.

## Remaining M12c gates

1. Exercise strict-debug Gate-H with one thread, and decide whether the
   configured restart tests need a separate strict-debug/one-thread run beyond
   the previously accepted v5 evidence.
2. Observe a successful GitHub Actions run on a clean Linux checkout for both
   profiles and thread counts.
3. Verify Docker build/run, or record and explicitly accept a concrete
   environment limitation.
4. Confirm compiler portability and document the supported toolchains and
   dependency versions.
5. Freeze a reproducibility manifest for the accepted build configuration,
   then update the roadmap status only after all acceptance criteria pass.
