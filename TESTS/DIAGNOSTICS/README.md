# Numerical diagnostics

This directory contains **independent numerical diagnostics and reference
cases**, not modules linked into the production IMEX-SfloW2D executable.
It was moved intact from `src/python_diagnostics/` during repository cleanup.

- [`HP_PCCU_1D/`](HP_PCCU_1D/README.md) contains reference Python
  implementations and comparison cases for one-dimensional numerical
  schemes, reconstruction and wet/dry transport.
- [`ARCHIVE/`](ARCHIVE/) retains earlier versions and the standalone
  Fortran Ritter diagnostic for historical reproducibility. Files are
  preserved without modifications to the original algorithms.

The Python diagnostics are exploratory/reference tools, not an alternative
production IMEX solver. Their output should not be interpreted as
scientific acceptance of the Fortran implementation without direct
comparison and verification.

The diagnostic files, including the archived Fortran example, are listed
in the root `Makefile.am` `EXTRA_DIST`. They are not linked into the
main Fortran binary or run automatically in the default test suite.
