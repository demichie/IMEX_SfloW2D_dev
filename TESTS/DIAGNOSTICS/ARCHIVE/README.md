# Historical numerical diagnostics

These files are preserved as standalone historical references rather
than active inputs to the production solver.

- `kt_v2_original_1d.py`: earlier Python reduced-model diagnostic.
- `kt_ritter_1d.f90`: standalone Fortran Ritter-flow diagnostic,
  formerly distributed from `src/python_diagnostics/`.

Their source bodies were not modified by the relocation. When comparing
past results, use the appropriate archived version and its historical
compiler/Python environment. Neither file is linked into the
IMEX-SfloW2D executable.
