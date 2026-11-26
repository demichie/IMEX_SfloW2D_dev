# Validation and scientific scope

## Recommended checks for a source revision

From GitHub **Actions**, choose the intended branch and run the following workflows explicitly:

1. **Static source and documentation audit:** Fortran style, Doxygen contracts, direct module dependencies and static invariants.
2. **Cell physics and numerical regressions:** state conversions, physical source laws, REAL/COMPLEX kernels and spatial operator checks.
3. **Solver integration and OpenMP validation:** optimized and strict/debug builds, representative solver cases, restart and OpenMP checks.

The build-only workflow checks compilation without substituting for the numerical regressions. The OpenMP benchmark is a performance diagnostic, not a correctness certification.

To reproduce a published result, record the Git commit SHA, compiler/build profile, `IMEX_SfloW2D.inp` and input raster checksums alongside output diagnostics.

## Acceptance limitations

The integrated two-layer model can enter non-hyperbolic density/shear regimes. Some historical N7 acceptance criteria are still marked pending revalidation against observables of the current solver. Neither an example run nor passing CI means that every model configuration is physically calibrated.

- [Two-layer characteristic analysis](../model/two-layer.md)
- [Historical acceptance evidence](https://github.com/demichie/IMEX_SfloW2D_dev/tree/master/docs/acceptance)
- [Current acceptance inventory](https://github.com/demichie/IMEX_SfloW2D_dev/blob/master/TESTS/ACCEPTANCE/acceptance_plan.json)
- [Current workflow catalogue](https://github.com/demichie/IMEX_SfloW2D_dev/blob/master/.github/README.md)
