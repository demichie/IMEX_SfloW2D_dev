# Spatial and temporal integration

## Finite-volume spatial operator

The current solver reconstructs cell states, evaluates high-resolution numerical fluxes and incorporates terrain geometry through the hydrostatic/HP-PCCU path. Single- and two-layer calculations use the same spatial framework and one shared active-cell domain.

## CFL policies

The `NUMERIC_PARAMETERS` namelist selects `CFL_POLICY='directional_min'` (default) or `CFL_POLICY='combined'`. For nonnegative directional speed bounds `lambda_x` and `lambda_y`:

- **directional_min:** `min(MAX_DT, CFL*dx/lambda_x, CFL*dy/lambda_y)`, omitting zero-speed directions.
- **combined:** `min(MAX_DT, CFL/(lambda_x/dx + lambda_y/dy))`, using only active directions.

In the two-layer model the total-depth characteristic envelope is **conditional on real gravitational characteristics**. It does not detect complex eigenvalues and cannot establish general hyperbolicity.

## Time integration and physical sources

Accepted conservative updates use one ordered Lie source executor. The source registry validates supported methods and per-layer scopes at initialization; physical modules evaluate their own laws. The layer-local state is not duplicated in a second permanent representation.

Physical-law modules keep uniform `EXPLICIT`, `SEMI_IMPLICIT`, `IMPLICIT` and `ANALYTIC` interfaces as appropriate, including intentional zero contributions or explicit unsupported-treatment errors.

## Scientific and implementation references

- [Two-layer characteristic limitations](two-layer.md)
- [Lie timestep and failure contract](https://github.com/demichie/IMEX_SfloW2D_dev/blob/master/src/numerics/README.md)
- [Canonical architecture](https://github.com/demichie/IMEX_SfloW2D_dev/blob/master/docs/TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md)
- [Numerical and OpenMP verification](../development/validation.md)
