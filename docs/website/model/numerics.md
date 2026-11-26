# Numerical methods

This page explains how IMEX-SfloW2D advances a simulation and how the
timestep is selected. The same spatial and temporal framework is used
with one or two layers.

## Finite-volume discretization

The computational area is represented by a Cartesian grid over the
specified terrain. In each cell, the model stores depth-integrated
conservative variables. The method reconstructs physical states on
cell faces and computes numerical fluxes between neighbouring cells.

The spatial discretization uses the **hydrostatic-path,
path-conservative central-upwind (HP-PCCU)** formulation. Its
reconstruction, hydrostatic pressure treatment and nonconservative
path contributions account for the influence of topography.

For two layers, the spatial operator also includes the
[coupled pressure terms](coupling.md). These terms are part of the
spatial equations rather than separate physical source operations.

Only cells and interfaces needed for the active flow and its
reconstruction neighbourhood are processed.

## CFL timestep selection

The `NUMERIC_PARAMETERS` group supports two CFL policies based on
nonnegative directional propagation-speed bounds
$\lambda_x$ and $\lambda_y$, cell sizes $\Delta x$, $\Delta y$,
and the configured CFL factor $C$.

For `CFL_POLICY='directional_min'` (the default):

$$
\Delta t_{\mathrm{CFL}}=
\min\left(\Delta t_{\max},
\frac{C\,\Delta x}{\lambda_x},
\frac{C\,\Delta y}{\lambda_y}\right).
$$

For `CFL_POLICY='combined'`:

$$
\Delta t_{\mathrm{CFL}}=
\min\left(\Delta t_{\max},
\frac{C}{\lambda_x/\Delta x+\lambda_y/\Delta y}\right).
$$

Zero-speed directions are omitted. Other timestep restrictions may
apply, including output scheduling and admissibility requirements.

!!! warning "CFL is not a hyperbolicity test"
    For two layers, the implemented total-depth speed estimate bounds
    *real* gravitational characteristics under specified conditions.
    It does not detect complex characteristics. Reducing the timestep
    cannot restore hyperbolicity to a physically non-hyperbolic state.
    See [characteristic limits](hyperbolicity.md).

## Spatial update followed by physical sources

Within a timestep, the current method applies the finite-volume
spatial update and then **ordered Lie-split source groups**.

The first source group, when configured as group 0, evaluates its
rates using the beginning-of-step physical state but applies the
limited source increment to the state after the spatial update.
Subsequent source groups act on the evolving cell-local candidate
in their configured order.

Physical laws may support explicit, semi-implicit, implicit
(Newton) or analytical integration, according to the individual
model. Unsupported combinations are rejected rather than silently
replaced by another method.

Only an admissible final state is accepted. If a candidate becomes
invalid and cannot be accepted under the implemented rules, the
simulation stops with diagnostics instead of silently clipping
the state or retrying with a modified physical model.

## What this means for a user

- **Select the CFL policy and numerical parameters** from an example
  suited to the problem; do not assume the same timestep is suitable
  for all material laws.
- **Choose physical sources deliberately.** Different source groups
  are not generally interchangeable because their order affects the
  evolution over a finite timestep.
- **Check convergence and budgets.** Compare simulations with refined
  grids and timesteps when quantitative predictions matter.
- **Interpret two-layer conditions separately.** Refer to the
  [two-layer formulation](two-layer.md), [interlayer processes](interlayer-processes.md)
  and [characteristic limits](hyperbolicity.md).

The optional [Fortran API](https://demichie.github.io/IMEX_SfloW2D_dev/doxygen/)
documents the underlying implementation for readers who need source details.
