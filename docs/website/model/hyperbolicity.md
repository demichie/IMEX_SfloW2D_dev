# Characteristic limits of the two-layer model

The pressure coupling between two mobile layers introduces wave
modes that are absent from independent single-layer simulations.
Depending on the local densities and relative velocities,
some gravitational characteristic speeds can become complex.

This page describes the restrictions relevant when interpreting
a two-layer result.

## The state considered by the analysis

Consider two **wet** layers over locally planar terrain, with
positive depth $h_1,h_2$, positive reduced gravities $g'_1,g'_2$,
constant terrain metric $G>0$, and locally frozen density and
composition. Choose one horizontal direction, with layer
velocities $U_1$ and $U_2$ along that direction.

Let

$$
a=Gg'_1 h_1,\qquad d=Gg'_2 h_2,\qquad
k=G^2\frac{\Gamma_2}{\rho_1}g'_2 h_1h_2.
$$

The four gravitational characteristic speeds satisfy

$$
\big((\lambda-U_1)^2-a\big)
\big((\lambda-U_2)^2-d\big)-k=0.
$$

This describes a local gravitational subsystem, **not** a full
eigenanalysis of wet/dry fronts, varying composition or the
discrete numerical flux.

## Density ordering when the two layers move together

For $U_1=U_2$ and both layers wet:

| Density condition | Consequence for the gravitational block |
| --- | --- |
| $\rho_1>\rho_2>\rho_a$ | Four distinct, real gravitational speeds at zero slip |
| $\rho_1=\rho_2>\rho_a$ | Repeated real speed with a defective eigenbasis; not strongly hyperbolic |
| $\rho_2>\rho_1>\rho_a$ | A complex-conjugate pair already at zero slip |

Thus an upper layer denser than the lower one violates this
gravitational hyperbolicity condition even if their velocities
are identical.

## Relative velocity and shear

Stable density ordering alone is not sufficient if the layers
move at different speeds. The transition between real and
complex roots depends on the depths, densities, terrain metric
and the relative normal velocity.

At a threshold, two roots of the characteristic polynomial
coincide. The exact condition is obtained by solving the
polynomial and its derivative together; it is **not** one
universal shear threshold.

In two dimensions, the relative velocity must be considered in
every horizontal direction. Real roots along one selected
direction do not establish hyperbolicity along all others.

## Why the CFL condition is not enough

The current method evaluates propagation-speed bounds for each
present layer using the total depth $H=h_1+h_2$. Under stable
density ordering and positive reduced gravity, this envelope
contains **the real gravitational speeds** of the locally
homogeneous wet/wet block.

However:

- It does **not** test whether the quartic has four real roots.
- It does **not** certify strong hyperbolicity.
- It does **not** cover every dry/wet, density-inverted or
  variable-composition configuration.
- A smaller timestep does not change complex eigenvalues of the
  continuous pressure system.

## Interpreting simulation output

For a quantitative two-layer study, assess local layer depths,
mixture densities and relative velocities from the simulated
state and, where necessary, evaluate the characteristic roots.
Do not infer physical validity from the absence of a solver error.

Changes in settling, elutriation or atmospheric entrainment can
alter the density contrast and shear, so the characteristic
conditions can evolve during a run.

See [Two-layer formulation](two-layer.md) and
[Interlayer processes](interlayer-processes.md) for the relevant
physical mechanisms.
