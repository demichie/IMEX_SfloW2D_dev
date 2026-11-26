# Coupled pressure and terrain geometry

This page describes the pressure interaction used in the two-layer
spatial equations. It is useful when comparing the model with other
depth-averaged two-layer formulations.

## Free surfaces and terrain metric

Let $B(x,y)$ be the terrain elevation, $h_1$ and $h_2$ the two
layer depths, and $\rho_1$, $\rho_2$ their local mixture densities.
The layer interface and top free surface are

$$
\eta_1=B+h_1,\qquad \eta_2=B+h_1+h_2.
$$

The terrain metric in the implemented pressure terms is

$$
G=\frac{1}{1+B_x^2+B_y^2}.
$$

For atmospheric density $\rho_a$ and gravitational acceleration $g$,
define the reduced gravity and pressure coefficient

$$
g'_i=g\left(1-\frac{\rho_a}{\rho_i}\right),\qquad
\Gamma_i=\rho_i g'_i .
$$

## Pressure coupling between layers

Define the upper-layer pressure load

$$
L=\Gamma_2 h_2.
$$

For each horizontal direction $j\in\{x,y\}$, the implemented
interlayer contributions are

$$
S_{\mathrm{inter},1,j}=-G h_1\,\partial_j L,
$$

$$
S_{\mathrm{inter},2,j}=-L\,\partial_j(Gh_1).
$$

The remaining upper-layer hydrostatic pressure contribution is
already part of the spatial pressure path; it must not be added
again. The pressure coupling is included in the finite-volume
spatial operator, **not** in the Lie groups for physical transfer.

Terrain curvature and slope are calculated from $B$, not
from the evolving interface $B+h_1$.

## What changes compared with a single layer?

The pressure gradients of L1 depend on the upper-layer pressure
load. Conversely, L2 responds to the lower-layer depth through
the height of the layer interface. A two-layer simulation is
therefore not equivalent to adding two independently calculated
single-layer solutions.

The numerical method reconstructs both layers on the shared
grid and evaluates the coupled pressure contributions within
the same spatial update.

## Validity of the characteristic approximation

The solver uses a total-depth characteristic-speed envelope for
the coupled spatial timestep. It encloses real gravitational
characteristics only under stated assumptions; it does not
guarantee that all eigenvalues are real.

The precise density and relative-velocity conditions are explained
in [Characteristic limits](hyperbolicity.md). For the physical
transfers that occur after the spatial update, see
[Interlayer processes](interlayer-processes.md).
