# Physical and numerical model

IMEX-SfloW2D simulates two-dimensional, depth-averaged geophysical flows over
a terrain surface. It combines finite-volume transport with a configurable
mixture description, constitutive rheology, and physical source processes.
Simulations can use one moving layer or two interacting layers.

The choice of rheology, exchanged materials, thermal processes and numerical
source treatment determines the physical model of a particular simulation.
**Two-layer coupling is an optional extension, not a prerequisite for the
other physical processes.**

## State variables and governing balances

For each flowing layer, the solver advances conservative quantities including
mixture mass, two horizontal momenta, thermal content and transported
constituent masses. With mixture density $\rho_m$, thickness $h$, horizontal
velocities $u,v$, mixture heat capacity $c_{p,m}$ and temperature $T$, the
leading conserved quantities per unit horizontal area are

$$
M=\rho_m h,\qquad Mu,\qquad Mv,\qquad Q_T=M c_{p,m}T.
$$

The thermal quantity is **thermal content**, not total mechanical energy.
Optional variables describe pore pressure, stochastic forcing and configured
additional constituents. Constitutive laws use the corresponding physical
states reconstructed from these conservative variables.

## Choose the physical processes

| Topic | What the model represents |
| --- | --- |
| [Rheological models](rheology.md) | Friction, yield stress, velocity-dependent and temperature/composition-dependent resistance |
| [Erosion, deposition and entrainment](mass-exchange.md) | Exchanges with terrain, atmosphere and carrier material; bed growth and erosion |
| [Thermal processes](thermal.md) | Radiative and convective cooling, temperature-dependent material behaviour, vaporization |
| [Pore pressure and other sources](additional-processes.md) | Darcy pore pressure, stochastic variations and external injections |
| [Lava and lahar applications](applications.md) | How multiple physical laws are combined for representative volcanic-flow settings |
| [Numerical methods](numerics.md) | Reconstruction, fluxes, source integration, CFL limits and admissibility |
| [Two-layer flows](two-layer.md) | Optional interaction of two moving mixtures through pressure, exchange and drag |

These processes are configurable, not necessarily active in every run.
Some laws require particular material assumptions or support only specific
source-integration treatments.

## One layer or two?

A **single-layer** run evolves one moving mixture, with the configured
rheology, heat transfer, sources and exchange with its environment.

A **two-layer** run evolves two moving mixtures over the same terrain,
with independent layer properties and additional pressure coupling and
interlayer interactions. It requires separate attention to density ordering,
relative velocity and the validity of the coupled characteristic speeds.
Read the [two-layer formulation](two-layer.md) and
[characteristic limits](hyperbolicity.md) for these assumptions.

## From equations to a simulation

The [example catalogue](../examples/index.md) provides complete cases with
specific parameter sets. Use it together with the physical-process pages to
identify which laws are active; a working example does not, by itself,
establish calibration or predictive validity for a different event.
