# Rheological models

A rheology determines the resistance experienced by a moving mixture. In
IMEX-SfloW2D, rheological momentum sources can depend on velocity, depth,
composition, temperature, local geometry and, for supported models, pore
pressure. The available laws are alternatives with different physical
assumptions, rather than one universally applicable friction formula.

## Available families

| Model | Main behaviour |
| --- | --- |
| **Voellmy** | Combines basal Coulomb resistance with a velocity-dependent turbulent contribution |
| **Plastic** | Yield-type opposition to motion |
| **Lava** | Temperature-dependent viscous drag and a yield contribution |
| **Lahar** | Fluid viscosity and yield resistance depending on temperature and material concentration, with turbulent resistance |
| **Empirical** | Depth- and speed-dependent drag |
| **Quadratic** | Quadratic velocity-dependent resistance |
| **Froude-dependent** | Friction coefficient varies with the Froude number |
| **Lucas** | Piecewise friction response governed by speed and depth |
| **Granular µ(I)** | Friction dependent on a regularized inertial number and effective normal stress |
| **µ(I) + Voellmy** | Combines granular basal resistance with turbulent drag |

## Examples of constitutive dependence

The Froude-dependent law uses a coefficient of the form

$$
\mu=\mu_\infty+(\mu_0-\mu_\infty)
\exp(-\mathrm{Fr}/\mathrm{Fr}_0),\qquad
\mathrm{Fr}=\frac{|\boldsymbol u|}{\sqrt{g'h}}.
$$

The lava law has a viscous term that depends exponentially on
temperature. The lahar law uses temperature- and
concentration-dependent fluid viscosity and yield contributions.
The granular µ(I) law depends on inertial number and effective normal
stress; pore pressure may therefore modify the resistance in
configurations where that coupling is enabled.

These are **descriptions of dependence**, not complete replacements for
the implemented dry-state branches, regularizations and parameter
definitions.

## Selecting a law

Choose the formulation on physical grounds: whether the mixture behaves
primarily as a granular mass, a yield-stress material, a temperature-sensitive
lava or a concentration-sensitive lahar. Also check whether the selected
law is compatible with the intended explicit, semi-implicit, implicit
or analytic source treatment. Availability differs by law and is
validated by the solver.

The [lava and lahar page](applications.md) shows how rheology fits
alongside other processes; the [example catalogue](../examples/index.md)
provides input settings and runnable cases.
