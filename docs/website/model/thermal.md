# Thermal processes

Temperature is part of the physical state and can influence both
material properties and source terms. IMEX-SfloW2D evolves
depth-integrated **thermal content**

$$
Q_T=M c_{p,m} T,
$$

where $M$ is mixture mass per unit area and $c_{p,m}$ its specific
heat capacity. This quantity is not total kinetic-plus-potential energy.

## Radiative and convective cooling

The configured local cooling source follows the implemented form

$$
S_{Q_T}=-r(T^4-T_{\mathrm{amb}}^4)
-c(T-T_{\mathrm{amb}}),
$$

with configured radiative and convective coefficients and ambient
temperature. Supported numerical treatments include explicit,
implicit and frozen-state analytic updates.

## Temperature-dependent rheology

The influence of temperature is not restricted to a cooling source.
The **lava rheology** incorporates temperature-dependent viscosity
and yield resistance. The **lahar rheology** uses temperature and
concentration when evaluating fluid resistance.

A thermal calculation and a temperature-dependent rheology are
different laws: activating one does not automatically determine the
parameterization of the other. See [rheological models](rheology.md).

## Thermal effects of material exchange

Erosion, deposition, carrier losses, atmospheric entrainment,
injections and vaporization may transport or modify thermal content
with the corresponding mass changes. The physical source models
apply their documented donor, substrate or ambient material
assumptions.

For the complete exchange context, see
[erosion and deposition](mass-exchange.md). Example settings are
available in the [catalogue](../examples/index.md).
