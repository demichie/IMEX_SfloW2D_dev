# Interlayer mass and momentum exchange

Two-layer simulations may transfer solids and carrier material
between L1 and L2. They may also exchange momentum through
interlayer drag. These processes are distinct from the
[coupled hydrostatic pressure](coupling.md), which acts within
the spatial transport calculation.

## Particle settling and elutriation

- **Settling** moves solid particles from L2 to L1.
- **Elutriation** moves solid particles from L1 to L2.

Proposed rates are evaluated for each configured solid class and
then limited according to the donor inventory, both layer
compositions and the accepted packing restrictions. Material
received in one source update is not immediately available for
a simultaneous counter-transfer.

L2 settling places material in the lower mobile layer; it is
**not** direct deposition into the static substrate. Basal
deposition is a separate process involving L1.

### Excess pore pressure in the dense layer

Only L1 carries excess pore pressure. The material received from dilute L2
has zero excess pore pressure, so settling adds mixture mass without adding
to the conservative L1 pressure load $Q_p=M_1p_e$. Existing L1 pressure is
therefore diluted. Conversely, elutriation removes $p_{e,1}$ times the
outgoing L1 mixture mass (solids and carrier); it does not transfer a
pore-pressure variable into L2. Simultaneous opposed streams are treated
by their gross fluxes, not the net mass change. No automatic lithostatic
pore-pressure seed is added to newly formed L1 deposits.

## Carrier transported with lifted solids

Elutriated particles may be accompanied by carrier gas.
The configuration parameter
`ELUTRIATION_SOLID_VOLUME_FRACTION` specifies the
fraction of solids **in the transferred stream**, not the
desired concentration of the resulting L2.

When a positive value $\alpha$ is prescribed, the requested
carrier volume accompanying lifted solids is

$$
\Delta V_c=\frac{1-\alpha}{\alpha}\,\Delta V_s.
$$

With value `-1`, the fraction is derived from the
current L1 donor composition. This is a *mean donor mixture*
assumption; it need not describe the exact material lifted
at the interface.

Two carrier policies are supported:

| Policy | Behaviour when carrier is insufficient |
| --- | --- |
| `LIMIT_SOLIDS` | Reduce lifted solids so the requested solid/carrier ratio can be maintained |
| `ALLOW_CONCENTRATED` | Permit a more concentrated lifted stream while respecting the donor inventory and packing limits |

Atmospheric air entrainment into L2 is a separate **external**
source that can dilute the upper mixture later in the evolution.

When optional Darcy gas loss is enabled in L1, the gas leaving L1 is
routed into L2 in a two-layer run, including when L2 was previously
dry. The transfer carries the donor's horizontal momentum, sensible
thermal content and gas-species composition. Darcy relaxation without
gas loss only reduces L1 excess pore pressure. In a single-layer run,
the existing atmospheric gas-loss behavior is preserved.

## Conservation and thermal content

A transfer between L1 and L2 produces opposite source terms
for the corresponding conservative species masses, total
mixture mass and transported momentum/thermal content.
Therefore the **combined two-layer mass** is unchanged by
interlayer exchange alone.

When carrier densities differ, equal transferred *mass* does
not imply equal carrier *volume* changes in the two layers.
Thermal transport uses the properties of donor material.
Combined layer budgets must additionally include external
entrainment, deposition, boundary fluxes and other enabled
sources.

## Interlayer drag

The implemented Chézy-type drag relaxes relative layer
velocity at fixed layer masses. It conserves the sum of
horizontal momenta and decreases resolved relative kinetic
energy.

!!! note "Energy approximation"
    Mechanical energy dissipated by interlayer drag is **not**
    converted into sensible heat in the current model.
    Combined mechanical-plus-thermal energy is therefore not
    conserved by this drag-only operation. This is an accepted
    physical approximation, not an output error.

## Source ordering

Physical transfer and drag are applied in configured Lie
source-group order, following the spatial finite-volume step.
Processes in the same group have a defined sequence.

For example, deposition can change the L1 inventory before
an elutriation process evaluates the lifted mixture composition.
Thus a different source ordering can lead to a different finite
timestep evolution, even with identical model coefficients.

## Limits and calibration

The numerical implementation enforces paired source budgets
and admissibility checks, but those checks do not establish that
an empirical transfer law is calibrated for a particular flow.
The gas co-lift assumption, simultaneous settling/elutriation,
compressible carrier thermodynamics and sensitivity to source
ordering require care in quantitative studies.

For the resulting changes in pressure-wave characteristics,
read [Characteristic limits](hyperbolicity.md). Case-specific
settings are documented in the [examples](../examples/index.md).
