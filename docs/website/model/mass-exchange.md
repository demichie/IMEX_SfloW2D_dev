# Erosion, deposition and environmental exchange

A flow can exchange material with the substrate or the surrounding
atmosphere. These processes alter the moving mixture's mass,
composition and thermal content. Some also modify the substrate
inventory or terrain elevation.

## Substrate erosion

Erosion transfers material from the erodible substrate into the moving
flow. Each eroded solid class contributes mass at a rate proportional
to its solid density and erosional volume flux. The empirical rate
depends on configured coefficients, flow conditions and available
erodible material.

Erosion is limited by the substrate inventory: the moving layer
cannot remove material that is not present.

## Deposition and settling

Deposition transfers material out of the moving flow into the bed.
The model evaluates settling for configured solid classes, with
limits imposed by the donor inventory and admissible mixture
composition. Carrier deposition is associated with the prescribed
bed porosity.

Deposition to the **static substrate** differs from transfer of
solids between two moving layers. In a two-layer simulation, settling
from the upper layer to the lower layer is an
[interlayer exchange](interlayer-processes.md), not direct bed
deposition.

## Other exchanges

| Process | Effect |
| --- | --- |
| Ambient-air entrainment | Adds atmospheric carrier mass and its thermal content |
| Liquid-carrier loss | Removes carrier mass and corresponding thermal content |
| Vaporization | Converts liquid to vapor under the implemented temperature-dependent law |
| Lava deposition | Removes mobile lava and contributes to bed growth |

The physical availability of these processes depends on the material
and selected input options. A moving layer's contact with the bed or
atmosphere can also depend on whether another layer is locally present.

## Conservation and terrain change

For each source process, mass and heat changes must be interpreted
together with substrate inventories, external inflows and outflows,
and transported constituents. A decrease in moving-layer mass may
represent deposition rather than a numerical loss.

The simulator combines proposed exchange rates with shared
availability and composition limits before accepting the source
increment. [Numerical methods](numerics.md) describes source-group
ordering; [examples](../examples/index.md) demonstrate the settings.
