# interlayer mass exchange

## Purpose and role in IMEX-SfloW2D

Paired conservative L1/L2 solid and carrier transfer without terrain or atmospheric exchange.

## Files in this directory

| File | Responsibility |
| --- | --- |
| [`interlayer_mass_exchange_parameters_2d.f90`](./interlayer_mass_exchange_parameters_2d.f90) | Pair-wide settling, elutriation and co-lift configuration |
| [`interlayer_settling_2d.f90`](./interlayer_settling_2d.f90) | L2-to-L1 settling sharing the basal deposition kernel |
| [`interlayer_elutriation_2d.f90`](./interlayer_elutriation_2d.f90) | L1-to-L2 shear/Rouse filtered elutriation |
| [`interlayer_mass_exchange_2d.f90`](./interlayer_mass_exchange_2d.f90) | Concurrent donor, carrier and packing limiter; conservative paired update |

## Physical equations

For accepted downward solid volume fluxes \(S_i\) (L2 to L1) and upward elutriation \(E_i\) (L1 to L2), both in m/s:

\[
\dot V_{s,i,1}=S_i-E_i,\qquad \dot V_{s,i,2}=-S_i+E_i.
\]

Mass transfers use species density \(\rho_{s,i}\), and carrier co-lift and downward carrier exchange enter the two layers with opposite signs. Each pairwise total-mass/species/momentum/thermal transfer conserves the sum of both layers. Upward LINEAR_USTAR elutriation includes the factor \(\epsilon\alpha_{s,i,1}u_*\phi A\), filtered by shear/Rouse and availability factors; SHEAR_POWER is a distinct branch. Downward settling delegates to the basal deposition-volume kernel using L2 properties. The coordinator limits both gross donors, carrier supply and final packing *simultaneously*: material received in the substep is not immediately available for withdrawal.

## Architectural and numerical contracts

The pair-wide configuration is owned once; no duplicate L1/L2 settings. LIMIT_SOLIDS and ALLOW_CONCENTRATED produce different co-lift constraints. Unsupported source treatments reject use.

Numerical thresholds, dry-state branches, regularization, and available methods are governed by the linked implementation. The equations above describe implemented physical terms; they are not a replacement for checking each kernel. See the [parent physics guide](../README.md) and [Fortran style](../../../docs/FORTRAN_STYLE.md).
