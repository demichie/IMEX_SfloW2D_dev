# mass exchange

## Purpose and role in IMEX-SfloW2D

Layer-local exchange with bed and atmosphere, assembled through one shared-inventory coordinator.

## Files in this directory

| File | Responsibility |
| --- | --- |
| [`mass_exchange_parameters_2d.f90`](./mass_exchange_parameters_2d.f90) | Layer-local process configuration |
| [`mass_exchange_2d.f90`](./mass_exchange_2d.f90) | Combined coordinator and shared resource bounds |
| [`mass_exchange_erosion_2d.f90`](./mass_exchange_erosion_2d.f90) | Substrate erosion |
| [`mass_exchange_deposition_2d.f90`](./mass_exchange_deposition_2d.f90) | Sediment deposition and settling kernel |
| [`mass_exchange_entrainment_2d.f90`](./mass_exchange_entrainment_2d.f90) | Ambient-air entrainment |
| [`mass_exchange_carrier_liquid_loss_2d.f90`](./mass_exchange_carrier_liquid_loss_2d.f90) | Carrier withdrawal |
| [`mass_exchange_vaporization_2d.f90`](./mass_exchange_vaporization_2d.f90) | Vaporization and heat transfer |
| [`mass_exchange_lava_deposition_2d.f90`](./mass_exchange_lava_deposition_2d.f90) | Lava deposition and bed growth |

## Physical equations

For a layer, conservative source contributions are additive: \(S_M=\sum_p R_M^{(p)}\) and \(S_{Q_T}=\sum_p R_{Q_T}^{(p)}\), with \(Q_T=M c_{p,m}T\). Positive mass rate enters the moving layer.

**Erosion:** solid volume influx \(E_i\) yields \(R_{s,i}=\rho_{s,i}E_i\); carrier influx \(E_c\) yields \(R_c=\rho_{c,\mathrm{bed}}E_c\), with thermal content evaluated at substrate temperature. The empirical rate depends on erosion coefficient, speed, depth and erodible inventory.

**Deposition:** accepted settling volume flux \(D_i\) gives \(R_{s,i}=-\rho_{s,i}D_i\); carrier deposition is determined by bed porosity. The settling kernel is shared with interlayer settling.

**Ambient-air entrainment:** \(R_M=\rho_{a,\mathrm{amb}} E_a\), \(R_{Q_T}=T_{\mathrm{amb}} c_{p,a}R_M\). **Liquid-carrier loss:** for \(L=\max(0,\mathrm{loss\_rate})\), \(R_M=-\rho_c L\) and \(R_{Q_T}=T c_{p,c}R_M\). **Lava deposition:** \(R_M=-\rho_l D_l\), \(R_{Q_T}=c_{p,l}T R_M\), proposed bed growth \(\dot B=D_l\). **Vaporization** uses a temperature-dependent liquid-to-vapor rate and corresponding gas thermal contribution; its thresholds and precise formula reside in the vaporization module.

## Architectural and numerical contracts

The coordinator enforces competing-withdrawal limits and packing on one common candidate. Laws do not commit simulation state. These process laws currently support explicit evaluation; unimplemented methods fail explicitly.

Numerical thresholds, dry-state branches, regularization, and available methods are governed by the linked implementation. The equations above describe implemented physical terms; they are not a replacement for checking each kernel. See the [parent physics guide](../README.md) and [Fortran style](../../../docs/FORTRAN_STYLE.md).
