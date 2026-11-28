# pore pressure

## Purpose and role in IMEX-SfloW2D

Darcy excess pore pressure in dense L1 and optional carrier-gas evolution.

## Files in this directory

| File | Responsibility |
| --- | --- |
| [`pore_pressure_parameters_2d.f90`](./pore_pressure_parameters_2d.f90) | Per-layer configuration and derived coefficient cache |
| [`pore_pressure_darcy_2d.f90`](./pore_pressure_darcy_2d.f90) | Darcy pressure/gas physics and analytic substep |
| [`pore_pressure_2d.f90`](./pore_pressure_2d.f90) | Source treatment dispatcher |

## Physical equations

The pore-pressure term has its own conservative component and can couple to interstitial carrier gas. The implemented *frozen-state analytic* solve returns a complete average conservative source:

$$
\overline{\boldsymbol S}_{\mathrm{Darcy}}=
\frac{\boldsymbol q_{\mathrm{after}}-\boldsymbol q_{\mathrm{before}}}{\Delta t},\qquad \Delta t>0.
$$

This rate represents the integrated Darcy pressure/gas update and is not an explicit Euler discretization. Pore pressure can reduce effective basal normal stress in the rheology law; the stress subtraction itself is owned by rheology. Detailed permeability and thermodynamic closures are implemented in the Darcy module.

In two-layer simulations L2 is dilute and cannot activate an excess
pore-pressure equation. When L1's Darcy step removes carrier gas, the Lie
executor credits the same mass, horizontal momentum, thermal energy and
gas-species masses to L2. The transfer can create an initially dry L2.
The analytic Darcy law remains unchanged; only the destination of its
accepted carrier loss differs from the single-layer atmospheric sink.
Settling into L1 does not seed excess pore pressure; elutriation removes
the pressure load associated with outgoing L1 mass.

## Architectural and numerical contracts

Canonical per-layer NAMELIST parameters and separately precomputed, read-only Pascal coefficient caches are distinguished. Unsupported method paths must fail rather than return a false zero rate.

Numerical thresholds, dry-state branches, regularization, and available methods are governed by the linked implementation. The equations above describe implemented physical terms; they are not a replacement for checking each kernel. See the [parent physics guide](../README.md) and [Fortran style](../../../docs/FORTRAN_STYLE.md).
