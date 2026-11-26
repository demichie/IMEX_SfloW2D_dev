# mass injection terms

## Purpose and role in IMEX-SfloW2D

Bottom circular/fissural, radial and boundary source material injection.

## Files in this directory

| File | Responsibility |
| --- | --- |
| [`mass_injection_coefficients_2d.f90`](./mass_injection_coefficients_2d.f90) | Source windows and injected material properties |
| [`mass_injection_2d.f90`](./mass_injection_2d.f90) | Conservative source dispatch |
| [`mass_bottom_circular_injection_2d.f90`](./mass_bottom_circular_injection_2d.f90) | Circular bottom injection |
| [`mass_bottom_fissural_injection_2d.f90`](./mass_bottom_fissural_injection_2d.f90) | Fissural bottom injection |
| [`mass_radial_injection_2d.f90`](./mass_radial_injection_2d.f90) | Radial lateral injection |
| [`mass_boundary_injection_2d.f90`](./mass_boundary_injection_2d.f90) | Boundary injection traces |

## Physical equations

For an injection activation coefficient \(a(t)\), interval-based source scheduling uses

\[
\bar a_{[t,t+\Delta t]}=\frac1{\Delta t}\int_t^{t+\Delta t}a(s)\,ds.
\]

A volumetric injection rate per horizontal area \(q_V\), from geometric source coverage and the selected time window, adds mass \(S_M=\rho_{\mathrm{inj}}q_V\). Constituent and thermal rates use the injected composition and temperature. Bottom circular, fissural and radial geometries compute different local \(q_V\); boundary injection instead supplies physical boundary traces and must not be treated as a cell-volume source.

## Architectural and numerical contracts

Material descriptors and time-window coefficients can be prepared once outside OpenMP cell loops. Injection is an external mass source, unlike antisymmetric interlayer exchange.

Numerical thresholds, dry-state branches, regularization, and available methods are governed by the linked implementation. The equations above describe implemented physical terms; they are not a replacement for checking each kernel. See the [parent physics guide](../README.md) and [Fortran style](../../../docs/FORTRAN_STYLE.md).
