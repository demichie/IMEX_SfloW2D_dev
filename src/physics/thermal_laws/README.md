# thermal laws

## Purpose and role in IMEX-SfloW2D

Cell-local radiative/convective cooling: explicit, implicit REAL/COMPLEX and frozen-state analytic treatments.

## Files in this directory

| File | Responsibility |
| --- | --- |
| [`thermal_parameters_2d.f90`](./thermal_parameters_2d.f90) | Thermal source configuration |
| [`thermal_laws_2d.f90`](./thermal_laws_2d.f90) | Model/method dispatcher |
| [`thermal_cooling_2d.f90`](./thermal_cooling_2d.f90) | Radiative-convective cooling and analytic flow |
| [`thermal_cooling_rate.inc`](./thermal_cooling_rate.inc) | Shared REAL/COMPLEX cooling expression |

## Physical equations

The exact shared thermal-cooling source in thermal_cooling_rate.inc is

$$
S_{Q_T}=-r(T^4-T_{\mathrm{amb}}^4)-c(T-T_{\mathrm{amb}}),
\qquad Q_T=M c_{p,m}T.
$$

Here $r,c$ are the configured radiative and convective coefficients; $Q_T$ is conservative *thermal content*, not total kinetic-plus-potential energy. The cooling-only analytic path integrates its supported frozen-composition/heat-capacity flow. The scalar-temperature Newton path avoids rebuilding the physical state where no other coupling requires it.

## Architectural and numerical contracts

REAL and COMPLEX branches use the same cooling-rate kernel, preserving complex-step differentiation. Unsupported selected methods must be rejected by the source catalog and dispatcher.

Numerical thresholds, dry-state branches, regularization, and available methods are governed by the linked implementation. The equations above describe implemented physical terms; they are not a replacement for checking each kernel. See the [parent physics guide](../README.md) and [Fortran style](../../../docs/FORTRAN_STYLE.md).
