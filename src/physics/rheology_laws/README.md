# rheology laws

## Purpose and role in IMEX-SfloW2D

Cell-local momentum resistance models, their shared REAL/COMPLEX mathematical kernels and frozen-coefficient analytic solvers.

## Files in this directory

| File | Responsibility |
| --- | --- |
| [`rheology_2d.f90`](./rheology_2d.f90) | Rheology dispatch and friction-only momentum limiter |
| [`rheology_analytic_utils_2d.f90`](./rheology_analytic_utils_2d.f90) | Exact frozen-law velocity flows |
| [`rheology_parameters_2d.f90`](./rheology_parameters_2d.f90) | Canonical per-layer rheology parameters |
| [`rheology_voellmy_2d.f90`](./rheology_voellmy_2d.f90) | Voellmy basal resistance |
| [`rheology_plastic_2d.f90`](./rheology_plastic_2d.f90) | Plastic yield |
| [`rheology_lava_2d.f90`](./rheology_lava_2d.f90) | Thermal lava viscous and yield drag |
| [`rheology_lahar_2d.f90`](./rheology_lahar_2d.f90) | Lahar temperature/composition drag |
| [`rheology_empirical_2d.f90`](./rheology_empirical_2d.f90) | Empirical depth-speed drag |
| [`rheology_quadratic_2d.f90`](./rheology_quadratic_2d.f90) | Quadratic drag |
| [`rheology_froude_2d.f90`](./rheology_froude_2d.f90) | Froude-dependent friction |
| [`rheology_lucas_2d.f90`](./rheology_lucas_2d.f90) | Lucas friction |
| [`rheology_mui_2d.f90`](./rheology_mui_2d.f90) | Regularized granular mu(I) |
| [`rheology_mui_voellmy_2d.f90`](./rheology_mui_voellmy_2d.f90) | Granular and Voellmy hybrid |
| [`rheology_empirical_drag.inc`](./rheology_empirical_drag.inc) | Shared empirical kernel |
| [`rheology_lahar_drag.inc`](./rheology_lahar_drag.inc) | Shared lahar kernel |
| [`rheology_lava_viscous.inc`](./rheology_lava_viscous.inc) | Shared lava kernel |
| [`rheology_plastic_yield.inc`](./rheology_plastic_yield.inc) | Shared plastic kernel |
| [`rheology_quadratic_drag.inc`](./rheology_quadratic_drag.inc) | Shared quadratic kernel |
| [`rheology_voellmy_turbulent.inc`](./rheology_voellmy_turbulent.inc) | Shared turbulent kernel |

## Physical equations

With horizontal velocity $\boldsymbol u_h$, mixture density $\rho_m$, thickness $h$, and moving-state speed $U=|\boldsymbol u_h|$, the model-specific source opposes motion.

| Law | Defining implemented relation |
| --- | --- |
| Voellmy | Coulomb basal stress plus turbulent $\boldsymbol S_t=-(\rho_m g'/\xi)|\boldsymbol u|^2\boldsymbol u_h/|\boldsymbol u|$. |
| Plastic | $\boldsymbol S=-\rho_m\tau\boldsymbol u_h/|\boldsymbol u|$ for moving states. |
| Lava | Viscous $\boldsymbol S_v=-\rho_m\gamma\boldsymbol u_h$, $\gamma=3\nu_{\mathrm{ref}}\exp[-b(T-T_{\mathrm{ref}})]/h_*$, plus yield $-\tau_0\boldsymbol u_h/|\boldsymbol u|$. |
| Lahar | Fluid viscosity $\eta_f(T,C)=\alpha_1(T)\exp(\beta_1 C)$ with temperature piecewise $\alpha_1$, turbulent resistance and yield $Y=\alpha_2[\exp(\beta_2 C)-1]$. |
| Empirical | $C=10^{-3}|\boldsymbol u|/(1+10h)$, $\boldsymbol S=-\rho_m C\boldsymbol u_h/U$. |
| Quadratic | $\boldsymbol S=-\rho_m f_q|\boldsymbol u|^2\boldsymbol u_h/U$. |
| Froude | $\mathrm{Fr}=|\boldsymbol u|/\sqrt{g'h}$, $\mu=\mu_\infty+(\mu_0-\mu_\infty)e^{-\mathrm{Fr}/\mathrm{Fr}_0}$. |
| Lucas | $\mu=\mu_0$ below threshold $U_w$ or for shallow $h<0.1$; otherwise $\mu_\infty+(\mu_0-\mu_\infty)U_w/U_f$. |
| mu(I) | $I=d\dot\gamma/\sqrt{\sigma'_n/\rho_p}$, $\mu(I)=(\mu_s I_0+\mu_2 I+\mu_{I,\infty}I^2)/(I_0+I)$, $\tau_b=\mu(I)\sigma'_n$. |
| mu(I)+Voellmy | Regularized granular basal stress plus turbulent Voellmy resistance. |

The exact dry-state, pore-pressure, curvature, and velocity regularizations must be read in each law. In isolated dissipative friction groups, the limiter restricts the momentum result to $\boldsymbol m_{\mathrm{after}}=\alpha\boldsymbol m_{\mathrm{before}}$, $0\le\alpha\le1$, preventing reversal or unwanted rotation.

## Architectural and numerical contracts

Model coefficients belong to the canonical layer configuration; the dispatch module chooses the law and supported explicit, semi-implicit, implicit REAL/COMPLEX or analytic path. Each physical-law module retains its uniform treatment entry points even when a contribution is identically zero or the treatment is rejected explicitly; catalog capabilities are a separate contract. Shared .inc files contain substantive formulas, not duplicate state. New models must register in source_catalog_2d.f90; consult docs/ADDING_A_RHEOLOGY_MODEL.md.

Numerical thresholds, dry-state branches, regularization, and available methods are governed by the linked implementation. The equations above describe implemented physical terms; they are not a replacement for checking each kernel. See the [parent physics guide](../README.md) and [Fortran style](../../../docs/FORTRAN_STYLE.md).
