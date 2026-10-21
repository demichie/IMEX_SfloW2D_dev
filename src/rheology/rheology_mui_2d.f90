!********************************************************************************
!> \brief Regularized \f$\mu(I)\f$ granular basal-friction model (model 11).
!> \details The inertial number is formed from particle diameter \f$d\f$,
!> shear rate \f$\dot\gamma=(5/2)|\boldsymbol{u}|/(h_{\rm eff}\sqrt{G})\f$
!> and effective normal stress
!> \f$\sigma'_n=\max(\epsilon,\sqrt{G}\max(0,\rho_m(g'+a_c)h_{\rm eff}
!> -\mathbf{1}_p p_e))\f$, where \f$h_{\rm eff}\f$ is depth regularized by
!> machine epsilon and \f$\mathbf{1}_p\f$ enables pore-pressure reduction:
!> \f[
!> I=\frac{d\dot\gamma}{\sqrt{\sigma'_n/\rho_p}},\qquad
!> \mu(I)=\frac{\mu_s I_0+\mu_2 I+\mu_{I,\infty}I^2}{I_0+I}.
!> \f]
!> The basal stress is
!> \f$\tau_b=\mu(I)\sigma'_n\f$ and the returned source is
!> \f$\boldsymbol{S}=-\tau_b\boldsymbol{u}_h/|\boldsymbol{u}|\f$.
!> \ref eval_mui_source evaluates the real-valued x/y source components,
!> choosing a Sauter diameter and average solid density for mixtures with
!> multiple solid classes.
!********************************************************************************

MODULE rheology_mui_2d

  USE parameters_2d, ONLY : wp, n_solid
  USE constitutive_parameters_2d, ONLY : diam_s, rho_s
  USE rheology_parameters_2d, ONLY : I_0, mu_2, mu_s, muI_inf
  USE state_conversion_2d, ONLY : sauter_diameter, average_density_solids

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_mui_source

CONTAINS

  !******************************************************************************
  !> \brief Evaluate model 11's regularized \f$\mu(I)\f$ source.
  !> \details Computes particle properties, effective normal stress, inertial
  !> number and regularized friction coefficient, then projects basal stress
  !> onto the horizontal velocity direction. Stationary states return zero.
  !> \param flow_depth Flow depth.
  !> \param velocity_x x component of horizontal velocity.
  !> \param velocity_y y component of horizontal velocity.
  !> \param speed Speed magnitude.
  !> \param mixture_density Mixture density.
  !> \param reduced_gravity Reduced gravity.
  !> \param gravity_coefficient Gravity multiplier.
  !> \param solid_fraction Solid volume fractions by class.
  !> \param curvature_acceleration Curvature contribution to vertical stress.
  !> \param use_pore_pressure Enable pore-pressure reduction.
  !> \param excess_pore_pressure Excess pore pressure.
  !> \param source_x Returned x momentum source.
  !> \param source_y Returned y momentum source.
  !******************************************************************************

  SUBROUTINE eval_mui_source( flow_depth, velocity_x, velocity_y, speed,     &
       mixture_density, reduced_gravity, gravity_coefficient, solid_fraction, &
       curvature_acceleration, use_pore_pressure, excess_pore_pressure,       &
       source_x, source_y )

    REAL(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y, speed
    REAL(wp), INTENT(IN) :: mixture_density, reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient
    REAL(wp), INTENT(IN) :: solid_fraction(n_solid)
    REAL(wp), INTENT(IN) :: curvature_acceleration
    LOGICAL, INTENT(IN) :: use_pore_pressure
    REAL(wp), INTENT(IN) :: excess_pore_pressure
    REAL(wp), INTENT(OUT) :: source_x, source_y

    REAL(wp) :: particle_diameter
    REAL(wp) :: particle_density
    REAL(wp) :: effective_depth
    REAL(wp) :: effective_vertical_stress
    REAL(wp) :: effective_normal_stress
    REAL(wp) :: shear_rate
    REAL(wp) :: inertial_number
    REAL(wp) :: friction_coefficient
    REAL(wp) :: friction_stress

    source_x = 0.0_wp
    source_y = 0.0_wp
    IF ( speed .LE. 0.0_wp ) RETURN

    IF ( n_solid .GT. 1 ) THEN
       particle_diameter = sauter_diameter(solid_fraction)
       particle_density = average_density_solids(solid_fraction)
    ELSE
       particle_diameter = diam_s(1)
       particle_density = rho_s(1)
    END IF

    effective_depth = flow_depth
    IF ( effective_depth .LT. EPSILON(1.0_wp) ) THEN
       effective_depth = EPSILON(1.0_wp)
    END IF

    effective_vertical_stress = mixture_density *                          &
         ( reduced_gravity + curvature_acceleration ) * effective_depth
    IF ( use_pore_pressure ) THEN
       effective_vertical_stress = effective_vertical_stress -              &
            excess_pore_pressure
    END IF
    effective_vertical_stress = MAX(effective_vertical_stress, 0.0_wp)

    effective_normal_stress = MAX(0.0_wp, effective_vertical_stress *       &
         SQRT(gravity_coefficient))
    IF ( effective_normal_stress .LT. EPSILON(1.0_wp) ) THEN
       effective_normal_stress = EPSILON(1.0_wp)
    END IF

    shear_rate = 5.0_wp/2.0_wp * speed /                                    &
         (effective_depth * SQRT(gravity_coefficient))
    inertial_number = particle_diameter * shear_rate /                      &
         SQRT(effective_normal_stress / particle_density)
    friction_coefficient = (mu_s * I_0 + mu_2 * inertial_number +            &
         muI_inf * inertial_number**2) / (I_0 + inertial_number)

    friction_stress = friction_coefficient * effective_normal_stress
    source_x = source_x - friction_stress * (velocity_x / speed)
    source_y = source_y - friction_stress * (velocity_y / speed)

  END SUBROUTINE eval_mui_source

END MODULE rheology_mui_2d
