!********************************************************************************
!> \brief Combined regularized \f$\mu(I)\f$ and Voellmy model (model 12).
!> \details The semi-implicit basal term uses the regularized granular law
!> \f[
!> I=\frac{d\dot\gamma}{\sqrt{\sigma'_n/\rho_p}},\quad
!> \mu(I)=\frac{\mu_s I_0+\mu_2 I+\mu_{I,\infty}I^2}{I_0+I},\quad
!> \tau_b=\mu(I)\sigma'_n,
!> \f]
!> with \f$\dot\gamma=(5/2)|\boldsymbol{u}|/(h_{\rm eff}\sqrt{G})\f$ and
!> \f$\sigma'_n=\max(\epsilon,\sqrt{G}\max(0,\rho_m(g'+a_c)h_{\rm eff}
!> -\mathbf{1}_p p_e))\f$. Here \f$h_{\rm eff}\f$ is depth regularized by
!> machine epsilon, and \f$\mathbf{1}_p\f$ enables pore-pressure reduction.
!> The basal source is
!> \f$\boldsymbol{S}_b=-\tau_b\boldsymbol{u}_h/|\boldsymbol{u}|\f$. The
!> implicit turbulent term is
!> \f$\boldsymbol{S}_t=-\rho_mg'|\boldsymbol{u}|^2
!> \boldsymbol{u}_h/(\xi|\boldsymbol{u}|)\f$. The module exposes
!> \ref eval_mui_voellmy_source for the regularized basal part and
!> \ref eval_mui_voellmy_implicit_source for the turbulent part. The granular
!> source uses Sauter particle diameter and average particle density when
!> multiple solid classes are present.
!> \ref solve_mui_voellmy_analytic_velocity integrates the combined
!> speed-dependent source with non-velocity state frozen explicitly.
!********************************************************************************

MODULE rheology_mui_voellmy_2d

  USE parameters_2d, ONLY : wp, n_solid
  USE constitutive_parameters_2d, ONLY : diam_s, rho_s
  USE rheology_parameters_2d, ONLY : I_0, mu_2, mu_s, muI_inf, xi
  USE state_conversion_2d, ONLY : sauter_diameter, average_density_solids
  USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs, &
       fatal_analytic_input_error, solve_exact_velocity_flow,                &
       EXACT_FLOW_MUI_VOELLMY

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_mui_voellmy_source, eval_mui_voellmy_implicit_source
  PUBLIC :: solve_mui_voellmy_analytic_velocity

CONTAINS

  !******************************************************************************
  !> \brief Integrate the frozen \f$\mu(I)\f$-Voellmy source in time.
  !> \details Holds particle properties, stress, depth, and speed ratio fixed
  !> while integrating the combined rational and quadratic decay law.
  !> \param integration_step Exact-flow integration interval.
  !> \param flow_depth Explicit flow depth.
  !> \param speed_ratio Explicit total-to-horizontal speed ratio.
  !> \param velocity_before Input horizontal velocity.
  !> \param mixture_density Explicit mixture density.
  !> \param reduced_gravity Explicit reduced gravity.
  !> \param gravity_coefficient Explicit gravity multiplier.
  !> \param solid_fraction Explicit solid volume fractions.
  !> \param curvature_acceleration Explicit curvature contribution.
  !> \param use_pore_pressure Enable explicit pore-pressure reduction.
  !> \param excess_pore_pressure Explicit excess pore pressure.
  !> \param use_stochastic_xi Enable stochastic \f$\xi\f$ adjustment.
  !> \param stochastic_xi_value Explicit stochastic \f$\xi\f$ increment.
  !> \param velocity_after Returned horizontal velocity.
  !******************************************************************************
  SUBROUTINE solve_mui_voellmy_analytic_velocity(integration_step,         &
       flow_depth, speed_ratio, velocity_before, mixture_density,             &
       reduced_gravity, gravity_coefficient, solid_fraction,                  &
       curvature_acceleration, use_pore_pressure, excess_pore_pressure,       &
       use_stochastic_xi, stochastic_xi_value, velocity_after)

    REAL(wp), INTENT(IN) :: integration_step, flow_depth, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2), mixture_density
    REAL(wp), INTENT(IN) :: reduced_gravity, gravity_coefficient
    REAL(wp), INTENT(IN) :: solid_fraction(n_solid)
    REAL(wp), INTENT(IN) :: curvature_acceleration, excess_pore_pressure
    LOGICAL, INTENT(IN) :: use_pore_pressure, use_stochastic_xi
    REAL(wp), INTENT(IN) :: stochastic_xi_value
    REAL(wp), INTENT(OUT) :: velocity_after(2)

    REAL(wp) :: speed_before, speed_after, particle_diameter, particle_density
    REAL(wp) :: effective_depth, vertical_stress, normal_stress
    REAL(wp) :: inertial_slope, basal_factor, sqrt_gravity
    REAL(wp) :: xi_effective, turbulent_drag

    velocity_after = velocity_before
    speed_before = SQRT(SUM(velocity_before**2))
    CALL validate_analytic_velocity_inputs(integration_step, flow_depth,    &
         speed_ratio, speed_before, 'mu(I)-Voellmy')
    IF (integration_step .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
    IF (mixture_density .LE. 0.0_wp .OR. gravity_coefficient .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive density or gravity multiplier in mu(I)-Voellmy analytic solver')
    END IF

    IF (n_solid .GT. 1) THEN
       particle_diameter = sauter_diameter(solid_fraction)
       particle_density = average_density_solids(solid_fraction)
    ELSE
       particle_diameter = diam_s(1)
       particle_density = rho_s(1)
    END IF
    IF (particle_diameter .LE. 0.0_wp .OR. particle_density .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive particle property in mu(I)-Voellmy analytic solver')
    END IF
    IF (I_0 .LT. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Negative I_0 in mu(I)-Voellmy analytic solver')
    END IF

    effective_depth = MAX(flow_depth, EPSILON(1.0_wp))
    vertical_stress = mixture_density *                                   &
         (reduced_gravity + curvature_acceleration) * effective_depth
    IF (use_pore_pressure) vertical_stress = vertical_stress -             &
         excess_pore_pressure
    vertical_stress = MAX(vertical_stress, 0.0_wp)
    normal_stress = MAX(EPSILON(1.0_wp),                                  &
         vertical_stress * SQRT(gravity_coefficient))

    sqrt_gravity = SQRT(gravity_coefficient)
    inertial_slope = 2.5_wp * particle_diameter * speed_ratio /             &
         (effective_depth * sqrt_gravity *                                  &
         SQRT(normal_stress/particle_density))
    basal_factor = normal_stress / (mixture_density*flow_depth*speed_ratio)

    xi_effective = xi
    IF (use_stochastic_xi) THEN
       xi_effective = xi + stochastic_xi_value
       xi_effective = MIN(1.0E5_wp, MAX(1.0_wp, xi_effective))
    ELSEIF (xi_effective .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive xi in mu(I)-Voellmy analytic solver')
    END IF
    turbulent_drag = reduced_gravity * speed_ratio / (xi_effective*flow_depth)
    CALL solve_exact_velocity_flow(EXACT_FLOW_MUI_VOELLMY,                  &
         [basal_factor*mu_s*I_0, basal_factor*mu_2*inertial_slope,         &
         basal_factor*muI_inf*inertial_slope**2, I_0, inertial_slope,       &
         turbulent_drag], integration_step, speed_before, speed_after)
    velocity_after = (speed_after / speed_before) * velocity_before

  END SUBROUTINE solve_mui_voellmy_analytic_velocity

  !******************************************************************************
  !> \brief Evaluate the implicit Voellmy turbulent source.
  !> \details Computes the velocity-squared stress with \f$\xi\f$ and returns
  !> its x/y components. Optional stochastic variation of \f$\xi\f$ is bounded
  !> to [1, 10^5]; a stationary state returns zero.
  !> \param velocity_x x component of horizontal velocity.
  !> \param velocity_y y component of horizontal velocity.
  !> \param mixture_density Mixture density.
  !> \param reduced_gravity Reduced gravity.
  !> \param speed Speed magnitude.
  !> \param speed_squared Squared speed.
  !> \param moving Whether the state is moving.
  !> \param use_stochastic_xi Enable stochastic \f$\xi\f$ variation.
  !> \param stochastic_value Stochastic increment to \f$\xi\f$.
  !> \param source_x Returned x momentum source.
  !> \param source_y Returned y momentum source.
  !******************************************************************************
  SUBROUTINE eval_mui_voellmy_implicit_source(velocity_x, velocity_y,        &
       mixture_density, reduced_gravity, speed, speed_squared, moving,       &
       use_stochastic_xi, stochastic_value, source_x, source_y)

    COMPLEX(wp), INTENT(IN) :: velocity_x, velocity_y
    COMPLEX(wp), INTENT(IN) :: mixture_density, reduced_gravity
    COMPLEX(wp), INTENT(IN) :: speed, speed_squared
    LOGICAL, INTENT(IN) :: moving, use_stochastic_xi
    REAL(wp), INTENT(IN) :: stochastic_value
    COMPLEX(wp), INTENT(OUT) :: source_x, source_y

    REAL(wp) :: xi_effective
    COMPLEX(wp) :: turbulent_stress

    source_x = CMPLX(0.0_wp, 0.0_wp, wp)
    source_y = CMPLX(0.0_wp, 0.0_wp, wp)
    IF (.NOT. moving) RETURN

    IF (use_stochastic_xi) THEN
       xi_effective = xi + stochastic_value
       IF (xi_effective .LT. 1.0_wp) xi_effective = 1.0_wp
       IF (xi_effective .GT. 1.0E5_wp) xi_effective = 1.0E5_wp
    ELSE
       xi_effective = xi
    END IF

    turbulent_stress = mixture_density * reduced_gravity / xi_effective * speed_squared
    source_x = -turbulent_stress * (velocity_x / speed)
    source_y = -turbulent_stress * (velocity_y / speed)

  END SUBROUTINE eval_mui_voellmy_implicit_source

  !******************************************************************************
  !> \brief Evaluate model 12's regularized \f$\mu(I)\f$ basal source.
  !> \details Computes the inertial number and effective normal stress,
  !> including optional pore-pressure reduction, then returns the basal
  !> friction source in x and y. Stationary states return zero.
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

  SUBROUTINE eval_mui_voellmy_source( flow_depth, velocity_x, velocity_y,    &
       speed, mixture_density, reduced_gravity, gravity_coefficient,         &
       solid_fraction, curvature_acceleration, use_pore_pressure,             &
       excess_pore_pressure, source_x, source_y )

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

  END SUBROUTINE eval_mui_voellmy_source

END MODULE rheology_mui_voellmy_2d
