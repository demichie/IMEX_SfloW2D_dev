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
!> \ref eval_rheology_mui_semi_implicit evaluates the real-valued x/y momentum rates,
!> choosing a Sauter diameter and average solid density for mixtures with
!> multiple solid classes. \ref solve_rheology_mui_analytic solves its
!> rational speed-dependent frozen-coefficient source in time.
!********************************************************************************
MODULE rheology_mui_2d
  USE diagnostics_2d, ONLY : fatal_error
  USE parameters_2d, ONLY : wp, n_solid
  USE constitutive_parameters_2d, ONLY : diam_s, rho_s
  USE rheology_parameters_2d, ONLY : I_0, mu_2, mu_s, muI_inf
  USE state_conversion_2d, ONLY : sauter_diameter, average_density_solids
  USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs, &
       fatal_analytic_input_error, solve_exact_velocity_flow,                &
       EXACT_FLOW_MUI
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_rheology_mui_semi_implicit, solve_rheology_mui_analytic
  ! Uniform source treatment entry points. Existing kernels remain unchanged.
  PUBLIC :: eval_rheology_mui_explicit
  PUBLIC :: eval_rheology_mui_implicit
CONTAINS
  !******************************************************************************
  !> \brief Integrate the regularized \f$\mu(I)\f$ source in time.
  !> \details Holds particle properties, effective stress, depth, and the
  !> speed ratio fixed while integrating the rational speed-decay law.
  !> \param[in] integration_step Exact-flow integration interval [s].
  !> \param[in] flow_depth Explicit flow depth [m].
  !> \param[in] speed_ratio Explicit total-to-horizontal speed ratio.
  !> \param[in] velocity_before Input horizontal velocity [m s^-1].
  !> \param[in] mixture_density Explicit mixture density [kg m^-3].
  !> \param[in] reduced_gravity Explicit reduced gravity [m s^-2].
  !> \param[in] gravity_coefficient Explicit gravity multiplier [dimensionless].
  !> \param[in] solid_fraction Explicit solid volume fractions by class.
  !> \param[in] curvature_acceleration Explicit curvature contribution [m s^-2].
  !> \param[in] use_pore_pressure Enable explicit pore-pressure reduction.
  !> \param[in] excess_pore_pressure Explicit excess pore pressure [Pa].
  !> \param[out] velocity_after Returned horizontal velocity [m s^-1].
  !******************************************************************************
  SUBROUTINE solve_rheology_mui_analytic(integration_step, flow_depth,     &
       speed_ratio, velocity_before, mixture_density, reduced_gravity,        &
       gravity_coefficient, solid_fraction, curvature_acceleration,           &
       use_pore_pressure, excess_pore_pressure, velocity_after)
    REAL(wp), INTENT(IN) :: integration_step, flow_depth, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2), mixture_density
    REAL(wp), INTENT(IN) :: reduced_gravity, gravity_coefficient
    REAL(wp), INTENT(IN) :: solid_fraction(n_solid)
    REAL(wp), INTENT(IN) :: curvature_acceleration, excess_pore_pressure
    LOGICAL, INTENT(IN) :: use_pore_pressure
    REAL(wp), INTENT(OUT) :: velocity_after(2)
    REAL(wp) :: speed_before, speed_after, particle_diameter, particle_density
    REAL(wp) :: effective_depth, vertical_stress, normal_stress
    REAL(wp) :: inertial_slope, basal_factor, sqrt_gravity
    velocity_after = velocity_before
    speed_before = SQRT(SUM(velocity_before**2))
    CALL validate_analytic_velocity_inputs(integration_step, flow_depth,    &
         speed_ratio, speed_before, 'mu(I)')
    IF (integration_step .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
    IF (mixture_density .LE. 0.0_wp .OR. gravity_coefficient .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive density or gravity multiplier in mu(I) analytic solver')
    END IF
    IF (n_solid .GT. 1) THEN
       particle_diameter = sauter_diameter(solid_fraction)
       particle_density = average_density_solids(solid_fraction)
    ELSE
       particle_diameter = diam_s(1)
       particle_density = rho_s(1)
    END IF
    IF (particle_diameter .LE. 0.0_wp .OR. particle_density .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive particle property in mu(I) analytic solver')
    END IF
    IF (I_0 .LT. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Negative I_0 in mu(I) analytic solver')
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
    CALL solve_exact_velocity_flow(EXACT_FLOW_MUI,                          &
         [basal_factor*mu_s*I_0, basal_factor*mu_2*inertial_slope,         &
         basal_factor*muI_inf*inertial_slope**2, I_0, inertial_slope],      &
         integration_step, speed_before, speed_after)
    velocity_after = (speed_after / speed_before) * velocity_before
  END SUBROUTINE solve_rheology_mui_analytic
  !******************************************************************************
  !> \brief Evaluate model 11's regularized \f$\mu(I)\f$ source.
  !> \details Computes particle properties, effective normal stress, inertial
  !> number and regularized friction coefficient, then projects basal stress
  !> onto the horizontal velocity direction. Stationary states return zero.
  !> \param[in] flow_depth Flow depth [m].
  !> \param[in] velocity_x x velocity [m s^-1].
  !> \param[in] velocity_y y velocity [m s^-1].
  !> \param[in] speed Speed magnitude [m s^-1].
  !> \param[in] mixture_density Mixture density [kg m^-3].
  !> \param[in] reduced_gravity Reduced gravity [m s^-2].
  !> \param[in] gravity_coefficient Gravity multiplier [dimensionless].
  !> \param[in] solid_fraction Solid volume fractions by class.
  !> \param[in] curvature_acceleration Curvature contribution [m s^-2].
  !> \param[in] use_pore_pressure Enable pore-pressure reduction.
  !> \param[in] excess_pore_pressure Excess pore pressure [Pa].
  !> \param[out] rate_x x momentum source rate.
  !> \param[out] rate_y y momentum source rate.
  !******************************************************************************
  SUBROUTINE eval_rheology_mui_semi_implicit( flow_depth, velocity_x, velocity_y, speed, &
       mixture_density, reduced_gravity, gravity_coefficient, solid_fraction, &
       curvature_acceleration, use_pore_pressure, excess_pore_pressure,       &
       rate_x, rate_y )
    REAL(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y, speed
    REAL(wp), INTENT(IN) :: mixture_density, reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient
    REAL(wp), INTENT(IN) :: solid_fraction(n_solid)
    REAL(wp), INTENT(IN) :: curvature_acceleration
    LOGICAL, INTENT(IN) :: use_pore_pressure
    REAL(wp), INTENT(IN) :: excess_pore_pressure
    REAL(wp), INTENT(OUT) :: rate_x, rate_y
    REAL(wp) :: particle_diameter
    REAL(wp) :: particle_density
    REAL(wp) :: effective_depth
    REAL(wp) :: effective_vertical_stress
    REAL(wp) :: effective_normal_stress
    REAL(wp) :: shear_rate
    REAL(wp) :: inertial_number
    REAL(wp) :: friction_coefficient
    REAL(wp) :: friction_stress
    rate_x = 0.0_wp
    rate_y = 0.0_wp
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
    rate_x = rate_x - friction_stress * (velocity_x / speed)
    rate_y = rate_y - friction_stress * (velocity_y / speed)
  END SUBROUTINE eval_rheology_mui_semi_implicit
  ! Unsupported physical treatment, rejected by the source catalog.
  SUBROUTINE eval_rheology_mui_explicit()
    CALL fatal_error('rheology_mui: explicit treatment is unsupported')
  END SUBROUTINE eval_rheology_mui_explicit
  ! Valid zero component of the configured Newton treatment.
  SUBROUTINE eval_rheology_mui_implicit(rate_x, rate_y)
    COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y
    rate_x = CMPLX(0.0_wp,0.0_wp,wp)
    rate_y = CMPLX(0.0_wp,0.0_wp,wp)
  END SUBROUTINE eval_rheology_mui_implicit
END MODULE rheology_mui_2d