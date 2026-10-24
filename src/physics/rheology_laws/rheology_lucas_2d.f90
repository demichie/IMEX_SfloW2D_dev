!********************************************************************************
!> \brief Lucas velocity-weakening Coulomb basal-friction model (model 10).
!> \details For friction speed \f$U_f\f$, the coefficient is constant below
!> the transition speed \f$U_w\f$ (and for shallow flow), and weakens above
!> that threshold:
!> \f[
!> \mu(U_f)=
!> \begin{cases}
!> \mu_0, & U_f<U_w\ \mathrm{or}\ h<0.1,\\
!> \mu_\infty+(\mu_0-\mu_\infty)U_w/U_f, & \mathrm{otherwise}.
!> \end{cases}
!> \f]
!> The source is \f$\boldsymbol{S}=-\rho_m\mu hG(g'+a_c)
!> \boldsymbol{u}_h/|\boldsymbol{u}_h|\f$. When stochastic mode is enabled,
!> its increment is added to the speed before evaluating the law. The module
!> exposes \ref eval_lucas_momentum_rate for the x/y momentum rates and
!> \ref solve_lucas_analytic_velocity for exact integration of the frozen
!> speed-decay law with non-velocity variables explicit.
!********************************************************************************

MODULE rheology_lucas_2d
  USE diagnostics_2d, ONLY : fatal_error

  USE parameters_2d, ONLY : wp
  USE rheology_parameters_2d, ONLY : mu_0, mu_inf, U_w
  USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs, &
       fatal_analytic_input_error, solve_exact_velocity_flow,                &
       EXACT_FLOW_LUCAS

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_lucas_momentum_rate
  PUBLIC :: solve_lucas_analytic_velocity

  ! Uniform source treatment entry points. Existing kernels remain unchanged.
  PUBLIC :: eval_rheology_lucas_explicit
  PUBLIC :: eval_rheology_lucas_semi_implicit
  INTERFACE eval_rheology_lucas_semi_implicit
    MODULE PROCEDURE eval_lucas_momentum_rate
  END INTERFACE eval_rheology_lucas_semi_implicit
  PUBLIC :: eval_rheology_lucas_implicit
  PUBLIC :: solve_rheology_lucas_analytic
  INTERFACE solve_rheology_lucas_analytic
    MODULE PROCEDURE solve_lucas_analytic_velocity
  END INTERFACE solve_rheology_lucas_analytic

CONTAINS

  !******************************************************************************
  !> \brief Integrate the Lucas velocity-weakening source exactly.
  !> \details Holds depth, gravity, curvature, stochastic speed shift, and
  !> speed ratio fixed while integrating the scalar speed-decay law.
  !> \param[in] integration_step Exact-flow integration interval [s].
  !> \param[in] flow_depth Explicit flow depth [m].
  !> \param[in] speed_ratio Explicit total-to-horizontal speed ratio.
  !> \param[in] velocity_before Input horizontal velocity [m s^-1].
  !> \param[in] reduced_gravity Explicit reduced gravity [m s^-2].
  !> \param[in] gravity_coefficient Explicit gravity multiplier [dimensionless].
  !> \param[in] curvature_term Explicit curvature acceleration [m s^-2].
  !> \param[in] stochastic Enable the explicit speed shift.
  !> \param[in] stochastic_value Explicit stochastic increment [m s^-1].
  !> \param[out] velocity_after Returned horizontal velocity [m s^-1].
  !******************************************************************************
  SUBROUTINE solve_lucas_analytic_velocity(integration_step, flow_depth,   &
       speed_ratio, velocity_before, reduced_gravity, gravity_coefficient,    &
       curvature_term, stochastic, stochastic_value, velocity_after)

    REAL(wp), INTENT(IN) :: integration_step, flow_depth, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2), reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient, curvature_term
    LOGICAL, INTENT(IN) :: stochastic
    REAL(wp), INTENT(IN) :: stochastic_value
    REAL(wp), INTENT(OUT) :: velocity_after(2)

    REAL(wp) :: speed_before, speed_after, shift, scale_factor, shallow

    velocity_after = velocity_before
    speed_before = SQRT(SUM(velocity_before**2))
    CALL validate_analytic_velocity_inputs(integration_step, flow_depth,    &
         speed_ratio, speed_before, 'Lucas')
    IF (integration_step .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
    IF (U_w .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive U_w in Lucas analytic solver')
    END IF

    scale_factor = gravity_coefficient * (reduced_gravity + curvature_term)
    IF (scale_factor .LE. 0.0_wp) RETURN
    shift = 0.0_wp
    IF (stochastic) shift = stochastic_value

    shallow = 0.0_wp
    IF (flow_depth .LT. 0.1_wp) shallow = 1.0_wp
    CALL solve_exact_velocity_flow(EXACT_FLOW_LUCAS,                        &
         [scale_factor, speed_ratio, shift, mu_0, mu_inf, U_w, shallow],    &
         integration_step, speed_before, speed_after)
    velocity_after = (speed_after / speed_before) * velocity_before

  END SUBROUTINE solve_lucas_analytic_velocity

  !******************************************************************************
  !> \brief Evaluate model 10's velocity-dependent Coulomb source.
  !> \details Applies the Lucas transition-speed law, including its shallow
  !> flow branch, and returns the horizontal momentum-source components. No
  !> source is produced at zero speed.
  !> \param[in] flow_depth Flow depth \f$h\f$ [m].
  !> \param[in] velocity_x x velocity [m s^-1].
  !> \param[in] velocity_y y velocity [m s^-1].
  !> \param[in] speed Speed used by the weakening law [m s^-1].
  !> \param[in] horizontal_speed Horizontal velocity magnitude [m s^-1].
  !> \param[in] mixture_density Mixture density \f$\rho_m\f$ [kg m^-3].
  !> \param[in] reduced_gravity Reduced gravity \f$g'\f$ [m s^-2].
  !> \param[in] gravity_coefficient Gravity multiplier \f$G\f$ [dimensionless].
  !> \param[in] curvature_term Curvature acceleration \f$a_c\f$ [m s^-2].
  !> \param[in] stochastic Enable the stochastic speed increment.
  !> \param[in] stochastic_value Stochastic speed increment [m s^-1].
  !> \param[out] rate_x x conservative momentum source rate [kg m^-1 s^-2].
  !> \param[out] rate_y y conservative momentum source rate [kg m^-1 s^-2].
  !******************************************************************************

  SUBROUTINE eval_lucas_momentum_rate( flow_depth, velocity_x, velocity_y, speed, &
       horizontal_speed, mixture_density, reduced_gravity, gravity_coefficient, &
       curvature_term, stochastic, stochastic_value, rate_x, rate_y )

    REAL(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y, speed
    REAL(wp), INTENT(IN) :: horizontal_speed, mixture_density, reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient, curvature_term
    LOGICAL, INTENT(IN) :: stochastic
    REAL(wp), INTENT(IN) :: stochastic_value
    REAL(wp), INTENT(OUT) :: rate_x, rate_y

    REAL(wp) :: friction_speed
    REAL(wp) :: friction_coefficient
    REAL(wp) :: friction_stress

    rate_x = 0.0_wp
    rate_y = 0.0_wp

    IF ( speed .GT. 0.0_wp ) THEN
       IF ( stochastic ) THEN
          friction_speed = speed + stochastic_value
       ELSE
          friction_speed = speed
       END IF

       IF ( friction_speed .LT. U_w ) THEN
          friction_coefficient = mu_0
       ELSE
          IF ( flow_depth .LT. 0.1_wp ) THEN
             friction_coefficient = mu_0
          ELSE
             friction_coefficient = mu_inf + (mu_0-mu_inf)/(friction_speed/U_w)
          END IF
       END IF

       friction_stress = mixture_density * ( friction_coefficient * flow_depth * &
            gravity_coefficient * ( reduced_gravity + curvature_term ) )
       rate_x = rate_x - friction_stress * velocity_x / horizontal_speed
       rate_y = rate_y - friction_stress * velocity_y / horizontal_speed
    END IF

  END SUBROUTINE eval_lucas_momentum_rate

  ! Unsupported physical treatment, rejected by the source catalog.
  SUBROUTINE eval_rheology_lucas_explicit()
    CALL fatal_error('rheology_lucas: explicit treatment is unsupported')
  END SUBROUTINE eval_rheology_lucas_explicit

  ! Unsupported physical treatment, rejected by the source catalog.
  SUBROUTINE eval_rheology_lucas_implicit()
    CALL fatal_error('rheology_lucas: semi-implicit treatment is unsupported')
  END SUBROUTINE eval_rheology_lucas_implicit

END MODULE rheology_lucas_2d
