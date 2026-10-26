!********************************************************************************
!> \brief Froude-dependent Coulomb basal-friction model (model 9).
!> \details The friction coefficient weakens exponentially with the Froude
!> number, \f$\mathrm{Fr}=|\boldsymbol{u}|/\sqrt{g' h}\f$:
!> \f[
!> \mu(\mathrm{Fr})=\mu_\infty+(\mu_0-\mu_\infty)
!>                 \exp(-\mathrm{Fr}/\mathrm{Fr}_0).
!> \f]
!> The horizontal source is
!> \f$\boldsymbol{S}=-\rho_m\mu h G(g'+a_c)
!> \boldsymbol{u}_h/|\boldsymbol{u}_h|\f$, where \f$G\f$ is the supplied
!> gravity coefficient and \f$a_c\f$ is the curvature acceleration. Optional
!> stochastic perturbations shift the Froude number for sufficiently deep
!> flow and clamp it at zero. The module provides
!> \ref eval_rheology_froude_semi_implicit to evaluate the two horizontal momentum rates.
!> \ref solve_rheology_froude_analytic integrates the frozen Froude-dependent
!> velocity source in time.
!********************************************************************************
MODULE rheology_froude_2d
  USE diagnostics_2d, ONLY : fatal_error
  USE parameters_2d, ONLY : wp
  USE rheology_parameters_2d, ONLY : mu_0, mu_inf, Fr_0
  USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs,  &
       fatal_analytic_input_error, solve_exact_velocity_flow,                &
       EXACT_FLOW_FROUDE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_rheology_froude_explicit
  PUBLIC :: eval_rheology_froude_semi_implicit
  PUBLIC :: eval_rheology_froude_implicit
  PUBLIC :: solve_rheology_froude_analytic
  ! Uniform source treatment entry points. Existing kernels remain unchanged.
CONTAINS

  SUBROUTINE eval_rheology_froude_explicit( &
       flow_depth, temperature, velocity_x, velocity_y, speed, horizontal_speed, &
       speed_squared, mixture_density, inverse_density, reduced_gravity, &
       gravity_coefficient, solid_fractions, curvature_term, curvature_xx, curvature_xy, &
       curvature_yy, use_pore_pressure, excess_pore_pressure, moving, use_stochastic, &
       stochastic_value, rate )
    REAL(wp), INTENT(IN) :: flow_depth, temperature, velocity_x, velocity_y
    REAL(wp), INTENT(IN) :: speed, horizontal_speed, speed_squared
    REAL(wp), INTENT(IN) :: mixture_density, inverse_density, reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient, solid_fractions(:)
    REAL(wp), INTENT(IN) :: curvature_xx, curvature_xy, curvature_yy
    REAL(wp), INTENT(IN) :: excess_pore_pressure, stochastic_value
    LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure, moving, use_stochastic
    REAL(wp), INTENT(OUT) :: rate(2)
    REAL(wp) :: sx, sy, curvature_accel
    COMPLEX(wp) :: ix, iy
    curvature_accel=0.0_wp
    IF (curvature_term) curvature_accel=curvature_xx*velocity_x**2 + &
         2.0_wp*curvature_xy*velocity_x*velocity_y+curvature_yy*velocity_y**2
    CALL eval_rheology_froude_semi_implicit( &
       flow_depth, velocity_x, velocity_y, speed, horizontal_speed, mixture_density, &
       reduced_gravity, gravity_coefficient, 0.0_wp, use_stochastic, stochastic_value, &
       sx, sy )
    CALL eval_rheology_froude_implicit( &
       ix, iy )
    rate=[sx,sy]+REAL([ix,iy],wp)
  END SUBROUTINE eval_rheology_froude_explicit

  SUBROUTINE eval_rheology_froude_semi_implicit( flow_depth, velocity_x, velocity_y, speed, &
       horizontal_speed, mixture_density, reduced_gravity, gravity_coefficient, &
       curvature_term, stochastic, stochastic_value, rate_x, rate_y )
    REAL(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y, speed
    REAL(wp), INTENT(IN) :: horizontal_speed, mixture_density, reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient, curvature_term
    LOGICAL, INTENT(IN) :: stochastic
    REAL(wp), INTENT(IN) :: stochastic_value
    REAL(wp), INTENT(OUT) :: rate_x, rate_y
    REAL(wp) :: froude_number
    REAL(wp) :: friction_coefficient
    REAL(wp) :: friction_stress
    rate_x = 0.0_wp
    rate_y = 0.0_wp
    IF ( speed .GT. 0.0_wp .AND. reduced_gravity .GT. 0.0_wp ) THEN
       froude_number = speed / SQRT(reduced_gravity * flow_depth)
       IF ( stochastic ) THEN
          IF ( flow_depth .GT. 0.01_wp ) froude_number = froude_number + stochastic_value
          IF ( froude_number .LT. 0.0_wp ) froude_number = 0.0_wp
       END IF
       friction_coefficient = mu_inf + (mu_0 - mu_inf) * EXP(-froude_number/Fr_0)
       friction_stress = mixture_density * ( friction_coefficient * flow_depth * &
            gravity_coefficient * ( reduced_gravity + curvature_term ) )
       rate_x = rate_x - friction_stress * velocity_x / horizontal_speed
       rate_y = rate_y - friction_stress * velocity_y / horizontal_speed
    END IF
  END SUBROUTINE eval_rheology_froude_semi_implicit

  SUBROUTINE eval_rheology_froude_implicit(rate_x, rate_y)
    COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y
    rate_x = CMPLX(0.0_wp,0.0_wp,wp)
    rate_y = CMPLX(0.0_wp,0.0_wp,wp)
  END SUBROUTINE eval_rheology_froude_implicit

  SUBROUTINE solve_rheology_froude_analytic(integration_step, flow_depth,  &
       speed_ratio, velocity_before, reduced_gravity, gravity_coefficient,   &
       curvature_term, stochastic, stochastic_value, velocity_after)
    REAL(wp), INTENT(IN) :: integration_step, flow_depth, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2), reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient, curvature_term
    LOGICAL, INTENT(IN) :: stochastic
    REAL(wp), INTENT(IN) :: stochastic_value
    REAL(wp), INTENT(OUT) :: velocity_after(2)
    REAL(wp) :: speed_before, speed_after, shift, base_drag
    REAL(wp) :: froude_slope, froude_value_slope
    velocity_after = velocity_before
    speed_before = SQRT(SUM(velocity_before**2))
    CALL validate_analytic_velocity_inputs(integration_step, flow_depth,    &
         speed_ratio, speed_before, 'Froude')
    IF (integration_step .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
    IF (flow_depth .LE. 0.0_wp .OR. reduced_gravity .LE. 0.0_wp .OR.       &
         speed_ratio .LE. 0.0_wp) RETURN
    IF (Fr_0 .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive Fr_0 in Froude analytic solver')
    END IF
    base_drag = gravity_coefficient * (reduced_gravity + curvature_term)
    IF (base_drag .LE. 0.0_wp) RETURN
    shift = 0.0_wp
    IF (stochastic .AND. flow_depth .GT. 0.01_wp) shift = stochastic_value
    froude_value_slope = speed_ratio / SQRT(reduced_gravity*flow_depth)
    froude_slope = froude_value_slope / Fr_0
    CALL solve_exact_velocity_flow(EXACT_FLOW_FROUDE,                       &
         [base_drag, froude_slope, shift/Fr_0, mu_0, mu_inf],                &
         integration_step, speed_before, speed_after)
    velocity_after = (speed_after / speed_before) * velocity_before
  END SUBROUTINE solve_rheology_froude_analytic

  !******************************************************************************
  !> \brief Integrate the Froude-dependent velocity source exactly.
  !> \details Holds depth, gravity, curvature, and the speed ratio fixed while
  !> integrating the frozen scalar speed-decay law.
  !> \param[in] integration_step Exact-flow integration interval [s].
  !> \param[in] flow_depth Explicit flow depth [m].
  !> \param[in] speed_ratio Explicit total-to-horizontal speed ratio.
  !> \param[in] velocity_before Input horizontal velocity [m s^-1].
  !> \param[in] reduced_gravity Explicit reduced gravity [m s^-2].
  !> \param[in] gravity_coefficient Explicit gravity multiplier [dimensionless].
  !> \param[in] curvature_term Explicit curvature acceleration [m s^-2].
  !> \param[in] stochastic Enable the explicit Froude shift.
  !> \param[in] stochastic_value Explicit stochastic increment [dimensionless].
  !> \param[out] velocity_after Returned horizontal velocity [m s^-1].
  !******************************************************************************
  !******************************************************************************
  !> \brief Evaluate model 9's Froude-dependent Coulomb source.
  !> \details Computes the Froude number and the exponentially weakening
  !> friction coefficient, then returns the x/y components of the basal
  !> friction source. A stochastic perturbation is applied only when enabled;
  !> the source is zero unless speed and reduced gravity are positive.
  !> \param[in] flow_depth Flow depth \f$h\f$ [m].
  !> \param[in] velocity_x x velocity [m s^-1].
  !> \param[in] velocity_y y velocity [m s^-1].
  !> \param[in] speed Magnitude used in the Froude number [m s^-1].
  !> \param[in] horizontal_speed Horizontal velocity magnitude [m s^-1].
  !> \param[in] mixture_density Mixture density \f$\rho_m\f$ [kg m^-3].
  !> \param[in] reduced_gravity Reduced gravity \f$g'\f$ [m s^-2].
  !> \param[in] gravity_coefficient Gravity multiplier \f$G\f$ [dimensionless].
  !> \param[in] curvature_term Curvature acceleration \f$a_c\f$ [m s^-2].
  !> \param[in] stochastic Enable the stochastic Froude perturbation.
  !> \param[in] stochastic_value Stochastic increment to the Froude number.
  !> \param[out] rate_x x conservative momentum source rate [kg m^-1 s^-2].
  !> \param[out] rate_y y conservative momentum source rate [kg m^-1 s^-2].
  !******************************************************************************
  ! Unsupported physical treatment, rejected by the source catalog.
  ! Valid zero component of the configured Newton treatment.
END MODULE rheology_froude_2d