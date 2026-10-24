!********************************************************************************
!> \brief Voellmy-Salm basal-friction model.
!> \details The model combines a Coulomb basal contribution and a
!> velocity-squared turbulent contribution. The implicit turbulent source is
!> \f[
!> \boldsymbol{S}_t=-\frac{\rho_m g'}{\xi}|\boldsymbol{u}|^2
!>                  \frac{\boldsymbol{u}_h}{|\boldsymbol{u}|}.
!> \f]
!> The semi-implicit Coulomb stress is
!> \f[
!> \tau_b=\max\!\left(0,\,
!> \mu\rho_m h(g'+a_c)\sqrt{G}
!> -\mathbf{1}_{p}\mu\sqrt{G}\,p_e\right),\qquad
!> \boldsymbol{S}_b=-\tau_b\boldsymbol{u}_h/|\boldsymbol{u}|,
!> \f]
!> where \f$a_c\f$ is the optional curvature acceleration and the pore
!> pressure subtraction is enabled only when requested. The module provides
!> \ref eval_voellmy_implicit_momentum_rate for the turbulent term and
!> \ref eval_voellmy_semi_implicit_momentum_rate for the Coulomb term; \ref mu and
!> \ref xi are their respective friction parameters. When enabled,
!> \f$a_c=c_{xx}u_x^2+2c_{xy}u_xu_y+c_{yy}u_y^2\f$.
!> \ref solve_voellmy_analytic_velocity exactly integrates the frozen
!> combined velocity source while holding non-velocity inputs fixed.
!********************************************************************************

MODULE rheology_voellmy_2d
  USE diagnostics_2d, ONLY : fatal_error

  USE parameters_2d, ONLY : wp
  USE rheology_parameters_2d, ONLY : xi
  USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs, &
       fatal_analytic_input_error, solve_exact_velocity_flow,                &
       EXACT_FLOW_QUADRATIC_POLYNOMIAL

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mu, eval_voellmy_implicit_momentum_rate
  PUBLIC :: eval_voellmy_semi_implicit_momentum_rate
  PUBLIC :: solve_voellmy_analytic_velocity

  REAL(wp) :: mu

  ! Uniform source treatment entry points. Existing kernels remain unchanged.
  PUBLIC :: eval_rheology_voellmy_explicit
  PUBLIC :: eval_rheology_voellmy_semi_implicit
  INTERFACE eval_rheology_voellmy_semi_implicit
    MODULE PROCEDURE eval_voellmy_semi_implicit_momentum_rate
  END INTERFACE eval_rheology_voellmy_semi_implicit
  PUBLIC :: eval_rheology_voellmy_implicit
  INTERFACE eval_rheology_voellmy_implicit
    MODULE PROCEDURE eval_voellmy_implicit_momentum_rate
  END INTERFACE eval_rheology_voellmy_implicit
  PUBLIC :: solve_rheology_voellmy_analytic
  INTERFACE solve_rheology_voellmy_analytic
    MODULE PROCEDURE solve_voellmy_analytic_velocity
  END INTERFACE solve_rheology_voellmy_analytic

CONTAINS

  !******************************************************************************
  !> \brief Integrate the Voellmy velocity source exactly.
  !> \details Combines the frozen Coulomb stress with quadratic turbulent
  !> drag. Depth, density, gravity, curvature, pore pressure, and the
  !> speed-to-horizontal-speed ratio remain fixed during integration.
  !> \param[in] integration_step Exact-flow integration interval [s].
  !> \param[in] flow_depth Explicit flow depth [m].
  !> \param[in] speed_ratio Explicit total-to-horizontal speed ratio.
  !> \param[in] velocity_before Input horizontal velocity [m s^-1].
  !> \param[in] mixture_density Explicit mixture density [kg m^-3].
  !> \param[in] reduced_gravity Explicit reduced gravity [m s^-2].
  !> \param[in] gravity_coefficient Explicit gravity multiplier [dimensionless].
  !> \param[in] curvature_acceleration Explicit curvature contribution [m s^-2].
  !> \param[in] use_pore_pressure Enable explicit pore-pressure reduction.
  !> \param[in] excess_pore_pressure Explicit excess pore pressure [Pa].
  !> \param[in] use_stochastic_xi Enable stochastic \f$\xi\f$ adjustment.
  !> \param[in] stochastic_xi_value Explicit stochastic \f$\xi\f$ increment.
  !> \param[out] velocity_after Returned horizontal velocity [m s^-1].
  !******************************************************************************
  SUBROUTINE solve_voellmy_analytic_velocity(integration_step, flow_depth, &
       speed_ratio, velocity_before, mixture_density, reduced_gravity,        &
       gravity_coefficient, curvature_acceleration, use_pore_pressure,       &
       excess_pore_pressure, use_stochastic_xi, stochastic_xi_value,          &
       velocity_after)

    REAL(wp), INTENT(IN) :: integration_step, flow_depth, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2), mixture_density
    REAL(wp), INTENT(IN) :: reduced_gravity, gravity_coefficient
    REAL(wp), INTENT(IN) :: curvature_acceleration, excess_pore_pressure
    LOGICAL, INTENT(IN) :: use_pore_pressure, use_stochastic_xi
    REAL(wp), INTENT(IN) :: stochastic_xi_value
    REAL(wp), INTENT(OUT) :: velocity_after(2)

    REAL(wp) :: speed_before, speed_after, xi_effective
    REAL(wp) :: coulomb_stress, constant_drag, quadratic_drag

    velocity_after = velocity_before
    speed_before = SQRT(SUM(velocity_before**2))
    CALL validate_analytic_velocity_inputs(integration_step, flow_depth,    &
         speed_ratio, speed_before, 'Voellmy')
    IF (integration_step .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
    IF (mixture_density .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive mixture density in Voellmy analytic solver')
    END IF
    IF (gravity_coefficient .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive gravity multiplier in Voellmy analytic solver')
    END IF

    xi_effective = xi
    IF (use_stochastic_xi) THEN
       xi_effective = xi + stochastic_xi_value
       xi_effective = MIN(1.0E5_wp, MAX(1.0_wp, xi_effective))
    ELSEIF (xi_effective .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive xi in Voellmy analytic solver')
    END IF

    coulomb_stress = mu * mixture_density * (reduced_gravity +              &
         curvature_acceleration) * flow_depth * SQRT(gravity_coefficient)
    coulomb_stress = MAX(0.0_wp, coulomb_stress)
    IF (use_pore_pressure) coulomb_stress = MAX(0.0_wp, coulomb_stress -    &
         mu * SQRT(gravity_coefficient) * excess_pore_pressure)
    constant_drag = coulomb_stress / (mixture_density * flow_depth * speed_ratio)
    quadratic_drag = reduced_gravity * speed_ratio / (xi_effective*flow_depth)

    CALL solve_exact_velocity_flow(EXACT_FLOW_QUADRATIC_POLYNOMIAL,         &
         [constant_drag, 0.0_wp, quadratic_drag], integration_step,          &
         speed_before, speed_after)
    velocity_after = (speed_after / speed_before) * velocity_before

  END SUBROUTINE solve_voellmy_analytic_velocity

  !******************************************************************************
  !> \brief Evaluate the implicit turbulent Voellmy source.
  !> \details Computes the velocity-squared stress using \f$\xi\f$ and returns
  !> its x/y components opposite to the velocity direction. Optional
  !> stochastic variation of \f$\xi\f$ is bounded to [1, 10^5].
  !> \param[in] velocity_x x velocity [m s^-1].
  !> \param[in] velocity_y y velocity [m s^-1].
  !> \param[in] mixture_density Mixture density \f$\rho_m\f$ [kg m^-3].
  !> \param[in] reduced_gravity Reduced gravity \f$g'\f$ [m s^-2].
  !> \param[in] speed Speed magnitude [m s^-1].
  !> \param[in] speed_squared Squared speed used in the turbulent stress [m^2 s^-2].
  !> \param[in] moving Whether the state is moving.
  !> \param[in] use_stochastic_xi Enable the stochastic \f$\xi\f$ variation.
  !> \param[in] stochastic_value Stochastic increment to \f$\xi\f$.
  !> \param[out] rate_x x conservative momentum source rate [kg m^-1 s^-2].
  !> \param[out] rate_y y conservative momentum source rate [kg m^-1 s^-2].
  !> \note Complex-valued state and rates preserve complex-step derivatives.
  !******************************************************************************
  SUBROUTINE eval_voellmy_implicit_momentum_rate(velocity_x, velocity_y, mixture_density, &
       reduced_gravity, speed, speed_squared, moving, use_stochastic_xi,          &
       stochastic_value, rate_x, rate_y)

    COMPLEX(wp), INTENT(IN) :: velocity_x, velocity_y
    COMPLEX(wp), INTENT(IN) :: mixture_density, reduced_gravity
    COMPLEX(wp), INTENT(IN) :: speed, speed_squared
    LOGICAL, INTENT(IN) :: moving, use_stochastic_xi
    REAL(wp), INTENT(IN) :: stochastic_value
    COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y

    REAL(wp) :: xi_effective
    COMPLEX(wp) :: turbulent_stress

    rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
    rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
    IF (.NOT. moving) RETURN

    IF (use_stochastic_xi) THEN
       xi_effective = xi + stochastic_value
       IF (xi_effective .LT. 1.0_wp) xi_effective = 1.0_wp
       IF (xi_effective .GT. 1.0E5_wp) xi_effective = 1.0E5_wp
    ELSE
       xi_effective = xi
    END IF

    turbulent_stress = mixture_density * reduced_gravity / xi_effective * speed_squared
    rate_x = -turbulent_stress * (velocity_x / speed)
    rate_y = -turbulent_stress * (velocity_y / speed)

  END SUBROUTINE eval_voellmy_implicit_momentum_rate

  !******************************************************************************
  !> \brief Evaluate the semi-implicit Coulomb Voellmy momentum rate.
  !> \details Builds the basal stress from flow depth, reduced gravity and
  !> optional curvature acceleration, clips negative stress to zero, and
  !> optionally reduces it by excess pore pressure. The returned x/y sources
  !> oppose the horizontal flow direction.
  !> \param[in] flow_depth Flow depth [m].
  !> \param[in] velocity_x x velocity [m s^-1].
  !> \param[in] velocity_y y velocity [m s^-1].
  !> \param[in] speed Speed magnitude [m s^-1].
  !> \param[in] mixture_density Mixture density [kg m^-3].
  !> \param[in] reduced_gravity Reduced gravity [m s^-2].
  !> \param[in] gravity_coefficient Gravity multiplier [dimensionless].
  !> \param[in] curvature_term Enable curvature acceleration.
  !> \param[in] curvature_xx xx curvature coefficient [m^-1].
  !> \param[in] curvature_xy xy curvature coefficient [m^-1].
  !> \param[in] curvature_yy yy curvature coefficient [m^-1].
  !> \param[in] use_pore_pressure Enable pore-pressure stress reduction.
  !> \param[in] excess_pore_pressure Excess pore pressure [Pa].
  !> \param[out] rate_x x momentum source rate.
  !> \param[out] rate_y y momentum source rate.
  !******************************************************************************
  SUBROUTINE eval_voellmy_semi_implicit_momentum_rate(flow_depth, velocity_x, velocity_y, &
       speed, mixture_density, reduced_gravity, gravity_coefficient,            &
       curvature_term, curvature_xx, curvature_xy, curvature_yy,                &
       use_pore_pressure, excess_pore_pressure, rate_x, rate_y)

    REAL(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y, speed
    REAL(wp), INTENT(IN) :: mixture_density, reduced_gravity, gravity_coefficient
    LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure
    REAL(wp), INTENT(IN) :: curvature_xx, curvature_xy, curvature_yy
    REAL(wp), INTENT(IN) :: excess_pore_pressure
    REAL(wp), INTENT(OUT) :: rate_x, rate_y

    REAL(wp) :: curvature_acceleration, friction_stress

    rate_x = 0.0_wp
    rate_y = 0.0_wp
    IF (speed .LE. 0.0_wp) RETURN

    curvature_acceleration = 0.0_wp
    IF (curvature_term) THEN
       curvature_acceleration = curvature_xx * velocity_x**2                &
            + 2.0_wp * curvature_xy * velocity_x * velocity_y                &
            + curvature_yy * velocity_y**2
    END IF

    friction_stress = mu * mixture_density *                                  &
         (reduced_gravity + curvature_acceleration) * flow_depth              &
         * SQRT(gravity_coefficient)
    friction_stress = MAX(0.0_wp, friction_stress)
    IF (use_pore_pressure) friction_stress = MAX(0.0_wp, friction_stress      &
         - mu * SQRT(gravity_coefficient) * excess_pore_pressure)

    rate_x = -friction_stress * velocity_x / speed
    rate_y = -friction_stress * velocity_y / speed

  END SUBROUTINE eval_voellmy_semi_implicit_momentum_rate

  ! Unsupported physical treatment, rejected by the source catalog.
  SUBROUTINE eval_rheology_voellmy_explicit()
    CALL fatal_error('rheology_voellmy: explicit treatment is unsupported')
  END SUBROUTINE eval_rheology_voellmy_explicit

END MODULE rheology_voellmy_2d
