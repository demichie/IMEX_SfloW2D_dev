!********************************************************************************
!> \brief Quadratic basal-friction model.
!> \details The implicit drag is quadratic in the supplied speed:
!> \f[
!> \boldsymbol{S}=-\rho_m f_q|\boldsymbol{u}|^2
!>                \frac{\boldsymbol{u}_h}{|\boldsymbol{u}_h|},
!> \f]
!> where \f$f_q\f$ is \ref friction_factor. The module provides
!> \ref eval_quadratic_implicit_momentum_rate, which returns zero for a stationary
!> state and otherwise evaluates the complex-valued x/y sources.
!> \ref solve_quadratic_analytic_velocity integrates the frozen quadratic
!> velocity source and reconstructs both components along their input direction.
!********************************************************************************

MODULE rheology_quadratic_2d

  USE parameters_2d, ONLY : wp
  USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs, &
       solve_exact_velocity_flow, EXACT_FLOW_QUADRATIC

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: friction_factor, eval_quadratic_implicit_momentum_rate
  PUBLIC :: solve_quadratic_analytic_velocity

  REAL(wp) :: friction_factor

CONTAINS

  !******************************************************************************
  !> \brief Integrate the quadratic-friction velocity source exactly.
  !> \details Holds depth and the speed-to-horizontal-speed ratio fixed and
  !> applies the reciprocal-speed exact solution of the quadratic drag law.
  !> \param[in] integration_step Exact-flow integration interval [s].
  !> \param[in] flow_depth Explicit flow depth [m].
  !> \param[in] speed_ratio Explicit total-to-horizontal speed ratio.
  !> \param[in] velocity_before Input horizontal velocity [m s^-1].
  !> \param[out] velocity_after Returned horizontal velocity [m s^-1].
  !******************************************************************************
  SUBROUTINE solve_quadratic_analytic_velocity(integration_step, flow_depth, &
       speed_ratio, velocity_before, velocity_after)

    REAL(wp), INTENT(IN) :: integration_step, flow_depth, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2)
    REAL(wp), INTENT(OUT) :: velocity_after(2)

    REAL(wp) :: speed_before, speed_after, quadratic_drag

    velocity_after = velocity_before
    speed_before = SQRT(SUM(velocity_before**2))
    CALL validate_analytic_velocity_inputs(integration_step, flow_depth,    &
         speed_ratio, speed_before, 'quadratic')
    IF (integration_step .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN

    quadratic_drag = friction_factor * speed_ratio**2 / flow_depth
    CALL solve_exact_velocity_flow(EXACT_FLOW_QUADRATIC, [quadratic_drag],  &
         integration_step, speed_before, speed_after)
    velocity_after = (speed_after / speed_before) * velocity_before

  END SUBROUTINE solve_quadratic_analytic_velocity

  !******************************************************************************
  !> \brief Evaluate the implicit quadratic-friction source.
  !> \details Applies the quadratic speed-dependent drag opposite to the
  !> horizontal velocity direction. Complex arithmetic retains sensitivities
  !> for the implicit Jacobian.
  !> \param[in] velocity_x x velocity [m s^-1].
  !> \param[in] velocity_y y velocity [m s^-1].
  !> \param[in] horizontal_speed Horizontal speed magnitude [m s^-1].
  !> \param[in] speed_squared Squared speed used by the quadratic law [m^2 s^-2].
  !> \param[in] mixture_density Mixture density [kg m^-3].
  !> \param[in] moving Whether the state is moving.
  !> \param[out] rate_x x conservative momentum source rate [kg m^-1 s^-2].
  !> \param[out] rate_y y conservative momentum source rate [kg m^-1 s^-2].
  !> \note Complex-valued state and rates preserve complex-step derivatives.
  !******************************************************************************
  SUBROUTINE eval_quadratic_implicit_momentum_rate(velocity_x, velocity_y,   &
       horizontal_speed, speed_squared, mixture_density, moving, rate_x, rate_y)

    COMPLEX(wp), INTENT(IN) :: velocity_x, velocity_y, horizontal_speed
    COMPLEX(wp), INTENT(IN) :: speed_squared, mixture_density
    LOGICAL, INTENT(IN) :: moving
    COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y

    rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
    rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
    IF (.NOT. moving) RETURN

    rate_x = -mixture_density * (velocity_x / horizontal_speed)            &
         * friction_factor * speed_squared
    rate_y = -mixture_density * (velocity_y / horizontal_speed)            &
         * friction_factor * speed_squared

  END SUBROUTINE eval_quadratic_implicit_momentum_rate

END MODULE rheology_quadratic_2d
