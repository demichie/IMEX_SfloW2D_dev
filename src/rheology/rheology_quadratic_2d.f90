!********************************************************************************
!> \brief Quadratic basal-friction model.
!> \details The implicit drag is quadratic in the supplied speed:
!> \f[
!> \boldsymbol{S}=-\rho_m f_q|\boldsymbol{u}|^2
!>                \frac{\boldsymbol{u}_h}{|\boldsymbol{u}_h|},
!> \f]
!> where \f$f_q\f$ is \ref friction_factor. The module provides
!> \ref eval_quadratic_implicit_source, which returns zero for a stationary
!> state and otherwise evaluates the complex-valued x/y sources.
!> \ref solve_quadratic_analytic_velocity solves the scalar quadratic DIRK
!> stage and reconstructs both velocity components along their input direction.
!********************************************************************************

MODULE rheology_quadratic_2d

  USE parameters_2d, ONLY : wp
  USE rheology_analytic_utils_2d, ONLY : solve_quadratic_in_interval,       &
       validate_analytic_velocity_inputs

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: friction_factor, eval_quadratic_implicit_source
  PUBLIC :: solve_quadratic_analytic_velocity

  REAL(wp) :: friction_factor

CONTAINS

  !******************************************************************************
  !> \brief Solve the quadratic-friction DIRK stage analytically for velocity.
  !> \details Holds depth and the speed-to-horizontal-speed ratio fixed. The
  !> stage speed is the largest root in [0, speed_before] of the quadratic
  !> \f$\lambda a r^2+r-r_0=0\f$ obtained from the model's speed-squared drag.
  !> \param integration_step Effective DIRK integration step (dt * a_ii).
  !> \param flow_depth Explicit flow depth.
  !> \param speed_ratio Explicit total-to-horizontal speed ratio.
  !> \param velocity_before Input horizontal velocity.
  !> \param velocity_after Returned horizontal velocity.
  !******************************************************************************
  SUBROUTINE solve_quadratic_analytic_velocity(integration_step, flow_depth, &
       speed_ratio, velocity_before, velocity_after)

    REAL(wp), INTENT(IN) :: integration_step, flow_depth, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2)
    REAL(wp), INTENT(OUT) :: velocity_after(2)

    REAL(wp) :: speed_before, speed_after, quadratic_drag
    LOGICAL :: root_found

    velocity_after = velocity_before
    speed_before = SQRT(SUM(velocity_before**2))
    CALL validate_analytic_velocity_inputs(integration_step, flow_depth,    &
         speed_ratio, speed_before, 'quadratic')
    IF (integration_step .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN

    quadratic_drag = friction_factor * speed_ratio**2 / flow_depth
    CALL solve_quadratic_in_interval(integration_step*quadratic_drag,      &
         1.0_wp, -speed_before, 0.0_wp, speed_before, speed_after, root_found)
    IF (.NOT. root_found) speed_after = 0.0_wp
    velocity_after = (speed_after / speed_before) * velocity_before

  END SUBROUTINE solve_quadratic_analytic_velocity

  !******************************************************************************
  !> \brief Evaluate the implicit quadratic-friction source.
  !> \details Applies the quadratic speed-dependent drag opposite to the
  !> horizontal velocity direction. Complex arithmetic retains sensitivities
  !> for the implicit Jacobian.
  !> \param velocity_x x component of horizontal velocity.
  !> \param velocity_y y component of horizontal velocity.
  !> \param speed Speed magnitude.
  !> \param horizontal_speed Horizontal speed magnitude.
  !> \param speed_squared Squared speed used by the quadratic law.
  !> \param mixture_density Mixture density.
  !> \param moving Whether the state is moving.
  !> \param source_x Returned x momentum source.
  !> \param source_y Returned y momentum source.
  !******************************************************************************
  SUBROUTINE eval_quadratic_implicit_source(velocity_x, velocity_y, speed,     &
       horizontal_speed, speed_squared, mixture_density, moving, source_x, source_y)

    COMPLEX(wp), INTENT(IN) :: velocity_x, velocity_y, speed, horizontal_speed
    COMPLEX(wp), INTENT(IN) :: speed_squared, mixture_density
    LOGICAL, INTENT(IN) :: moving
    COMPLEX(wp), INTENT(OUT) :: source_x, source_y

    source_x = CMPLX(0.0_wp, 0.0_wp, wp)
    source_y = CMPLX(0.0_wp, 0.0_wp, wp)
    IF (.NOT. moving) RETURN

    source_x = -mixture_density * (velocity_x / horizontal_speed)            &
         * friction_factor * speed_squared
    source_y = -mixture_density * (velocity_y / horizontal_speed)            &
         * friction_factor * speed_squared

  END SUBROUTINE eval_quadratic_implicit_source

END MODULE rheology_quadratic_2d
