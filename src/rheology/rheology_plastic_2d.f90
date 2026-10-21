!********************************************************************************
!> \brief Constant-strength plastic basal-friction model.
!> \details The model applies a speed-independent resisting strength
!> \f$\tau\f$ opposite to the horizontal flow direction:
!> \f[
!> \boldsymbol{S}=-\rho_m\tau\,\boldsymbol{u}_h/|\boldsymbol{u}|.
!> \f]
!> The supplied speed is the normalization used by the implementation. The
!> source is zero when the caller indicates that the flow is not moving.
!> \ref eval_plastic_implicit_source evaluates the complex-valued x/y sources
!> for the implicit solve, while \ref solve_plastic_analytic_velocity solves
!> its algebraic DIRK velocity stage with the other state held explicit.
!> \ref tau stores the model strength.
!********************************************************************************

MODULE rheology_plastic_2d

  USE parameters_2d, ONLY : wp
  USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tau, eval_plastic_implicit_source
  PUBLIC :: solve_plastic_analytic_velocity

  REAL(wp) :: tau

CONTAINS

  !******************************************************************************
  !> \brief Solve a plastic-friction DIRK stage analytically for velocity.
  !> \details Holds depth and the speed-to-horizontal-speed ratio fixed. The
  !> speed is reduced by the constant friction acceleration until it reaches
  !> zero; the horizontal velocity direction is preserved:
  !> \f$r^+=\max(0,r_0-\lambda\tau/(h\beta))\f$.
  !> \param[in] integration_step Effective DIRK integration step (dt * a_ii).
  !> \param[in] flow_depth Explicit flow depth.
  !> \param[in] speed_ratio Explicit ratio of total to horizontal speed.
  !> \param[in] velocity_before Horizontal velocity before the rheology stage.
  !> \param[out] velocity_after Horizontal velocity after the rheology stage.
  !******************************************************************************
  SUBROUTINE solve_plastic_analytic_velocity(integration_step, flow_depth,  &
       speed_ratio, velocity_before, velocity_after)

    REAL(wp), INTENT(IN) :: integration_step, flow_depth, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2)
    REAL(wp), INTENT(OUT) :: velocity_after(2)

    REAL(wp) :: speed_before, speed_after

    velocity_after = velocity_before
    speed_before = SQRT(SUM(velocity_before**2))
    CALL validate_analytic_velocity_inputs(integration_step, flow_depth,    &
         speed_ratio, speed_before, 'plastic')
    IF (integration_step .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN

    speed_after = MAX(0.0_wp, speed_before - integration_step * tau /       &
         (flow_depth * speed_ratio))
    velocity_after = (speed_after / speed_before) * velocity_before

  END SUBROUTINE solve_plastic_analytic_velocity

  !******************************************************************************
  !> \brief Evaluate the implicit plastic-friction source.
  !> \details Returns \f$-\rho_m\tau(u_x,u_y)/|\boldsymbol{u}|\f$ for a
  !> moving state and zero otherwise. Complex arithmetic preserves derivatives
  !> with respect to the state in the implicit Jacobian.
  !> \param velocity_x x component of horizontal velocity.
  !> \param velocity_y y component of horizontal velocity.
  !> \param speed Speed used to normalize the direction.
  !> \param mixture_density Mixture density \f$\rho_m\f$.
  !> \param moving Whether the state is moving.
  !> \param source_x Returned x momentum source.
  !> \param source_y Returned y momentum source.
  !******************************************************************************
  SUBROUTINE eval_plastic_implicit_source(velocity_x, velocity_y, speed,       &
       mixture_density, moving, source_x, source_y)

    COMPLEX(wp), INTENT(IN) :: velocity_x, velocity_y, speed, mixture_density
    LOGICAL, INTENT(IN) :: moving
    COMPLEX(wp), INTENT(OUT) :: source_x, source_y

    source_x = CMPLX(0.0_wp, 0.0_wp, wp)
    source_y = CMPLX(0.0_wp, 0.0_wp, wp)
    IF (.NOT. moving) RETURN

    source_x = -mixture_density * tau * (velocity_x / speed)
    source_y = -mixture_density * tau * (velocity_y / speed)

  END SUBROUTINE eval_plastic_implicit_source

END MODULE rheology_plastic_2d
