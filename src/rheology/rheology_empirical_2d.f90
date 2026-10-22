!********************************************************************************
!> \brief Empirical depth- and speed-dependent drag model.
!> \details The implemented friction coefficient is
!> \f$C(h,|\boldsymbol{u}|)=10^{-3}|\boldsymbol{u}|/(1+10h)\f$. The
!> horizontal source is
!> \f$\boldsymbol{S}=-\rho_m C\boldsymbol{u}_h/|\boldsymbol{u}_h|\f$.
!> \ref eval_empirical_implicit_momentum_rate evaluates its complex-valued x/y
!> components for a moving state; the empirical constants are fixed in the
!> routine. \ref solve_empirical_analytic_velocity integrates its frozen
!> linear-drag velocity source exactly.
!********************************************************************************

MODULE rheology_empirical_2d

  USE parameters_2d, ONLY : wp
  USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs,  &
       solve_exact_velocity_flow, EXACT_FLOW_LINEAR

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_empirical_implicit_momentum_rate
  PUBLIC :: solve_empirical_analytic_velocity

CONTAINS

  !******************************************************************************
  !> \brief Integrate the empirical-drag velocity source exactly.
  !> \details Holds depth and the speed-to-horizontal-speed ratio fixed. Since
  !> the drag acceleration is proportional to speed, the input velocity is
  !> contracted exponentially.
  !> \param integration_step Exact-flow interval, normally dt * a_ii.
  !> \param flow_depth Explicit flow depth.
  !> \param speed_ratio Explicit total-to-horizontal speed ratio.
  !> \param velocity_before Input horizontal velocity.
  !> \param velocity_after Returned horizontal velocity.
  !******************************************************************************
  SUBROUTINE solve_empirical_analytic_velocity(integration_step, flow_depth, &
       speed_ratio, velocity_before, velocity_after)

    REAL(wp), INTENT(IN) :: integration_step, flow_depth, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2)
    REAL(wp), INTENT(OUT) :: velocity_after(2)

    REAL(wp) :: linear_drag, speed_before, speed_after

    velocity_after = velocity_before
    speed_before = SQRT(SUM(velocity_before**2))
    CALL validate_analytic_velocity_inputs(integration_step, flow_depth,    &
         speed_ratio, speed_before, 'empirical')
    IF (integration_step .EQ. 0.0_wp) RETURN
    IF (speed_before .EQ. 0.0_wp) RETURN

    linear_drag = 1.0E-3_wp * speed_ratio /                               &
         (flow_depth * (1.0_wp + 10.0_wp*flow_depth))
    CALL solve_exact_velocity_flow(EXACT_FLOW_LINEAR, [linear_drag],        &
         integration_step, speed_before, speed_after)
    velocity_after = (speed_after / speed_before) * velocity_before

  END SUBROUTINE solve_empirical_analytic_velocity

  !******************************************************************************
  !> \brief Evaluate the implicit empirical drag source.
  !> \details Forms the depth-dependent coefficient proportional to speed and
  !> returns the x/y momentum sources opposite to horizontal flow.
  !> \param flow_depth Flow depth.
  !> \param velocity_x x component of horizontal velocity.
  !> \param velocity_y y component of horizontal velocity.
  !> \param horizontal_speed Horizontal speed magnitude.
  !> \param speed Speed magnitude used by the empirical coefficient.
  !> \param mixture_density Mixture density.
  !> \param moving Whether the state is moving.
  !> \param rate_x Returned x momentum rate.
  !> \param rate_y Returned y momentum rate.
  !******************************************************************************
  SUBROUTINE eval_empirical_implicit_momentum_rate(flow_depth, velocity_x, velocity_y, &
       horizontal_speed, speed, mixture_density, moving, rate_x, rate_y)

    COMPLEX(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y
    COMPLEX(wp), INTENT(IN) :: horizontal_speed, speed, mixture_density
    LOGICAL, INTENT(IN) :: moving
    COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y

    COMPLEX(wp) :: friction_coefficient

    rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
    rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
    IF (.NOT. moving) RETURN

    friction_coefficient = 1.0E-3_wp / (1.0_wp + 10.0_wp * flow_depth) * speed
    rate_x = -mixture_density * friction_coefficient * (velocity_x / horizontal_speed)
    rate_y = -mixture_density * friction_coefficient * (velocity_y / horizontal_speed)

  END SUBROUTINE eval_empirical_implicit_momentum_rate

END MODULE rheology_empirical_2d
