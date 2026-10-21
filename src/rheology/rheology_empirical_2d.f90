!********************************************************************************
!> \brief Empirical depth- and speed-dependent drag model.
!> \details The implemented friction coefficient is
!> \f$C(h,|\boldsymbol{u}|)=10^{-3}|\boldsymbol{u}|/(1+10h)\f$. The
!> horizontal source is
!> \f$\boldsymbol{S}=-\rho_m C\boldsymbol{u}_h/|\boldsymbol{u}_h|\f$.
!> \ref eval_empirical_implicit_source evaluates its complex-valued x/y
!> components for a moving state; the empirical constants are fixed in the
!> routine. \ref solve_empirical_analytic_velocity solves its algebraic DIRK
!> velocity stage as a linear contraction with depth and speed ratio explicit.
!********************************************************************************

MODULE rheology_empirical_2d

  USE parameters_2d, ONLY : wp
  USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_empirical_implicit_source
  PUBLIC :: solve_empirical_analytic_velocity

CONTAINS

  !******************************************************************************
  !> \brief Solve the empirical-drag DIRK stage analytically for velocity.
  !> \details Holds depth and the speed-to-horizontal-speed ratio fixed. Since
  !> the drag acceleration is proportional to speed, the stage is a scalar
  !> linear contraction of the input velocity,
  !> \f$r^+=r_0/(1+\lambda a)\f$.
  !> \param integration_step Effective DIRK integration step (dt * a_ii).
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

    REAL(wp) :: linear_drag

    velocity_after = velocity_before
    CALL validate_analytic_velocity_inputs(integration_step, flow_depth,    &
         speed_ratio, SQRT(SUM(velocity_before**2)), 'empirical')
    IF (integration_step .EQ. 0.0_wp) RETURN
    IF (SUM(velocity_before**2) .EQ. 0.0_wp) RETURN

    linear_drag = 1.0E-3_wp * speed_ratio /                               &
         (flow_depth * (1.0_wp + 10.0_wp*flow_depth))
    velocity_after = velocity_before / (1.0_wp + integration_step*linear_drag)

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
  !> \param source_x Returned x momentum source.
  !> \param source_y Returned y momentum source.
  !******************************************************************************
  SUBROUTINE eval_empirical_implicit_source(flow_depth, velocity_x, velocity_y, &
       horizontal_speed, speed, mixture_density, moving, source_x, source_y)

    COMPLEX(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y
    COMPLEX(wp), INTENT(IN) :: horizontal_speed, speed, mixture_density
    LOGICAL, INTENT(IN) :: moving
    COMPLEX(wp), INTENT(OUT) :: source_x, source_y

    COMPLEX(wp) :: friction_coefficient

    source_x = CMPLX(0.0_wp, 0.0_wp, wp)
    source_y = CMPLX(0.0_wp, 0.0_wp, wp)
    IF (.NOT. moving) RETURN

    friction_coefficient = 1.0E-3_wp / (1.0_wp + 10.0_wp * flow_depth) * speed
    source_x = -mixture_density * friction_coefficient * (velocity_x / horizontal_speed)
    source_y = -mixture_density * friction_coefficient * (velocity_y / horizontal_speed)

  END SUBROUTINE eval_empirical_implicit_source

END MODULE rheology_empirical_2d
