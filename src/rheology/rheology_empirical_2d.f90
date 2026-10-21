!********************************************************************************
!> \brief Empirical depth- and speed-dependent drag model.
!> \details The implemented friction coefficient is
!> \f$C(h,|\boldsymbol{u}|)=10^{-3}|\boldsymbol{u}|/(1+10h)\f$. The
!> horizontal source is
!> \f$\boldsymbol{S}=-\rho_m C\boldsymbol{u}_h/|\boldsymbol{u}_h|\f$.
!> \ref eval_empirical_implicit_source evaluates its complex-valued x/y
!> components for a moving state; the empirical constants are fixed in the
!> routine.
!********************************************************************************

MODULE rheology_empirical_2d

  USE parameters_2d, ONLY : wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_empirical_implicit_source

CONTAINS

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
