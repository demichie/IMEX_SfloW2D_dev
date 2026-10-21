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
!********************************************************************************

MODULE rheology_quadratic_2d

  USE parameters_2d, ONLY : wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: friction_factor, eval_quadratic_implicit_source

  REAL(wp) :: friction_factor

CONTAINS

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
