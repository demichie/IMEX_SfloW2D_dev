!********************************************************************************
!> \brief Quadratic basal-friction rheology model terms.
!********************************************************************************

MODULE rheology_quadratic_2d

  USE parameters_2d, ONLY : wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: friction_factor, eval_quadratic_implicit_source

  REAL(wp) :: friction_factor

CONTAINS

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
