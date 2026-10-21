!********************************************************************************
!> \brief Plastic rheology model terms.
!********************************************************************************

MODULE rheology_plastic_2d

  USE parameters_2d, ONLY : wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tau, eval_plastic_implicit_source

  REAL(wp) :: tau

CONTAINS

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
