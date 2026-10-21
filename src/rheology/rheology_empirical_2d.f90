!********************************************************************************
!> \brief Empirical rheology model terms.
!********************************************************************************

MODULE rheology_empirical_2d

  USE parameters_2d, ONLY : wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_empirical_implicit_source

CONTAINS

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
