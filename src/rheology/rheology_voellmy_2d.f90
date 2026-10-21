!********************************************************************************
!> \brief Voellmy-Salm rheology model terms.
!********************************************************************************

MODULE rheology_voellmy_2d

  USE parameters_2d, ONLY : wp
  USE rheology_parameters_2d, ONLY : xi

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mu, eval_voellmy_implicit_source, eval_voellmy_semi_impl_source

  REAL(wp) :: mu

CONTAINS

  SUBROUTINE eval_voellmy_implicit_source(velocity_x, velocity_y, mixture_density, &
       reduced_gravity, speed, speed_squared, moving, use_stochastic_xi,          &
       stochastic_value, source_x, source_y)

    COMPLEX(wp), INTENT(IN) :: velocity_x, velocity_y
    COMPLEX(wp), INTENT(IN) :: mixture_density, reduced_gravity
    COMPLEX(wp), INTENT(IN) :: speed, speed_squared
    LOGICAL, INTENT(IN) :: moving, use_stochastic_xi
    REAL(wp), INTENT(IN) :: stochastic_value
    COMPLEX(wp), INTENT(OUT) :: source_x, source_y

    REAL(wp) :: xi_effective
    COMPLEX(wp) :: turbulent_stress

    source_x = CMPLX(0.0_wp, 0.0_wp, wp)
    source_y = CMPLX(0.0_wp, 0.0_wp, wp)
    IF (.NOT. moving) RETURN

    IF (use_stochastic_xi) THEN
       xi_effective = xi + stochastic_value
       IF (xi_effective .LT. 1.0_wp) xi_effective = 1.0_wp
       IF (xi_effective .GT. 1.0E5_wp) xi_effective = 1.0E5_wp
    ELSE
       xi_effective = xi
    END IF

    turbulent_stress = mixture_density * reduced_gravity / xi_effective * speed_squared
    source_x = -turbulent_stress * (velocity_x / speed)
    source_y = -turbulent_stress * (velocity_y / speed)

  END SUBROUTINE eval_voellmy_implicit_source

  SUBROUTINE eval_voellmy_semi_impl_source(flow_depth, velocity_x, velocity_y, &
       speed, mixture_density, reduced_gravity, gravity_coefficient,            &
       curvature_term, curvature_xx, curvature_xy, curvature_yy,                &
       use_pore_pressure, excess_pore_pressure, source_x, source_y)

    REAL(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y, speed
    REAL(wp), INTENT(IN) :: mixture_density, reduced_gravity, gravity_coefficient
    LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure
    REAL(wp), INTENT(IN) :: curvature_xx, curvature_xy, curvature_yy
    REAL(wp), INTENT(IN) :: excess_pore_pressure
    REAL(wp), INTENT(OUT) :: source_x, source_y

    REAL(wp) :: curvature_acceleration, friction_stress

    source_x = 0.0_wp
    source_y = 0.0_wp
    IF (speed .LE. 0.0_wp) RETURN

    curvature_acceleration = 0.0_wp
    IF (curvature_term) THEN
       curvature_acceleration = curvature_xx * velocity_x**2                &
            + 2.0_wp * curvature_xy * velocity_x * velocity_y                &
            + curvature_yy * velocity_y**2
    END IF

    friction_stress = mu * mixture_density *                                  &
         (reduced_gravity + curvature_acceleration) * flow_depth              &
         * SQRT(gravity_coefficient)
    friction_stress = MAX(0.0_wp, friction_stress)
    IF (use_pore_pressure) friction_stress = MAX(0.0_wp, friction_stress      &
         - mu * SQRT(gravity_coefficient) * excess_pore_pressure)

    source_x = -friction_stress * velocity_x / speed
    source_y = -friction_stress * velocity_y / speed

  END SUBROUTINE eval_voellmy_semi_impl_source

END MODULE rheology_voellmy_2d
