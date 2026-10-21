!********************************************************************************
!> \brief Lahar rheology model terms and viscosity relation.
!********************************************************************************

MODULE rheology_lahar_2d

  USE parameters_2d, ONLY : wp, sp, four_thirds, neg_four_thirds
  USE constitutive_parameters_2d, ONLY : grav, inv_grav
  USE rheology_parameters_2d, ONLY : T_ref

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: alpha2, beta2, alpha1_coeff, beta1, Kappa, n_td, n_td2
  PUBLIC :: eval_lahar_fluid_viscosity, eval_lahar_implicit_source
  PUBLIC :: eval_lahar_semi_impl_source, eval_lahar_deposition_rate

  REAL(wp) :: alpha2
  REAL(wp) :: beta2
  REAL(wp) :: alpha1_coeff
  REAL(wp) :: beta1
  REAL(wp) :: Kappa
  REAL(wp) :: n_td
  REAL(wp) :: n_td2

  INTERFACE eval_lahar_fluid_viscosity
     MODULE PROCEDURE eval_lahar_fluid_viscosity_real
     MODULE PROCEDURE eval_lahar_fluid_viscosity_complex
  END INTERFACE eval_lahar_fluid_viscosity

CONTAINS

  FUNCTION eval_lahar_fluid_viscosity_real(temperature, solid_fraction_sum)  &
       RESULT(fluid_viscosity)

    REAL(wp), INTENT(IN) :: temperature, solid_fraction_sum
    REAL(wp) :: fluid_viscosity
    REAL(wp) :: temperature_celsius, exponent_a, exponent_b, alpha1

    temperature_celsius = temperature - 273.15_wp
    IF (temperature_celsius .LT. 20.0_wp) THEN
       exponent_a = 1301.0_wp / (998.333_wp + 8.1855_wp *                    &
            (temperature_celsius - 20.0_wp) + 0.00585_wp *                   &
            (temperature_celsius - 20.0_wp)**2) - 1.30223_wp
       alpha1 = alpha1_coeff * 1.0E-3_wp * 10.0_wp**exponent_a
    ELSE
       exponent_b = (1.3272_wp * (20.0_wp - temperature_celsius)             &
            - 0.001053_wp * (temperature_celsius - 20.0_wp)**2)              &
            / (temperature_celsius + 105.0_wp)
       alpha1 = alpha1_coeff * 1.002E-3_wp * 10.0_wp**exponent_b
    END IF

    fluid_viscosity = alpha1 * EXP(beta1 * solid_fraction_sum)

  END FUNCTION eval_lahar_fluid_viscosity_real

  FUNCTION eval_lahar_fluid_viscosity_complex(temperature, solid_fraction_sum) &
       RESULT(fluid_viscosity)

    COMPLEX(wp), INTENT(IN) :: temperature, solid_fraction_sum
    COMPLEX(wp) :: fluid_viscosity
    COMPLEX(wp) :: temperature_celsius, exponent_a, exponent_b, alpha1

    temperature_celsius = temperature - 273.15_wp
    IF (REAL(temperature_celsius, wp) .LT. 20.0_wp) THEN
       exponent_a = 1301.0_wp / (998.333_wp + 8.1855_wp *                    &
            (temperature_celsius - 20.0_wp) + 0.00585_wp *                   &
            (temperature_celsius - 20.0_wp)**2) - 1.30223_wp
       alpha1 = alpha1_coeff * 1.0E-3_wp * 10.0_wp**exponent_a
    ELSE
       exponent_b = (1.3272_wp * (20.0_wp - temperature_celsius)             &
            - 0.001053_wp * (temperature_celsius - 20.0_wp)**2)              &
            / (temperature_celsius + 105.0_wp)
       alpha1 = alpha1_coeff * 1.002E-3_wp * 10.0_wp**exponent_b
    END IF

    fluid_viscosity = alpha1 * EXP(beta1 * solid_fraction_sum)

  END FUNCTION eval_lahar_fluid_viscosity_complex

  PURE FUNCTION eval_lahar_deposition_rate(solid_fraction, settling_velocity, &
       flow_depth, time_step) RESULT(deposition_rate)

    REAL(wp), INTENT(IN) :: solid_fraction, settling_velocity
    REAL(wp), INTENT(IN) :: flow_depth, time_step
    REAL(wp) :: deposition_rate

    deposition_rate = MIN(solid_fraction * settling_velocity,                 &
         flow_depth * solid_fraction / time_step)

  END FUNCTION eval_lahar_deposition_rate

  SUBROUTINE eval_lahar_implicit_source(flow_depth, velocity_x, velocity_y,    &
       speed, horizontal_speed, temperature, mixture_density, inverse_density, &
       solid_fractions, moving, source_x, source_y)

    COMPLEX(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y
    COMPLEX(wp), INTENT(IN) :: speed, horizontal_speed, temperature
    COMPLEX(wp), INTENT(IN) :: mixture_density, inverse_density
    COMPLEX(wp), INTENT(IN) :: solid_fractions(:)
    LOGICAL, INTENT(IN) :: moving
    COMPLEX(wp), INTENT(OUT) :: source_x, source_y

    COMPLEX(wp) :: fluid_viscosity, inverse_depth, viscous_slope
    COMPLEX(wp) :: turbulent_slope, total_slope, friction_term
    REAL(wp) :: depth_threshold

    source_x = CMPLX(0.0_wp, 0.0_wp, wp)
    source_y = CMPLX(0.0_wp, 0.0_wp, wp)
    IF (wp .EQ. sp) THEN
       depth_threshold = 1.0E-10_wp
    ELSE
       depth_threshold = 1.0E-20_wp
    END IF

    fluid_viscosity = eval_lahar_fluid_viscosity(temperature, SUM(solid_fractions))
    IF (REAL(flow_depth, wp) .GT. depth_threshold) THEN
       inverse_depth = 1.0_wp / flow_depth
       viscous_slope = Kappa * fluid_viscosity * speed * 0.125_wp *          &
            inverse_density * inv_grav * inverse_depth**2
       turbulent_slope = n_td2 * speed**2 * inverse_depth**four_thirds
    ELSE
       viscous_slope = Kappa * fluid_viscosity * speed /                    &
            (8.0_wp * mixture_density * grav * depth_threshold**2)
       turbulent_slope = n_td2 * speed**2 * depth_threshold**neg_four_thirds
    END IF

    total_slope = viscous_slope + turbulent_slope
    IF (.NOT. moving) RETURN

    friction_term = grav * mixture_density * flow_depth * total_slope / horizontal_speed
    source_x = -velocity_x * friction_term
    source_y = -velocity_y * friction_term

  END SUBROUTINE eval_lahar_implicit_source

  SUBROUTINE eval_lahar_semi_impl_source(flow_depth, velocity_x, velocity_y, &
       horizontal_speed, mixture_density, solid_fraction_sum, moving,          &
       source_x, source_y)

    REAL(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y
    REAL(wp), INTENT(IN) :: horizontal_speed, mixture_density, solid_fraction_sum
    LOGICAL, INTENT(IN) :: moving
    REAL(wp), INTENT(OUT) :: source_x, source_y

    REAL(wp) :: yield_strength, yield_slope, friction_term
    REAL(wp), PARAMETER :: depth_threshold = 1.0E-20_wp

    source_x = 0.0_wp
    source_y = 0.0_wp
    IF (.NOT. moving) RETURN

    yield_strength = alpha2 * (EXP(beta2 * solid_fraction_sum) - 1.0_wp)
    IF (flow_depth .GT. depth_threshold) THEN
       yield_slope = yield_strength / (grav * mixture_density * flow_depth)
    ELSE
       yield_slope = yield_strength / (grav * mixture_density * depth_threshold)
    END IF
    friction_term = grav * mixture_density * flow_depth * yield_slope
    source_x = -friction_term * velocity_x / horizontal_speed
    source_y = -friction_term * velocity_y / horizontal_speed

  END SUBROUTINE eval_lahar_semi_impl_source

END MODULE rheology_lahar_2d
