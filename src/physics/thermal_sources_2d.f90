!********************************************************************************
!> \brief General thermal source terms for the depth-averaged equations.
!> \details Provides an independent exact scalar ODE flow for lava cooling;
!>          callers decide whether and how to use it.
!********************************************************************************

MODULE thermal_sources_2d

  USE parameters_2d, ONLY : wp
  USE diagnostics_2d, ONLY : fatal_error
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_inlet_thermal_energy_source
  PUBLIC :: eval_lava_cooling_source
  PUBLIC :: advance_lava_cooling_temperature
  PUBLIC :: eval_mass_exchange_thermal_energy
  PUBLIC :: eval_vaporization_sources

CONTAINS

  PURE FUNCTION eval_inlet_thermal_energy_source(time_coefficient, injection_rate, &
       mixture_density, mixture_heat_capacity, source_temperature) RESULT(source)

    REAL(wp), INTENT(IN) :: time_coefficient
    REAL(wp), INTENT(IN) :: injection_rate
    REAL(wp), INTENT(IN) :: mixture_density
    REAL(wp), INTENT(IN) :: mixture_heat_capacity
    REAL(wp), INTENT(IN) :: source_temperature
    REAL(wp) :: source

    source = time_coefficient * injection_rate * mixture_density *             &
         mixture_heat_capacity * source_temperature

  END FUNCTION eval_inlet_thermal_energy_source

  PURE FUNCTION eval_lava_cooling_source(temperature, radiative_coefficient,   &
       convective_coefficient, ambient_temperature) RESULT(source)

    COMPLEX(wp), INTENT(IN) :: temperature
    REAL(wp), INTENT(IN) :: radiative_coefficient
    REAL(wp), INTENT(IN) :: convective_coefficient
    REAL(wp), INTENT(IN) :: ambient_temperature
    COMPLEX(wp) :: source

    source = -radiative_coefficient * (temperature**4 - ambient_temperature**4) &
         - convective_coefficient * (temperature - ambient_temperature)

  END FUNCTION eval_lava_cooling_source

  !******************************************************************************
  !> \brief Advance the lava cooling ODE exactly over a time interval.
  !> \details Separates the autonomous temperature equation and evaluates its
  !>          analytic time-temperature primitive. A bisection inversion finds
  !>          the temperature after the requested DIRK integration step. The
  !>          areal thermal capacity and cooling coefficients are held fixed.
  !> \param integration_step Effective DIRK integration step (dt * a_ii).
  !******************************************************************************
  SUBROUTINE advance_lava_cooling_temperature(temperature, integration_step,    &
       areal_heat_capacity, radiative_coefficient, convective_coefficient,     &
       ambient_temperature, temperature_after)

    REAL(wp), INTENT(IN) :: temperature, integration_step, areal_heat_capacity
    REAL(wp), INTENT(IN) :: radiative_coefficient, convective_coefficient
    REAL(wp), INTENT(IN) :: ambient_temperature
    REAL(wp), INTENT(OUT) :: temperature_after

    REAL(wp) :: beta, cubic_root, lower, upper, midpoint
    REAL(wp) :: initial_primitive, target_primitive, midpoint_primitive
    REAL(wp) :: elapsed_scale, tolerance
    INTEGER :: iteration

    IF (.NOT. ieee_is_finite(temperature) .OR.                              &
         .NOT. ieee_is_finite(integration_step) .OR.                        &
         .NOT. ieee_is_finite(areal_heat_capacity) .OR.                    &
         .NOT. ieee_is_finite(radiative_coefficient) .OR.                  &
         .NOT. ieee_is_finite(convective_coefficient) .OR.                 &
         .NOT. ieee_is_finite(ambient_temperature)) THEN
       CALL fatal_error('Nonfinite input to exact lava cooling integrator')
    END IF
    IF (temperature .LE. 0.0_wp .OR. ambient_temperature .LE. 0.0_wp) THEN
       CALL fatal_error('Nonpositive temperature in exact lava cooling integrator')
    END IF
    IF (integration_step .LT. 0.0_wp .OR. areal_heat_capacity .LE. 0.0_wp) THEN
       CALL fatal_error('Invalid time or heat capacity in exact lava cooling integrator')
    END IF
    IF (radiative_coefficient .LT. 0.0_wp .OR.                            &
         convective_coefficient .LT. 0.0_wp) THEN
       CALL fatal_error('Negative cooling coefficient in exact lava cooling integrator')
    END IF

    temperature_after = temperature
    IF (integration_step .EQ. 0.0_wp .OR.                                   &
         (radiative_coefficient .EQ. 0.0_wp .AND.                          &
         convective_coefficient .EQ. 0.0_wp) .OR.                          &
         temperature .EQ. ambient_temperature) RETURN

    IF (radiative_coefficient .EQ. 0.0_wp) THEN
       temperature_after = ambient_temperature + (temperature -             &
            ambient_temperature) * EXP(-convective_coefficient *          &
            integration_step / areal_heat_capacity)
       RETURN
    END IF

    beta = convective_coefficient / (radiative_coefficient *               &
         ambient_temperature**3)
    IF (.NOT. ieee_is_finite(beta)) THEN
       CALL fatal_error('Cooling coefficient ratio overflow in lava integrator')
    END IF
    CALL lava_cooling_cubic_root(beta, cubic_root)

    initial_primitive = lava_cooling_time_primitive(                        &
         temperature / ambient_temperature, beta, cubic_root)
    elapsed_scale = radiative_coefficient * ambient_temperature**3 *        &
         integration_step / areal_heat_capacity
    IF (.NOT. ieee_is_finite(elapsed_scale)) THEN
       temperature_after = ambient_temperature
       RETURN
    END IF
    target_primitive = initial_primitive - elapsed_scale

    IF (temperature .GT. ambient_temperature) THEN
       lower = 1.0_wp
       upper = temperature / ambient_temperature
    ELSE
       lower = temperature / ambient_temperature
       upper = 1.0_wp
    END IF

    tolerance = 8.0_wp * EPSILON(1.0_wp) *                                  &
         MAX(1.0_wp, ABS(lower), ABS(upper))
    DO iteration = 1, 256
       midpoint = lower + 0.5_wp * (upper - lower)
       IF (midpoint .EQ. lower .OR. midpoint .EQ. upper .OR.                &
            ABS(upper - lower) .LE. tolerance) EXIT
       IF (midpoint .EQ. 1.0_wp) THEN
          IF (temperature .GT. ambient_temperature) THEN
             lower = midpoint
          ELSE
             upper = midpoint
          END IF
          CYCLE
       END IF

       midpoint_primitive = lava_cooling_time_primitive(midpoint, beta,       &
            cubic_root)
       IF (temperature .GT. ambient_temperature) THEN
          IF (midpoint_primitive .LT. target_primitive) THEN
             lower = midpoint
          ELSE
             upper = midpoint
          END IF
       ELSE
          IF (midpoint_primitive .GT. target_primitive) THEN
             lower = midpoint
          ELSE
             upper = midpoint
          END IF
       END IF
    END DO
    IF (iteration .GT. 256) THEN
       CALL fatal_error('Exact lava cooling inversion failed to converge')
    END IF

    midpoint = lower + 0.5_wp * (upper - lower)
    temperature_after = ambient_temperature * midpoint
    IF (temperature .GT. ambient_temperature) THEN
       temperature_after = MAX(ambient_temperature,                        &
            MIN(temperature, temperature_after))
    ELSE
       temperature_after = MIN(ambient_temperature,                        &
            MAX(temperature, temperature_after))
    END IF

  END SUBROUTINE advance_lava_cooling_temperature

  !******************************************************************************
  !> \brief Find the unique real root below -1 of the normalized cooling cubic.
  !******************************************************************************
  SUBROUTINE lava_cooling_cubic_root(beta, root)

    REAL(wp), INTENT(IN) :: beta
    REAL(wp), INTENT(OUT) :: root

    REAL(wp) :: lower, upper, midpoint, scale, polynomial
    INTEGER :: iteration

    IF (beta .EQ. 0.0_wp) THEN
       root = -1.0_wp
       RETURN
    END IF

    scale = MAX(1.0_wp, (1.0_wp + beta)**(1.0_wp / 3.0_wp))
    lower = -2.0_wp
    upper = -1.0_wp / scale
    DO iteration = 1, 256
       midpoint = lower + 0.5_wp * (upper - lower)
       IF (midpoint .EQ. lower .OR. midpoint .EQ. upper) EXIT
       polynomial = midpoint**3 + midpoint**2/scale + midpoint/scale**2 +  &
            1.0_wp
       IF (polynomial .LT. 0.0_wp) THEN
          lower = midpoint
       ELSE
          upper = midpoint
       END IF
    END DO
    IF (iteration .GT. 256) THEN
       CALL fatal_error('Cooling cubic root calculation failed to converge')
    END IF
    root = scale * (lower + 0.5_wp * (upper - lower))

  END SUBROUTINE lava_cooling_cubic_root

  !******************************************************************************
  !> \brief Evaluate the analytic primitive for the normalized cooling law.
  !******************************************************************************
  PURE FUNCTION lava_cooling_time_primitive(x, beta, cubic_root) RESULT(value)

    REAL(wp), INTENT(IN) :: x, beta, cubic_root
    REAL(wp) :: value

    REAL(wp) :: quadratic_linear, quadratic_constant, discriminant
    REAL(wp) :: partial_a, partial_b, partial_c, partial_d, residual

    quadratic_linear = 1.0_wp + cubic_root
    quadratic_constant = 1.0_wp + cubic_root + cubic_root**2
    discriminant = 4.0_wp * quadratic_constant - quadratic_linear**2

    partial_a = 1.0_wp / (4.0_wp + beta)
    partial_b = 1.0_wp / ((cubic_root - 1.0_wp) *                        &
         (3.0_wp * cubic_root**2 + 2.0_wp * cubic_root + 1.0_wp))
    partial_c = -partial_a - partial_b
    partial_d = -partial_a * (cubic_root + 2.0_wp) -                     &
         partial_b * (2.0_wp * cubic_root + 1.0_wp)
    residual = partial_d - 0.5_wp * partial_c * quadratic_linear

    value = partial_a * LOG(ABS(x - 1.0_wp)) +                            &
         partial_b * LOG(ABS(x - cubic_root)) +                           &
         0.5_wp * partial_c * LOG(x**2 + quadratic_linear*x +             &
         quadratic_constant) + 2.0_wp * residual / SQRT(discriminant) *   &
         ATAN((2.0_wp*x + quadratic_linear) / SQRT(discriminant))

  END FUNCTION lava_cooling_time_primitive

  PURE FUNCTION eval_mass_exchange_thermal_energy(flow_temperature,            &
       solid_density, solid_heat_capacity, deposition_rate, erosion_rate,      &
       carrier_density, carrier_heat_capacity, carrier_loss_rate,              &
       substrate_temperature, substrate_carrier_density,                       &
       carrier_erosion_rate, ambient_temperature, air_heat_capacity,            &
       ambient_air_density, air_entrainment_rate) RESULT(source)

    REAL(wp), INTENT(IN) :: flow_temperature
    REAL(wp), INTENT(IN) :: solid_density(:)
    REAL(wp), INTENT(IN) :: solid_heat_capacity(:)
    REAL(wp), INTENT(IN) :: deposition_rate(:)
    REAL(wp), INTENT(IN) :: erosion_rate(:)
    REAL(wp), INTENT(IN) :: carrier_density
    REAL(wp), INTENT(IN) :: carrier_heat_capacity
    REAL(wp), INTENT(IN) :: carrier_loss_rate
    REAL(wp), INTENT(IN) :: substrate_temperature
    REAL(wp), INTENT(IN) :: substrate_carrier_density
    REAL(wp), INTENT(IN) :: carrier_erosion_rate
    REAL(wp), INTENT(IN) :: ambient_temperature
    REAL(wp), INTENT(IN) :: air_heat_capacity
    REAL(wp), INTENT(IN) :: ambient_air_density
    REAL(wp), INTENT(IN) :: air_entrainment_rate
    REAL(wp) :: source

    source = -flow_temperature * ( SUM( solid_density * solid_heat_capacity *  &
         deposition_rate ) + carrier_density * carrier_heat_capacity *         &
         carrier_loss_rate ) + substrate_temperature * ( SUM( solid_density *  &
         solid_heat_capacity * erosion_rate ) + substrate_carrier_density *    &
         carrier_heat_capacity * carrier_erosion_rate ) + ambient_temperature &
         * air_heat_capacity * ambient_air_density * air_entrainment_rate

  END FUNCTION eval_mass_exchange_thermal_energy

  PURE SUBROUTINE eval_vaporization_sources(flow_temperature, solid_density,    &
       solid_heat_capacity, deposition_rate, steam_fraction, gas_heat_capacity,&
       vapor_mass_rate, vapor_energy_source)

    REAL(wp), INTENT(IN) :: flow_temperature
    REAL(wp), INTENT(IN) :: solid_density(:)
    REAL(wp), INTENT(IN) :: solid_heat_capacity(:)
    REAL(wp), INTENT(IN) :: deposition_rate(:)
    REAL(wp), INTENT(IN) :: steam_fraction
    REAL(wp), INTENT(IN) :: gas_heat_capacity
    REAL(wp), INTENT(OUT) :: vapor_mass_rate
    REAL(wp), INTENT(OUT) :: vapor_energy_source

    REAL(wp), PARAMETER :: liquid_temperature = 290.0_wp
    REAL(wp), PARAMETER :: boiling_temperature = 373.15_wp
    REAL(wp), PARAMETER :: latent_heat = 2264705.0_wp
    REAL(wp), PARAMETER :: liquid_heat_capacity = 4184.0_wp

    vapor_mass_rate = steam_fraction * MAX(0.0_wp, flow_temperature -         &
         boiling_temperature) * SUM(solid_density * solid_heat_capacity *      &
         deposition_rate) / (liquid_heat_capacity * (boiling_temperature -    &
         liquid_temperature) + latent_heat)
    vapor_energy_source = vapor_mass_rate * gas_heat_capacity *                &
         boiling_temperature

  END SUBROUTINE eval_vaporization_sources

END MODULE thermal_sources_2d
