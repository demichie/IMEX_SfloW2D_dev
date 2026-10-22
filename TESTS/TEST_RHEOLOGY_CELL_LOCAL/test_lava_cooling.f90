PROGRAM test_lava_cooling

  USE parameters_2d, ONLY : wp
  USE thermal_sources_2d, ONLY : eval_lava_cooling_source,                   &
       advance_lava_cooling_temperature, solve_lava_analytic_temperature,     &
       solve_lava_cooling_stage_temperature

  IMPLICIT NONE

  REAL(wp), PARAMETER :: ambient = 300.0_wp
  REAL(wp), PARAMETER :: initial_temperature = 1200.0_wp
  REAL(wp), PARAMETER :: areal_heat_capacity = 1.0E6_wp
  REAL(wp), PARAMETER :: radiative_coefficient = 5.103E-8_wp
  REAL(wp), PARAMETER :: convective_coefficient = 9.0_wp
  REAL(wp) :: after_first, after_second, after_full, after_short
  REAL(wp) :: expected_derivative, measured_derivative, expected_temperature
  REAL(wp) :: stage_temperature, stage_residual, tau

  CALL advance_lava_cooling_temperature(initial_temperature, 10.0_wp,       &
       areal_heat_capacity, radiative_coefficient,                           &
       convective_coefficient, ambient, after_first)
  CALL advance_lava_cooling_temperature(after_first, 15.0_wp,                &
       areal_heat_capacity, radiative_coefficient,                           &
       convective_coefficient, ambient, after_second)
  CALL advance_lava_cooling_temperature(initial_temperature, 25.0_wp,        &
       areal_heat_capacity, radiative_coefficient,                           &
       convective_coefficient, ambient, after_full)
  CALL assert_close(after_second, after_full, 2.0E-11_wp,                    &
       'combined radiative-convective flow violates the semigroup property')
  CALL solve_lava_analytic_temperature(initial_temperature, 25.0_wp,         &
       areal_heat_capacity, radiative_coefficient, convective_coefficient,    &
       ambient, expected_temperature)
  CALL assert_close(after_full, expected_temperature, 0.0_wp,                &
       'legacy and canonical exact-flow names return different temperatures')

  CALL advance_lava_cooling_temperature(initial_temperature, 10.0_wp,       &
       areal_heat_capacity, radiative_coefficient, 0.0_wp, ambient,          &
       after_first)
  CALL advance_lava_cooling_temperature(after_first, 15.0_wp,                &
       areal_heat_capacity, radiative_coefficient, 0.0_wp, ambient,          &
       after_second)
  CALL advance_lava_cooling_temperature(initial_temperature, 25.0_wp,        &
       areal_heat_capacity, radiative_coefficient, 0.0_wp, ambient,          &
       after_full)
  CALL assert_close(after_second, after_full, 2.0E-11_wp,                    &
       'radiation-only flow violates the semigroup property')

  CALL advance_lava_cooling_temperature(initial_temperature, 0.01_wp,        &
       areal_heat_capacity, radiative_coefficient,                           &
       convective_coefficient, ambient, after_short)
  expected_derivative = REAL(eval_lava_cooling_source(                       &
       CMPLX(initial_temperature, 0.0_wp, wp), radiative_coefficient,         &
       convective_coefficient, ambient), wp) / areal_heat_capacity
  measured_derivative = (after_short - initial_temperature) / 0.01_wp
  CALL assert_close(measured_derivative, expected_derivative, 2.0E-5_wp,     &
       'combined cooling flow does not match the cooling ODE derivative')

  CALL advance_lava_cooling_temperature(initial_temperature, 25.0_wp,        &
       areal_heat_capacity, 0.0_wp, convective_coefficient, ambient,         &
       after_full)
  expected_temperature = ambient + (initial_temperature - ambient) *         &
       EXP(-convective_coefficient * 25.0_wp / areal_heat_capacity)
  CALL assert_close(after_full, expected_temperature, 2.0E-14_wp,            &
       'convection-only cooling does not match its exact exponential solution')

  CALL advance_lava_cooling_temperature(ambient - 100.0_wp, 25.0_wp,        &
       areal_heat_capacity, radiative_coefficient,                           &
       convective_coefficient, ambient, after_full)
  IF (after_full .LE. ambient - 100.0_wp .OR. after_full .GE. ambient) THEN
     ERROR STOP 'sub-ambient cooling flow did not approach ambient temperature'
  END IF

  CALL advance_lava_cooling_temperature(initial_temperature, 10.0_wp,       &
       areal_heat_capacity, 0.0_wp, 0.0_wp, ambient, after_full)
  CALL assert_close(after_full, initial_temperature, 0.0_wp,                  &
       'zero cooling coefficients changed the temperature')

  tau = 500.0_wp
  CALL solve_lava_cooling_stage_temperature(initial_temperature, tau,        &
       areal_heat_capacity, radiative_coefficient,                           &
       convective_coefficient, ambient, stage_temperature)
  stage_residual = (areal_heat_capacity * (stage_temperature -               &
       initial_temperature) - tau * REAL(eval_lava_cooling_source(            &
       CMPLX(stage_temperature,0.0_wp,wp), radiative_coefficient,             &
       convective_coefficient, ambient),wp)) /                               &
       (areal_heat_capacity * initial_temperature)
  CALL assert_close(stage_residual, 0.0_wp, 2.0E-12_wp,                      &
       'combined DIRK cooling stage does not satisfy its algebraic equation')

  CALL solve_lava_cooling_stage_temperature(initial_temperature, tau,        &
       areal_heat_capacity, 0.0_wp, convective_coefficient, ambient,         &
       stage_temperature)
  expected_temperature = (initial_temperature + tau *                      &
       convective_coefficient * ambient / areal_heat_capacity) /             &
       (1.0_wp + tau * convective_coefficient / areal_heat_capacity)
  CALL assert_close(stage_temperature, expected_temperature, 2.0E-14_wp,     &
       'convection-only DIRK cooling stage is not algebraically exact')

  tau = 1.0E12_wp
  CALL solve_lava_cooling_stage_temperature(initial_temperature, tau,        &
       areal_heat_capacity, radiative_coefficient, 0.0_wp, ambient,         &
       stage_temperature)
  stage_residual = (areal_heat_capacity * (stage_temperature -               &
       initial_temperature) - tau * REAL(eval_lava_cooling_source(            &
       CMPLX(stage_temperature,0.0_wp,wp), radiative_coefficient, 0.0_wp,     &
       ambient),wp)) / (areal_heat_capacity * initial_temperature)
  CALL assert_close(stage_residual, 0.0_wp, 1.0E-9_wp,                       &
       'radiation-dominated quartic root does not satisfy the DIRK equation')

  WRITE(*,*) 'All exact lava cooling checks passed.'

CONTAINS

  SUBROUTINE assert_close(actual, expected, relative_tolerance, message)

    REAL(wp), INTENT(IN) :: actual, expected, relative_tolerance
    CHARACTER(LEN=*), INTENT(IN) :: message
    REAL(wp) :: tolerance

    tolerance = relative_tolerance * MAX(1.0_wp, ABS(expected))
    IF (ABS(actual - expected) .GT. tolerance) THEN
       WRITE(*,*) 'FAIL: ', TRIM(message)
       WRITE(*,*) 'actual, expected, tolerance:', actual, expected, tolerance
       ERROR STOP 'exact lava cooling test failed'
    END IF

  END SUBROUTINE assert_close

END PROGRAM test_lava_cooling
