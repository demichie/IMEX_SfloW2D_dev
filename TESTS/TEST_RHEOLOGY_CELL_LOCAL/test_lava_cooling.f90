PROGRAM test_lava_cooling

  USE parameters_2d, ONLY : wp
  USE thermal_cooling_2d, ONLY : solve_thermal_cooling_analytic, &
       eval_thermal_cooling_implicit

  IMPLICIT NONE

  REAL(wp), PARAMETER :: ambient = 300.0_wp
  REAL(wp), PARAMETER :: initial_temperature = 1200.0_wp
  REAL(wp), PARAMETER :: areal_heat_capacity = 1.0E6_wp
  REAL(wp), PARAMETER :: radiative_coefficient = 5.103E-8_wp
  REAL(wp), PARAMETER :: convective_coefficient = 9.0_wp
  REAL(wp) :: after_first, after_second, after_full, after_short
  REAL(wp) :: expected_derivative, measured_derivative, expected_temperature

  CALL test_cooling_analytic(initial_temperature, 0.0_wp,        &
       areal_heat_capacity, radiative_coefficient, convective_coefficient,   &
       ambient, after_full)
  CALL assert_close(after_full, initial_temperature, 0.0_wp,               &
       'zero-duration cooling changed the temperature')

  CALL test_cooling_analytic(initial_temperature, 10.0_wp,       &
       areal_heat_capacity, radiative_coefficient,                           &
       convective_coefficient, ambient, after_first)
  CALL test_cooling_analytic(after_first, 15.0_wp,                &
       areal_heat_capacity, radiative_coefficient,                           &
       convective_coefficient, ambient, after_second)
  CALL test_cooling_analytic(initial_temperature, 25.0_wp,        &
       areal_heat_capacity, radiative_coefficient,                           &
       convective_coefficient, ambient, after_full)
  CALL assert_close(after_second, after_full, 2.0E-11_wp,                    &
       'combined radiative-convective flow violates the semigroup property')
  CALL test_cooling_analytic(initial_temperature, 25.0_wp,         &
       areal_heat_capacity, radiative_coefficient, convective_coefficient,    &
       ambient, expected_temperature)
  CALL assert_close(after_full, expected_temperature, 0.0_wp,                &
       'legacy and canonical exact-flow names return different temperatures')

  CALL test_cooling_analytic(initial_temperature, 10.0_wp,       &
       areal_heat_capacity, radiative_coefficient, 0.0_wp, ambient,          &
       after_first)
  CALL test_cooling_analytic(after_first, 15.0_wp,                &
       areal_heat_capacity, radiative_coefficient, 0.0_wp, ambient,          &
       after_second)
  CALL test_cooling_analytic(initial_temperature, 25.0_wp,        &
       areal_heat_capacity, radiative_coefficient, 0.0_wp, ambient,          &
       after_full)
  CALL assert_close(after_second, after_full, 2.0E-11_wp,                    &
       'radiation-only flow violates the semigroup property')

  CALL test_cooling_analytic(initial_temperature, 0.01_wp,        &
       areal_heat_capacity, radiative_coefficient,                           &
       convective_coefficient, ambient, after_short)
  expected_derivative = REAL(test_cooling_implicit(                  &
       CMPLX(initial_temperature, 0.0_wp, wp), radiative_coefficient,         &
       convective_coefficient, ambient), wp) / areal_heat_capacity
  measured_derivative = (after_short - initial_temperature) / 0.01_wp
  CALL assert_close(measured_derivative, expected_derivative, 2.0E-5_wp,     &
       'combined cooling flow does not match the cooling ODE derivative')

  CALL test_cooling_analytic(initial_temperature, 25.0_wp,        &
       areal_heat_capacity, 0.0_wp, convective_coefficient, ambient,         &
       after_full)
  expected_temperature = ambient + (initial_temperature - ambient) *         &
       EXP(-convective_coefficient * 25.0_wp / areal_heat_capacity)
  CALL assert_close(after_full, expected_temperature, 2.0E-14_wp,            &
       'convection-only cooling does not match its exact exponential solution')

  CALL test_cooling_analytic(ambient - 100.0_wp, 25.0_wp,        &
       areal_heat_capacity, radiative_coefficient,                           &
       convective_coefficient, ambient, after_full)
  IF (after_full .LE. ambient - 100.0_wp .OR. after_full .GE. ambient) THEN
     ERROR STOP 'sub-ambient cooling flow did not approach ambient temperature'
  END IF

  CALL test_cooling_analytic(ambient, 25.0_wp,                   &
       areal_heat_capacity, radiative_coefficient, convective_coefficient,   &
       ambient, after_full)
  CALL assert_close(after_full, ambient, 0.0_wp,                             &
       'ambient equilibrium was not preserved')

  CALL test_cooling_analytic(initial_temperature, 1.0E9_wp,      &
       areal_heat_capacity, 0.0_wp, convective_coefficient, ambient,        &
       after_full)
  CALL assert_close(after_full, ambient, 2.0E-14_wp,                         &
       'long-time convection limit did not reach ambient temperature')

  CALL test_cooling_analytic(initial_temperature, 10.0_wp,       &
       areal_heat_capacity, 0.0_wp, 0.0_wp, ambient, after_full)
  CALL assert_close(after_full, initial_temperature, 0.0_wp,                  &
       'zero cooling coefficients changed the temperature')

  WRITE(*,*) 'All exact lava cooling checks passed.'

CONTAINS
  SUBROUTINE test_cooling_analytic(temperature,step,capacity,radiative,convective,ambient,result)
    USE state_conversion_2d, ONLY : physical_state_type
    REAL(wp), INTENT(IN) :: temperature,step,capacity,radiative,convective,ambient
    REAL(wp), INTENT(OUT) :: result
    TYPE(physical_state_type) :: state
    state%T=temperature
    CALL solve_thermal_cooling_analytic(state,step,capacity,radiative,convective,ambient,result)
  END SUBROUTINE test_cooling_analytic
  FUNCTION test_cooling_implicit(temperature,radiative,convective,ambient) RESULT(result)
    USE state_conversion_2d, ONLY : complex_physical_state_type
    COMPLEX(wp), INTENT(IN) :: temperature
    REAL(wp), INTENT(IN) :: radiative,convective,ambient
    COMPLEX(wp) :: result
    TYPE(complex_physical_state_type) :: state
    state%T=temperature
    CALL eval_thermal_cooling_implicit(state,radiative,convective,ambient,result)
  END FUNCTION test_cooling_implicit


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
