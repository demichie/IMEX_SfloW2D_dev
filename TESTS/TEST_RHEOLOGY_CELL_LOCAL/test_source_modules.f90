PROGRAM test_source_modules

  USE parameters_2d, ONLY : wp
  USE pore_pressure_2d, ONLY : eval_pore_pressure_inlet_source,                &
       eval_pore_pressure_mass_exchange_source
  USE thermal_sources_2d, ONLY : eval_inlet_thermal_energy_source,             &
       eval_lava_cooling_source, eval_mass_exchange_thermal_energy,            &
       eval_vaporization_sources

  IMPLICIT NONE

  REAL(wp) :: actual
  REAL(wp) :: expected
  REAL(wp) :: expected_derivative
  REAL(wp) :: temperature_perturbation
  COMPLEX(wp) :: complex_actual
  REAL(wp) :: vapor_mass_rate
  REAL(wp) :: vapor_energy_source
  REAL(wp) :: solid_density(2)
  REAL(wp) :: solid_heat_capacity(2)
  REAL(wp) :: deposition_rate(2)
  REAL(wp) :: erosion_rate(2)

  solid_density = [2500.0_wp, 2700.0_wp]
  solid_heat_capacity = [800.0_wp, 900.0_wp]
  deposition_rate = [0.01_wp, 0.02_wp]
  erosion_rate = [0.003_wp, 0.004_wp]

  actual = eval_inlet_thermal_energy_source(0.5_wp, 0.02_wp, 1200.0_wp,       &
       800.0_wp, 500.0_wp)
  expected = 0.5_wp * 0.02_wp * 1200.0_wp * 800.0_wp * 500.0_wp
  CALL assert_close(actual, expected, 'inlet thermal energy')

  temperature_perturbation = 1.0E-20_wp
  complex_actual = eval_lava_cooling_source(                              &
       CMPLX(900.0_wp, temperature_perturbation, wp), 3.0E-8_wp,           &
       4.0_wp, 300.0_wp)
  expected = -3.0E-8_wp * (900.0_wp**4 - 300.0_wp**4) -                    &
       4.0_wp * (900.0_wp - 300.0_wp)
  expected_derivative = -4.0_wp * 3.0E-8_wp * 900.0_wp**3 - 4.0_wp
  CALL assert_close(REAL(complex_actual,wp), expected, 'lava cooling value')
  CALL assert_close(AIMAG(complex_actual) / temperature_perturbation,        &
       expected_derivative,                                                 &
       'lava cooling complex-step derivative')

  actual = eval_mass_exchange_thermal_energy(650.0_wp, solid_density,          &
       solid_heat_capacity, deposition_rate, erosion_rate, 1.2_wp, 1005.0_wp,  &
       0.004_wp, 300.0_wp, 1.1_wp, 0.005_wp, 290.0_wp, 1005.0_wp, 1.18_wp,    &
       0.006_wp)
  expected = -650.0_wp * (SUM(solid_density * solid_heat_capacity *           &
       deposition_rate) + 1.2_wp * 1005.0_wp * 0.004_wp) + 300.0_wp *         &
       (SUM(solid_density * solid_heat_capacity * erosion_rate) +              &
       1.1_wp * 1005.0_wp * 0.005_wp) + 290.0_wp * 1005.0_wp * 1.18_wp *      &
       0.006_wp
  CALL assert_close(actual, expected, 'mass-exchange thermal energy')

  CALL eval_vaporization_sources(400.0_wp, solid_density, solid_heat_capacity, &
       deposition_rate, 0.114_wp, 2010.0_wp, vapor_mass_rate,                  &
       vapor_energy_source)
  expected = 0.114_wp * (400.0_wp - 373.15_wp) *                              &
       SUM(solid_density * solid_heat_capacity * deposition_rate) /            &
       (4184.0_wp * (373.15_wp - 290.0_wp) + 2264705.0_wp)
  CALL assert_close(vapor_mass_rate, expected, 'vaporization mass rate')
  CALL assert_close(vapor_energy_source, expected * 2010.0_wp * 373.15_wp,    &
       'vaporization energy')

  actual = eval_pore_pressure_inlet_source(0.4_wp, 1.2_wp, 0.03_wp, 0.02_wp,  &
       900.0_wp, 0.7_wp, 5000.0_wp)
  expected = 0.4_wp * (1.2_wp * 0.03_wp * 0.02_wp * 900.0_wp * 0.7_wp +       &
       5000.0_wp * 0.02_wp * 900.0_wp)
  CALL assert_close(actual, expected, 'pore-pressure inlet source')

  actual = eval_pore_pressure_mass_exchange_source(-2.5_wp, 1200.0_wp)
  CALL assert_close(actual, -3000.0_wp, 'pore-pressure mass-exchange source')

CONTAINS

  SUBROUTINE assert_close(value, reference, label)

    REAL(wp), INTENT(IN) :: value
    REAL(wp), INTENT(IN) :: reference
    CHARACTER(LEN=*), INTENT(IN) :: label

    IF (ABS(value-reference) .GT. 1.0E-13_wp * MAX(1.0_wp,ABS(reference))) THEN
       WRITE(*,*) 'FAIL: ', label, value, reference
       ERROR STOP 1
    END IF

  END SUBROUTINE assert_close

END PROGRAM test_source_modules
