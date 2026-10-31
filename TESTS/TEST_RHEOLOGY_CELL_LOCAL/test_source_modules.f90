PROGRAM test_source_modules

  USE parameters_2d, ONLY : wp
  USE rheology_lahar_2d, ONLY : alpha1_coeff, beta1, eval_lahar_fluid_viscosity
  USE state_2d, ONLY : correct_friction_momentum
  USE mass_exchange_vaporization_2d, ONLY : eval_vaporization_mass_rate, &
       eval_vaporization_energy_rate
  USE mass_exchange_deposition_2d, ONLY : eval_hindered_settling_deposition_rate, &
       eval_lahar_deposition_rate, eval_settling_fluid_properties, eval_solid_deposition_rate
  USE mass_exchange_erosion_2d, ONLY : eval_erosion_rates
  USE mass_exchange_entrainment_2d, ONLY : eval_air_entrainment_rate
  USE thermal_cooling_2d, ONLY : eval_thermal_cooling_implicit

  IMPLICIT NONE

  REAL(wp) :: actual
  REAL(wp) :: expected
  REAL(wp) :: expected_derivative
  REAL(wp) :: temperature_perturbation
  REAL(wp) :: temperature_points(2), finite_difference_step
  REAL(wp) :: finite_difference_derivative, complex_derivative
  COMPLEX(wp) :: complex_actual
  REAL(wp) :: vapor_mass_rate
  REAL(wp) :: vapor_energy_rate
  REAL(wp) :: solid_density(2)
  REAL(wp) :: solid_heat_capacity(2)
  REAL(wp) :: deposition_rate(2)
  REAL(wp) :: erosion_rate(2)
  REAL(wp) :: momentum_before_friction(2)
  REAL(wp) :: momentum_after_friction(2)
  REAL(wp) :: solid_erosion_rate(2), continuous_erosion_rate
  REAL(wp) :: carrier_density, inverse_kinematic_viscosity, expected_viscosity
  INTEGER :: i

  solid_density = [2500.0_wp, 2700.0_wp]
  solid_heat_capacity = [800.0_wp, 900.0_wp]
  deposition_rate = [0.01_wp, 0.02_wp]
  erosion_rate = [0.003_wp, 0.004_wp]

  actual = eval_lahar_deposition_rate(0.2_wp, 0.5_wp, 2.0_wp, 0.5_wp)
  CALL assert_close(actual, 0.1_wp, 'lahar settling-limited deposition')
  actual = eval_lahar_deposition_rate(0.2_wp, 0.5_wp, 0.02_wp, 0.5_wp)
  CALL assert_close(actual, 0.008_wp, 'lahar available-mass deposition cap')
  actual = eval_hindered_settling_deposition_rate(0.2_wp, 0.4_wp, 1.0_wp,   &
       1.0_wp, 0.5_wp, 1.0_wp, 1.0_wp)
  CALL assert_close(actual, 0.04_wp, 'hindered-settling deposition')
  actual = eval_hindered_settling_deposition_rate(0.2_wp, 0.4_wp, 0.1_wp,  &
       1.0_wp, 0.5_wp, 1.0_wp, 1.0_wp)
  CALL assert_close(actual, 0.02_wp, 'hindered-settling available-mass cap')
  actual = eval_solid_deposition_rate('lahar', 0.2_wp, 0.5_wp,    &
       2.0_wp, 0.5_wp, 0.5_wp, 0.6_wp, 4.65_wp)
  CALL assert_close(actual, 0.1_wp, 'lahar deposition dispatch')
  actual = eval_solid_deposition_rate('none', 0.2_wp, 0.4_wp,     &
       0.1_wp, 1.0_wp, 0.5_wp, 1.0_wp, 1.0_wp)
  CALL assert_close(actual, 0.02_wp, 'hindered-settling deposition dispatch')

  alpha1_coeff = 1.0_wp
  beta1 = 2.0_wp
  expected_viscosity = eval_lahar_fluid_viscosity(293.15_wp, 0.3_wp)
  CALL eval_settling_fluid_properties('lahar', 293.15_wp,         &
       0.3_wp, 1800.0_wp, 1000.0_wp, 1.0E-6_wp, carrier_density,            &
       inverse_kinematic_viscosity)
  CALL assert_close(carrier_density, 1800.0_wp, 'lahar settling density')
  CALL assert_close(inverse_kinematic_viscosity,                            &
       1800.0_wp/expected_viscosity, 'lahar settling viscosity dispatch')
  temperature_points=[280.0_wp,320.0_wp]
  temperature_perturbation=1.0E-20_wp
  finite_difference_step=1.0E-4_wp
  DO i=1,SIZE(temperature_points)
     expected_viscosity=eval_lahar_fluid_viscosity(temperature_points(i),0.3_wp)
     complex_actual=eval_lahar_fluid_viscosity( &
          CMPLX(temperature_points(i),temperature_perturbation,wp), &
          CMPLX(0.3_wp,0.0_wp,wp))
     CALL assert_close(REAL(complex_actual,wp),expected_viscosity, &
          'lahar real/complex viscosity value')
     complex_derivative=AIMAG(complex_actual)/temperature_perturbation
     finite_difference_derivative=( &
          eval_lahar_fluid_viscosity(temperature_points(i)+finite_difference_step,0.3_wp)- &
          eval_lahar_fluid_viscosity(temperature_points(i)-finite_difference_step,0.3_wp))/ &
          (2.0_wp*finite_difference_step)
     IF (ABS(complex_derivative-finite_difference_derivative) .GT. &
          1.0E-6_wp*MAX(1.0_wp,ABS(finite_difference_derivative))) THEN
        WRITE(*,*) 'FAIL: lahar complex-step derivative',i, &
             complex_derivative,finite_difference_derivative
        ERROR STOP 1
     END IF
  END DO
  CALL eval_settling_fluid_properties('none', 293.15_wp,           &
       0.3_wp, 1800.0_wp, 1000.0_wp, 2.0E-6_wp, carrier_density,            &
       inverse_kinematic_viscosity)
  CALL assert_close(carrier_density, 1000.0_wp, 'carrier settling density')
  CALL assert_close(inverse_kinematic_viscosity, 5.0E5_wp,                  &
       'carrier settling viscosity')
  CALL eval_erosion_rates(2.0_wp, 3.0_wp, 0.2_wp, 0.5_wp, 0.4_wp,          &
       [0.25_wp, 0.75_wp], [100.0_wp, 100.0_wp], 0.5_wp, solid_erosion_rate, &
       continuous_erosion_rate)
  CALL assert_close(solid_erosion_rate(1), 0.36_wp, 'first solid erosion rate')
  CALL assert_close(solid_erosion_rate(2), 1.08_wp, 'second solid erosion rate')
  CALL assert_close(continuous_erosion_rate, 0.96_wp,                        &
       'carrier erosion rate')
  CALL eval_erosion_rates(2.0_wp, 3.0_wp, 0.2_wp, 0.5_wp, 0.4_wp,          &
       [0.25_wp, 0.75_wp], [0.1_wp, 0.2_wp], 0.5_wp, solid_erosion_rate,    &
       continuous_erosion_rate)
  CALL assert_close(solid_erosion_rate(1), 0.2_wp, 'first erosion inventory cap')
  CALL assert_close(solid_erosion_rate(2), 0.4_wp, 'second erosion inventory cap')
  CALL assert_close(continuous_erosion_rate, 0.4_wp,                        &
       'capped carrier erosion rate')
  ! Carrier loss and packing are covered by the D18.MX production and limiter tests.
  actual = eval_air_entrainment_rate(0.25_wp, 25.0_wp)
  expected = 0.075_wp / SQRT(1.0_wp + 718.0_wp * 0.25_wp**2.4_wp) * 5.0_wp
  CALL assert_close(actual, expected, 'air entrainment rate')

  momentum_before_friction = [3.0_wp, 4.0_wp]
  momentum_after_friction = [1.8_wp, 2.4_wp]
  CALL correct_friction_momentum(momentum_before_friction, momentum_after_friction)
  CALL assert_close(momentum_after_friction(1), 1.8_wp,                    &
       'friction correction preserves x direction')
  CALL assert_close(momentum_after_friction(2), 2.4_wp,                    &
       'friction correction preserves y direction')

  momentum_after_friction = [6.0_wp, 8.0_wp]
  CALL correct_friction_momentum(momentum_before_friction, momentum_after_friction)
  CALL assert_close(momentum_after_friction(1), 3.0_wp,                    &
       'friction correction caps momentum increase x')
  CALL assert_close(momentum_after_friction(2), 4.0_wp,                    &
       'friction correction caps momentum increase y')

  momentum_after_friction = [-1.0_wp, 3.0_wp]
  CALL correct_friction_momentum(momentum_before_friction, momentum_after_friction)
  CALL assert_close(momentum_after_friction(1), 1.08_wp,                   &
       'friction correction projects x without rotation')
  CALL assert_close(momentum_after_friction(2), 1.44_wp,                   &
       'friction correction projects y without rotation')

  momentum_after_friction = [-3.0_wp, -4.0_wp]
  CALL correct_friction_momentum(momentum_before_friction, momentum_after_friction)
  CALL assert_close(momentum_after_friction(1), 0.0_wp,                    &
       'friction correction stops reversed x momentum')
  CALL assert_close(momentum_after_friction(2), 0.0_wp,                    &
       'friction correction stops reversed y momentum')

  momentum_before_friction = 0.0_wp
  momentum_after_friction = [2.0_wp, 3.0_wp]
  CALL correct_friction_momentum(momentum_before_friction, momentum_after_friction)
  CALL assert_close(momentum_after_friction(1), 0.0_wp,                    &
       'friction correction keeps zero x momentum at rest')
  CALL assert_close(momentum_after_friction(2), 0.0_wp,                    &
       'friction correction keeps zero y momentum at rest')

  temperature_perturbation = 1.0E-20_wp
  complex_actual = test_cooling_implicit(     &
       CMPLX(900.0_wp, temperature_perturbation, wp), 3.0E-8_wp,           &
       4.0_wp, 300.0_wp)
  expected = -3.0E-8_wp * (900.0_wp**4 - 300.0_wp**4) -                    &
       4.0_wp * (900.0_wp - 300.0_wp)
  expected_derivative = -4.0_wp * 3.0E-8_wp * 900.0_wp**3 - 4.0_wp
  CALL assert_close(REAL(complex_actual,wp), expected, 'lava cooling value')
  CALL assert_close(AIMAG(complex_actual) / temperature_perturbation,        &
       expected_derivative,                                                 &
       'lava cooling complex-step derivative')

  vapor_mass_rate = eval_vaporization_mass_rate(400.0_wp, solid_density,    &
       solid_heat_capacity, deposition_rate, 0.114_wp)
  vapor_energy_rate = eval_vaporization_energy_rate(vapor_mass_rate,        &
       2010.0_wp)
  expected = 0.114_wp * (400.0_wp - 373.15_wp) *                              &
       SUM(solid_density * solid_heat_capacity * deposition_rate) /            &
       (4184.0_wp * (373.15_wp - 290.0_wp) + 2264705.0_wp)
  CALL assert_close(vapor_mass_rate, expected, 'vaporization mass rate')
  CALL assert_close(vapor_energy_rate, expected * 2010.0_wp * 373.15_wp,    &
       'vaporization energy')
  actual = eval_vaporization_mass_rate(373.15_wp, solid_density,            &
       solid_heat_capacity, deposition_rate, 0.114_wp)
  CALL assert_close(actual, 0.0_wp, 'vaporization threshold mass rate')
  actual = eval_vaporization_energy_rate(actual, 2010.0_wp)
  CALL assert_close(actual, 0.0_wp, 'vaporization threshold energy rate')


CONTAINS
  FUNCTION test_cooling_implicit(temperature,radiative,convective,ambient) RESULT(result)
    USE state_conversion_2d, ONLY : complex_physical_state_type
    COMPLEX(wp), INTENT(IN) :: temperature
    REAL(wp), INTENT(IN) :: radiative,convective,ambient
    COMPLEX(wp) :: result
    TYPE(complex_physical_state_type) :: state
    state%T=temperature
    CALL eval_thermal_cooling_implicit(state,radiative,convective,ambient,result)
  END FUNCTION test_cooling_implicit


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
