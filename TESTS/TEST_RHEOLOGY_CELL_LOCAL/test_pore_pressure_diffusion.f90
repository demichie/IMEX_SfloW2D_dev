PROGRAM test_pore_pressure_diffusion

  USE parameters_2d, ONLY : wp, n_solid, gas_flag, sutherland_flag
  USE constitutive_parameters_2d, ONLY : f_inhibit_mode, hydraulic_permeability, &
       kin_visc_c, sp_heat_a, sp_gas_const_a, pres, dynamic_permeability_flag
  USE geometry_2d, ONLY : pi_g
  USE pore_pressure_2d, ONLY : eval_pore_pressure_dissipation_implicit,       &
       advance_pore_pressure_diffusion, solve_pore_pressure_dissipation_analytic, &
       solve_frozen_darcy_mass_pressure

  IMPLICIT NONE

  REAL(wp), PARAMETER :: flow_depth = 2.0_wp
  REAL(wp), PARAMETER :: temperature = 350.0_wp
  REAL(wp), PARAMETER :: mixture_density = 1200.0_wp
  REAL(wp), PARAMETER :: initial_pressure = 5000.0_wp
  REAL(wp), PARAMETER :: integration_step = 0.4_wp
  REAL(wp), ALLOCATABLE :: solid_fractions(:)
  REAL(wp) :: after_first, after_second, after_full, after_short
  REAL(wp) :: decay_rate, expected_pressure, measured_derivative
  COMPLEX(wp) :: source
  REAL(wp) :: m_after,p_after,qp_after,gas_lost
  REAL(wp) :: initial_mass,loss_coefficient,analytic_loss
  LOGICAL :: admissible

  n_solid = 2
  ALLOCATE(solid_fractions(n_solid))
  gas_flag = .FALSE.
  sutherland_flag = .FALSE.
  dynamic_permeability_flag = .FALSE.
  f_inhibit_mode = 'OFF'
  hydraulic_permeability = 1.0E-10_wp
  kin_visc_c = 1.5E-5_wp
  sp_heat_a = 1005.0_wp
  sp_gas_const_a = 287.0_wp
  pres = 101325.0_wp
  pi_g = 4.0_wp * ATAN(1.0_wp)
  solid_fractions = 0.2_wp / REAL(n_solid,wp)

  CALL eval_pore_pressure_dissipation_implicit(                                  &
       CMPLX(flow_depth,0.0_wp,wp), CMPLX(temperature,0.0_wp,wp),            &
       CMPLX(mixture_density,0.0_wp,wp),                                     &
       CMPLX(solid_fractions,0.0_wp,wp),                                     &
       CMPLX(initial_pressure,0.0_wp,wp), source)
  decay_rate = -REAL(source,wp) /                                            &
       (mixture_density * flow_depth * initial_pressure)
  expected_pressure = initial_pressure * EXP(-decay_rate * integration_step)

  CALL advance_pore_pressure_diffusion(initial_pressure, integration_step,   &
       flow_depth, temperature, solid_fractions, after_full)
  CALL assert_close(after_full, expected_pressure, 2.0E-14_wp,               &
       'diffusion integrator disagrees with the source-derived decay rate')
  CALL solve_pore_pressure_dissipation_analytic(initial_pressure, integration_step,       &
       flow_depth, temperature, solid_fractions, after_first)
  CALL assert_close(after_first, after_full, 0.0_wp,                          &
       'legacy and canonical exact-flow names return different pressures')

  CALL advance_pore_pressure_diffusion(initial_pressure, 0.15_wp,            &
       flow_depth, temperature, solid_fractions, after_first)
  CALL advance_pore_pressure_diffusion(after_first, 0.25_wp,                  &
       flow_depth, temperature, solid_fractions, after_second)
  CALL assert_close(after_second, after_full, 2.0E-14_wp,                    &
       'frozen-state pore-pressure flow violates the semigroup property')

  CALL advance_pore_pressure_diffusion(initial_pressure, 1.0E-5_wp,          &
       flow_depth, temperature, solid_fractions, after_short)
  measured_derivative = (after_short - initial_pressure) / 1.0E-5_wp
  CALL assert_close(measured_derivative, REAL(source,wp) /                   &
       (mixture_density * flow_depth), 2.0E-5_wp,                            &
       'diffusion integrator initial derivative disagrees with source')

  CALL advance_pore_pressure_diffusion(initial_pressure, 0.0_wp,            &
       flow_depth, temperature, solid_fractions, after_full)
  CALL assert_close(after_full, initial_pressure, 0.0_wp,                    &
       'zero-duration diffusion changed the pressure')

  CALL advance_pore_pressure_diffusion(0.0_wp, integration_step,            &
       flow_depth, temperature, solid_fractions, after_full)
  CALL assert_close(after_full, 0.0_wp, 0.0_wp,                             &
       'zero excess pressure was not preserved')

  CALL advance_pore_pressure_diffusion(-initial_pressure, integration_step,  &
       flow_depth, temperature, solid_fractions, after_full)
  CALL assert_close(after_full, -initial_pressure, 0.0_wp,                    &
       'negative excess pressure should remain unchanged')

  CALL advance_pore_pressure_diffusion(initial_pressure, integration_step,    &
       0.0_wp, temperature, solid_fractions, after_full)
  CALL assert_close(after_full, initial_pressure, 0.0_wp,                    &
       'dry state should remain unchanged')

  solid_fractions = 0.0_wp
  solid_fractions(1) = 1.0_wp
  CALL advance_pore_pressure_diffusion(initial_pressure, integration_step,    &
       flow_depth, temperature, solid_fractions, after_full)
  CALL assert_close(after_full, initial_pressure, 0.0_wp,                    &
       'state without pore volume should remain unchanged')

  solid_fractions = 0.2_wp / REAL(n_solid,wp)
  CALL advance_pore_pressure_diffusion(initial_pressure, HUGE(1.0_wp),       &
       flow_depth, temperature, solid_fractions, after_full)
  CALL assert_close(after_full, 0.0_wp, 0.0_wp,                              &
       'long-time diffusion limit did not decay excess pressure to zero')

  hydraulic_permeability = 0.0_wp
  CALL advance_pore_pressure_diffusion(initial_pressure, integration_step,   &
       flow_depth, temperature, solid_fractions, after_full)
  CALL assert_close(after_full, initial_pressure, 0.0_wp,                    &
       'zero permeability changed excess pore pressure')


  ! Exact coupled frozen-coefficient map: p decays exponentially, and
  ! the integral of p controls the carrier-gas mass withdrawn.
  initial_mass=4000.0_wp
  loss_coefficient=0.03_wp
  CALL solve_frozen_darcy_mass_pressure(initial_mass,initial_pressure, &
       integration_step,decay_rate,loss_coefficient,.TRUE., &
       m_after,p_after,qp_after,gas_lost,admissible)
  IF (.NOT.admissible) ERROR STOP 'coupled Darcy analytic state rejected'
  analytic_loss=loss_coefficient*initial_pressure* &
       (1.0_wp-EXP(-decay_rate*integration_step))/decay_rate
  CALL assert_close(gas_lost,analytic_loss,1.0e-12_wp, &
       'coupled Darcy gas loss differs from exact pressure integral')
  CALL assert_close(m_after,initial_mass-gas_lost,1.0e-12_wp, &
       'coupled Darcy total mass balance')
  CALL assert_close(p_after,expected_pressure,2.0e-14_wp, &
       'coupled Darcy pressure must follow the same exponential')
  CALL assert_close(qp_after,m_after*p_after,2.0e-14_wp, &
       'conservative pore-pressure product is inconsistent')

  CALL solve_frozen_darcy_mass_pressure(initial_mass,initial_pressure, &
       integration_step,decay_rate,loss_coefficient,.FALSE., &
       m_after,p_after,qp_after,gas_lost,admissible)
  IF (.NOT.admissible) ERROR STOP 'disabled Darcy gas loss rejected'
  CALL assert_close(m_after,initial_mass,0.0_wp, &
       'gas_loss_flag off changed mass')
  CALL assert_close(gas_lost,0.0_wp,0.0_wp, &
       'gas_loss_flag off removed gas')
  CALL assert_close(p_after,expected_pressure,2.0e-14_wp, &
       'gas_loss_flag off disabled pressure dissipation')

  CALL solve_frozen_darcy_mass_pressure(initial_mass,initial_pressure, &
       integration_step,0.0_wp,loss_coefficient,.TRUE., &
       m_after,p_after,qp_after,gas_lost,admissible)
  IF (.NOT.admissible) ERROR STOP 'zero diffusion coupled Darcy state rejected'
  CALL assert_close(gas_lost,loss_coefficient*initial_pressure*integration_step, &
       2.0e-14_wp,'zero diffusion limit of gas loss')
  CALL assert_close(p_after,initial_pressure,0.0_wp, &
       'zero diffusion changed the pressure')

  CALL solve_frozen_darcy_mass_pressure(1.0_wp,initial_pressure, &
       integration_step,0.0_wp,loss_coefficient,.TRUE., &
       m_after,p_after,qp_after,gas_lost,admissible)
  IF (admissible) ERROR STOP 'Darcy accepted loss exceeding available mass'

  WRITE(*,*) 'All exact pore-pressure diffusion and coupled Darcy checks passed.'

CONTAINS

  SUBROUTINE assert_close(actual, expected, relative_tolerance, message)

    REAL(wp), INTENT(IN) :: actual, expected, relative_tolerance
    CHARACTER(LEN=*), INTENT(IN) :: message
    REAL(wp) :: tolerance

    tolerance = relative_tolerance * MAX(1.0_wp, ABS(expected))
    IF (ABS(actual - expected) .GT. tolerance) THEN
       WRITE(*,*) 'FAIL: ', TRIM(message)
       WRITE(*,*) 'actual, expected, tolerance:', actual, expected, tolerance
       ERROR STOP 'exact pore-pressure diffusion test failed'
    END IF

  END SUBROUTINE assert_close

END PROGRAM test_pore_pressure_diffusion
