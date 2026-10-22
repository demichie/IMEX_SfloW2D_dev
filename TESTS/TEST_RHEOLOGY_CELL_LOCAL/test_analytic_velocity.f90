PROGRAM test_analytic_velocity

  USE parameters_2d, ONLY : wp, n_solid
  USE constitutive_parameters_2d, ONLY : grav, inv_grav, diam_s, rho_s
  USE rheology_parameters_2d, ONLY : I_0, mu_2, mu_s, muI_inf, mu_0, mu_inf, &
       Fr_0, U_w, xi, T_ref
  USE rheology_analytic_utils_2d, ONLY : lambert_w0, solve_quadratic_in_interval, &
       solve_cubic_in_interval
  USE rheology_voellmy_2d, ONLY : voellmy_mu => mu,                          &
       eval_voellmy_implicit_source, eval_voellmy_semi_impl_source,          &
       solve_voellmy_analytic_velocity
  USE rheology_plastic_2d, ONLY : tau, eval_plastic_implicit_source,         &
       solve_plastic_analytic_velocity
  USE rheology_lava_2d, ONLY : nu_ref, visc_par, tau0,                       &
       eval_lava_implicit_source, eval_lava_semi_impl_source,                &
       solve_lava_analytic_velocity
  USE rheology_lahar_2d, ONLY : alpha2, beta2, alpha1_coeff, beta1, Kappa,   &
       n_td, n_td2, eval_lahar_implicit_source, eval_lahar_semi_impl_source, &
       solve_lahar_analytic_velocity
  USE rheology_empirical_2d, ONLY : eval_empirical_implicit_source,          &
       solve_empirical_analytic_velocity
  USE rheology_quadratic_2d, ONLY : friction_factor,                        &
       eval_quadratic_implicit_source, solve_quadratic_analytic_velocity
  USE rheology_froude_2d, ONLY : eval_froude_source, solve_froude_analytic_velocity
  USE rheology_lucas_2d, ONLY : eval_lucas_source, solve_lucas_analytic_velocity
  USE rheology_mui_2d, ONLY : eval_mui_source, solve_mui_analytic_velocity
  USE rheology_mui_voellmy_2d, ONLY : eval_mui_voellmy_source,              &
       eval_mui_voellmy_implicit_source, solve_mui_voellmy_analytic_velocity

  IMPLICIT NONE

  REAL(wp), PARAMETER :: integration_step = 0.01_wp
  REAL(wp), PARAMETER :: flow_depth = 1.8_wp
  REAL(wp), PARAMETER :: mixture_density = 1800.0_wp
  REAL(wp), PARAMETER :: reduced_gravity = 5.5_wp
  REAL(wp), PARAMETER :: gravity_coefficient = 1.2_wp
  REAL(wp), PARAMETER :: speed_ratio = 1.15_wp
  REAL(wp), PARAMETER :: temperature = 293.15_wp
  REAL(wp) :: velocity_before(2), velocity_after(2), solid_fraction(1)
  REAL(wp) :: speed_after, horizontal_speed, source_x, source_y
  REAL(wp) :: source(2), lambert_value, expected_speed, linear_drag, yield_drag
  REAL(wp) :: root
  REAL(wp) :: stochastic_value
  COMPLEX(wp) :: complex_source_x, complex_source_y
  LOGICAL :: lambert_success, root_found

  n_solid = 1
  ALLOCATE(diam_s(1), rho_s(1))
  diam_s(1) = 0.01_wp
  rho_s(1) = 2500.0_wp
  grav = 9.81_wp
  inv_grav = 1.0_wp / grav
  solid_fraction = 0.3_wp
  velocity_before = [3.0_wp, 4.0_wp]
  stochastic_value = 0.0_wp

  CALL solve_quadratic_in_interval(1.0_wp, 0.0_wp, -4.0_wp, 0.0_wp, 3.0_wp, &
       root, root_found)
  IF (.NOT. root_found) ERROR STOP 'Quadratic root selection failed'
  CALL assert_close(root, 2.0_wp, 'quadratic largest root')
  CALL solve_cubic_in_interval(1.0_wp, 0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp,     &
       2.0_wp, root, root_found)
  IF (.NOT. root_found) ERROR STOP 'Cubic three-real-root selection failed'
  CALL assert_close(root, 1.0_wp, 'cubic largest admissible root')
  CALL solve_cubic_in_interval(1.0_wp, 0.0_wp, 0.0_wp, -8.0_wp, 0.0_wp,     &
       3.0_wp, root, root_found)
  IF (.NOT. root_found) ERROR STOP 'Cubic one-real-root selection failed'
  CALL assert_close(root, 2.0_wp, 'cubic single real root')

  tau = 1.0_wp
  CALL solve_plastic_analytic_velocity(integration_step, flow_depth,       &
       speed_ratio, velocity_before, velocity_after)
  expected_speed = MAX(0.0_wp, SQRT(SUM(velocity_before**2)) -              &
       integration_step*tau/(flow_depth*speed_ratio))
  CALL assert_close(SQRT(SUM(velocity_after**2)), expected_speed,           &
       'plastic exact constant-drag flow')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  CALL eval_plastic_implicit_source(CMPLX(velocity_after(1),0.0_wp,wp),     &
       CMPLX(velocity_after(2),0.0_wp,wp), CMPLX(speed_after,0.0_wp,wp),     &
       CMPLX(mixture_density,0.0_wp,wp), .TRUE., complex_source_x,          &
       complex_source_y)
  CALL assert_stage(velocity_before, velocity_after, integration_step,     &
       mixture_density*flow_depth, [REAL(complex_source_x,wp),              &
       REAL(complex_source_y,wp)], 'plastic')

  nu_ref = 0.2_wp
  visc_par = 0.01_wp
  tau0 = 100.0_wp
  T_ref = 300.0_wp
  CALL solve_lava_analytic_velocity(integration_step, flow_depth,           &
       temperature, mixture_density, speed_ratio, velocity_before, velocity_after)
  linear_drag = 3.0_wp*nu_ref/flow_depth**2 *                               &
       EXP(-visc_par*(temperature-T_ref))
  yield_drag = tau0/(mixture_density*flow_depth*speed_ratio)
  expected_speed = MAX(0.0_wp, (SQRT(SUM(velocity_before**2)) +             &
       yield_drag/linear_drag)*EXP(-linear_drag*integration_step) -         &
       yield_drag/linear_drag)
  CALL assert_close(SQRT(SUM(velocity_after**2)), expected_speed,           &
       'lava exact affine-drag flow')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  CALL eval_lava_implicit_source(CMPLX(flow_depth,0.0_wp,wp),                &
       CMPLX(temperature,0.0_wp,wp), CMPLX(mixture_density,0.0_wp,wp),       &
       CMPLX(velocity_after(1),0.0_wp,wp), CMPLX(velocity_after(2),0.0_wp,wp), &
       CMPLX(speed_after,0.0_wp,wp), .TRUE., complex_source_x, complex_source_y)
  CALL eval_lava_semi_impl_source(velocity_after(1), velocity_after(2),      &
       speed_after, source_x, source_y)
  source = [REAL(complex_source_x,wp)+source_x,                              &
       REAL(complex_source_y,wp)+source_y]
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, source, 'lava')

  alpha2 = 100.0_wp
  beta2 = 2.0_wp
  alpha1_coeff = 1.0_wp
  beta1 = 2.0_wp
  Kappa = 1.0_wp
  n_td = 0.1_wp
  n_td2 = 0.01_wp
  CALL solve_lahar_analytic_velocity(integration_step, flow_depth,          &
       temperature, SUM(solid_fraction), mixture_density, speed_ratio,       &
       velocity_before, velocity_after)
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  CALL eval_lahar_implicit_source(CMPLX(flow_depth,0.0_wp,wp),               &
       CMPLX(velocity_after(1),0.0_wp,wp), CMPLX(velocity_after(2),0.0_wp,wp), &
       CMPLX(speed_after,0.0_wp,wp), CMPLX(horizontal_speed,0.0_wp,wp),      &
       CMPLX(temperature,0.0_wp,wp), CMPLX(mixture_density,0.0_wp,wp),      &
       CMPLX(1.0_wp/mixture_density,0.0_wp,wp),                             &
       CMPLX(solid_fraction,0.0_wp,wp), .TRUE., complex_source_x,            &
       complex_source_y)
  CALL eval_lahar_semi_impl_source(flow_depth, velocity_after(1),             &
       velocity_after(2), horizontal_speed, mixture_density,                 &
       SUM(solid_fraction), .TRUE., source_x, source_y)
  source = [REAL(complex_source_x,wp)+source_x,                              &
       REAL(complex_source_y,wp)+source_y]
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, source, 'lahar')

  CALL solve_empirical_analytic_velocity(integration_step, flow_depth,      &
       speed_ratio, velocity_before, velocity_after)
  linear_drag = 1.0E-3_wp*speed_ratio/                                     &
       (flow_depth*(1.0_wp+10.0_wp*flow_depth))
  expected_speed = SQRT(SUM(velocity_before**2))*                           &
       EXP(-linear_drag*integration_step)
  CALL assert_close(SQRT(SUM(velocity_after**2)), expected_speed,           &
       'empirical exact linear-drag flow')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  CALL eval_empirical_implicit_source(CMPLX(flow_depth,0.0_wp,wp),           &
       CMPLX(velocity_after(1),0.0_wp,wp), CMPLX(velocity_after(2),0.0_wp,wp), &
       CMPLX(horizontal_speed,0.0_wp,wp), CMPLX(speed_after,0.0_wp,wp),     &
       CMPLX(mixture_density,0.0_wp,wp), .TRUE., complex_source_x,           &
       complex_source_y)
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, [REAL(complex_source_x,wp),               &
       REAL(complex_source_y,wp)], 'empirical')

  friction_factor = 0.1_wp
  CALL solve_quadratic_analytic_velocity(integration_step, flow_depth,      &
       speed_ratio, velocity_before, velocity_after)
  linear_drag = friction_factor*speed_ratio**2/flow_depth
  expected_speed = SQRT(SUM(velocity_before**2))/                           &
       (1.0_wp+linear_drag*SQRT(SUM(velocity_before**2))*integration_step)
  CALL assert_close(SQRT(SUM(velocity_after**2)), expected_speed,           &
       'quadratic exact reciprocal-speed flow')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  CALL eval_quadratic_implicit_source(                                      &
       CMPLX(velocity_after(1),0.0_wp,wp), CMPLX(velocity_after(2),0.0_wp,wp), &
       CMPLX(speed_after,0.0_wp,wp), CMPLX(horizontal_speed,0.0_wp,wp),     &
       CMPLX(speed_after**2,0.0_wp,wp), CMPLX(mixture_density,0.0_wp,wp),  &
       .TRUE., complex_source_x, complex_source_y)
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, [REAL(complex_source_x,wp),               &
       REAL(complex_source_y,wp)], 'quadratic')

  voellmy_mu = 0.3_wp
  xi = 500.0_wp
  CALL solve_voellmy_analytic_velocity(integration_step, flow_depth,       &
       speed_ratio, velocity_before, mixture_density, reduced_gravity,        &
       gravity_coefficient, 0.0_wp, .FALSE., 0.0_wp, .FALSE., 0.0_wp,        &
       velocity_after)
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  CALL eval_voellmy_implicit_source(                                        &
       CMPLX(velocity_after(1),0.0_wp,wp), CMPLX(velocity_after(2),0.0_wp,wp), &
       CMPLX(mixture_density,0.0_wp,wp), CMPLX(reduced_gravity,0.0_wp,wp),  &
       CMPLX(speed_after,0.0_wp,wp), CMPLX(speed_after**2,0.0_wp,wp),       &
       .TRUE., .FALSE., 0.0_wp, complex_source_x, complex_source_y)
  CALL eval_voellmy_semi_impl_source(flow_depth, velocity_after(1),          &
       velocity_after(2), speed_after, mixture_density, reduced_gravity,    &
       gravity_coefficient, .FALSE., 0.0_wp, 0.0_wp, 0.0_wp, .FALSE.,        &
       0.0_wp, source_x, source_y)
  source = [REAL(complex_source_x,wp)+source_x,                              &
       REAL(complex_source_y,wp)+source_y]
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, source, 'voellmy')

  CALL lambert_w0(1.0_wp, lambert_value, lambert_success)
  IF (.NOT. lambert_success) ERROR STOP 'Lambert W evaluation failed'
  CALL assert_close(lambert_value*EXP(lambert_value), 1.0_wp, 'Lambert W identity')
  CALL lambert_w0(-0.2_wp, lambert_value, lambert_success)
  IF (.NOT. lambert_success) ERROR STOP 'Lambert W negative-branch evaluation failed'
  CALL assert_close(lambert_value*EXP(lambert_value), -0.2_wp,                &
       'Lambert W negative identity')

  mu_0 = 0.5_wp
  mu_inf = 0.2_wp
  Fr_0 = 0.8_wp
  CALL solve_froude_analytic_velocity(integration_step, flow_depth,        &
       speed_ratio, velocity_before, reduced_gravity, gravity_coefficient,   &
       0.0_wp, .FALSE., stochastic_value, velocity_after)
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  CALL eval_froude_source(flow_depth, velocity_after(1), velocity_after(2), &
       speed_after, horizontal_speed, mixture_density, reduced_gravity,      &
       gravity_coefficient, 0.0_wp, .FALSE., stochastic_value,              &
       source_x, source_y)
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, [source_x, source_y], 'Froude')

  stochastic_value = -2.0_wp
  CALL solve_froude_analytic_velocity(integration_step, flow_depth,        &
       speed_ratio, velocity_before, reduced_gravity, gravity_coefficient,   &
       0.0_wp, .TRUE., stochastic_value, velocity_after)
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  CALL eval_froude_source(flow_depth, velocity_after(1), velocity_after(2), &
       speed_after, horizontal_speed, mixture_density, reduced_gravity,      &
       gravity_coefficient, 0.0_wp, .TRUE., stochastic_value,               &
       source_x, source_y)
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, [source_x, source_y], 'Froude clamp')
  stochastic_value = 0.0_wp

  mu_0 = 0.0_wp
  mu_inf = 0.0_wp
  CALL solve_froude_analytic_velocity(integration_step, flow_depth,        &
       speed_ratio, velocity_before, reduced_gravity, gravity_coefficient,   &
       0.0_wp, .FALSE., stochastic_value, velocity_after)
  CALL assert_close(velocity_after(1), velocity_before(1), 'Froude zero-drag x')
  CALL assert_close(velocity_after(2), velocity_before(2), 'Froude zero-drag y')
  mu_0 = 0.5_wp
  mu_inf = 0.2_wp

  U_w = 2.0_wp
  CALL solve_lucas_analytic_velocity(integration_step, flow_depth,          &
       speed_ratio, velocity_before, reduced_gravity, gravity_coefficient,    &
       0.0_wp, .FALSE., stochastic_value, velocity_after)
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  CALL eval_lucas_source(flow_depth, velocity_after(1), velocity_after(2),  &
       speed_after, horizontal_speed, mixture_density, reduced_gravity,      &
       gravity_coefficient, 0.0_wp, .FALSE., stochastic_value,               &
       source_x, source_y)
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, [source_x, source_y], 'Lucas')

  stochastic_value = -10.0_wp
  CALL solve_lucas_analytic_velocity(integration_step, flow_depth,          &
       speed_ratio, velocity_before, reduced_gravity, gravity_coefficient,    &
       0.0_wp, .TRUE., stochastic_value, velocity_after)
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  CALL eval_lucas_source(flow_depth, velocity_after(1), velocity_after(2),  &
       speed_after, horizontal_speed, mixture_density, reduced_gravity,      &
       gravity_coefficient, 0.0_wp, .TRUE., stochastic_value,               &
       source_x, source_y)
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, [source_x, source_y], 'Lucas stochastic branch')
  stochastic_value = 0.0_wp

  I_0 = 0.3_wp
  mu_s = 0.4_wp
  mu_2 = 0.7_wp
  muI_inf = 0.8_wp
  CALL solve_mui_analytic_velocity(integration_step, flow_depth,            &
       speed_ratio, velocity_before, mixture_density, reduced_gravity,        &
       gravity_coefficient, solid_fraction, 0.0_wp, .FALSE., 0.0_wp,         &
       velocity_after)
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  CALL eval_mui_source(flow_depth, velocity_after(1), velocity_after(2),    &
       speed_after, mixture_density, reduced_gravity, gravity_coefficient,    &
       solid_fraction, 0.0_wp, .FALSE., 0.0_wp, source_x, source_y)
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, [source_x, source_y], 'mu(I)')

  CALL solve_mui_analytic_velocity(integration_step, flow_depth,            &
       speed_ratio, velocity_before, mixture_density, reduced_gravity,        &
       gravity_coefficient, solid_fraction, 0.0_wp, .TRUE., 1000.0_wp,       &
       velocity_after)
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  CALL eval_mui_source(flow_depth, velocity_after(1), velocity_after(2),    &
       speed_after, mixture_density, reduced_gravity, gravity_coefficient,    &
       solid_fraction, 0.0_wp, .TRUE., 1000.0_wp, source_x, source_y)
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, [source_x, source_y], 'mu(I) pore pressure')

  CALL solve_mui_voellmy_analytic_velocity(integration_step, flow_depth,   &
       speed_ratio, velocity_before, mixture_density, reduced_gravity,       &
       gravity_coefficient, solid_fraction, 0.0_wp, .FALSE., 0.0_wp,         &
       .FALSE., 0.0_wp, velocity_after)
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  CALL eval_mui_voellmy_source(flow_depth, velocity_after(1),                &
       velocity_after(2), speed_after, mixture_density, reduced_gravity,     &
       gravity_coefficient, solid_fraction, 0.0_wp, .FALSE., 0.0_wp,         &
       source_x, source_y)
  CALL eval_mui_voellmy_implicit_source(                                    &
       CMPLX(velocity_after(1),0.0_wp,wp), CMPLX(velocity_after(2),0.0_wp,wp), &
       CMPLX(mixture_density,0.0_wp,wp), CMPLX(reduced_gravity,0.0_wp,wp),  &
       CMPLX(speed_after,0.0_wp,wp), CMPLX(speed_after**2,0.0_wp,wp),       &
       .TRUE., .FALSE., 0.0_wp, complex_source_x, complex_source_y)
  source = [source_x+REAL(complex_source_x,wp),                              &
       source_y+REAL(complex_source_y,wp)]
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, source, 'mu(I)-Voellmy')

  CALL solve_mui_voellmy_analytic_velocity(integration_step, flow_depth,   &
       speed_ratio, velocity_before, mixture_density, reduced_gravity,       &
       gravity_coefficient, solid_fraction, 0.0_wp, .TRUE., 500.0_wp,        &
       .TRUE., 100.0_wp, velocity_after)
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  CALL eval_mui_voellmy_source(flow_depth, velocity_after(1),                &
       velocity_after(2), speed_after, mixture_density, reduced_gravity,     &
       gravity_coefficient, solid_fraction, 0.0_wp, .TRUE., 500.0_wp,        &
       source_x, source_y)
  CALL eval_mui_voellmy_implicit_source(                                    &
       CMPLX(velocity_after(1),0.0_wp,wp), CMPLX(velocity_after(2),0.0_wp,wp), &
       CMPLX(mixture_density,0.0_wp,wp), CMPLX(reduced_gravity,0.0_wp,wp),  &
       CMPLX(speed_after,0.0_wp,wp), CMPLX(speed_after**2,0.0_wp,wp),       &
       .TRUE., .TRUE., 100.0_wp, complex_source_x, complex_source_y)
  source = [source_x+REAL(complex_source_x,wp),                              &
       source_y+REAL(complex_source_y,wp)]
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, source, 'mu(I)-Voellmy stochastic pore pressure')

  WRITE(*,*) 'All analytic rheology exact-flow checks passed.'

CONTAINS

  SUBROUTINE assert_stage(before, after, coefficient, mass_depth, source, label)

    REAL(wp), INTENT(IN) :: before(2), after(2), coefficient, mass_depth
    REAL(wp), INTENT(IN) :: source(2)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp) :: tolerance, source_power, initial_speed, final_speed
    REAL(wp) :: average_decay_rate, final_decay_rate, rate_tolerance

    tolerance = 2.0E-10_wp * MAX(1.0_wp, MAXVAL(ABS(before)),              &
         ABS(coefficient)*MAXVAL(ABS(source))/mass_depth)
    initial_speed = SQRT(SUM(before**2))
    final_speed = SQRT(SUM(after**2))
    source_power = SUM(after*source)
    IF (source_power .GT. tolerance*mass_depth) THEN
       WRITE(*,*) 'FAIL: exact rheology source accelerates the flow ', label
       ERROR STOP 1
    END IF
    IF (SUM(after*before) .LT. 0.0_wp) THEN
       WRITE(*,*) 'FAIL: analytic stage reversed direction ', label
       ERROR STOP 1
    END IF
    IF (final_speed .GT. initial_speed + tolerance) THEN
       WRITE(*,*) 'FAIL: analytic friction increased speed ', label
       ERROR STOP 1
    END IF
    IF (final_speed .GT. tolerance) THEN
       average_decay_rate = (initial_speed-final_speed)/coefficient
       final_decay_rate = SQRT(SUM(source**2))/mass_depth
       rate_tolerance = 5.0E-2_wp*MAX(average_decay_rate, final_decay_rate, &
            TINY(1.0_wp))
       IF (ABS(average_decay_rate-final_decay_rate) .GT. rate_tolerance) THEN
          WRITE(*,*) 'FAIL: exact-flow decay disagrees with terminal source ', &
               label, average_decay_rate, final_decay_rate
          ERROR STOP 1
       END IF
    END IF

  END SUBROUTINE assert_stage

  SUBROUTINE assert_close(value, reference, label)

    REAL(wp), INTENT(IN) :: value, reference
    CHARACTER(LEN=*), INTENT(IN) :: label

    IF (ABS(value-reference) .GT. 2.0E-12_wp*MAX(1.0_wp,ABS(reference))) THEN
       WRITE(*,*) 'FAIL: ', label, value, reference
       ERROR STOP 1
    END IF

  END SUBROUTINE assert_close

END PROGRAM test_analytic_velocity
