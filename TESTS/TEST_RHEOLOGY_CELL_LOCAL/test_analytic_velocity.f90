PROGRAM test_analytic_velocity

  USE parameters_2d, ONLY : wp, n_solid
  USE constitutive_parameters_2d, ONLY : grav, inv_grav, diam_s, rho_s
  USE rheology_parameters_2d, ONLY : I_0, mu_2, mu_s, muI_inf, mu_0, mu_inf, &
       Fr_0, U_w, xi, T_ref
  USE state_conversion_2d, ONLY : physical_state_type
  USE rheology_analytic_utils_2d, ONLY : lambert_w0, solve_quadratic_in_interval, &
       solve_cubic_in_interval, solve_exact_velocity_flow,                    &
       EXACT_FLOW_CONSTANT, EXACT_FLOW_LINEAR, EXACT_FLOW_AFFINE,              &
       EXACT_FLOW_QUADRATIC, EXACT_FLOW_QUADRATIC_POLYNOMIAL,                  &
       EXACT_FLOW_FROUDE, EXACT_FLOW_LUCAS, EXACT_FLOW_MUI,                    &
       EXACT_FLOW_MUI_VOELLMY
  USE rheology_voellmy_2d, ONLY : voellmy_mu => mu,                          &
       eval_rheology_voellmy_implicit,                                  &
       eval_rheology_voellmy_semi_implicit,                             &
       eval_rheology_voellmy_analytic
  USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
  USE rheology_plastic_2d, ONLY : tau, eval_rheology_plastic_implicit,  &
       eval_rheology_plastic_analytic
  USE rheology_lava_2d, ONLY : nu_ref, visc_par, tau0,                       &
       eval_rheology_lava_implicit, eval_rheology_lava_semi_implicit, &
       eval_rheology_lava_analytic
  USE rheology_lahar_2d, ONLY : alpha2, beta2, alpha1_coeff, beta1, Kappa,   &
       n_td, n_td2, eval_rheology_lahar_implicit,                       &
       eval_rheology_lahar_semi_implicit,                               &
       eval_rheology_lahar_analytic
  USE rheology_empirical_2d, ONLY : eval_rheology_empirical_implicit,   &
       eval_rheology_empirical_analytic
  USE rheology_quadratic_2d, ONLY : friction_factor,                        &
       eval_rheology_quadratic_implicit, eval_rheology_quadratic_analytic
  USE rheology_froude_2d, ONLY : eval_rheology_froude_semi_implicit,                  &
       eval_rheology_froude_analytic
  USE rheology_lucas_2d, ONLY : eval_rheology_lucas_semi_implicit,                    &
       eval_rheology_lucas_analytic
  USE rheology_mui_2d, ONLY : eval_rheology_mui_semi_implicit, eval_rheology_mui_analytic
  USE rheology_mui_voellmy_2d, ONLY : eval_rheology_mui_voellmy_semi_implicit,        &
       eval_rheology_mui_voellmy_implicit, eval_rheology_mui_voellmy_analytic
  USE rheology_2d, ONLY : eval_rheology_analytic

  IMPLICIT NONE

  TYPE(physical_state_type) :: rheo_state
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
  BLOCK
    TYPE(physical_state_type) :: ps
    ps%h = flow_depth
    ps%u = velocity_before(1)
    ps%v = velocity_before(2)
    CALL eval_rheology_plastic_analytic(integration_step, ps, speed_ratio, velocity_after)
  END BLOCK
  CALL assert_dispatch_matches('rheology_plastic', .FALSE., 0.0_wp,       &
       .FALSE., 0.0_wp, velocity_after, 'plastic dispatch')
  expected_speed = MAX(0.0_wp, SQRT(SUM(velocity_before**2)) -              &
       integration_step*tau/(flow_depth*speed_ratio))
  CALL assert_close(SQRT(SUM(velocity_after**2)), expected_speed,           &
       'plastic exact constant-drag flow')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  BLOCK
    TYPE(complex_physical_state_type) :: cs
    cs%u = CMPLX(velocity_after(1),0.0_wp,wp)
    cs%v = CMPLX(velocity_after(2),0.0_wp,wp)
    cs%rho_m = CMPLX(mixture_density,0.0_wp,wp)
    cs%mod_vel=CMPLX(speed_after,0.0_wp,wp)
    cs%mod_vel2=cs%mod_vel**2
    cs%mod_vel_hor=CMPLX(horizontal_speed,0.0_wp,wp)
    CALL eval_rheology_plastic_implicit(cs, complex_source_x, complex_source_y)
  END BLOCK
  CALL assert_stage(velocity_before, velocity_after, integration_step,     &
       mixture_density*flow_depth, [REAL(complex_source_x,wp),              &
       REAL(complex_source_y,wp)], 'plastic')

  nu_ref = 0.2_wp
  visc_par = 0.01_wp
  tau0 = 100.0_wp
  T_ref = 300.0_wp
  BLOCK
    TYPE(physical_state_type) :: ps
    ps%h=flow_depth
    ps%T=temperature
    ps%rho_m=mixture_density
    ps%u=velocity_before(1)
    ps%v=velocity_before(2)
    CALL eval_rheology_lava_analytic(integration_step, ps, speed_ratio, velocity_after)
  END BLOCK
  CALL assert_dispatch_matches('rheology_lava', .FALSE., 0.0_wp,           &
       .FALSE., 0.0_wp, velocity_after, 'lava dispatch')
  linear_drag = 3.0_wp*nu_ref/flow_depth**2 *                               &
       EXP(-visc_par*(temperature-T_ref))
  yield_drag = tau0/(mixture_density*flow_depth*speed_ratio)
  expected_speed = MAX(0.0_wp, (SQRT(SUM(velocity_before**2)) +             &
       yield_drag/linear_drag)*EXP(-linear_drag*integration_step) -         &
       yield_drag/linear_drag)
  CALL assert_close(SQRT(SUM(velocity_after**2)), expected_speed,           &
       'lava exact affine-drag flow')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  BLOCK
    TYPE(complex_physical_state_type) :: cs
    cs%h=CMPLX(flow_depth,0.0_wp,wp)
    cs%T=CMPLX(temperature,0.0_wp,wp)
    cs%rho_m=CMPLX(mixture_density,0.0_wp,wp)
    cs%u=CMPLX(velocity_after(1),0.0_wp,wp)
    cs%v=CMPLX(velocity_after(2),0.0_wp,wp)
    cs%mod_vel=CMPLX(speed_after,0.0_wp,wp)
    cs%mod_vel2=cs%mod_vel**2
    cs%mod_vel_hor=CMPLX(horizontal_speed,0.0_wp,wp)
    CALL eval_rheology_lava_implicit(cs, complex_source_x, complex_source_y)
  END BLOCK
  BLOCK
    TYPE(physical_state_type) :: ps
    ps%u=velocity_after(1)
    ps%v=velocity_after(2)
    ps%mod_vel=speed_after
    ps%mod_vel2=ps%mod_vel**2
    ps%mod_vel_hor=horizontal_speed
    CALL eval_rheology_lava_semi_implicit(ps, source_x, source_y)
  END BLOCK
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
  BLOCK
    TYPE(physical_state_type) :: ps
    ps%h=flow_depth
    ps%T=temperature
    ps%u=velocity_before(1)
    ps%v=velocity_before(2)
    ps%rho_m=mixture_density
    ps%alphas=solid_fraction
    CALL eval_rheology_lahar_analytic(integration_step, ps, speed_ratio, velocity_after)
  END BLOCK
  CALL assert_dispatch_matches('rheology_lahar', .FALSE., 0.0_wp,          &
       .FALSE., 0.0_wp, velocity_after, 'lahar dispatch')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  BLOCK
    TYPE(complex_physical_state_type) :: cs
    cs%h=CMPLX(flow_depth,0.0_wp,wp)
    cs%T=CMPLX(temperature,0.0_wp,wp)
    cs%u=CMPLX(velocity_after(1),0.0_wp,wp)
    cs%v=CMPLX(velocity_after(2),0.0_wp,wp)
    cs%rho_m=CMPLX(mixture_density,0.0_wp,wp)
    cs%inv_rhom=CMPLX(1.0_wp/mixture_density,0.0_wp,wp)
    cs%alphas=CMPLX(solid_fraction,0.0_wp,wp)
    cs%mod_vel=CMPLX(speed_after,0.0_wp,wp)
    cs%mod_vel2=cs%mod_vel**2
    cs%mod_vel_hor=CMPLX(horizontal_speed,0.0_wp,wp)
    CALL eval_rheology_lahar_implicit(cs, complex_source_x, complex_source_y)
  END BLOCK
  BLOCK
    TYPE(physical_state_type) :: ps
    ps%h=flow_depth
    ps%u=velocity_after(1)
    ps%v=velocity_after(2)
    ps%rho_m=mixture_density
    ps%alphas=solid_fraction
    ps%mod_vel=speed_after
    ps%mod_vel2=ps%mod_vel**2
    ps%mod_vel_hor=horizontal_speed
    CALL eval_rheology_lahar_semi_implicit(ps, source_x, source_y)
  END BLOCK
  source = [REAL(complex_source_x,wp)+source_x,                              &
       REAL(complex_source_y,wp)+source_y]
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, source, 'lahar')

  BLOCK
    TYPE(physical_state_type) :: ps
    ps%h=flow_depth
    ps%T=temperature
    ps%rho_m=mixture_density
    ps%u=velocity_before(1)
    ps%v=velocity_before(2)
    CALL eval_rheology_empirical_analytic(integration_step, ps, speed_ratio, velocity_after)
  END BLOCK
  CALL assert_dispatch_matches('rheology_empirical', .FALSE., 0.0_wp,      &
       .FALSE., 0.0_wp, velocity_after, 'empirical dispatch')
  linear_drag = 1.0E-3_wp*speed_ratio/                                     &
       (flow_depth*(1.0_wp+10.0_wp*flow_depth))
  expected_speed = SQRT(SUM(velocity_before**2))*                           &
       EXP(-linear_drag*integration_step)
  CALL assert_close(SQRT(SUM(velocity_after**2)), expected_speed,           &
       'empirical exact linear-drag flow')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  BLOCK
    TYPE(complex_physical_state_type) :: cs
    cs%h=CMPLX(flow_depth,0.0_wp,wp)
    cs%T=CMPLX(temperature,0.0_wp,wp)
    cs%rho_m=CMPLX(mixture_density,0.0_wp,wp)
    cs%u=CMPLX(velocity_after(1),0.0_wp,wp)
    cs%v=CMPLX(velocity_after(2),0.0_wp,wp)
    cs%mod_vel=CMPLX(speed_after,0.0_wp,wp)
    cs%mod_vel2=cs%mod_vel**2
    cs%mod_vel_hor=CMPLX(horizontal_speed,0.0_wp,wp)
    CALL eval_rheology_empirical_implicit(cs, complex_source_x, complex_source_y)
  END BLOCK
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, [REAL(complex_source_x,wp),               &
       REAL(complex_source_y,wp)], 'empirical')

  friction_factor = 0.1_wp
  BLOCK
    TYPE(physical_state_type) :: ps
    ps%h=flow_depth
    ps%u=velocity_before(1)
    ps%v=velocity_before(2)
    CALL eval_rheology_quadratic_analytic(integration_step, ps, speed_ratio, velocity_after)
  END BLOCK
  CALL assert_dispatch_matches('rheology_quadratic', .FALSE., 0.0_wp,      &
       .FALSE., 0.0_wp, velocity_after, 'quadratic dispatch')
  linear_drag = friction_factor*speed_ratio**2/flow_depth
  expected_speed = SQRT(SUM(velocity_before**2))/                           &
       (1.0_wp+linear_drag*SQRT(SUM(velocity_before**2))*integration_step)
  CALL assert_close(SQRT(SUM(velocity_after**2)), expected_speed,           &
       'quadratic exact reciprocal-speed flow')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  BLOCK
    TYPE(complex_physical_state_type) :: cs
    cs%u=CMPLX(velocity_after(1),0.0_wp,wp)
    cs%v=CMPLX(velocity_after(2),0.0_wp,wp)
    cs%rho_m=CMPLX(mixture_density,0.0_wp,wp)
    cs%mod_vel=CMPLX(speed_after,0.0_wp,wp)
    cs%mod_vel2=cs%mod_vel**2
    cs%mod_vel_hor=CMPLX(horizontal_speed,0.0_wp,wp)
    CALL eval_rheology_quadratic_implicit(cs, complex_source_x, complex_source_y)
  END BLOCK
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, [REAL(complex_source_x,wp),               &
       REAL(complex_source_y,wp)], 'quadratic')

  voellmy_mu = 0.3_wp
  xi = 500.0_wp
    BLOCK
    TYPE(physical_state_type) :: rs
    rs%h = flow_depth
    rs%u = velocity_before(1)
    rs%v = velocity_before(2)
    rs%rho_m = mixture_density
    rs%red_grav = reduced_gravity
    rs%exc_pore_pres = 0.0_wp
    CALL eval_rheology_voellmy_analytic(integration_step, rs, speed_ratio, gravity_coefficient, &
         0.0_wp, .FALSE., .FALSE., 0.0_wp, velocity_after)
  END BLOCK
  CALL assert_dispatch_matches('rheology_voellmy', .FALSE., 0.0_wp,        &
       .FALSE., 0.0_wp, velocity_after, 'Voellmy dispatch')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  BLOCK
    TYPE(complex_physical_state_type) :: cs
    cs%u=CMPLX(velocity_after(1),0.0_wp,wp)
    cs%v=CMPLX(velocity_after(2),0.0_wp,wp)
    cs%rho_m=CMPLX(mixture_density,0.0_wp,wp)
    cs%red_grav=CMPLX(reduced_gravity,0.0_wp,wp)
    cs%mod_vel=CMPLX(speed_after,0.0_wp,wp)
    cs%mod_vel2=CMPLX(speed_after**2,0.0_wp,wp)
    CALL eval_rheology_voellmy_implicit(cs, .FALSE., 0.0_wp, &
         complex_source_x, complex_source_y)
  END BLOCK
  BLOCK
    TYPE(physical_state_type) :: rs
    rs%h=flow_depth
    rs%u=velocity_after(1)
    rs%v=velocity_after(2)
    rs%rho_m=mixture_density
    rs%red_grav=reduced_gravity
    rs%exc_pore_pres=0.0_wp
    rs%mod_vel=speed_after
    CALL eval_rheology_voellmy_semi_implicit(rs, gravity_coefficient, .FALSE., source_x, source_y)
  END BLOCK
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
  rheo_state%h=flow_depth
  rheo_state%u=velocity_before(1)
  rheo_state%v=velocity_before(2)
  rheo_state%rho_m=mixture_density
  rheo_state%red_grav=reduced_gravity
  CALL eval_rheology_froude_analytic(integration_step,rheo_state,speed_ratio,gravity_coefficient,0.0_wp, &
       .FALSE.,stochastic_value,velocity_after)
  CALL assert_dispatch_matches('rheology_froude', .FALSE., 0.0_wp,          &
       .FALSE., stochastic_value, velocity_after, 'Froude dispatch')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  rheo_state%h=flow_depth
  rheo_state%u=velocity_after(1)
  rheo_state%v=velocity_after(2)
  rheo_state%rho_m=mixture_density
  rheo_state%red_grav=reduced_gravity
      rheo_state%mod_vel=speed_after
    rheo_state%mod_vel_hor=horizontal_speed
CALL eval_rheology_froude_semi_implicit(rheo_state, &
             gravity_coefficient, &
             .FALSE., &
             stochastic_value, &
             source_x, &
             source_y)
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, [source_x, source_y], 'Froude')

  stochastic_value = -2.0_wp
  rheo_state%h=flow_depth
  rheo_state%u=velocity_before(1)
  rheo_state%v=velocity_before(2)
  rheo_state%rho_m=mixture_density
  rheo_state%red_grav=reduced_gravity
  CALL eval_rheology_froude_analytic(integration_step,rheo_state,speed_ratio,gravity_coefficient,0.0_wp, &
       .TRUE.,stochastic_value,velocity_after)
  CALL assert_dispatch_matches('rheology_froude', .FALSE., 0.0_wp,         &
       .TRUE., stochastic_value, velocity_after, 'Froude stochastic dispatch')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  rheo_state%h=flow_depth
  rheo_state%u=velocity_after(1)
  rheo_state%v=velocity_after(2)
  rheo_state%rho_m=mixture_density
  rheo_state%red_grav=reduced_gravity
      rheo_state%mod_vel=speed_after
    rheo_state%mod_vel_hor=horizontal_speed
CALL eval_rheology_froude_semi_implicit(rheo_state, &
             gravity_coefficient, &
             .TRUE., &
             stochastic_value, &
             source_x, &
             source_y)
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, [source_x, source_y], 'Froude clamp')
  stochastic_value = 0.0_wp

  mu_0 = 0.0_wp
  mu_inf = 0.0_wp
  rheo_state%h=flow_depth
  rheo_state%u=velocity_before(1)
  rheo_state%v=velocity_before(2)
  rheo_state%rho_m=mixture_density
  rheo_state%red_grav=reduced_gravity
  CALL eval_rheology_froude_analytic(integration_step,rheo_state,speed_ratio,gravity_coefficient,0.0_wp, &
       .FALSE.,stochastic_value,velocity_after)
  CALL assert_dispatch_matches('rheology_froude', .FALSE., 0.0_wp,         &
       .FALSE., stochastic_value, velocity_after, 'Froude zero-drag dispatch')
  CALL assert_close(velocity_after(1), velocity_before(1), 'Froude zero-drag x')
  CALL assert_close(velocity_after(2), velocity_before(2), 'Froude zero-drag y')
  mu_0 = 0.5_wp
  mu_inf = 0.2_wp

  U_w = 2.0_wp
  rheo_state%h=flow_depth
  rheo_state%u=velocity_before(1)
  rheo_state%v=velocity_before(2)
  rheo_state%rho_m=mixture_density
  rheo_state%red_grav=reduced_gravity
  CALL eval_rheology_lucas_analytic(integration_step,rheo_state,speed_ratio,gravity_coefficient,0.0_wp, &
       .FALSE.,stochastic_value,velocity_after)
  CALL assert_dispatch_matches('rheology_lucas', .FALSE., 0.0_wp,          &
       .FALSE., stochastic_value, velocity_after, 'Lucas dispatch')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  rheo_state%h=flow_depth
  rheo_state%u=velocity_after(1)
  rheo_state%v=velocity_after(2)
  rheo_state%rho_m=mixture_density
  rheo_state%red_grav=reduced_gravity
      rheo_state%mod_vel=speed_after
    rheo_state%mod_vel_hor=horizontal_speed
CALL eval_rheology_lucas_semi_implicit(rheo_state, &
             gravity_coefficient, &
             .FALSE., &
             stochastic_value, &
             source_x, &
             source_y)
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, [source_x, source_y], 'Lucas')

  stochastic_value = -10.0_wp
  rheo_state%h=flow_depth
  rheo_state%u=velocity_before(1)
  rheo_state%v=velocity_before(2)
  rheo_state%rho_m=mixture_density
  rheo_state%red_grav=reduced_gravity
  CALL eval_rheology_lucas_analytic(integration_step,rheo_state,speed_ratio,gravity_coefficient,0.0_wp, &
       .TRUE.,stochastic_value,velocity_after)
  CALL assert_dispatch_matches('rheology_lucas', .FALSE., 0.0_wp,          &
       .TRUE., stochastic_value, velocity_after, 'Lucas stochastic dispatch')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  rheo_state%h=flow_depth
  rheo_state%u=velocity_after(1)
  rheo_state%v=velocity_after(2)
  rheo_state%rho_m=mixture_density
  rheo_state%red_grav=reduced_gravity
      rheo_state%mod_vel=speed_after
    rheo_state%mod_vel_hor=horizontal_speed
CALL eval_rheology_lucas_semi_implicit(rheo_state, &
             gravity_coefficient, &
             .TRUE., &
             stochastic_value, &
             source_x, &
             source_y)
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, [source_x, source_y], 'Lucas stochastic branch')
  stochastic_value = 0.0_wp

  I_0 = 0.3_wp
  mu_s = 0.4_wp
  mu_2 = 0.7_wp
  muI_inf = 0.8_wp
  BLOCK
    TYPE(physical_state_type) :: ps
    ps%h=flow_depth
    ps%u=velocity_before(1)
    ps%v=velocity_before(2)
    ps%rho_m=mixture_density
    ps%red_grav=reduced_gravity
    ps%alphas=solid_fraction
    ps%exc_pore_pres=0.0_wp
    CALL eval_rheology_mui_analytic(integration_step, ps, speed_ratio, &
         gravity_coefficient, 0.0_wp, .FALSE., velocity_after)
  END BLOCK
  CALL assert_dispatch_matches('rheology_mui', .FALSE., 0.0_wp, &
       .FALSE., 0.0_wp, velocity_after, 'mu(I) dispatch')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  BLOCK
    TYPE(physical_state_type) :: ps
    ps%h=flow_depth
    ps%u=velocity_after(1)
    ps%v=velocity_after(2)
    ps%rho_m=mixture_density
    ps%red_grav=reduced_gravity
    ps%alphas=solid_fraction
    ps%exc_pore_pres=0.0_wp
        ps%mod_vel=speed_after
    ps%mod_vel_hor=horizontal_speed
CALL eval_rheology_mui_semi_implicit(ps, &
             gravity_coefficient, &
             .FALSE., &
             source_x, &
             source_y)
  END BLOCK
  CALL assert_stage(velocity_before, velocity_after, integration_step, &
       mixture_density*flow_depth, [source_x, source_y], 'mu(I)')

  BLOCK
    TYPE(physical_state_type) :: ps
    ps%h=flow_depth
    ps%u=velocity_before(1)
    ps%v=velocity_before(2)
    ps%rho_m=mixture_density
    ps%red_grav=reduced_gravity
    ps%alphas=solid_fraction
    ps%exc_pore_pres=1000.0_wp
    CALL eval_rheology_mui_analytic(integration_step, ps, speed_ratio, gravity_coefficient, &
         0.0_wp, .TRUE., velocity_after)
  END BLOCK
  CALL assert_dispatch_matches('rheology_mui', .TRUE., 1000.0_wp,           &
       .FALSE., 0.0_wp, velocity_after, 'mu(I) pore-pressure dispatch')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  BLOCK
    TYPE(physical_state_type) :: ps
    ps%h=flow_depth
    ps%u=velocity_after(1)
    ps%v=velocity_after(2)
    ps%rho_m=mixture_density
    ps%red_grav=reduced_gravity
    ps%alphas=solid_fraction
    ps%exc_pore_pres=1000.0_wp
        ps%mod_vel=speed_after
    ps%mod_vel_hor=horizontal_speed
CALL eval_rheology_mui_semi_implicit(ps, &
             gravity_coefficient, &
             .TRUE., &
             source_x, &
             source_y)
  END BLOCK
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, [source_x, source_y], 'mu(I) pore pressure')

  xi = 500.0_wp
  BLOCK
    TYPE(physical_state_type) :: ps
    ps%h=flow_depth
    ps%u=velocity_before(1)
    ps%v=velocity_before(2)
    ps%rho_m=mixture_density
    ps%red_grav=reduced_gravity
    ps%alphas=solid_fraction
    ps%exc_pore_pres=0.0_wp
    CALL eval_rheology_mui_voellmy_analytic(integration_step, ps, speed_ratio, gravity_coefficient, &
         0.0_wp, .FALSE., .FALSE., 0.0_wp, velocity_after)
  END BLOCK
  CALL assert_dispatch_matches('rheology_mui_voellmy', .FALSE., 0.0_wp,     &
       .FALSE., 0.0_wp, velocity_after, 'mu(I)-Voellmy dispatch')
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  horizontal_speed = SQRT(SUM(velocity_after**2))
  BLOCK
    TYPE(physical_state_type) :: ps
    ps%h=flow_depth
    ps%u=velocity_after(1)
    ps%v=velocity_after(2)
    ps%rho_m=mixture_density
    ps%red_grav=reduced_gravity
    ps%alphas=solid_fraction
    ps%exc_pore_pres=0.0_wp
        ps%mod_vel=speed_after
    ps%mod_vel_hor=horizontal_speed
CALL eval_rheology_mui_voellmy_semi_implicit(ps, &
             gravity_coefficient, &
             .FALSE., &
             source_x, &
             source_y)
  END BLOCK
  BLOCK
    TYPE(complex_physical_state_type) :: cs
    cs%u=CMPLX(velocity_after(1),0.0_wp,wp)
    cs%v=CMPLX(velocity_after(2),0.0_wp,wp)
    cs%rho_m=CMPLX(mixture_density,0.0_wp,wp)
    cs%red_grav=CMPLX(reduced_gravity,0.0_wp,wp)
        cs%mod_vel=CMPLX(speed_after,0.0_wp,wp)
    cs%mod_vel2=cs%mod_vel**2
    cs%mod_vel_hor=CMPLX(horizontal_speed,0.0_wp,wp)
CALL eval_rheology_mui_voellmy_implicit(cs, &
             .FALSE., &
             0.0_wp, &
             complex_source_x, &
             complex_source_y)
  END BLOCK
  source = [source_x+REAL(complex_source_x,wp),                              &
       source_y+REAL(complex_source_y,wp)]
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, source, 'mu(I)-Voellmy')

  BLOCK
    TYPE(physical_state_type) :: ps
    ps%h=flow_depth
    ps%u=velocity_before(1)
    ps%v=velocity_before(2)
    ps%rho_m=mixture_density
    ps%red_grav=reduced_gravity
    ps%alphas=solid_fraction
    ps%exc_pore_pres=500.0_wp
    CALL eval_rheology_mui_voellmy_analytic(integration_step, ps, speed_ratio, gravity_coefficient, &
         0.0_wp, .TRUE., .TRUE., 100.0_wp, velocity_after)
  END BLOCK
  speed_after = speed_ratio * SQRT(SUM(velocity_after**2))
  BLOCK
    TYPE(physical_state_type) :: ps
    ps%h=flow_depth
    ps%u=velocity_after(1)
    ps%v=velocity_after(2)
    ps%rho_m=mixture_density
    ps%red_grav=reduced_gravity
    ps%alphas=solid_fraction
    ps%exc_pore_pres=500.0_wp
        ps%mod_vel=speed_after
    ps%mod_vel_hor=horizontal_speed
CALL eval_rheology_mui_voellmy_semi_implicit(ps, &
             gravity_coefficient, &
             .TRUE., &
             source_x, &
             source_y)
  END BLOCK
  BLOCK
    TYPE(complex_physical_state_type) :: cs
    cs%u=CMPLX(velocity_after(1),0.0_wp,wp)
    cs%v=CMPLX(velocity_after(2),0.0_wp,wp)
    cs%rho_m=CMPLX(mixture_density,0.0_wp,wp)
    cs%red_grav=CMPLX(reduced_gravity,0.0_wp,wp)
        cs%mod_vel=CMPLX(speed_after,0.0_wp,wp)
    cs%mod_vel2=cs%mod_vel**2
    cs%mod_vel_hor=CMPLX(horizontal_speed,0.0_wp,wp)
CALL eval_rheology_mui_voellmy_implicit(cs, &
             .TRUE., &
             100.0_wp, &
             complex_source_x, &
             complex_source_y)
  END BLOCK
  CALL assert_dispatch_matches('rheology_mui_voellmy', .TRUE., 500.0_wp,   &
       .TRUE., 100.0_wp, velocity_after,                                   &
       'mu(I)-Voellmy stochastic pore-pressure dispatch')
  source = [source_x+REAL(complex_source_x,wp),                              &
       source_y+REAL(complex_source_y,wp)]
  CALL assert_stage(velocity_before, velocity_after, integration_step,      &
       mixture_density*flow_depth, source, 'mu(I)-Voellmy stochastic pore pressure')

  CALL assert_flow_semigroup(EXACT_FLOW_CONSTANT, [0.6_wp], 0.5_wp,       &
       'constant drag')
  CALL assert_flow_semigroup(EXACT_FLOW_LINEAR, [0.2_wp], 0.5_wp,         &
       'linear drag')
  CALL assert_flow_semigroup(EXACT_FLOW_AFFINE, [0.5_wp, 0.2_wp], 0.5_wp, &
       'affine drag')
  CALL assert_flow_semigroup(EXACT_FLOW_QUADRATIC, [0.08_wp], 0.5_wp,     &
       'quadratic drag')
  CALL assert_flow_semigroup(EXACT_FLOW_QUADRATIC_POLYNOMIAL,              &
       [0.2_wp, 0.1_wp, 0.02_wp], 0.5_wp, 'polynomial drag')
  CALL assert_flow_semigroup(EXACT_FLOW_FROUDE,                            &
       [0.5_wp, 0.2_wp, -0.86_wp, 0.5_wp, 0.2_wp], 1.0_wp,                &
       'Froude threshold crossing')
  CALL assert_flow_semigroup(EXACT_FLOW_LUCAS,                             &
       [0.5_wp, 1.0_wp, 0.0_wp, 0.5_wp, 0.2_wp, 4.3_wp, 0.0_wp],          &
       1.0_wp, 'Lucas threshold crossing')
  CALL assert_flow_semigroup(EXACT_FLOW_MUI,                               &
       [0.04_wp, 0.01_wp, 0.002_wp, 0.3_wp, 0.08_wp], 0.5_wp,             &
       'mu(I) drag')
  CALL assert_flow_semigroup(EXACT_FLOW_MUI_VOELLMY,                       &
       [0.04_wp, 0.01_wp, 0.002_wp, 0.3_wp, 0.08_wp, 0.02_wp], 0.5_wp,    &
       'mu(I)-Voellmy drag')
  CALL solve_exact_velocity_flow(EXACT_FLOW_CONSTANT, [2.0_wp], 1.0_wp,   &
       0.5_wp, speed_after)
  CALL assert_close(speed_after, 0.0_wp, 'constant-drag stopping limit')
  CALL solve_exact_velocity_flow(EXACT_FLOW_MUI,                          &
       [0.04_wp, 0.01_wp, 0.002_wp, 0.3_wp, 0.08_wp], 0.0_wp,             &
       4.5_wp, speed_after)
  CALL assert_close(speed_after, 4.5_wp, 'zero-interval identity')
  CALL assert_flow_limits(EXACT_FLOW_CONSTANT, [0.6_wp], [0.0_wp], 'constant')
  CALL assert_flow_limits(EXACT_FLOW_LINEAR, [0.2_wp], [0.0_wp], 'linear')
  CALL assert_flow_limits(EXACT_FLOW_AFFINE, [0.5_wp, 0.2_wp],            &
       [0.0_wp, 0.0_wp], 'affine')
  CALL assert_flow_limits(EXACT_FLOW_QUADRATIC, [0.08_wp], [0.0_wp],      &
       'quadratic')
  CALL assert_flow_limits(EXACT_FLOW_QUADRATIC_POLYNOMIAL,                &
       [0.2_wp, 0.1_wp, 0.02_wp], [0.0_wp, 0.0_wp, 0.0_wp], 'polynomial')
  CALL assert_flow_limits(EXACT_FLOW_FROUDE,                              &
       [0.5_wp, 0.2_wp, -0.86_wp, 0.5_wp, 0.2_wp],                        &
       [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp], 'Froude')
  CALL assert_flow_limits(EXACT_FLOW_LUCAS,                               &
       [0.5_wp, 1.0_wp, 0.0_wp, 0.5_wp, 0.2_wp, 4.3_wp, 0.0_wp],          &
       [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp], 'Lucas')
  CALL assert_flow_limits(EXACT_FLOW_MUI,                                 &
       [0.04_wp, 0.01_wp, 0.002_wp, 0.3_wp, 0.08_wp],                     &
       [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp], 'mu(I)')
  CALL assert_flow_limits(EXACT_FLOW_MUI_VOELLMY,                         &
       [0.04_wp, 0.01_wp, 0.002_wp, 0.3_wp, 0.08_wp, 0.02_wp],            &
       [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp], 'mu(I)-Voellmy')

  CALL solve_exact_velocity_flow(EXACT_FLOW_QUADRATIC_POLYNOMIAL,         &
       [0.0_wp, 0.2_wp, 0.0_wp], 0.7_wp, 4.5_wp, speed_after)
  CALL assert_close(speed_after, 4.5_wp*EXP(-0.2_wp*0.7_wp),               &
       'polynomial linear-only limit')
  CALL solve_exact_velocity_flow(EXACT_FLOW_QUADRATIC_POLYNOMIAL,         &
       [0.0_wp, 0.0_wp, 0.02_wp], 0.7_wp, 4.5_wp, speed_after)
  CALL assert_close(speed_after, 4.5_wp/(1.0_wp+0.02_wp*4.5_wp*0.7_wp),   &
       'polynomial quadratic-only limit')
  CALL solve_exact_velocity_flow(EXACT_FLOW_QUADRATIC_POLYNOMIAL,         &
       [0.0_wp, 0.2_wp, 0.02_wp], 0.7_wp, 4.5_wp, speed_after)
  expected_speed = 4.5_wp*EXP(-0.2_wp*0.7_wp) /                           &
       (1.0_wp+(0.02_wp*4.5_wp/0.2_wp)*(1.0_wp-EXP(-0.2_wp*0.7_wp)))
  CALL assert_close(speed_after, expected_speed, 'polynomial zero-constant limit')
  CALL solve_exact_velocity_flow(EXACT_FLOW_QUADRATIC_POLYNOMIAL,         &
       [0.5_wp, 0.2_wp, 0.0_wp], 0.7_wp, 4.5_wp, speed_after)
  expected_speed = (4.5_wp+0.5_wp/0.2_wp)*EXP(-0.2_wp*0.7_wp)-0.5_wp/0.2_wp
  CALL assert_close(speed_after, expected_speed, 'polynomial zero-quadratic limit')
  CALL solve_exact_velocity_flow(EXACT_FLOW_QUADRATIC_POLYNOMIAL,         &
       [1.0_wp, 2.0_wp, 1.0_wp], 0.2_wp, 4.5_wp, speed_after)
  expected_speed = 5.5_wp/(1.0_wp+5.5_wp*0.2_wp)-1.0_wp
  CALL assert_close(speed_after, expected_speed, 'polynomial repeated-root limit')
  CALL assert_flow_semigroup(EXACT_FLOW_QUADRATIC_POLYNOMIAL,              &
       [1.0_wp, 3.0_wp, 1.0_wp], 0.1_wp,                                  &
       'polynomial negative-discriminant branch')
  CALL solve_exact_velocity_flow(EXACT_FLOW_QUADRATIC_POLYNOMIAL,         &
       [1.0_wp, 0.2_wp, 0.02_wp], 10.0_wp, 4.5_wp, speed_after)
  CALL assert_close(speed_after, 0.0_wp, 'polynomial finite stopping limit')

  CALL assert_dispatch_special_cases()

  WRITE(*,*) 'All analytic rheology exact-flow checks passed.'

CONTAINS

  SUBROUTINE assert_dispatch_matches(model_id, use_pore_pressure,           &
       excess_pore_pressure, use_stochastic, stochastic, expected, label)

    CHARACTER(*), INTENT(IN) :: model_id
    LOGICAL, INTENT(IN) :: use_pore_pressure, use_stochastic
    REAL(wp), INTENT(IN) :: excess_pore_pressure, stochastic, expected(2)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp) :: actual(2)
    TYPE(physical_state_type) :: physical_state
    LOGICAL :: model_supported

    physical_state%h=flow_depth
    physical_state%T=temperature
    physical_state%rho_m=mixture_density
    physical_state%u=velocity_before(1)
    physical_state%v=velocity_before(2)
    physical_state%red_grav=reduced_gravity
    physical_state%alphas=solid_fraction
    physical_state%exc_pore_pres=excess_pore_pressure
    CALL eval_rheology_analytic([model_id],integration_step,physical_state, &
         speed_ratio,gravity_coefficient,0.0_wp,use_pore_pressure, &
         use_stochastic,stochastic,actual,model_supported, &
         0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp)
    IF (.NOT. model_supported) THEN
       WRITE(*,*) 'FAIL: analytic dispatcher rejected supported model ', label
       ERROR STOP 1
    END IF
    CALL assert_close(actual(1), expected(1), label//' x')
    CALL assert_close(actual(2), expected(2), label//' y')

  END SUBROUTINE assert_dispatch_matches

  SUBROUTINE assert_dispatch_special_cases()

    REAL(wp) :: actual(2)
    TYPE(physical_state_type) :: physical_state
    LOGICAL :: model_supported

    CALL assert_dispatch_matches('none', .FALSE., 0.0_wp,        &
         .FALSE., 0.0_wp, velocity_before, 'none no-op')
    physical_state%h=flow_depth
    physical_state%T=temperature
    physical_state%rho_m=mixture_density
    physical_state%u=velocity_before(1)
    physical_state%v=velocity_before(2)
    physical_state%red_grav=reduced_gravity
    physical_state%alphas=solid_fraction
    physical_state%exc_pore_pres=0.0_wp
    CALL eval_rheology_analytic(['rheology_unsupported'],integration_step,physical_state, &
         speed_ratio,gravity_coefficient,0.0_wp, .FALSE., &
         .FALSE.,0.0_wp,actual,model_supported, &
         0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp)
    IF (model_supported) ERROR STOP 'Unsupported rheology name was accepted'
    CALL assert_close(actual(1), velocity_before(1), 'unsupported model x')
    CALL assert_close(actual(2), velocity_before(2), 'unsupported model y')

  END SUBROUTINE assert_dispatch_special_cases

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

  SUBROUTINE assert_flow_semigroup(flow_model, coefficients, interval, label)

    INTEGER, INTENT(IN) :: flow_model
    REAL(wp), INTENT(IN) :: coefficients(:), interval
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), PARAMETER :: initial_speed = 4.5_wp
    REAL(wp) :: full_step, first_half, two_half_steps, reference_step, tolerance

    CALL solve_exact_velocity_flow(flow_model, coefficients, interval,     &
         initial_speed, full_step)
    CALL solve_exact_velocity_flow(flow_model, coefficients, 0.5_wp*interval, &
         initial_speed, first_half)
    CALL solve_exact_velocity_flow(flow_model, coefficients, 0.5_wp*interval, &
         first_half, two_half_steps)
    reference_step = reference_flow_solution(flow_model, coefficients,      &
         initial_speed, interval)

    tolerance = 2.0E-9_wp*MAX(1.0_wp, ABS(full_step))
    IF (ABS(full_step-two_half_steps) .GT. tolerance) THEN
       WRITE(*,*) 'FAIL: exact-flow semigroup property ', label,            &
            full_step, two_half_steps
       ERROR STOP 1
    END IF
    tolerance = 2.0E-8_wp*MAX(1.0_wp, ABS(reference_step))
    IF (ABS(full_step-reference_step) .GT. tolerance) THEN
       WRITE(*,*) 'FAIL: exact-flow endpoint disagrees with independent RK4 ', &
            label, full_step, reference_step
       ERROR STOP 1
    END IF
    IF (full_step .LT. 0.0_wp .OR. full_step .GT. initial_speed) THEN
       WRITE(*,*) 'FAIL: exact-flow semigroup left dissipative interval ',  &
            label, full_step
       ERROR STOP 1
    END IF

  END SUBROUTINE assert_flow_semigroup

  SUBROUTINE assert_flow_limits(flow_model, coefficients, zero_coefficients, label)

    INTEGER, INTENT(IN) :: flow_model
    REAL(wp), INTENT(IN) :: coefficients(:), zero_coefficients(:)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), PARAMETER :: initial_speed = 4.5_wp
    REAL(wp) :: result

    CALL solve_exact_velocity_flow(flow_model, coefficients, 0.0_wp,       &
         initial_speed, result)
    CALL assert_close(result, initial_speed, TRIM(label)//' zero-time limit')
    CALL solve_exact_velocity_flow(flow_model, coefficients, 0.5_wp,       &
         0.0_wp, result)
    CALL assert_close(result, 0.0_wp, TRIM(label)//' zero-speed limit')
    CALL solve_exact_velocity_flow(flow_model, zero_coefficients, 0.5_wp,  &
         initial_speed, result)
    CALL assert_close(result, initial_speed, TRIM(label)//' zero-drag limit')

  END SUBROUTINE assert_flow_limits

  FUNCTION reference_flow_solution(flow_model, coefficients, initial_speed, &
       interval) RESULT(speed)

    INTEGER, INTENT(IN) :: flow_model
    REAL(wp), INTENT(IN) :: coefficients(:), initial_speed, interval
    INTEGER, PARAMETER :: reference_steps = 32768
    REAL(wp) :: speed, step_size, k1, k2, k3, k4
    INTEGER :: i

    speed = initial_speed
    step_size = interval/REAL(reference_steps,wp)
    DO i = 1, reference_steps
       k1 = -reference_flow_rate(flow_model, coefficients, speed)
       k2 = -reference_flow_rate(flow_model, coefficients,                   &
            speed+0.5_wp*step_size*k1)
       k3 = -reference_flow_rate(flow_model, coefficients,                   &
            speed+0.5_wp*step_size*k2)
       k4 = -reference_flow_rate(flow_model, coefficients, speed+step_size*k3)
       speed = speed + step_size*(k1+2.0_wp*k2+2.0_wp*k3+k4)/6.0_wp
       IF (speed .LE. 0.0_wp) THEN
          speed = 0.0_wp
          EXIT
       END IF
    END DO

  END FUNCTION reference_flow_solution

  PURE FUNCTION reference_flow_rate(flow_model, coefficients, speed)       &
       RESULT(rate)

    INTEGER, INTENT(IN) :: flow_model
    REAL(wp), INTENT(IN) :: coefficients(:), speed
    REAL(wp) :: rate, friction_speed, denominator

    SELECT CASE (flow_model)
    CASE (EXACT_FLOW_CONSTANT)
       rate = coefficients(1)
    CASE (EXACT_FLOW_LINEAR)
       rate = coefficients(1)*speed
    CASE (EXACT_FLOW_AFFINE)
       rate = coefficients(1)+coefficients(2)*speed
    CASE (EXACT_FLOW_QUADRATIC)
       rate = coefficients(1)*speed**2
    CASE (EXACT_FLOW_QUADRATIC_POLYNOMIAL)
       rate = coefficients(1)+coefficients(2)*speed+                       &
            coefficients(3)*speed**2
    CASE (EXACT_FLOW_FROUDE)
       IF (coefficients(2)*speed+coefficients(3) .LE. 0.0_wp) THEN
          rate = coefficients(1)*coefficients(4)
       ELSE
          rate = coefficients(1)*(coefficients(5)+(coefficients(4)-        &
               coefficients(5))*EXP(-coefficients(2)*speed-coefficients(3)))
       END IF
    CASE (EXACT_FLOW_LUCAS)
       friction_speed = coefficients(2)*speed+coefficients(3)
       IF (coefficients(7) .GT. 0.5_wp .OR.                               &
            friction_speed .LT. coefficients(6)) THEN
          rate = coefficients(1)*coefficients(4)
       ELSEIF (friction_speed .GT. 0.0_wp) THEN
          rate = coefficients(1)*(coefficients(5)+(coefficients(4)-        &
               coefficients(5))*coefficients(6)/friction_speed)
       ELSE
          rate = coefficients(1)*coefficients(4)
       END IF
    CASE (EXACT_FLOW_MUI)
       denominator = coefficients(4)+coefficients(5)*speed
       rate = (coefficients(1)+coefficients(2)*speed+                      &
            coefficients(3)*speed**2)/denominator
    CASE (EXACT_FLOW_MUI_VOELLMY)
       denominator = coefficients(4)+coefficients(5)*speed
       rate = (coefficients(1)+coefficients(2)*speed+                      &
            coefficients(3)*speed**2)/denominator+coefficients(6)*speed**2
    CASE DEFAULT
       rate = 0.0_wp
    END SELECT

  END FUNCTION reference_flow_rate

  SUBROUTINE assert_close(value, reference, label)

    REAL(wp), INTENT(IN) :: value, reference
    CHARACTER(LEN=*), INTENT(IN) :: label

    IF (ABS(value-reference) .GT. 2.0E-12_wp*MAX(1.0_wp,ABS(reference))) THEN
       WRITE(*,*) 'FAIL: ', label, value, reference
       ERROR STOP 1
    END IF

  END SUBROUTINE assert_close

END PROGRAM test_analytic_velocity
