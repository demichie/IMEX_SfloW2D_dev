PROGRAM test_grouped_implicit_stage

  USE parameters_2d
  USE constitutive_parameters_2d
  USE rheology_parameters_2d, ONLY : T_ref
  USE rheology_lava_2d, ONLY : nu_ref, visc_par, tau0
  USE rheology_voellmy_2d, ONLY : mu, xi, solve_rheology_voellmy_analytic
  USE geometry_2d, ONLY : pi_g
  USE equation_metadata_2d, ONLY : equation_partition_type,                 &
       IM_GROUP_METHOD_ANALYTIC, IM_GROUP_METHOD_NEWTON
  USE source_terms_2d, ONLY : eval_source_implicit, eval_source_semi_implicit, &
       eval_source_analytic, analytic_term_context_type
  USE source_registry_2d, ONLY : source_registry_type, source_group_type
  USE source_configuration_2d, ONLY : initialize_configured_sources
  USE source_input_2d, ONLY : rheology_method, rheology_cadence, rheology_group, &
       thermal_method, thermal_cadence, thermal_group, &
       pore_pressure_method, pore_pressure_cadence, pore_pressure_group
  USE source_catalog_2d, ONLY : SOURCE_METHOD_NEWTON, SOURCE_METHOD_ANALYTIC, &
       SOURCE_CADENCE_STAGE, SOURCE_CADENCE_STEP
  USE state_2d, ONLY : correct_friction_momentum
  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver,               &
       finalize_nonlinear_solver, solve_imex_stage, solve_rk_step
  USE state_conversion_2d, ONLY : qc_to_qp, r_phys_var
  USE pore_pressure_laws_2d, ONLY : solve_pore_pressure_dissipation_analytic
  USE thermal_laws_2d, ONLY : solve_thermal_cooling_analytic

  IMPLICIT NONE

  TYPE(equation_partition_type) :: partition
  TYPE(source_registry_type) :: source_registry
  TYPE(equation_partition_type) :: multi_partition
  TYPE(source_registry_type) :: multi_registry
  TYPE(source_group_type) :: filter_group
  TYPE(source_group_type) :: filter_groups(3)
  LOGICAL :: implicit_mask(7), converged, line_search_failed
  INTEGER :: group_map(7), iterations, linear_info, i, equation
  REAL(wp) :: q_old(7), q_stage(7), known_residual(7), source(7), si_source(7)
  REAL(wp) :: filtered_implicit_source(7), filtered_semi_implicit_source(7)
  REAL(wp) :: q_stage_check(7), q_expected(7), effective_stage_term(7)
  REAL(wp) :: qp_known(9), qp_expected(9), qp_actual(9)
  REAL(wp) :: dynamic_pressure, stage_rate_error
  REAL(wp) :: areal_heat_capacity, exact_temperature, exact_energy
  REAL(wp) :: pressure_known, expected_pressure
  REAL(wp) :: expected_h, expected_u, expected_v, expected_density
  REAL(wp) :: expected_alphas(2), expected_alphag(0)
  REAL(wp) :: expected_temperature, expected_red_grav, expected_p_dyn
  REAL(wp) :: expected_Zs, expected_pore_pressure
  REAL(wp) :: speed_ratio, analytic_velocity(2)
  REAL(wp) :: residual, scale, dt_step, a_diag
  COMPLEX(wp) :: complex_q(7), complex_source(7)
  CHARACTER(LEN=256) :: source_error

  n_solid = 2
  n_add_gas = 0
  n_stoch_vars = 0
  n_pore_vars = 1
  n_vars = 7
  n_eqns = 7
  n_RK = 2
  idx_totMassEqn = 1
  idx_uEqn = 2
  idx_vEqn = 3
  idx_engyEqn = 4
  idx_solid_first = 5
  idx_solid_last = 6
  idx_solidEqn_first = 5
  idx_solidEqn_last = 6
  idx_add_gas_first = 7
  idx_add_gas_last = 6
  idx_pore = 7
  idx_poreEqn = 7
  idx_u = n_vars + 1
  idx_v = n_vars + 2
  idx_stoch = 0
  rheology_flag = .TRUE.
  rheology_model='voellmy'
  ! D17 decouples cooling from rheology: enable it in this thermal fixture.
  thermal_cooling_flag=.TRUE.
  pore_pressure_flag = .TRUE.
  gas_flag = .FALSE.
  liquid_flag = .TRUE.
  lava_flag = .FALSE.
  gas_loss_flag = .FALSE.
  liquid_vaporization_flag = .FALSE.
  erosion_coeff = 0.0_wp
  settling_flag = .FALSE.
  coeff_porosity = 0.0_wp
  entrainment_flag = .FALSE.
  sutherland_flag = .FALSE.
  dynamic_permeability_flag = .FALSE.
  stochastic_flag = .FALSE.
  stoch_transport_flag = .FALSE.
  f_inhibit_mode = 'OFF'
  slope_correction_flag = .FALSE.
  curvature_term_flag = .FALSE.
  stochastic_flag = .FALSE.

  ALLOCATE(rho_s(2), inv_rho_s(2), sp_heat_s(2), diam_s(2), sphericity_s(2))
  ALLOCATE(sp_heat_g(0), sp_gas_const_g(0))
  rho_s = [2500.0_wp, 1800.0_wp]
  inv_rho_s = 1.0_wp / rho_s
  sp_heat_s = [900.0_wp, 1100.0_wp]
  diam_s = [1.0E-4_wp, 4.0E-4_wp]
  sphericity_s = [0.8_wp, 0.6_wp]
  rho_l = 1000.0_wp
  inv_rho_l = 1.0_wp / rho_l
  sp_heat_l = 4180.0_wp
  grav = 9.81_wp
  pi_g = ACOS(-1.0_wp)
  pres = 101325.0_wp
  T_ambient = 300.0_wp
  rho_a_amb = pres / (287.0_wp * T_ambient)
  sp_heat_a = 1005.0_wp
  sp_gas_const_a = 287.0_wp
  kin_visc_c = 1.5E-5_wp
  hydraulic_permeability = 1.0E-10_wp
  mu = 0.3_wp
  xi = 500.0_wp
  nu_ref = 0.02_wp
  visc_par = 0.003_wp
  tau0 = 0.0_wp
  T_ref = 900.0_wp
  radiative_term_coeff = 5.103E-8_wp
  convective_term_coeff = 9.0_wp
  T_env = 300.0_wp
  eps_sing = 1.0E-8_wp
  eps_sing4 = eps_sing**4
  implicit_mask = [.FALSE., .TRUE., .TRUE., .TRUE., .FALSE., .FALSE., .TRUE.]
  group_map = [0, 3, 3, 1, 0, 0, 2]
  CALL partition%initialize(implicit_mask, group_map,                       &
       [IM_GROUP_METHOD_ANALYTIC, IM_GROUP_METHOD_ANALYTIC,                 &
       IM_GROUP_METHOD_NEWTON])
  rheology_method='NEWTON'; rheology_cadence='STAGE'; rheology_group=3
  thermal_method='ANALYTIC'; thermal_cadence='STAGE'; thermal_group=1
  pore_pressure_method='ANALYTIC'; pore_pressure_cadence='STAGE'; pore_pressure_group=2
  CALL initialize_configured_sources(source_registry,partition,source_error)
  IF (LEN_TRIM(source_error)>0) ERROR STOP TRIM(source_error)
  IF (SIZE(source_registry%stage_groups())/=4) &
       ERROR STOP 'default configuration did not build the expected stage groups'
  filter_group%id=1
  filter_group%cadence=SOURCE_CADENCE_STAGE
  filter_group%method=SOURCE_METHOD_NEWTON

  filter_group%source_names=[CHARACTER(LEN=64):: &
       'rheology_voellmy', 'thermal_cooling', 'pore_pressure_dissipation']
  filter_group%equation_ids=[2,3,4,7]
  filter_groups(1)=filter_group
  filter_groups(1)%source_names=[CHARACTER(LEN=64)::'rheology_voellmy']
  filter_groups(1)%equation_ids=[2,3]
  filter_groups(2)=filter_group
  filter_groups(2)%source_names=[CHARACTER(LEN=64)::'thermal_cooling']
  filter_groups(2)%equation_ids=[4]
  filter_groups(3)=filter_group
  filter_groups(3)%source_names=[CHARACTER(LEN=64)::'pore_pressure_dissipation']
  filter_groups(3)%equation_ids=[7]
  q_old = [1200.0_wp, 2400.0_wp, -1200.0_wp, 4.0E9_wp,                     &
       180.0_wp, 120.0_wp, 6.0E6_wp]
  q_stage = q_old
  known_residual = 0.0_wp
  dt_step = 0.2_wp
  a_diag = 0.75_wp
  CALL initialize_nonlinear_solver
  CALL solve_rk_step(partition, source_registry%stage_groups(),              &
       q_stage, q_old, dt_step, a_diag,           &
       known_residual, 0.0_wp, 0.0_wp, 0.0_wp, iterations, converged,      &
       linear_info, line_search_failed)
  IF (.NOT. converged .OR. linear_info .NE. 0 .OR. line_search_failed) THEN
     WRITE(*,*) 'FAIL: grouped solve status', converged, linear_info,       &
          line_search_failed, iterations
     ERROR STOP 'grouped implicit stage did not converge'
  END IF

  CALL eval_source_implicit(0.0_wp, 0.0_wp, 0.0_wp, r_qj=q_stage,            &
       r_nh_term_impl=source)
  CALL eval_source_implicit(0.0_wp,0.0_wp,0.0_wp,r_qj=q_stage,                &
       r_nh_term_impl=filtered_implicit_source,group=filter_groups(1))
  IF (MAXVAL(ABS(filtered_implicit_source(2:3)-source(2:3))) .GT.           &
       1.0E-12_wp*MAX(1.0_wp,MAXVAL(ABS(source(2:3)))) .OR.                  &
       ANY(filtered_implicit_source([1,4,5,6,7])/=0.0_wp)) THEN
     ERROR STOP 'rheology solver selection leaked unrelated implicit terms'
  END IF
  CALL eval_source_implicit(0.0_wp,0.0_wp,0.0_wp,r_qj=q_stage,                &
       r_nh_term_impl=filtered_implicit_source,group=filter_groups(2))
  IF (ABS(filtered_implicit_source(4)-source(4)) .GT.                       &
       1.0E-12_wp*MAX(1.0_wp,ABS(source(4))) .OR.                            &
       ANY(filtered_implicit_source([1,2,3,5,6,7])/=0.0_wp)) THEN
     ERROR STOP 'thermal solver selection leaked unrelated implicit terms'
  END IF
  CALL eval_source_implicit(0.0_wp,0.0_wp,0.0_wp,r_qj=q_stage,                &
       r_nh_term_impl=filtered_implicit_source,group=filter_groups(3))
  IF (ABS(filtered_implicit_source(7)-source(7)) .GT.                       &
       1.0E-12_wp*MAX(1.0_wp,ABS(source(7))) .OR.                            &
       ANY(filtered_implicit_source([1,2,3,4,5,6])/=0.0_wp)) THEN
     ERROR STOP 'pore-pressure solver selection leaked unrelated implicit terms'
  END IF

  CALL eval_source_implicit(0.0_wp,0.0_wp,0.0_wp,r_qj=q_stage,                &
       r_nh_term_impl=filtered_implicit_source,group=filter_groups(1))
  complex_q=CMPLX(q_stage,0.0_wp,wp)
  complex_q(2)=CMPLX(q_stage(2),1.0E-12_wp,wp)
  CALL eval_source_implicit(0.0_wp,0.0_wp,0.0_wp,c_qj=complex_q,             &
       c_nh_term_impl=complex_source,group=filter_groups(1))
  IF (MAXVAL(ABS(REAL(complex_source(2:3),wp)-                              &
       filtered_implicit_source(2:3))) .GT.                                  &
       1.0E-12_wp*MAX(1.0_wp,MAXVAL(ABS(source(2:3))))) THEN
     ERROR STOP 'complex rheology solver selection changed the real rate'
  END IF
  IF (MAXVAL(ABS(AIMAG(complex_source(2:3))))==0.0_wp .OR.                  &
       ANY(complex_source([1,4,5,6,7])/=(0.0_wp,0.0_wp))) THEN
     ERROR STOP 'complex rheology solver selection lost derivatives or leaked terms'
  END IF
  CALL eval_source_implicit(0.0_wp,0.0_wp,0.0_wp,r_qj=q_stage,                &
       r_nh_term_impl=filtered_implicit_source,group=filter_groups(2))
  complex_q=CMPLX(q_stage,0.0_wp,wp)
  complex_q(4)=CMPLX(q_stage(4),1.0E-12_wp,wp)
  CALL eval_source_implicit(0.0_wp,0.0_wp,0.0_wp,c_qj=complex_q,             &
       c_nh_term_impl=complex_source,group=filter_groups(2))
  IF (ABS(REAL(complex_source(4),wp)-filtered_implicit_source(4)) .GT.      &
       1.0E-12_wp*MAX(1.0_wp,ABS(source(4))) .OR.                           &
       ABS(AIMAG(complex_source(4)))==0.0_wp .OR.                           &
       ANY(complex_source([1,2,3,5,6,7])/=(0.0_wp,0.0_wp))) THEN
     ERROR STOP 'complex thermal solver selection failed'
  END IF
  CALL eval_source_implicit(0.0_wp,0.0_wp,0.0_wp,r_qj=q_stage,                &
       r_nh_term_impl=filtered_implicit_source,group=filter_groups(3))
  complex_q=CMPLX(q_stage,0.0_wp,wp)
  complex_q(7)=CMPLX(q_stage(7),1.0E-12_wp,wp)
  CALL eval_source_implicit(0.0_wp,0.0_wp,0.0_wp,c_qj=complex_q,             &
       c_nh_term_impl=complex_source,group=filter_groups(3))
  IF (ABS(REAL(complex_source(7),wp)-filtered_implicit_source(7)) .GT.      &
       1.0E-12_wp*MAX(1.0_wp,ABS(source(7))) .OR.                           &
       ABS(AIMAG(complex_source(7)))==0.0_wp .OR.                           &
       ANY(complex_source([1,2,3,4,5,6])/=(0.0_wp,0.0_wp))) THEN
     ERROR STOP 'complex pore-pressure solver selection failed'
  END IF

  CALL qc_to_qp(q_stage,qp_known,dynamic_pressure)
  CALL eval_source_semi_implicit(0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,       &
       grav,qp_known,si_source,0.0_wp)
  CALL eval_source_semi_implicit(0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,       &
       grav,qp_known,filtered_semi_implicit_source,0.0_wp,                 &
       group=filter_groups(1))
  IF (MAXVAL(ABS(filtered_semi_implicit_source-si_source)) .GT.            &
       1.0E-12_wp*MAX(1.0_wp,MAXVAL(ABS(si_source)))) THEN
     ERROR STOP 'semi-implicit rheology solver selection changed the rate'
  END IF
  CALL eval_source_semi_implicit(0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,       &
       grav,qp_known,filtered_semi_implicit_source,0.0_wp,                 &
       group=filter_groups(2))
  IF (ANY(filtered_semi_implicit_source/=0.0_wp)) THEN
     ERROR STOP 'thermal solver selection leaked a semi-implicit rate'
  END IF
  CALL eval_source_semi_implicit(0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,       &
       grav,qp_known,filtered_semi_implicit_source,0.0_wp,                 &
       group=filter_groups(3))
  IF (ANY(filtered_semi_implicit_source/=0.0_wp)) THEN
     ERROR STOP 'pore-pressure solver selection leaked a semi-implicit rate'
  END IF
  DO i = 1, SIZE(partition%implicit_map)
     IF (partition%implicit_map(i) .EQ. 4 .OR.                             &
          partition%implicit_map(i) .EQ. 7) CYCLE
     residual = q_stage(partition%implicit_map(i)) -                       &
          q_old(partition%implicit_map(i)) - dt_step*a_diag*                &
          source(partition%implicit_map(i))
     scale = MAX(1.0_wp, ABS(q_old(partition%implicit_map(i))))
     IF (ABS(residual) / scale .GT. 5.0E-11_wp) THEN
        WRITE(*,*) 'FAIL: stage residual', partition%implicit_map(i),       &
             residual/scale
        ERROR STOP 'grouped implicit stage residual is too large'
     END IF
  END DO

  IF (q_stage(4) .GE. q_old(4) .OR. q_stage(7) .GE. q_old(7)) THEN
     ERROR STOP 'analytic cooling or pore-pressure group did not relax'
  END IF
  IF (q_stage(1) .NE. q_old(1) .OR.                                    &
       ANY(q_stage(5:6) .NE. q_old(5:6))) THEN
     ERROR STOP 'group solves modified equations outside their partitions'
  END IF

  q_stage_check = q_old
  q_stage_check(2:3) = 0.0_wp
  CALL qc_to_qp(q_stage_check, qp_known, dynamic_pressure)
  CALL solve_imex_stage(partition, source_registry%stage_groups(),           &
       q_stage_check, q_stage_check,              &
       dt_step, a_diag, known_residual, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,    &
       0.0_wp, grav, 0.0_wp, .FALSE., .FALSE., q_stage,                     &
       effective_stage_term, iterations, converged, linear_info,             &
       line_search_failed)
  IF (.NOT. converged .OR. linear_info .NE. 0 .OR. line_search_failed) THEN
     ERROR STOP 'IMEX stage coordinator did not converge'
  END IF
  stage_rate_error = MAXVAL(ABS(q_stage - q_stage_check -                  &
       dt_step*a_diag*effective_stage_term(1:7)))
  scale = MAX(1.0_wp, MAXVAL(ABS(q_stage)))
  IF (stage_rate_error / scale .GT. 5.0E-14_wp) THEN
     WRITE(*,*) 'FAIL: effective IM stage rate', stage_rate_error / scale
     ERROR STOP 'effective IM rate does not match the net stage increment'
  END IF

  CALL initialize_configured_sources(source_registry,partition,source_error)
  IF (LEN_TRIM(source_error)>0) ERROR STOP TRIM(source_error)
  IF (partition%IM_group_method(1) .NE. IM_GROUP_METHOD_ANALYTIC .OR.      &
       partition%IM_group_method(2) .NE. IM_GROUP_METHOD_ANALYTIC .OR.      &
       partition%IM_group_method(3) .NE. IM_GROUP_METHOD_NEWTON) THEN
     ERROR STOP 'equation groups do not select only analytic or Newton methods'
  END IF

  CALL qc_to_qp(q_stage_check, qp_known, dynamic_pressure)
  areal_heat_capacity = q_stage_check(4) / qp_known(4)
  CALL solve_thermal_cooling_analytic(qp_known(4), dt_step*a_diag,         &
       areal_heat_capacity, radiative_term_coeff, convective_term_coeff,    &
       T_env, exact_temperature)
  exact_energy = areal_heat_capacity * exact_temperature
  q_expected = q_stage_check
  q_expected(4) = exact_energy
  CALL r_phys_var(q_expected, expected_h, expected_u, expected_v,           &
       expected_alphas, expected_density, expected_temperature,              &
       expected_alphag, expected_red_grav, expected_p_dyn,   &
       expected_Zs, expected_pore_pressure)
  pressure_known = q_expected(7) / q_expected(1)
  CALL solve_pore_pressure_dissipation_analytic(pressure_known, dt_step*a_diag,          &
       expected_h, expected_temperature, expected_alphas, expected_pressure)
  CALL solve_imex_stage(partition, source_registry%stage_groups(),            &
       q_stage_check, q_stage_check,               &
       dt_step, a_diag, known_residual, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,    &
       0.0_wp, grav, 0.0_wp, .FALSE., .FALSE., q_stage,                     &
       effective_stage_term, iterations, converged, linear_info,             &
       line_search_failed)
  IF (.NOT. converged .OR. linear_info .NE. 0 .OR. line_search_failed) THEN
     ERROR STOP 'analytic DIRK lava stage did not converge'
  END IF
  scale = MAX(1.0_wp, ABS(exact_energy))
  IF (ABS(q_stage(4) - exact_energy) / scale .GT. 5.0E-14_wp) THEN
     ERROR STOP 'lava analytic group did not update the energy conservatively'
  END IF
  scale = MAX(1.0_wp, ABS(q_expected(1) * expected_pressure))
  IF (ABS(q_stage(7) - q_expected(1) * expected_pressure) / scale .GT.      &
       5.0E-14_wp) THEN
     WRITE(*,*) 'pore exact-flow result, expected:', q_stage(7),             &
          q_expected(1)*expected_pressure
     ERROR STOP 'pore-pressure analytic group did not use the updated state'
  END IF
  IF (q_stage(1) .NE. q_stage_check(1) .OR.                               &
       ANY(q_stage(5:6) .NE. q_stage_check(5:6))) THEN
     ERROR STOP 'lava analytic group modified unrelated conservative variables'
  END IF
  stage_rate_error = MAXVAL(ABS(q_stage - q_stage_check -                  &
       dt_step*a_diag*effective_stage_term(1:7)))
  scale = MAX(1.0_wp, MAXVAL(ABS(q_stage)))
  IF (stage_rate_error / scale .GT. 5.0E-14_wp) THEN
     ERROR STOP 'analytic effective rate does not match its net stage increment'
  END IF

  ! A combined ANALYTIC group is forbidden: analytical maps are
  ! individual physical operators, even for disjoint scalar equations.
  thermal_method='ANALYTIC'; thermal_group=1
  pore_pressure_method='ANALYTIC'; pore_pressure_group=1
  rheology_method='ANALYTIC'; rheology_group=2
  CALL initialize_configured_sources(source_registry,partition,source_error)
  IF (LEN_TRIM(source_error)==0) THEN
     ERROR STOP 'shared analytic stage group should be rejected'
  END IF

  ! Separate ordered groups preserve the sequential exact-flow checks.
  pore_pressure_group=2
  rheology_group=3
  CALL initialize_configured_sources(source_registry,partition,source_error)
  IF (LEN_TRIM(source_error)>0) ERROR STOP TRIM(source_error)
  IF (partition%n_IM_groups/=3 .OR. &
       partition%IM_group(4)/=1 .OR. partition%IM_group(7)/=2 .OR. &
       partition%IM_group(2)/=3 .OR. &
       ANY(partition%IM_group_method(1:3)/=IM_GROUP_METHOD_ANALYTIC)) THEN
     ERROR STOP 'separate analytic stage groups were not preserved'
  END IF
  CALL solve_imex_stage(partition, source_registry%stage_groups(), &
       q_stage_check, q_stage_check, dt_step, a_diag, known_residual, &
       0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, grav, 0.0_wp, &
       .FALSE., .FALSE., q_stage, effective_stage_term, iterations, &
       converged, linear_info, line_search_failed)
  IF (.NOT.converged .OR. linear_info/=0 .OR. line_search_failed) &
       ERROR STOP 'separate analytic stage groups did not converge'
  IF (ABS(q_stage(4)-exact_energy)/MAX(1.0_wp,ABS(exact_energy))> &
       5.0E-14_wp .OR. &
       ABS(q_stage(7)-q_expected(1)*expected_pressure)/ &
       MAX(1.0_wp,ABS(q_expected(1)*expected_pressure))>5.0E-14_wp) THEN
     ERROR STOP 'separate analytic groups produced inconsistent stage'
  END IF

  thermal_method='NEWTON'; thermal_group=2
  pore_pressure_method='ANALYTIC'; pore_pressure_group=1
  rheology_method='ANALYTIC'; rheology_group=3
  CALL initialize_configured_sources(source_registry,partition,source_error)
  IF (LEN_TRIM(source_error)>0) ERROR STOP TRIM(source_error)
  IF (partition%n_IM_groups .NE. 3 .OR.                               &
       partition%IM_group(7) .NE. 1 .OR.                               &
       partition%IM_group(4) .NE. 2 .OR.                               &
       partition%IM_group(2) .NE. 3 .OR.                               &
       partition%IM_group(3) .NE. 3 .OR.                               &
       partition%IM_group_method(1) .NE. IM_GROUP_METHOD_ANALYTIC .OR. &
       partition%IM_group_method(2) .NE. IM_GROUP_METHOD_NEWTON .OR.   &
       partition%IM_group_method(3) .NE. IM_GROUP_METHOD_ANALYTIC) THEN
     ERROR STOP 'configured IMEX group order or methods were not preserved'
  END IF

  thermal_method='NEWTON'; thermal_group=1
  pore_pressure_method='NEWTON'; pore_pressure_group=1
  rheology_method='NEWTON'; rheology_group=1
  CALL initialize_configured_sources(source_registry,partition,source_error)
  IF (LEN_TRIM(source_error)>0) ERROR STOP TRIM(source_error)
  IF (partition%n_IM_groups .NE. 1 .OR.                               &
       partition%n_IM_group_equations(1) .NE. 4 .OR.                   &
       partition%IM_group_method(1) .NE. IM_GROUP_METHOD_NEWTON) THEN
     ERROR STOP 'all implicit quantities could not be assigned to one Newton group'
  END IF

  CALL qc_to_qp(q_stage_check, qp_known, dynamic_pressure)
  CALL solve_imex_stage(partition, source_registry%stage_groups(),           &
       q_stage_check, q_stage_check,              &
       dt_step, a_diag, known_residual, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,    &
       0.0_wp, grav, 0.0_wp, .FALSE., .FALSE., q_stage,                     &
       effective_stage_term, iterations, converged, linear_info,             &
       line_search_failed)
  IF (.NOT. converged .OR. linear_info .NE. 0 .OR. line_search_failed) THEN
     ERROR STOP 'single Newton group did not converge'
  END IF
  CALL eval_source_implicit(0.0_wp, 0.0_wp, 0.0_wp, r_qj=q_stage,            &
       r_nh_term_impl=source)
  DO i = 1, SIZE(partition%implicit_map)
     equation = partition%implicit_map(i)
     residual = q_stage(equation) - q_stage_check(equation) -               &
          dt_step*a_diag*source(equation)
     scale = MAX(1.0_wp, ABS(q_stage_check(equation)))
     IF (ABS(residual) / scale .GT. 2.0E-5_wp) THEN
        WRITE(*,*) 'FAIL: single Newton group residual', equation, residual/scale
        ERROR STOP 'single Newton group does not satisfy the joint stage equation'
     END IF
  END DO

  ! A large yield impulse must be constrained before Newton, not only after it.
  tau0 = 20000.0_wp
  q_stage_check = q_old
  CALL qc_to_qp(q_stage_check, qp_known, dynamic_pressure)
  CALL eval_source_semi_implicit(0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,   &
       grav, qp_known, si_source, 0.0_wp)
  q_expected = q_stage_check
  q_expected(2:3) = q_expected(2:3) +                                     &
       dt_step*a_diag*si_source(2:3)
  CALL correct_friction_momentum(q_stage_check(2:3), q_expected(2:3))
  CALL solve_imex_stage(partition, source_registry%stage_groups(),           &
       q_stage_check, q_stage_check,              &
       dt_step, a_diag, known_residual, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,    &
       0.0_wp, grav, 0.0_wp, .FALSE., .FALSE., q_stage,                     &
       effective_stage_term, iterations, converged, linear_info,             &
       line_search_failed)
  IF (.NOT. converged .OR. linear_info .NE. 0 .OR. line_search_failed) THEN
     ERROR STOP 'large semi-implicit yield stage did not converge'
  END IF
  scale = MAX(1.0_wp, MAXVAL(ABS(q_expected(2:3))))
  IF (MAXVAL(ABS(q_stage(2:3) - q_expected(2:3))) / scale .GT. 2.0E-12_wp) THEN
     ERROR STOP 'semi-implicit yield friction reversed the momentum'
  END IF
  tau0 = 0.0_wp

  q_stage_check = q_old
  CALL qc_to_qp(q_stage_check, qp_known, dynamic_pressure)
  CALL eval_source_semi_implicit(0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,   &
       grav, qp_known, si_source, 0.0_wp)
  CALL solve_imex_stage(partition, source_registry%stage_groups(),          &
       q_stage_check, q_stage_check,             &
       dt_step, a_diag, known_residual, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,   &
       0.0_wp, grav, 0.0_wp, .FALSE., .FALSE., q_stage,                    &
       effective_stage_term, iterations, converged, linear_info,            &
       line_search_failed)
  IF (.NOT. converged .OR. linear_info .NE. 0 .OR. line_search_failed) THEN
     ERROR STOP 'Newton stage with SI momentum source did not converge'
  END IF
  CALL eval_source_implicit(0.0_wp, 0.0_wp, 0.0_wp, r_qj=q_stage,            &
       r_nh_term_impl=source)
  scale = MAX(1.0_wp, MAXVAL(ABS(q_stage(2:3))))
  IF (MAXVAL(ABS(q_stage(2:3) - q_stage_check(2:3) -                       &
       dt_step*a_diag*(si_source(2:3)+source(2:3)))) / scale .GT. 2.0E-5_wp) THEN
     ERROR STOP 'Newton momentum group did not combine SI and implicit sources'
  END IF

  rheology_method='ANALYTIC'; rheology_group=1
  thermal_method='ANALYTIC'; thermal_group=2
  pore_pressure_method='ANALYTIC'; pore_pressure_group=3
  CALL initialize_configured_sources(source_registry,partition,source_error)
  IF (LEN_TRIM(source_error)>0) ERROR STOP TRIM(source_error)
  q_stage_check = q_old
  CALL qc_to_qp(q_stage_check, qp_known, dynamic_pressure)
  CALL r_phys_var(q_stage_check, expected_h, expected_u, expected_v,      &
       expected_alphas, expected_density, expected_temperature,           &
       expected_alphag, expected_red_grav, expected_p_dyn,&
       expected_Zs, expected_pore_pressure)
  speed_ratio = 1.0_wp
  CALL solve_rheology_voellmy_analytic(dt_step*a_diag, expected_h,         &
       speed_ratio, [expected_u, expected_v], expected_density,             &
       expected_red_grav, grav, 0.0_wp, pore_pressure_flag,               &
       expected_pore_pressure, .FALSE., 0.0_wp, analytic_velocity)
  CALL solve_imex_stage(partition, source_registry%stage_groups(),          &
       q_stage_check, q_stage_check,             &
       dt_step, a_diag, known_residual, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,   &
       0.0_wp, grav, 0.0_wp, .FALSE., .FALSE., q_stage,                    &
       effective_stage_term, iterations, converged, linear_info,            &
       line_search_failed)
  IF (.NOT. converged .OR. linear_info .NE. 0 .OR. line_search_failed) THEN
     ERROR STOP 'analytic rheology group did not complete the stage'
  END IF
  scale = MAX(1.0_wp, MAXVAL(ABS(q_stage(2:3))))
  IF (MAXVAL(ABS(q_stage(2:3) - q_stage_check(1)*analytic_velocity)) /     &
       scale .GT. 5.0E-14_wp) THEN
     ERROR STOP 'analytic rheology group did not replace the SI+implicit solve'
  END IF

  CALL test_ordered_multi_solver_group
  CALL test_multi_source_newton_step
  CALL test_lava_newton_stage_cooling_analytic_step

  CALL partition%finalize
  CALL source_registry%finalize
  CALL multi_partition%finalize
  CALL multi_registry%finalize
  CALL finalize_nonlinear_solver
  WRITE(*,*) 'PASS: grouped DIRK and effective IMEX stage-rate contracts'

CONTAINS
  ! D18: lava momentum solved during DIRK, cooling once after the RK step.
  SUBROUTINE test_lava_newton_stage_cooling_analytic_step
    TYPE(source_group_type), ALLOCATABLE :: stage_groups(:),step_groups(:)
    TYPE(analytic_term_context_type) :: context
    REAL(wp) :: state_before(7),state_stage(7),state_cooled(7)
    REAL(wp) :: stage_source(7),heat_capacity,reference_temperature
    INTEGER :: iters,info,status,g
    LOGICAL :: success,failed,found_cooling
    CHARACTER(LEN=256) :: error

    rheology_method='NEWTON';rheology_cadence='STAGE';rheology_group=1
    thermal_method='ANALYTIC';thermal_cadence='STEP';thermal_group=2
    pore_pressure_method='ANALYTIC';pore_pressure_cadence='STAGE';pore_pressure_group=3
    CALL initialize_configured_sources(source_registry,partition,error)
    IF (LEN_TRIM(error)>0) ERROR STOP TRIM(error)
    stage_groups=source_registry%stage_groups()
    step_groups=source_registry%step_groups()
    IF (partition%n_IM_groups/=2) ERROR STOP 'D18: cooling entered DIRK partition'
    IF (ANY(partition%IM_group(4:4)/=0)) ERROR STOP 'D18: energy in DIRK partition'
    IF (partition%IM_group_method(1)/=IM_GROUP_METHOD_NEWTON) &
         ERROR STOP 'D18: lava stage no longer Newton'
    found_cooling=.FALSE.
    DO g=1,SIZE(step_groups)
      IF (.NOT.ANY(step_groups(g)%source_names=='thermal_cooling')) CYCLE
      IF (step_groups(g)%method/=SOURCE_METHOD_ANALYTIC .OR. &
           step_groups(g)%cadence/=SOURCE_CADENCE_STEP) &
           ERROR STOP 'D18: incorrect cooling STEP triplet'
      found_cooling=.TRUE.
    END DO
    IF (.NOT.found_cooling) ERROR STOP 'D18: cooling missing from STEP schedule'
    state_before=q_old
    state_stage=state_before
    CALL solve_imex_stage(partition,stage_groups,state_before,state_before, &
         dt_step,a_diag,known_residual,0.0_wp,0.0_wp,0.0_wp,0.0_wp, &
         0.0_wp,grav,0.0_wp,.FALSE.,.FALSE.,state_stage,stage_source, &
         iters,success,info,failed)
    IF (.NOT.success .OR. info/=0 .OR. failed) &
         ERROR STOP 'D18: lava Newton/STAGE did not converge'
    IF (MAXVAL(ABS(state_stage(2:3)-state_before(2:3)))<=0.0_wp) &
         ERROR STOP 'D18: Newton/STAGE did not update lava momentum'
    IF (state_stage(4)/=state_before(4)) &
         ERROR STOP 'D18: STEP cooling erroneously applied during DIRK'
    CALL qc_to_qp(state_stage,qp_known,dynamic_pressure)
    heat_capacity=state_stage(4)/qp_known(4)
    CALL solve_thermal_cooling_analytic(qp_known(4),dt_step,heat_capacity, &
         radiative_term_coeff,convective_term_coeff,T_env,reference_temperature)
    DO g=1,SIZE(step_groups)
      IF (.NOT.ANY(step_groups(g)%source_names=='thermal_cooling')) CYCLE
      context=analytic_term_context_type()
      CALL eval_source_analytic(step_groups(g),state_stage,context,dt_step, &
           state_cooled,status)
      IF (status/=0) ERROR STOP 'D18: analytic cooling STEP failed'
      EXIT
    END DO
    IF (ABS(state_cooled(4)-heat_capacity*reference_temperature)> &
         5.0E-14_wp*MAX(1.0_wp,ABS(state_cooled(4)))) &
         ERROR STOP 'D18: cooling STEP does not match exact flow'
    IF (ANY(state_cooled([1,2,3,5,6,7])/=state_stage([1,2,3,5,6,7]))) &
         ERROR STOP 'D18: cooling STEP modified unrelated state'
    WRITE(*,*) 'PASS D18 lava Newton/STAGE and analytic cooling/STEP'
  END SUBROUTINE test_lava_newton_stage_cooling_analytic_step

  ! D18: STEP Newton advances all selected equations in one joint solve.
  SUBROUTINE test_multi_source_newton_step
    TYPE(equation_partition_type) :: step_partition
    TYPE(source_registry_type) :: step_registry
    TYPE(source_group_type) :: step_group
    TYPE(source_group_type), ALLOCATABLE :: groups(:)
    LOGICAL :: mask(7),success,failed
    INTEGER :: memberships(7),iters,info,target
    REAL(wp) :: before(7),after(7),rates(7),normalized
    CHARACTER(LEN=256) :: error

    mask=.FALSE.;mask([2,3,4,7])=.TRUE.
    memberships=0;memberships([2,3,4,7])=1
    CALL step_partition%initialize(mask,memberships,[IM_GROUP_METHOD_NEWTON])
    step_group=filter_group
    step_group%id=1
    step_group%source_names=[CHARACTER(LEN=64):: &
         'rheology_voellmy','thermal_cooling','pore_pressure_dissipation']
    step_group%equation_ids=[2,3,4,7]
    step_group%cadence=SOURCE_CADENCE_STEP
    CALL step_registry%initialize(n_vars,idx_pore,[source_group_type::], &
         [step_group],error)
    IF (LEN_TRIM(error)>0) ERROR STOP TRIM(error)
    groups=step_registry%step_groups()
    IF (SIZE(groups)/=1) ERROR STOP 'D18: STEP Newton group split'
    before=q_old
    after=before
    CALL solve_rk_step(step_partition,groups,after,before,dt_step,1.0_wp, &
         known_residual,0.0_wp,0.0_wp,0.0_wp,iters,success,info,failed)
    IF (.NOT.success .OR. failed .OR. info/=0) &
         ERROR STOP 'D18: joint NEWTON/STEP failed'
    IF (ANY(after([1,5,6])/=before([1,5,6]))) &
         ERROR STOP 'D18: joint NEWTON/STEP changed explicit equations'
    CALL eval_source_implicit(0.0_wp,0.0_wp,0.0_wp,r_qj=after, &
         r_nh_term_impl=rates,group=groups(1))
    DO target=1,4
      SELECT CASE(target)
      CASE(1);equation=2
      CASE(2);equation=3
      CASE(3);equation=4
      CASE(4);equation=7
      END SELECT
      normalized=ABS(after(equation)-before(equation)-dt_step*rates(equation))/ &
           MAX(1.0_wp,ABS(before(equation)))
      IF (normalized>2.0E-5_wp) THEN
        WRITE(*,*) 'FAIL NEWTON/STEP multi-source residual',equation,normalized
        ERROR STOP 'D18: joint NEWTON/STEP residual'
      END IF
    END DO
    CALL step_registry%finalize
    CALL step_partition%finalize
    WRITE(*,*) 'PASS D18 multi-source coupled NEWTON/STEP'
  END SUBROUTINE test_multi_source_newton_step

  SUBROUTINE test_ordered_multi_solver_group
    TYPE(source_group_type), ALLOCATABLE :: groups(:)
    LOGICAL :: multi_mask(7), multi_converged, multi_line_search_failed
    INTEGER :: multi_map(7), multi_iterations, multi_linear_info, target
    REAL(wp) :: multi_state(7), multi_stage(7), multi_effective_rate(7)
    REAL(wp) :: joint_rate(7), normalized_error
    CHARACTER(LEN=256) :: error

    multi_mask=.FALSE.
    multi_mask([2,3,4,7])=.TRUE.
    multi_map=0
    multi_map([2,3,4,7])=1
    CALL multi_partition%initialize(multi_mask,multi_map, &
         [IM_GROUP_METHOD_NEWTON])
    ! One joint Newton group owns all coupled velocity, temperature and
    ! pore-pressure equations; no ordered internal solver schedule exists.
    CALL multi_registry%initialize(n_vars,idx_pore,[filter_group], &
         [source_group_type::],error)
    IF (LEN_TRIM(error)>0) ERROR STOP TRIM(error)
    groups=multi_registry%stage_groups()
    multi_state=q_old
    CALL solve_imex_stage(multi_partition,groups,multi_state,q_old, &
         dt_step,a_diag,known_residual,0.0_wp,0.0_wp,0.0_wp,0.0_wp, &
         0.0_wp,grav,0.0_wp,.FALSE.,.FALSE.,multi_stage, &
         multi_effective_rate,multi_iterations,multi_converged, &
         multi_linear_info,multi_line_search_failed)
    IF (.NOT.multi_converged .OR. multi_linear_info/=0 .OR. &
         multi_line_search_failed) THEN
       ERROR STOP 'joint Newton group did not converge'
    END IF
    ! The stage coordinator includes the semi-implicit momentum increment
    ! before the joint Newton solve. The physical implicit law alone is not
    ! the complete stage rate when lava yield/friction is active.
    CALL qc_to_qp(q_old,qp_known,dynamic_pressure)
    CALL eval_source_semi_implicit(0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp, &
         grav,qp_known,si_source,0.0_wp,group=groups(1))
    CALL eval_source_implicit(0.0_wp,0.0_wp,0.0_wp,r_qj=multi_stage, &
         r_nh_term_impl=joint_rate,group=groups(1))
    DO target=1,4
       SELECT CASE(target)
       CASE(1); equation=2
       CASE(2); equation=3
       CASE(3); equation=4
       CASE(4); equation=7
       END SELECT
       ! Independently check the complete rate, including the SI momentum
       ! contribution, against the stage increment. Numerical projection
       ! of frictional momentum is covered by the separate correction tests.
       normalized_error=ABS(multi_stage(equation)-q_old(equation)- &
            dt_step*a_diag*(joint_rate(equation)+si_source(equation)))/ &
            MAX(1.0_wp,ABS(q_old(equation)))
       IF (normalized_error>2.0E-5_wp) THEN
          WRITE(*,*) 'FAIL: joint residual',equation,normalized_error, &
               joint_rate(equation),si_source(equation)
          ERROR STOP 'joint Newton group residual'
       END IF
    END DO
  END SUBROUTINE test_ordered_multi_solver_group



END PROGRAM test_grouped_implicit_stage
