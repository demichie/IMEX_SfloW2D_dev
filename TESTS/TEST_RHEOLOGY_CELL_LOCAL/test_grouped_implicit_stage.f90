PROGRAM test_grouped_implicit_stage

  USE parameters_2d
  USE constitutive_parameters_2d
  USE rheology_parameters_2d, ONLY : T_ref
  USE rheology_lava_2d, ONLY : nu_ref, visc_par, tau0
  USE geometry_2d, ONLY : pi_g
  USE equation_metadata_2d, ONLY : equation_partition_type,                 &
       IM_GROUP_METHOD_ANALYTIC, IM_GROUP_METHOD_NEWTON
  USE equation_terms_2d, ONLY : eval_implicit_terms,                     &
       eval_semi_implicit_terms, init_problem_param
  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver,               &
       finalize_nonlinear_solver, solve_imex_stage, solve_rk_step
  USE state_conversion_2d, ONLY : qc_to_qp, r_phys_var
  USE pore_pressure_laws_2d, ONLY : solve_analytic_pore_pressure
  USE thermal_laws_2d, ONLY : solve_lava_analytic_temperature
  USE rheology_lava_2d, ONLY : solve_lava_analytic_velocity

  IMPLICIT NONE

  TYPE(equation_partition_type) :: partition
  LOGICAL :: implicit_mask(7), converged, line_search_failed
  INTEGER :: group_map(7), iterations, linear_info, i, equation
  REAL(wp) :: q_old(7), q_stage(7), known_residual(7), source(7), si_source(7)
  REAL(wp) :: q_stage_check(7), q_expected(7), effective_stage_term(7)
  REAL(wp) :: qp_known(9), qp_expected(9), qp_actual(9)
  REAL(wp) :: dynamic_pressure, stage_rate_error
  REAL(wp) :: areal_heat_capacity, exact_temperature, exact_energy
  REAL(wp) :: pressure_known, expected_pressure
  REAL(wp) :: expected_h, expected_u, expected_v, expected_density
  REAL(wp) :: expected_alphas(2), expected_alphag(0), expected_alphal
  REAL(wp) :: expected_temperature, expected_red_grav, expected_p_dyn
  REAL(wp) :: expected_Zs, expected_pore_pressure
  REAL(wp) :: speed_ratio, analytic_velocity(2)
  REAL(wp) :: residual, scale, dt_step, a_diag

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
  rheology_model_id = rheology_id_lava
  pore_pressure_flag = .TRUE.
  gas_flag = .FALSE.
  liquid_flag = .TRUE.
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

  q_old = [1200.0_wp, 2400.0_wp, -1200.0_wp, 4.0E9_wp,                     &
       180.0_wp, 120.0_wp, 6.0E6_wp]
  q_stage = q_old
  known_residual = 0.0_wp
  dt_step = 0.2_wp
  a_diag = 0.75_wp
  CALL initialize_nonlinear_solver
  CALL solve_rk_step(partition, q_stage, q_old, dt_step, a_diag,            &
       known_residual, 0.0_wp, 0.0_wp, 0.0_wp, iterations, converged,      &
       linear_info, line_search_failed)
  IF (.NOT. converged .OR. linear_info .NE. 0 .OR. line_search_failed) THEN
     WRITE(*,*) 'FAIL: grouped solve status', converged, linear_info,       &
          line_search_failed, iterations
     ERROR STOP 'grouped implicit stage did not converge'
  END IF

  CALL eval_implicit_terms(0.0_wp, 0.0_wp, 0.0_wp, r_qj=q_stage,            &
       r_nh_term_impl=source)
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
  CALL solve_imex_stage(partition, q_stage_check, q_stage_check, qp_known,  &
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

  CALL init_problem_param(partition)
  IF (partition%IM_group_method(1) .NE. IM_GROUP_METHOD_ANALYTIC .OR.      &
       partition%IM_group_method(2) .NE. IM_GROUP_METHOD_ANALYTIC .OR.      &
       partition%IM_group_method(3) .NE. IM_GROUP_METHOD_NEWTON) THEN
     ERROR STOP 'equation groups do not select only analytic or Newton methods'
  END IF

  CALL qc_to_qp(q_stage_check, qp_known, dynamic_pressure)
  areal_heat_capacity = q_stage_check(4) / qp_known(4)
  CALL solve_lava_analytic_temperature(qp_known(4), dt_step*a_diag,         &
       areal_heat_capacity, radiative_term_coeff, convective_term_coeff,    &
       T_env, exact_temperature)
  exact_energy = areal_heat_capacity * exact_temperature
  q_expected = q_stage_check
  q_expected(4) = exact_energy
  CALL r_phys_var(q_expected, expected_h, expected_u, expected_v,           &
       expected_alphas, expected_density, expected_temperature,              &
       expected_alphal, expected_alphag, expected_red_grav, expected_p_dyn,   &
       expected_Zs, expected_pore_pressure)
  pressure_known = q_expected(7) / q_expected(1)
  CALL solve_analytic_pore_pressure(pressure_known, dt_step*a_diag,          &
       expected_h, expected_temperature, expected_alphas, expected_pressure)
  CALL solve_imex_stage(partition, q_stage_check, q_stage_check, qp_known,   &
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

  imex_energy_group = 1
  imex_pore_pressure_group = 1
  imex_velocity_group = 2
  imex_group_method(1:2) = [CHARACTER(LEN=16) :: 'ANALYTIC', 'ANALYTIC']
  CALL init_problem_param(partition)
  IF (partition%n_IM_groups .NE. 2 .OR.                                  &
       partition%n_IM_group_equations(1) .NE. 2 .OR.                      &
       partition%IM_group(4) .NE. 1 .OR.                                 &
       partition%IM_group(7) .NE. 1 .OR.                                 &
       partition%IM_group(2) .NE. 2 .OR.                                 &
       partition%IM_group_method(1) .NE. IM_GROUP_METHOD_ANALYTIC .OR.   &
       partition%IM_group_method(2) .NE. IM_GROUP_METHOD_ANALYTIC) THEN
     ERROR STOP 'energy and pore pressure could not share an analytic group'
  END IF
  CALL solve_imex_stage(partition, q_stage_check, q_stage_check, qp_known, &
       dt_step, a_diag, known_residual, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,   &
       0.0_wp, grav, 0.0_wp, .FALSE., .FALSE., q_stage,                    &
       effective_stage_term, iterations, converged, linear_info,            &
       line_search_failed)
  IF (.NOT. converged .OR. linear_info .NE. 0 .OR. line_search_failed) THEN
     ERROR STOP 'mixed analytic group stage did not converge'
  END IF
  IF (ABS(q_stage(4) - exact_energy) / MAX(1.0_wp, ABS(exact_energy)) .GT. &
       5.0E-14_wp .OR.                                                   &
       ABS(q_stage(7) - q_expected(1)*expected_pressure) /                &
       MAX(1.0_wp, ABS(q_expected(1)*expected_pressure)) .GT. 5.0E-14_wp) THEN
     ERROR STOP 'shared analytic group produced an inconsistent stage'
  END IF

  imex_energy_group = 2
  imex_pore_pressure_group = 1
  imex_velocity_group = 3
  imex_group_method = [CHARACTER(LEN=16) :: 'ANALYTIC', 'NEWTON', 'ANALYTIC']
  CALL init_problem_param(partition)
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

  imex_energy_group = 1
  imex_pore_pressure_group = 1
  imex_velocity_group = 1
  imex_group_method(1) = 'newton'
  CALL init_problem_param(partition)
  IF (partition%n_IM_groups .NE. 1 .OR.                               &
       partition%n_IM_group_equations(1) .NE. 4 .OR.                   &
       partition%IM_group_method(1) .NE. IM_GROUP_METHOD_NEWTON) THEN
     ERROR STOP 'all implicit quantities could not be assigned to one Newton group'
  END IF

  CALL qc_to_qp(q_stage_check, qp_known, dynamic_pressure)
  CALL solve_imex_stage(partition, q_stage_check, q_stage_check, qp_known,  &
       dt_step, a_diag, known_residual, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,    &
       0.0_wp, grav, 0.0_wp, .FALSE., .FALSE., q_stage,                     &
       effective_stage_term, iterations, converged, linear_info,             &
       line_search_failed)
  IF (.NOT. converged .OR. linear_info .NE. 0 .OR. line_search_failed) THEN
     ERROR STOP 'single Newton group did not converge'
  END IF
  CALL eval_implicit_terms(0.0_wp, 0.0_wp, 0.0_wp, r_qj=q_stage,            &
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

  q_stage_check = q_old
  CALL qc_to_qp(q_stage_check, qp_known, dynamic_pressure)
  CALL eval_semi_implicit_terms(0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,   &
       grav, qp_known, si_source, 0.0_wp)
  CALL solve_imex_stage(partition, q_stage_check, q_stage_check, qp_known, &
       dt_step, a_diag, known_residual, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,   &
       0.0_wp, grav, 0.0_wp, .FALSE., .FALSE., q_stage,                    &
       effective_stage_term, iterations, converged, linear_info,            &
       line_search_failed)
  IF (.NOT. converged .OR. linear_info .NE. 0 .OR. line_search_failed) THEN
     ERROR STOP 'Newton stage with SI momentum source did not converge'
  END IF
  CALL eval_implicit_terms(0.0_wp, 0.0_wp, 0.0_wp, r_qj=q_stage,            &
       r_nh_term_impl=source)
  scale = MAX(1.0_wp, MAXVAL(ABS(q_stage(2:3))))
  IF (MAXVAL(ABS(q_stage(2:3) - q_stage_check(2:3) -                       &
       dt_step*a_diag*(si_source(2:3)+source(2:3)))) / scale .GT. 2.0E-5_wp) THEN
     ERROR STOP 'Newton momentum group did not combine SI and implicit sources'
  END IF

  imex_velocity_group = 1
  imex_energy_group = 2
  imex_pore_pressure_group = 3
  imex_group_method = [CHARACTER(LEN=16) ::                             &
       'ANALYTIC', 'ANALYTIC', 'ANALYTIC']
  CALL init_problem_param(partition)
  q_stage_check = q_old
  CALL qc_to_qp(q_stage_check, qp_known, dynamic_pressure)
  CALL r_phys_var(q_stage_check, expected_h, expected_u, expected_v,      &
       expected_alphas, expected_density, expected_temperature,           &
       expected_alphal, expected_alphag, expected_red_grav, expected_p_dyn,&
       expected_Zs, expected_pore_pressure)
  speed_ratio = 1.0_wp
  CALL solve_lava_analytic_velocity(dt_step*a_diag, expected_h,            &
       expected_temperature, expected_density, speed_ratio,                &
       [expected_u, expected_v], analytic_velocity)
  CALL solve_imex_stage(partition, q_stage_check, q_stage_check, qp_known, &
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

  CALL partition%finalize
  CALL finalize_nonlinear_solver
  WRITE(*,*) 'PASS: grouped DIRK and effective IMEX stage-rate contracts'

END PROGRAM test_grouped_implicit_stage
