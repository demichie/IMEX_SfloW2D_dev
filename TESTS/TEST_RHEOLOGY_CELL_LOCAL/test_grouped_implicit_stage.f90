PROGRAM test_grouped_implicit_stage

  USE parameters_2d
  USE constitutive_parameters_2d
  USE rheology_parameters_2d, ONLY : T_ref
  USE rheology_lava_2d, ONLY : nu_ref, visc_par, tau0
  USE geometry_2d, ONLY : pi_g
  USE equation_metadata_2d, ONLY : equation_partition_type,                 &
       IM_GROUP_METHOD_ANALYTIC, IM_GROUP_METHOD_NEWTON
  USE equation_terms_2d, ONLY : eval_implicit_terms
  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver,               &
       finalize_nonlinear_solver, solve_imex_stage, solve_rk_step
  USE state_conversion_2d, ONLY : qc_to_qp

  IMPLICIT NONE

  TYPE(equation_partition_type) :: partition
  LOGICAL :: implicit_mask(7), converged, line_search_failed
  INTEGER :: group_map(7), iterations, linear_info, i
  REAL(wp) :: q_old(7), q_stage(7), known_residual(7), source(7)
  REAL(wp) :: q_stage_check(7), effective_stage_term(7), qp_known(9)
  REAL(wp) :: dynamic_pressure, stage_rate_error
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

  CALL partition%finalize
  CALL finalize_nonlinear_solver
  WRITE(*,*) 'PASS: grouped DIRK solves and effective IMEX stage rate'

END PROGRAM test_grouped_implicit_stage
