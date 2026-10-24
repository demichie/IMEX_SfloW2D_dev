! D17.20: full production nonlinear_solver_2d stage sweep, 2x2 and 3x3.
PROGRAM test_d17_production_newton
  USE parameters_2d
  USE constitutive_parameters_2d
  USE equation_metadata_2d, ONLY : equation_partition_type, IM_GROUP_METHOD_NEWTON
  USE source_registry_2d, ONLY : source_group_type,source_registry_type
  USE source_configuration_2d, ONLY : initialize_configured_sources, &
       initialize_native_source_configuration
  USE source_input_2d
  USE source_catalog_2d, ONLY : SOURCE_METHOD_NEWTON,SOURCE_CADENCE_STAGE, &
       SOURCE_PHASE_RK_IMPLICIT
  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver,solve_rk_step,solve_imex_stage
  USE source_terms_2d, ONLY : eval_source_implicit
  USE rheology_lava_2d, ONLY : nu_ref,visc_par,tau0
  USE rheology_parameters_2d, ONLY : T_ref
  IMPLICIT NONE
  TYPE(equation_partition_type) :: partition
  TYPE(source_group_type) :: group(1)
  TYPE(source_registry_type) :: legacy_registry,native_registry
  TYPE(equation_partition_type) :: native_partition
  REAL(wp) :: legacy_q(4),native_q(4)
  CHARACTER(LEN=256) :: config_error
  INTEGER :: legacy_iterations,native_iterations,legacy_info,native_info
  LOGICAL :: legacy_converged,native_converged,legacy_search,native_search
  REAL(wp) :: q(4),old(4),known(4),dt,expected(2),rate_real(4)
  REAL(wp) :: im_term(4),raw_si(2),si_ratio,alpha
  COMPLEX(wp) :: complex_state(4),rate_complex(4)
  REAL(wp) :: cross_derivative(2)
  INTEGER :: iterations,info,k
  LOGICAL :: converged,search_failed
  n_vars=4;n_eqns=4;n_solid=0;n_add_gas=0;n_stoch_vars=0;n_pore_vars=0
  n_RK=2
  idx_totMassEqn=1;idx_uEqn=2;idx_vEqn=3;idx_engyEqn=4
  idx_solid_first=5;idx_solid_last=4
  idx_solidEqn_first=5;idx_solidEqn_last=4
  idx_add_gas_first=5;idx_add_gas_last=4
  idx_pore=0;idx_poreEqn=0;idx_stoch=0
  idx_u=5;idx_v=6
  rheology_flag=.TRUE.;rheology_model='newmodel'
  pore_pressure_flag=.FALSE.;thermal_cooling_flag=.FALSE.
  slope_correction_flag=.FALSE.;curvature_term_flag=.FALSE.
  stochastic_flag=.FALSE.;stoch_transport_flag=.FALSE.
  liquid_flag=.TRUE.;gas_flag=.FALSE.
  rho_l=1000.0_wp;inv_rho_l=1.0_wp/rho_l;sp_heat_l=4000.0_wp
  sp_heat_a=1005.0_wp;sp_gas_const_a=287.0_wp
  grav=9.81_wp;rho_a_amb=1.2_wp
  eps_sing=1.0E-8_wp;eps_sing4=eps_sing**4
  tol_abs=1.0E-10_wp;tol_rel=1.0E-12_wp;max_nl_iter=30
  ALLOCATE(rho_s(0),inv_rho_s(0),sp_heat_s(0),diam_s(0),sphericity_s(0))
  ALLOCATE(sp_heat_g(0),sp_gas_const_g(0))
  CALL initialize_nonlinear_solver
  DO k=1,2
    IF(k==1) THEN
      CALL partition%initialize([.FALSE.,.TRUE.,.TRUE.,.FALSE.], &
           [0,1,1,0],[IM_GROUP_METHOD_NEWTON])
      group(1)%equation_ids=[2,3]
      group(1)%source_names=[CHARACTER(LEN=64)::'rheology_newmodel']
    ELSE
      CALL partition%initialize([.FALSE.,.TRUE.,.TRUE.,.TRUE.], &
           [0,1,1,1],[IM_GROUP_METHOD_NEWTON])
      group(1)%equation_ids=[2,3,4]
      group(1)%source_names=[CHARACTER(LEN=64)::'rheology_newmodel','thermal_cooling']
    END IF
    group(1)%id=1;group(1)%cadence=SOURCE_CADENCE_STAGE
    group(1)%method=SOURCE_METHOD_NEWTON
    group(1)%phase=SOURCE_PHASE_RK_IMPLICIT
    group(1)%execution_order=1
    old=[1000.0_wp,3000.0_wp,-4000.0_wp,1.2E9_wp]
    q=old;known=0.0_wp;dt=0.2_wp
    ! A zero cooling coefficient keeps temperature invariant while the
    ! three-equation group still enters the production 3x3 pivoted solver.
    radiative_term_coeff=0.0_wp;convective_term_coeff=0.0_wp
    T_env=300.0_wp
    CALL solve_rk_step(partition,group,q,old,dt,0.75_wp,known, &
         0.0_wp,0.0_wp,0.0_wp,iterations,converged,info,search_failed)
    IF (.NOT.converged .OR. info/=0 .OR. search_failed) THEN
      PRINT *, 'FAIL production Newton status',k,iterations,info
      ERROR STOP 1
    END IF
    expected=old(2:3)/(1.0_wp+0.1_wp*dt*0.75_wp)
    CALL check(q(2:3),expected,'Newton momentum exact update')
    IF (ABS(q(1)-old(1))>1.0E-12_wp .OR. &
         ABS(q(4)-old(4))>1.0E-5_wp) &
         ERROR STOP 'non-source state changed'
    IF (iterations<1) ERROR STOP 'Newton not exercised'
    PRINT *, 'PASS D17 production Newton group dimension',SIZE(group(1)%equation_ids)
  END DO
  ! Physically coupled case: lava drag depends on temperature and thermal
  ! cooling participates in the SAME 3x3 Newton block.
  rheology_model='lava';thermal_cooling_flag=.TRUE.
  nu_ref=0.02_wp;visc_par=0.003_wp;tau0=0.0_wp;T_ref=900.0_wp
  T_env=280.0_wp;convective_term_coeff=10.0_wp
  radiative_term_coeff=0.0_wp
  CALL partition%initialize([.FALSE.,.TRUE.,.TRUE.,.TRUE.], &
       [0,1,1,1],[IM_GROUP_METHOD_NEWTON])
  group(1)%equation_ids=[2,3,4]
  group(1)%source_names=[CHARACTER(LEN=64)::'rheology_lava','thermal_cooling']
  old=[1000.0_wp,3000.0_wp,-4000.0_wp,1.2E9_wp]
  q=old;known=0.0_wp;dt=0.2_wp
  complex_state=CMPLX(old,0.0_wp,wp)
  complex_state(4)=complex_state(4)+CMPLX(0.0_wp,1.0E-10_wp,wp)
  CALL eval_source_implicit(0.0_wp,0.0_wp,0.0_wp, &
       c_qj=complex_state,c_nh_term_impl=rate_complex,group=group(1))
  cross_derivative=AIMAG(rate_complex(2:3))/1.0E-10_wp
  IF (MAXVAL(ABS(cross_derivative))<=1.0E-13_wp) &
       ERROR STOP 'lava momentum-temperature coupling absent'
  CALL solve_rk_step(partition,group,q,old,dt,0.75_wp,known, &
       0.0_wp,0.0_wp,0.0_wp,iterations,converged,info,search_failed)
  IF (.NOT.converged .OR. info/=0 .OR. search_failed) THEN
    PRINT *, 'FAIL coupled production Newton 3x3 status',iterations,info
    ERROR STOP 1
  END IF
  CALL eval_source_implicit(0.0_wp,0.0_wp,0.0_wp, &
       r_qj=q,r_nh_term_impl=rate_real,group=group(1))
  IF (ANY(ABS(q(2:4)-old(2:4)-dt*0.75_wp*rate_real(2:4))> &
       1.0E-8_wp*MAX(1.0_wp,MAXVAL(ABS(old(2:4)))))) &
       ERROR STOP 'coupled 3x3 stage residual nonzero'
  IF (ABS(q(4)-old(4))<=1.0E-8_wp) &
       ERROR STOP 'coupled thermal equation did not evolve'
  PRINT *, 'PASS D17 physically coupled Newton 3x3 lava/cooling'
  ! D17.10.1: run the production IMEX stage with the semi-implicit
  ! lava yield stress strong enough to reverse the raw momentum.
  ! The post-SI projection must stop it before the Newton solve.
  thermal_cooling_flag=.FALSE.;convective_term_coeff=0.0_wp
  CALL partition%initialize([.FALSE.,.TRUE.,.TRUE.,.FALSE.], &
       [0,1,1,0],[IM_GROUP_METHOD_NEWTON])
  group(1)%equation_ids=[2,3]
  group(1)%source_names=[CHARACTER(LEN=64)::'rheology_lava']
  tau0=50000.0_wp
  old=[1000.0_wp,3000.0_wp,-4000.0_wp,1.2E9_wp]
  q=old;known=0.0_wp;dt=0.2_wp
  si_ratio=SQRT(SUM((old(2:3)/old(1))**2))
  raw_si=-tau0*(old(2:3)/old(1))/si_ratio
  IF (DOT_PRODUCT(old(2:3),old(2:3)+dt*0.75_wp*raw_si)>=0.0_wp) &
       ERROR STOP 'fixture does not induce SI reversal'
  CALL solve_imex_stage(partition,group,old,old,dt,0.75_wp,known, &
       0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp, &
       .FALSE.,.FALSE.,q,im_term,iterations,converged,info,search_failed)
  IF (.NOT.converged .OR. info/=0 .OR. search_failed) &
       ERROR STOP 'post-SI corrected production solve failed'
  IF (MAXVAL(ABS(q(2:3)))>1.0E-6_wp) &
       ERROR STOP 'post-SI projection failed to prevent reversed momentum'
  IF (ABS(q(1)-old(1))>1.0E-12_wp .OR. ABS(q(4)-old(4))>1.0E-5_wp) &
       ERROR STOP 'post-SI projection changed unrelated state'
  PRINT *, 'PASS D17 production post-SI friction reversal guard'
  ! D17.10.1: the post-Newton projection must also preserve the
  ! original momentum ray for the strictly dissipative linear law.
  rheology_model='newmodel'
  group(1)%source_names=[CHARACTER(LEN=64)::'rheology_newmodel']
  old=[1000.0_wp,3000.0_wp,-4000.0_wp,1.2E9_wp]
  q=old;dt=0.2_wp
  CALL solve_imex_stage(partition,group,old,old,dt,0.75_wp,known, &
       0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp, &
       .FALSE.,.FALSE.,q,im_term,iterations,converged,info,search_failed)
  IF (.NOT.converged .OR. info/=0 .OR. search_failed) &
       ERROR STOP 'post-Newton production solve failed'
  alpha=DOT_PRODUCT(q(2:3),old(2:3))/SUM(old(2:3)**2)
  IF (alpha<0.0_wp .OR. alpha>1.0_wp) &
       ERROR STOP 'post-Newton projection did not preserve dissipative scale'
  CALL check(q(2:3),alpha*old(2:3),'post-Newton preserved momentum ray')
  expected=old(2:3)/(1.0_wp+0.1_wp*dt*0.75_wp)
  CALL check(q(2:3),expected,'post-Newton correct implicit drag')
  PRINT *, 'PASS D17 production post-Newton friction ray guard'
  ! Quantitative legacy/native configuration equivalence: both paths
  ! construct their own validated source groups, then call the real solver.
  rheology_flag=.TRUE.;rheology_model='newmodel'
  thermal_cooling_flag=.FALSE.;pore_pressure_flag=.FALSE.
  curvature_term_flag=.FALSE.;bottom_radial_source_flag=.FALSE.
  bottom_fissural_source_flag=.FALSE.;radial_source_flag=.FALSE.
  lateral_source_flag=.FALSE.;gas_flag=.FALSE.
  gas_loss_flag=.FALSE.;liquid_vaporization_flag=.FALSE.
  liquid_flag=.FALSE.;settling_flag=.FALSE.;entrainment_flag=.FALSE.
  erosion_coeff=0.0_wp
  imex_velocity_group=1;imex_energy_group=2;imex_pore_pressure_group=3
  imex_group_method=[CHARACTER(LEN=16)::'NEWTON','NEWTON','NEWTON']
  CALL partition%initialize([.FALSE.,.TRUE.,.TRUE.,.FALSE.], &
       [0,1,1,0],[IM_GROUP_METHOD_NEWTON])
  CALL initialize_configured_sources(legacy_registry,partition,config_error)
  IF (LEN_TRIM(config_error)>0) THEN
    PRINT *, TRIM(config_error)
    ERROR STOP 'legacy numeric configuration rejected'
  END IF
  CALL clear_source_input()
  native_source_configuration_present=.TRUE.
  source_group_count=1
  source_group_id(1)=27
  source_group_order(1)=1
  source_group_cadence(1)='STAGE'
  source_group_method(1)='NEWTON'
  source_group_term_count(1)=1
  source_group_terms(1,1)='RHEOLOGY_NEWMODEL'
  source_group_target_count(1)=1
  source_group_targets(1,1)='VELOCITY'
  CALL initialize_native_source_configuration(native_registry,native_partition,config_error)
  IF (LEN_TRIM(config_error)>0) THEN
    PRINT *, TRIM(config_error)
    ERROR STOP 'native numeric configuration rejected'
  END IF
  old=[1000.0_wp,3000.0_wp,-4000.0_wp,1.2E9_wp]
  legacy_q=old;native_q=old;known=0.0_wp;dt=0.2_wp
  CALL solve_rk_step(partition,legacy_registry%stage_groups(),legacy_q,old, &
       dt,0.75_wp,known,0.0_wp,0.0_wp,0.0_wp, &
       legacy_iterations,legacy_converged,legacy_info,legacy_search)
  CALL solve_rk_step(native_partition,native_registry%stage_groups(),native_q,old, &
       dt,0.75_wp,known,0.0_wp,0.0_wp,0.0_wp, &
       native_iterations,native_converged,native_info,native_search)
  IF (.NOT.legacy_converged .OR. .NOT.native_converged .OR. &
       legacy_info/=0 .OR. native_info/=0 .OR. &
       legacy_search .OR. native_search) &
       ERROR STOP 'legacy/native numeric Newton failed'
  IF (legacy_iterations/=native_iterations) &
       ERROR STOP 'legacy/native Newton iteration counts differ'
  IF (ANY(ABS(legacy_q-native_q)>1.0E-10_wp)) &
       ERROR STOP 'legacy/native numerical stage diverged'
  expected=old(2:3)/(1.0_wp+0.1_wp*dt*0.75_wp)
  CALL check(legacy_q(2:3),expected,'legacy/native exact stage reference')
  IF (ABS(legacy_q(1)-old(1))>1.0E-12_wp .OR. &
       ABS(legacy_q(4)-old(4))>1.0E-5_wp) &
       ERROR STOP 'legacy/native non-source equations changed'
  CALL clear_source_input()
  PRINT *, 'PASS D17 legacy/native quantitative Newton stage equivalence'
  ! Stage-local diagonal coefficient regression, including zero diagonal
  ! explicit stages and representative DIRK coefficients used by RK schemes.
  CALL partition%initialize([.FALSE.,.TRUE.,.TRUE.,.FALSE.], &
       [0,1,1,0],[IM_GROUP_METHOD_NEWTON])
  group(1)%equation_ids=[2,3]
  group(1)%source_names=[CHARACTER(LEN=64)::'rheology_newmodel']
  DO k=1,4
    SELECT CASE(k)
    CASE(1)
      dt=0.02_wp;alpha=0.25_wp
    CASE(2)
      dt=0.2_wp;alpha=0.5_wp
    CASE(3)
      dt=0.4_wp;alpha=0.75_wp
    CASE(4)
      dt=0.8_wp;alpha=1.0_wp
    END SELECT
    old=[1000.0_wp,3000.0_wp,-4000.0_wp,1.2E9_wp]
    q=old;known=0.0_wp
    CALL solve_rk_step(partition,group,q,old,dt,alpha,known, &
         0.0_wp,0.0_wp,0.0_wp,iterations,converged,info,search_failed)
    IF (.NOT.converged .OR. info/=0 .OR. search_failed) &
         ERROR STOP 'stage diagonal Newton solve failed'
    expected=old(2:3)/(1.0_wp+0.1_wp*dt*alpha)
    CALL check(q(2:3),expected,'DIRK diagonal exact stage update')
  END DO
  PRINT *, 'PASS D17 four stage-diagonal Newton coefficients'
  ! Verify the production stage residual with nonzero known forcing
  ! for n_RK=2,3,4. This exercises stage composition, not the spatial tableau.
  DO k=2,4
    n_RK=k
    old=[1000.0_wp,3000.0_wp,-4000.0_wp,1.2E9_wp]
    dt=0.2_wp
    known=0.0_wp
    known(2:3)=0.02_wp*old(2:3)
    q=old
    CALL solve_rk_step(partition,group,q,old,dt,0.75_wp,known, &
         0.0_wp,0.0_wp,0.0_wp,iterations,converged,info,search_failed)
    IF (.NOT.converged .OR. info/=0 .OR. search_failed) &
         ERROR STOP 'known forcing RK stage failed'
    expected=(old(2:3)-dt*known(2:3))/(1.0_wp+0.1_wp*dt*0.75_wp)
    CALL check(q(2:3),expected,'RK stage with known forcing')
  END DO
  PRINT *, 'PASS D17 known forcing for n_RK=2,3,4'
CONTAINS
  SUBROUTINE check(actual,reference,label)
    REAL(wp),INTENT(IN)::actual(2),reference(2)
    CHARACTER(*),INTENT(IN)::label
    IF(ANY(ABS(actual-reference)>1.0E-7_wp*MAX(1.0_wp,MAXVAL(ABS(reference))))) THEN
       PRINT *, 'FAIL ',label,actual,reference
       ERROR STOP 1
    END IF
  END SUBROUTINE check
END PROGRAM test_d17_production_newton
