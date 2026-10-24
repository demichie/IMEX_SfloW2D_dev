! D17.20: full production nonlinear_solver_2d stage sweep, 2x2 and 3x3.
PROGRAM test_d17_production_newton
  USE parameters_2d
  USE constitutive_parameters_2d
  USE equation_metadata_2d, ONLY : equation_partition_type, IM_GROUP_METHOD_NEWTON
  USE source_registry_2d, ONLY : source_group_type
  USE source_catalog_2d, ONLY : SOURCE_METHOD_NEWTON,SOURCE_CADENCE_STAGE, &
       SOURCE_PHASE_RK_IMPLICIT
  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver,solve_rk_step
  IMPLICIT NONE
  TYPE(equation_partition_type) :: partition
  TYPE(source_group_type) :: group(1)
  REAL(wp) :: q(4),old(4),known(4),dt,expected(2)
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
