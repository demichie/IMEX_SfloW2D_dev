! Lie backward-Euler regression: full conservative source-rate footprint.
PROGRAM test_d17_production_newton
  USE thermal_parameters_2d, ONLY : thermal_l1
  USE parameters_2d
  USE constitutive_parameters_2d
  USE rheology_parameters_2d, ONLY : rheology_l1
  USE source_registry_2d, ONLY : source_group_type
  USE source_catalog_2d, ONLY : SOURCE_METHOD_NEWTON
  USE nonlinear_solver_2d, ONLY : solve_backward_euler_group
  USE source_terms_2d, ONLY : eval_source_implicit
  IMPLICIT NONE
  TYPE(source_group_type) :: group
  REAL(wp) :: q(4),old(4),dt,rate(4),residual(4)
  COMPLEX(wp) :: rate_complex(4)
  INTEGER :: iterations,info
  LOGICAL :: converged,search_failed
  n_vars=4;n_eqns=4;n_solid=0;n_add_gas=0;n_stoch_vars=0;n_pore_vars=0
  idx_totMassEqn=1;idx_uEqn=2;idx_vEqn=3;idx_engyEqn=4
  idx_solid_first=5;idx_solid_last=4
  idx_solidEqn_first=5;idx_solidEqn_last=4
  idx_add_gas_first=5;idx_add_gas_last=4
  idx_pore=0;idx_poreEqn=0;idx_stoch=0;idx_u=2;idx_v=3
  rheology_flag=.TRUE.;rheology_l1%rheology_model='newmodel'
  pore_pressure_flag=.FALSE.;thermal_l1%thermal_cooling_flag=.FALSE.
  slope_correction_flag=.FALSE.;curvature_term_flag=.FALSE.
  stochastic_flag=.FALSE.;stoch_transport_flag=.FALSE.
  liquid_flag=.TRUE.;gas_flag=.FALSE.;lava_flag=.FALSE.
  rho_l=1000.0_wp;inv_rho_l=1.0_wp/rho_l;sp_heat_l=4000.0_wp
  sp_heat_a=1005.0_wp;sp_gas_const_a=287.0_wp
  grav=9.81_wp;rho_a_amb=1.2_wp
  eps_sing=1.0E-8_wp;eps_sing4=eps_sing**4
  ALLOCATE(rho_s(0),inv_rho_s(0),sp_heat_s(0),diam_s(0),sphericity_s(0))
  ALLOCATE(sp_heat_g(0),sp_gas_const_g(0))
  group%id=1;group%method=SOURCE_METHOD_NEWTON
  group%equation_ids=[2,3]
  group%source_names=[CHARACTER(LEN=64)::'rheology_newmodel']
  old=[1000.0_wp,3000.0_wp,-4000.0_wp,1.2E9_wp]
  q=old;dt=0.2_wp
  CALL solve_backward_euler_group(group,q,dt,0.0_wp,0.0_wp,0.0_wp,0.0_wp, &
       iterations,converged,info,search_failed,thermal=thermal_l1)
  IF (.NOT.converged .OR. info/=0 .OR. search_failed) &
       ERROR STOP 'Lie Newton nonlinear solve failed'
  IF (iterations<1) ERROR STOP 'Newton iterations not recorded'
  IF (q(1)/=old(1) .OR. q(4)/=old(4)) &
       ERROR STOP 'Inactive conservative equations changed'
  CALL eval_source_implicit(0.0_wp,0.0_wp,0.0_wp, &
       c_qj=CMPLX(q,0.0_wp,wp),c_nh_term_impl=rate_complex,group=group,thermal=thermal_l1)
  rate=REAL(rate_complex,wp)
  residual=q-old-dt*rate
  IF (MAXVAL(ABS(residual(group%equation_ids)))>1.0E-6_wp) &
       ERROR STOP 'Lie Newton backward Euler residual'
  PRINT *, 'PASS Lie Newton conservative footprint and backward Euler residual'
END PROGRAM test_d17_production_newton
