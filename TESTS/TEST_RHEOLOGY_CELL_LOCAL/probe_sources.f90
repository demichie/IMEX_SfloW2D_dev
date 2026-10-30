!> \brief Isolated production-source probe for D02 numerical characterization.
!> \details One process per fixture isolates floating-point traps and shared writes.
PROGRAM probe_sources
  USE omp_lib, ONLY : omp_get_num_threads
  USE parameters_2d
  USE constitutive_parameters_2d
  USE rheology_parameters_2d, ONLY : mu_s, mu_2, I_0, muI_inf,               &
       mu_0, mu_inf, Fr_0, U_w, xi, T_ref
  USE rheology_voellmy_2d, ONLY : mu
  USE rheology_plastic_2d, ONLY : tau
  USE rheology_lava_2d, ONLY : nu_ref, visc_par, tau0
  USE rheology_lahar_2d, ONLY : alpha1_coeff, beta1, alpha2, beta2, Kappa, n_td2
  USE rheology_quadratic_2d, ONLY : friction_factor
  USE state_conversion_2d, ONLY : qc_to_qp
  USE geometry_2d, ONLY : pi_g
  USE source_terms_2d, ONLY : eval_source_implicit, eval_source_semi_implicit
  USE mass_exchange_2d, ONLY : eval_mass_exchange_explicit
  USE equation_metadata_2d, ONLY : equation_partition_type, IM_GROUP_METHOD_NEWTON
  USE source_registry_2d, ONLY : source_registry_type
  USE source_configuration_2d, ONLY : initialize_configured_sources
  USE source_input_2d, ONLY : rheology_method, rheology_cadence, rheology_group, &
       thermal_method, thermal_cadence, thermal_group, &
       pore_pressure_method, pore_pressure_cadence, pore_pressure_group
  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver, solve_rk_step
  IMPLICIT NONE
  INTEGER :: operation, solids, column, j, level, iterations, linear_info, iz
  LOGICAL :: derivatives, converged, search_failed, one_sided_pressure
  TYPE(equation_partition_type) :: partition
  TYPE(source_registry_type) :: source_registry
  REAL(wp) :: input_q(7), bx, by, bxx, bxy, byy, z, dt, diagonal, pressure
  REAL(wp) :: q(7), qp(9), real_source(7), si(7), plus(7), minus(7), scale(7)
  REAL(wp) :: cs(7,7), normalized_cs(7,7), fd(7,7,4), steps(4), delta, env_before
  REAL(wp) :: stage(7), residual(7), erosion(2), deposition(2), exchange(7)
  REAL(wp) :: carrier_erosion, carrier_loss, topo, substrate(2)
  TYPE(mass_exchange_rates_type) :: total_rates
  REAL(wp) :: test_z(4), expected_parallel_source(7,4), parallel_source(7)
  REAL(wp) :: parallel_error, local_error
  INTEGER :: iteration, actual_threads
  COMPLEX(wp) :: cq(7), complex_source(7), base_complex(7)
  LOGICAL :: mask(7)
  CHARACTER(LEN=256) :: registry_error
  NAMELIST /fixture/ operation, solids, input_q, bx, by, bxx, bxy, byy, z, &
       dt, diagonal, derivatives, one_sided_pressure, rheology_flag, rheology_model, gas_flag, &
       pore_pressure_flag, dynamic_permeability_flag, sutherland_flag,    &
       stochastic_flag, slope_correction_flag, curvature_term_flag,      &
       f_inhibit_mode, mu, xi, tau, nu_ref, visc_par, T_ref, tau0,        &
       alpha1_coeff, beta1, alpha2, beta2, Kappa, n_td2, friction_factor, &
       mu_0, mu_inf, Fr_0, U_w, mu_s, mu_2, muI_inf, I_0, T_env,        &
       radiative_term_coeff, convective_term_coeff, hydraulic_permeability, &
       maximum_solid_packing, alpha_trans, N_inh, gas_loss_flag

  ! Explicit fixture defaults; no production initializer or parser supplies
  ! derived values to the independent Python oracle.
  operation=0; solids=2; derivatives=.TRUE.; one_sided_pressure=.FALSE.; input_q=0.0_wp
  bx=0.0_wp; by=0.0_wp; bxx=0.0_wp; bxy=0.0_wp; byy=0.0_wp; z=0.0_wp
  dt=0.05_wp; diagonal=1.0_wp
  rheology_flag=.TRUE.; rheology_model='quadratic'; gas_flag=.FALSE.
  pore_pressure_flag=.FALSE.; dynamic_permeability_flag=.FALSE.
  sutherland_flag=.FALSE.; stochastic_flag=.FALSE.
  slope_correction_flag=.FALSE.; curvature_term_flag=.FALSE.
  mu=0.3_wp; xi=500.0_wp; tau=0.4_wp; nu_ref=0.02_wp
  visc_par=0.003_wp; T_ref=900.0_wp; tau0=20.0_wp
  alpha1_coeff=1.0_wp; beta1=2.0_wp; alpha2=5.0_wp; beta2=3.0_wp
  Kappa=24.0_wp; n_td2=0.01_wp; friction_factor=0.03_wp
  mu_0=0.2_wp; mu_inf=0.6_wp; Fr_0=1.3_wp; U_w=2.0_wp
  mu_s=0.4_wp; mu_2=0.7_wp; muI_inf=0.005_wp; I_0=0.3_wp
  T_env=300.0_wp; radiative_term_coeff=3.0E-8_wp; convective_term_coeff=4.0_wp
  hydraulic_permeability=1.0E-12_wp; maximum_solid_packing=0.64_wp
  alpha_trans=0.2_wp; N_inh=2; f_inhibit_mode='OFF'; gas_loss_flag=.FALSE.
  READ(*,NML=fixture)
  ! The independent oracle applies radiative/convective cooling only for lava.
  ! D17 makes the cooling switch independent of the rheology model.
  thermal_cooling_flag=(TRIM(rheology_model)=='lava')
  IF (solids<1 .OR. solids>2) ERROR STOP 'fixture solids must be one or two'
  n_solid=solids; n_add_gas=0; n_stoch_vars=0
  n_pore_vars=MERGE(1,0,pore_pressure_flag)
  n_vars=4+solids+n_pore_vars; n_eqns=n_vars; n_RK=2
  kt_trace_reconstruction='froude_blend'
  IF (TRIM(rheology_model)=='lava') kt_trace_reconstruction='h_x'
  idx_h=1; idx_hu=2; idx_hv=3; idx_T=4
  idx_totMassEqn=1; idx_uEqn=2; idx_vEqn=3; idx_engyEqn=4
  idx_solid_first=5; idx_solid_last=4+solids
  idx_solidEqn_first=5; idx_solidEqn_last=4+solids
  idx_add_gas_first=5+solids; idx_add_gas_last=4+solids
  idx_addGasEqn_first=idx_add_gas_first; idx_addGasEqn_last=idx_add_gas_last
  idx_pore=n_vars; idx_poreEqn=n_vars; idx_stoch=0; idx_stochEqn=0
  idx_u=n_vars+1; idx_v=n_vars+2
  lava_flag=(TRIM(rheology_model)=='lava')
  liquid_flag=.NOT.gas_flag .AND. .NOT.lava_flag
  IF (lava_flag) gas_flag=.FALSE.
  stoch_transport_flag=.FALSE.
  ALLOCATE(rho_s(solids),inv_rho_s(solids),sp_heat_s(solids),diam_s(solids),sphericity_s(solids))
  ALLOCATE(sp_heat_g(0),sp_gas_const_g(0))
  rho_s=2500.0_wp; sp_heat_s=900.0_wp; diam_s=1.0E-4_wp; sphericity_s=0.8_wp
  IF (solids==2) THEN
     rho_s(2)=1800.0_wp; sp_heat_s(2)=1100.0_wp
     diam_s(2)=4.0E-4_wp; sphericity_s(2)=0.6_wp
  END IF
  inv_rho_s=1.0_wp/rho_s
  rho_l=1000.0_wp; inv_rho_l=1.0_wp/rho_l; sp_heat_l=4180.0_wp
  grav=9.81_wp; inv_grav=1.0_wp/grav
  pi_g=ACOS(-1.0_wp)
  pres=101300.0_wp; inv_pres=1.0_wp/pres
  sp_heat_a=1005.0_wp; sp_gas_const_a=287.05_wp; T_ambient=300.0_wp
  rho_a_amb=pres/(sp_gas_const_a*T_ambient)
  rho_c_sub=rho_a_amb
  kin_visc_c=1.5E-5_wp; kin_visc_a=kin_visc_c
  Tref_Suth=273.15_wp; muRef_Suth=1.716E-5_wp; S_mu=110.4_wp
  eps_sing=1.0E-8_wp; eps_sing4=eps_sing**4
  verbose_level=-1
  settling_flag=.FALSE.; entrainment_flag=.FALSE.; erosion_coeff=0.0_wp
  loss_rate=0.0_wp; alphastot_min=0.0_wp; erodible_porosity=0.4_wp
  coeff_porosity=erodible_porosity/(1.0_wp-erodible_porosity)
  erodible_deposit_flag=.FALSE.; liquid_vaporization_flag=.FALSE.
  T_erodible=300.0_wp
  CALL precompute_pascal_coefficient(N_inh)
  q=input_q; env_before=T_env
  IF (operation==1) THEN
     ! Hazard probes deliberately avoid unrelated SI or derivative calls.
     CALL implicit(q,real_source)
     WRITE(*,'(A,7ES26.17E3)') 'REAL ',real_source(1:n_vars)
     STOP
  END IF
  IF (operation==7) THEN
     cq=CMPLX(q,0.0_wp,wp)
     CALL eval_source_implicit(bx,by,z,c_qj=cq(1:n_vars),c_nh_term_impl=complex_source(1:n_vars))
     WRITE(*,'(A,7ES26.17E3)') 'COMPLEX_REAL ',REAL(complex_source(1:n_vars),wp)
     STOP
  END IF
  CALL qc_to_qp(q(1:n_vars),qp(1:n_vars+2),pressure)
  IF (operation==6) THEN
     CALL eval_source_semi_implicit(bx,by,bxx,bxy,byy,1.0_wp,qp(1:n_vars+2),si(1:n_vars),z)
     WRITE(*,'(A,7ES26.17E3)') 'SI ',si(1:n_vars)
     STOP
  END IF
  IF (operation==3) THEN
     substrate=0.0_wp
     CALL eval_mass_exchange_explicit(qp(1:n_vars+2),0,bx,by, &
          substrate(1:solids),dt,total_rates,exchange(1:n_vars))
     carrier_loss=total_rates%carrier_withdrawal_volume_rate
     WRITE(*,'(A,ES26.17E3)') 'LOSS ',carrier_loss
     WRITE(*,'(A,7ES26.17E3)') 'EXCHANGE ',exchange(1:n_vars)
     STOP
  END IF
  IF (operation==9) THEN
     IF ((TRIM(rheology_model)/='voellmy' .AND. TRIM(rheology_model)/='mui_voellmy') .OR. .NOT.stochastic_flag) &
          ERROR STOP 'parallel workspace probe requires stochastic Voellmy rheology'
     test_z=[-1000.0_wp,0.0_wp,500.0_wp,100000.0_wp]
     expected_parallel_source=0.0_wp
     DO j=1,SIZE(test_z)
        CALL eval_source_implicit(bx,by,test_z(j),r_qj=q(1:n_vars), &
             r_nh_term_impl=expected_parallel_source(1:n_vars,j))
     END DO
     parallel_error=0.0_wp; actual_threads=0
     !$OMP PARALLEL DO DEFAULT(NONE) SHARED(test_z,expected_parallel_source,q,bx,by,n_vars) &
     !$OMP& PRIVATE(iz,parallel_source,local_error) REDUCTION(MAX:parallel_error,actual_threads) &
     !$OMP& NUM_THREADS(4) SCHEDULE(STATIC)
     DO iteration=1,100000
        iz=MOD(iteration-1,SIZE(test_z))+1
        parallel_source=0.0_wp
        CALL eval_source_implicit(bx,by,test_z(iz),r_qj=q(1:n_vars), &
             r_nh_term_impl=parallel_source(1:n_vars))
        local_error=MAXVAL(ABS(parallel_source(1:n_vars)-expected_parallel_source(1:n_vars,iz)))
        parallel_error=MAX(parallel_error,local_error)
        actual_threads=MAX(actual_threads,omp_get_num_threads())
     END DO
     !$OMP END PARALLEL DO
     WRITE(*,'(A,ES26.17E3)') 'PARALLEL_XI ',parallel_error
     WRITE(*,'(A,I0)') 'PARALLEL_THREADS ',actual_threads
  END IF
  real_source=0.0_wp; si=0.0_wp; cs=0.0_wp; normalized_cs=0.0_wp; fd=0.0_wp
  CALL implicit(q,real_source)
  cq=CMPLX(q,0.0_wp,wp); base_complex=0.0_wp
  CALL eval_source_implicit(bx,by,z,c_qj=cq(1:n_vars),c_nh_term_impl=base_complex(1:n_vars))
  CALL eval_source_semi_implicit(bx,by,bxx,bxy,byy,                  &
       1.0_wp/(1.0_wp+bx**2+by**2),qp(1:n_vars+2),si(1:n_vars),z)
  scale=MAX(ABS(q),1.0_wp)
  steps=[1.0E-5_wp,1.0E-6_wp,1.0E-7_wp,1.0E-8_wp]
  IF (derivatives) THEN
     DO column=1,n_vars
        cq=CMPLX(q,0.0_wp,wp); delta=1.0E-25_wp*scale(column)
        cq(column)=cq(column)+CMPLX(0.0_wp,delta,wp)
        CALL eval_source_implicit(bx,by,z,c_qj=cq(1:n_vars),c_nh_term_impl=complex_source(1:n_vars))
        cs(1:n_vars,column)=AIMAG(complex_source(1:n_vars))/delta
        cq=CMPLX(q,0.0_wp,wp)
        delta=n_vars*EPSILON(1.0_wp)
        cq(column)=cq(column)+CMPLX(0.0_wp,delta*MAX(ABS(q(column)),1.0E-3_wp),wp)
        CALL eval_source_implicit(bx,by,z,c_qj=cq(1:n_vars),c_nh_term_impl=complex_source(1:n_vars))
        normalized_cs(1:n_vars,column)=AIMAG(complex_source(1:n_vars))/delta
        DO level=1,4
           delta=steps(level)*scale(column)
           plus=q; minus=q; plus(column)=plus(column)+delta; minus(column)=minus(column)-delta
           CALL implicit(plus,stage); CALL implicit(minus,residual)
           IF (one_sided_pressure .AND. column==n_vars) THEN
              fd(1:n_vars,column,level)=(stage(1:n_vars)-real_source(1:n_vars))/delta
           ELSE
              fd(1:n_vars,column,level)=(stage(1:n_vars)-residual(1:n_vars))/(2.0_wp*delta)
           END IF
        END DO
     END DO
  END IF
  WRITE(*,'(A,I0)') 'N ',n_vars
  WRITE(*,'(A,7ES26.17E3)') 'REAL ',real_source(1:n_vars)
  WRITE(*,'(A,7ES26.17E3)') 'COMPLEX_REAL ',REAL(base_complex(1:n_vars),wp)
  WRITE(*,'(A,7ES26.17E3)') 'COMPLEX_IMAG ',AIMAG(base_complex(1:n_vars))
  WRITE(*,'(A,7ES26.17E3)') 'SI ',si(1:n_vars)
  DO j=1,n_vars
     WRITE(*,'(A,I0,1X,7ES26.17E3)') 'CS',j,cs(1:n_vars,j)
     WRITE(*,'(A,I0,1X,7ES26.17E3)') 'NCS',j,normalized_cs(1:n_vars,j)
     DO level=1,4
        WRITE(*,'(A,I0,A,I0,1X,7ES26.17E3)') 'FD',j,'_',level,fd(1:n_vars,j,level)
     END DO
  END DO
  WRITE(*,'(A,2ES26.17E3)') 'GLOBALS ',env_before,T_env
  IF (operation==2) THEN
     mask=.FALSE.; mask(2:3)=.TRUE.
     IF (TRIM(rheology_model)=='lava') mask(4)=.TRUE.
     IF (pore_pressure_flag) mask(n_vars)=.TRUE.
     CALL partition%initialize(mask(1:n_vars))
     CALL initialize_nonlinear_solver
     rheology_method='NEWTON'; rheology_cadence='STAGE'; rheology_group=1
     thermal_method='NEWTON'; thermal_cadence='STAGE'; thermal_group=1
     pore_pressure_method='NEWTON'; pore_pressure_cadence='STAGE'; pore_pressure_group=1
     CALL initialize_configured_sources(source_registry,partition,registry_error)
     IF (LEN_TRIM(registry_error)>0) ERROR STOP TRIM(registry_error)
     stage=q
     CALL solve_rk_step(partition,source_registry%stage_groups(), &
          stage(1:n_vars),q(1:n_vars),dt,diagonal, &
          0.0_wp*q(1:n_vars),bx,by,z,iterations,converged,linear_info,search_failed)
     CALL implicit(stage,residual)
     residual(1:n_vars)=stage(1:n_vars)-q(1:n_vars)-dt*diagonal*residual(1:n_vars)
     WRITE(*,'(A,7ES26.17E3)') 'STAGE ',stage(1:n_vars)
     WRITE(*,'(A,7ES26.17E3)') 'RESIDUAL ',residual(1:n_vars)
     WRITE(*,'(A,4I8)') 'STATUS ',MERGE(1,0,converged),iterations,linear_info,MERGE(1,0,search_failed)
     CALL source_registry%finalize
  END IF
CONTAINS
  !> \brief Evaluate the REAL production source without changing caller storage.
  !> \param[in] state Conservative fixture state, padded to seven entries.
  !> \param[out] source Production implicit source, padded with zeros.
  SUBROUTINE implicit(state,source)
    REAL(wp), INTENT(IN) :: state(7)
    REAL(wp), INTENT(OUT) :: source(7)
    source=0.0_wp
    CALL eval_source_implicit(bx,by,z,r_qj=state(1:n_vars),r_nh_term_impl=source(1:n_vars))
  END SUBROUTINE implicit
END PROGRAM probe_sources
