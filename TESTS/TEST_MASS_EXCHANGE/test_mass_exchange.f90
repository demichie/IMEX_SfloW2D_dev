!> \brief N8 production exchange budgets, limited nodal updates and true checkpoint round trips.
PROGRAM test_mass_exchange

  USE parameters_2d
  USE constitutive_parameters_2d
  USE geometry_2d
  USE domain_2d, ONLY : domain_type
  USE state_2d, ONLY : state_type
  USE model_layout_2d, ONLY : model_layout_type
  USE runtime_2d, ONLY : runtime_state_type
  USE stochastic_module, ONLY : stochastic_workspace_type, sym_noise, noise_pow_val
  USE state_conversion_2d, ONLY : qc_to_qp, mixt_var
  USE mass_exchange_deposition_2d, ONLY : settling_velocity
  USE mass_exchange_entrainment_2d, ONLY : eval_air_entrainment_rate
  USE mass_exchange_laws_2d, ONLY : eval_mass_exchange_terms
  USE step_source_executor_2d, ONLY : apply_mass_exchange_step => execute_step_source_groups
  USE equation_metadata_2d, ONLY : equation_partition_type
  USE source_registry_2d, ONLY : source_registry_type, source_group_type
  USE source_terms_2d, ONLY : eval_cell_source_explicit
  USE source_catalog_2d, ONLY : SOURCE_CADENCE_STEP, &
       SOURCE_CADENCE_STAGE, SOURCE_METHOD_EXPLICIT
  USE source_configuration_2d, ONLY : initialize_configured_sources
  USE source_input_2d, ONLY : thermal_method,thermal_cadence,thermal_group, &
       rheology_method,rheology_cadence,rheology_group
  USE rheology_lava_2d, ONLY : nu_ref,visc_par,tau0
  USE rheology_parameters_2d, ONLY : T_ref
  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver, &
       finalize_nonlinear_solver,solve_imex_stage
  USE thermal_laws_2d, ONLY : solve_thermal_cooling_analytic
  USE inpout_2d, ONLY : write_restart_file, read_restart_file, output_idx
  USE omp_lib
  USE, INTRINSIC :: iso_fortran_env, ONLY : int64
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite

  IMPLICIT NONE

  INTEGER, PARAMETER :: nx=8, ny=7, ns=2
  INTEGER :: nq=6
  REAL(wp), PARAMETER :: tol=2048.0_wp*EPSILON(1.0_wp)
  TYPE(state_type) :: state
  TYPE(domain_type) :: domain
  TYPE(model_layout_type) :: layout
  TYPE(runtime_state_type) :: runtime
  TYPE(stochastic_workspace_type) :: stochastic
  TYPE(equation_partition_type) :: equation_partition
  TYPE(source_registry_type) :: source_registry
  CHARACTER(LEN=64) :: selected, argument
  INTEGER :: actual_threads, requested_threads, status, snapshot_unit, extension_unit

  CALL get_command_argument(1, argument)
  READ(argument,*) requested_threads
  actual_threads = 0
  !$OMP PARALLEL
  !$OMP SINGLE
  actual_threads = omp_get_num_threads()
  !$OMP END SINGLE
  !$OMP END PARALLEL
  IF (actual_threads /= requested_threads) ERROR STOP 'wrong actual OpenMP team'
  WRITE(*,*) 'Actual OpenMP team:', actual_threads
  CALL get_environment_variable('N8_CASE', selected, STATUS=status)
  IF (status /= 0) selected = 'all'
  CALL initialize_fixture
  CALL test_lava_deposition_transaction
  CALL test_integrated_lava_deposition_cooling_restart
  CALL test_liquid_loss_without_solids
  CALL test_empty_step_registry
  CALL test_analytic_cooling_step_executor
  CALL test_explicit_multi_source_groups
  CALL test_combined_step_sources
  CALL test_official_cooling_restart
  OPEN(NEWUNIT=snapshot_unit, FILE='snapshot.bin', ACCESS='STREAM',      &
       FORM='UNFORMATTED', STATUS='REPLACE')
  OPEN(NEWUNIT=extension_unit, FILE='flat_evaluations.bin', ACCESS='STREAM', &
       FORM='UNFORMATTED', STATUS='REPLACE')

  SELECT CASE(TRIM(selected))
  CASE('all')
     CALL run_case('deposition')
     CALL run_case('erosion')
     CALL run_case('combined')
     CALL run_case('masked')
     CALL run_case('masked_fissural')
     CALL run_case('cutoff')
     CALL run_case('loss_only')
     CALL run_case('loss_reserve')
     CALL run_case('thermal')
     CALL run_restart
  CASE('deposition','erosion','combined','masked','masked_fissural','cutoff','loss_only','loss_reserve','thermal')
     CALL run_case(TRIM(selected))
  CASE('restart')
     CALL run_restart
  CASE DEFAULT
     ERROR STOP 'unknown N8_CASE'
  END SELECT

  CLOSE(snapshot_unit)
  IF (TRIM(selected)=='all') CALL run_gas_cases
  CLOSE(extension_unit)
  CALL release_topography_workspace
  CALL release_topography_workspace
  CALL stochastic%finalize
  CALL domain%finalize
  CALL state%finalize
  CALL layout%finalize
  WRITE(*,*) 'PASS: production inventory limits, bed volume, geometry and restart'

CONTAINS

  !> \brief Allocate a rectangular two-solid liquid fixture and the real restart dependencies.
  SUBROUTINE initialize_fixture
    INTEGER :: j,k,l
    comp_cells_x=nx; comp_cells_y=ny
    comp_interfaces_x=nx+1; comp_interfaces_y=ny+1
    dx=0.75_wp; dy=1.25_wp; cell_size=dx
    n_solid=ns; n_add_gas=0; n_stoch_vars=0; n_pore_vars=0
    kt_trace_reconstruction = 'froude_blend'
    n_vars=nq; n_eqns=nq; n_RK=2
    n_thickness_levels=1; n_dyn_pres_levels=1
    idx_h=1; idx_hu=2; idx_hv=3; idx_T=4
    idx_solid_first=5; idx_solid_last=6
    idx_solidEqn_first=5; idx_solidEqn_last=6
    idx_add_gas_first=7; idx_add_gas_last=6
    idx_stoch=0; idx_pore=0; idx_u=7; idx_v=8
    gas_flag=.FALSE.; liquid_flag=.TRUE.; rheology_flag=.FALSE.
    rheology_model='none'; stochastic_flag=.FALSE.; stoch_transport_flag=.FALSE.
    pore_pressure_flag=.FALSE.; gas_loss_flag=.FALSE.
    liquid_vaporization_flag=.FALSE.; sutherland_flag=.FALSE.
    entrainment_flag=.FALSE.; slope_correction_flag=.TRUE.
    topo_change_flag=.TRUE.; radial_source_flag=.FALSE.
    bottom_fissural_source_flag=.FALSE.; bottom_radial_source_flag=.FALSE.
    grav=9.81_wp; T_ambient=300.0_wp; rho_a_amb=1.2_wp
    rho_l=1000.0_wp; inv_rho_l=1.0_wp/rho_l; sp_heat_l=4180.0_wp
    sp_heat_a=998.0_wp; sp_gas_const_a=287.05_wp
    pres=101300.0_wp; inv_pres=1.0_wp/pres
    rho_c_sub=rho_l; kin_visc_c=1.0E-6_wp; hydraulic_permeability=0.0_wp
    eps_sing=1.0E-8_wp; eps_sing4=eps_sing**4
    maximum_solid_packing=0.6_wp
    ALLOCATE(rho_s(ns), inv_rho_s(ns), sp_heat_s(ns), diam_s(ns), erodible_fract(ns))
    ALLOCATE(sp_heat_g(0), sp_gas_const_g(0), loss_rate)
    rho_s=[2000.0_wp,3000.0_wp]; inv_rho_s=1.0_wp/rho_s
    sp_heat_s=[900.0_wp,1100.0_wp]; diam_s=[1.0E-4_wp,2.0E-4_wp]
    erodible_fract=[0.6_wp,0.4_wp]
    ALLOCATE(B_vertex(nx+1,ny+1), B_cent(nx,ny), B_face_x(nx+1,ny), B_face_y(nx,ny+1))
    ALLOCATE(B_prime_x_geom(nx,ny), B_prime_y_geom(nx,ny))
    ALLOCATE(B_second_xx_geom(nx,ny), B_second_xy_geom(nx,ny), B_second_yy_geom(nx,ny))
    ALLOCATE(grav_coeff(nx,ny), grav_coeff_stag_x(nx+1,ny), grav_coeff_stag_y(nx,ny+1))
    ALLOCATE(deposit(nx,ny,ns), erosion(nx,ny,ns), erodible(ns,nx,ny),lava_deposit(nx,ny))
    ALLOCATE(B_zone(nx,ny), cell_source_fractions(nx,ny))
    ! The shared explicit source adapter also accesses injection geometry.
    ! This mass-exchange fixture has no fissures or lateral injection arcs.
    ALLOCATE(cell_fissure_fractions(nx,ny,0))
    ALLOCATE(cell_arc_perim(nx,ny),cell_arc_n_x(nx,ny),cell_arc_n_y(nx,ny))
    cell_fissure_fractions=0.0_wp
    cell_arc_perim=0.0_wp
    cell_arc_n_x=0.0_wp
    cell_arc_n_y=0.0_wp
    CALL layout%initialize(1,nq,nq,ns,0,0,0,.FALSE.)
    CALL state%initialize(layout)
    CALL domain%initialize
    ! All cells except one wet inactive corner are on the production workset.
    ! Two dry cells remain on the list as a realistic solve/halo boundary.
    l=0
    DO k=1,ny
       DO j=1,nx
          IF (j==nx .AND. k==ny) CYCLE
          l=l+1; domain%j_cent(l)=j; domain%k_cent(l)=k
       END DO
    END DO
    domain%solve_cells=l
    sym_noise=0.0_wp; noise_pow_val=1.0_wp
    CALL stochastic%initialize
    output_idx=0
    runtime%t_output=0.0_wp; runtime%t_runout=0.0_wp; runtime%t_probes=0.0_wp
  END SUBROUTINE initialize_fixture

  !> \brief Build the production source schedule from the current fixture settings.
  SUBROUTINE configure_source_registry
    CHARACTER(LEN=256) :: error
    CALL initialize_configured_sources(source_registry,equation_partition,error)
    IF (LEN_TRIM(error)>0) ERROR STOP 'failed to configure mass-exchange source registry'
  END SUBROUTINE configure_source_registry

  !> \\brief Exercise prescribed liquid loss with zero transported solid fraction.
  !> \\details The limiter must preserve carrier availability, velocity and
  !> \brief Verify prescribed lava solidification as an immobile STEP deposit.
  !> \details Check mass-energy transfer, nodal geometry and official restart.
  SUBROUTINE test_lava_deposition_transaction
    USE mass_exchange_lava_deposition_2d, ONLY : lava_deposition_flag,lava_deposition_model, &
         lava_deposition_rate
    REAL(wp), PARAMETER :: delta=0.1_wp, rate=0.01_wp
    REAL(wp) :: mass_before,energy_before,bed_before(nx+1,ny+1)
    REAL(wp) :: checkpoint(nx,ny),post_restart(nx,ny)
    REAL(wp) :: expected_mass,expected_energy
    INTEGER :: j,k
    CALL reset_fixture
    rheology_model='lava'
    lava_deposition_flag=.TRUE.
    lava_deposition_model='CONSTANT'
    lava_deposition_rate=rate
    DO k=1,ny
       DO j=1,nx
          CALL set_cell(j,k,0.5_wp,[0.0_wp,0.0_wp],900.0_wp)
       END DO
    END DO
    state%q(:,5,4)=0.0_wp
    mass_before=state%q(1,2,2)
    energy_before=state%q(4,2,2)
    bed_before=B_vertex
    CALL configure_source_registry
    CALL apply_mass_exchange_step(state%q,state%qp,delta,domain,source_registry)
    expected_mass=mass_before-rho_l*rate*delta
    expected_energy=energy_before-rho_l*sp_heat_l*900.0_wp*rate*delta
    IF (ABS(state%q(1,2,2)-expected_mass)>1.0E-10_wp*mass_before) &
         ERROR STOP 'D18 lava deposit mass imbalance'
    IF (ABS(state%q(4,2,2)-expected_energy)>1.0E-10_wp*energy_before) &
         ERROR STOP 'D18 lava deposit thermal imbalance'
    IF (ABS(lava_deposit(2,2)-rate*delta)>1.0E-12_wp) &
         ERROR STOP 'D18 lava deposit inventory not updated'
    IF (ALL(B_vertex==bed_before)) ERROR STOP 'D18 lava deposit did not raise topography'
    IF (ANY(deposit/=0.0_wp) .OR. ANY(erosion/=0.0_wp)) &
         ERROR STOP 'D18 lava deposit altered granular inventories'
    checkpoint=lava_deposit
    runtime%t=delta;runtime%dt=delta
    CALL write_restart_file('lava_deposit_checkpoint.bin',runtime,stochastic,state,domain)
    lava_deposit=-999.0_wp
    CALL read_restart_file('lava_deposit_checkpoint.bin',runtime,stochastic,state,domain)
    post_restart=lava_deposit
    CALL same_bits('D18 lava deposit restart',RESHAPE(post_restart,[SIZE(post_restart)]), &
         RESHAPE(checkpoint,[SIZE(checkpoint)]))
    CALL source_registry%finalize
    CALL equation_partition%finalize
    lava_deposition_flag=.FALSE.;lava_deposition_rate=0.0_wp
    rheology_model='none'
    WRITE(*,*) 'PASS D18 immobile lava deposition STEP and production restart'
  END SUBROUTINE test_lava_deposition_transaction

  !> \brief Full D18 lava stage/deposition/cooling/restart integration gate.
  !> \details Production NEWTON stage, two ordered STEP groups, physical
  !>          viscosity feedback and bitwise continuation from a real checkpoint.
  SUBROUTINE test_integrated_lava_deposition_cooling_restart
    USE mass_exchange_lava_deposition_2d, ONLY : lava_deposition_flag,lava_deposition_model, &
         lava_deposition_rate,lava_deposition_visc_ref,lava_deposition_exponent
    USE source_input_2d, ONLY : mass_exchange_group
    USE source_catalog_2d, ONLY : SOURCE_METHOD_EXPLICIT, SOURCE_METHOD_ANALYTIC
    TYPE(source_group_type), ALLOCATABLE :: groups(:)
    REAL(wp) :: initial_mass, temp_after_one, initial_temp
    REAL(wp) :: before_stage(nq),after_stage(nq),im_rate(nq),residual(nq)
    REAL(wp) :: before_deposition,first_deposition,second_deposition
    REAL(wp) :: saved_q(nq,nx,ny),full_q(nq,nx,ny)
    REAL(wp) :: saved_b(nx+1,ny+1),full_b(nx+1,ny+1)
    REAL(wp) :: saved_lava(nx,ny),full_lava(nx,ny)
    REAL(wp) :: saved_deposit(nx,ny,ns),full_deposit(nx,ny,ns)
    REAL(wp) :: saved_erosion(nx,ny,ns),full_erosion(nx,ny,ns)
    REAL(wp) :: saved_erodible(ns,nx,ny),full_erodible(ns,nx,ny)
    REAL(wp) :: pdyn
    REAL(wp), PARAMETER :: delta=0.1_wp
    INTEGER :: l,j,k,iterations,linear_info,mass_group,cooling_group,i
    LOGICAL :: converged,line_search_failed
    CHARACTER(LEN=256) :: error

    CALL reset_fixture
    rheology_flag=.TRUE.;rheology_model='lava'
    rheology_method='NEWTON';rheology_cadence='STAGE';rheology_group=0
    thermal_cooling_flag=.TRUE.
    thermal_method='ANALYTIC';thermal_cadence='STEP';thermal_group=1
    mass_exchange_group=0
    lava_deposition_flag=.TRUE.
    lava_deposition_model='VISCOSITY'
    lava_deposition_rate=0.02_wp
    lava_deposition_visc_ref=0.03_wp
    lava_deposition_exponent=1.0_wp
    nu_ref=0.02_wp;visc_par=0.003_wp;tau0=0.0_wp;T_ref=900.0_wp
    radiative_term_coeff=5.103E-8_wp
    convective_term_coeff=9.0_wp
    T_env=300.0_wp
    DO k=1,ny
       DO j=1,nx
          CALL set_cell(j,k,0.5_wp,[0.0_wp,0.0_wp],900.0_wp)
       END DO
    END DO
    state%q(:,5,4)=0.0_wp
    CALL configure_source_registry
    IF (equation_partition%n_IM_groups/=1) &
         ERROR STOP 'D18 integrated lava stage partition'
    groups=source_registry%step_groups()
    IF (SIZE(groups)/=2) ERROR STOP 'D18 integrated lava STEP group count'
    mass_group=0;cooling_group=0
    DO i=1,SIZE(groups)
       IF (ANY(groups(i)%source_names=='mass_exchange_lava_deposition')) mass_group=i
       IF (ANY(groups(i)%source_names=='thermal_cooling')) cooling_group=i
    END DO
    IF (mass_group/=1 .OR. cooling_group/=2) ERROR STOP 'D18 lava STEP group ordering'
    IF (groups(mass_group)%method/=SOURCE_METHOD_EXPLICIT .OR. &
         groups(cooling_group)%method/=SOURCE_METHOD_ANALYTIC) &
         ERROR STOP 'D18 integrated lava STEP methods'

    CALL initialize_nonlinear_solver
    before_stage=state%q(:,2,2)
    initial_mass=before_stage(1)
    residual=0.0_wp
    CALL solve_imex_stage(equation_partition,source_registry%stage_groups(), &
         before_stage,before_stage,delta,0.75_wp,residual, &
         B_prime_x_geom(2,2),B_prime_y_geom(2,2), &
         B_second_xx_geom(2,2),B_second_xy_geom(2,2), &
         B_second_yy_geom(2,2),grav_coeff(2,2),0.0_wp, &
         .FALSE.,.FALSE.,after_stage,im_rate,iterations,converged, &
         linear_info,line_search_failed)
    IF (.NOT.converged .OR. linear_info/=0 .OR. line_search_failed) &
         ERROR STOP 'D18 integrated lava Newton stage failed'
    IF (MAXVAL(ABS(after_stage(2:3)-before_stage(2:3)))<=0.0_wp) &
         ERROR STOP 'D18 integrated lava stage inactive'
    IF (after_stage(4)/=before_stage(4)) ERROR STOP 'D18 stage performed STEP cooling'
    state%q(:,2,2)=after_stage
    CALL finalize_nonlinear_solver
    DO l=1,domain%solve_cells
       j=domain%j_cent(l);k=domain%k_cent(l)
       IF (state%q(1,j,k)<=0.0_wp) CYCLE
       CALL qc_to_qp(state%q(:,j,k),state%qp(:,j,k),pdyn)
    END DO
    initial_temp=state%qp(4,2,2)
    CALL apply_mass_exchange_step(state%q,state%qp,delta,domain,source_registry)
    first_deposition=lava_deposit(2,2)
    IF (first_deposition<=0.0_wp) ERROR STOP 'D18 first lava STEP deposition inactive'
    IF (ABS(initial_mass-(state%q(1,2,2)+rho_l*first_deposition))> &
         1.0E-10_wp*initial_mass) ERROR STOP 'D18 combined mobile+immobile mass'
    temp_after_one=state%qp(4,2,2)
    IF (temp_after_one>=initial_temp) ERROR STOP 'D18 STEP cooling did not lower temperature'
    IF (ALL(lava_deposit<=0.0_wp)) ERROR STOP 'D18 lava deposit inventory missing'
    saved_q=state%q;saved_b=B_vertex;saved_lava=lava_deposit
    saved_deposit=deposit;saved_erosion=erosion;saved_erodible=erodible
    runtime%t=delta;runtime%dt=delta
    CALL write_restart_file('integrated_lava_checkpoint.bin',runtime,stochastic,state,domain)

    CALL apply_mass_exchange_step(state%q,state%qp,delta,domain,source_registry)
    second_deposition=lava_deposit(2,2)-first_deposition
    IF (second_deposition<=first_deposition) &
         ERROR STOP 'D18 cooling did not increase viscosity-controlled deposition'
    full_q=state%q;full_b=B_vertex;full_lava=lava_deposit
    full_deposit=deposit;full_erosion=erosion;full_erodible=erodible

    state%q=-999.0_wp;state%qp=-999.0_wp;lava_deposit=-999.0_wp
    CALL read_restart_file('integrated_lava_checkpoint.bin',runtime,stochastic,state,domain)
    CALL same_bits('integrated lava restart q',RESHAPE(state%q,[SIZE(saved_q)]), &
         RESHAPE(saved_q,[SIZE(saved_q)]))
    CALL same_bits('integrated lava restart bed',RESHAPE(B_vertex,[SIZE(saved_b)]), &
         RESHAPE(saved_b,[SIZE(saved_b)]))
    CALL same_bits('integrated lava restart inventory',RESHAPE(lava_deposit,[SIZE(saved_lava)]), &
         RESHAPE(saved_lava,[SIZE(saved_lava)]))
    CALL source_registry%finalize
    CALL equation_partition%finalize
    CALL configure_source_registry
    state%qp=0.0_wp;state%qp(4,:,:)=T_ambient
    DO l=1,domain%solve_cells
       j=domain%j_cent(l);k=domain%k_cent(l)
       IF (state%q(1,j,k)<=0.0_wp) CYCLE
       CALL qc_to_qp(state%q(:,j,k),state%qp(:,j,k),pdyn)
    END DO
    CALL apply_mass_exchange_step(state%q,state%qp,delta,domain,source_registry)
    CALL same_bits('integrated resumed lava q',RESHAPE(state%q,[SIZE(full_q)]), &
         RESHAPE(full_q,[SIZE(full_q)]))
    CALL same_bits('integrated resumed lava bed',RESHAPE(B_vertex,[SIZE(full_b)]), &
         RESHAPE(full_b,[SIZE(full_b)]))
    CALL same_bits('integrated resumed lava inventory',RESHAPE(lava_deposit,[SIZE(full_lava)]), &
         RESHAPE(full_lava,[SIZE(full_lava)]))
    CALL same_bits('integrated resumed granular deposit',RESHAPE(deposit,[SIZE(full_deposit)]), &
         RESHAPE(full_deposit,[SIZE(full_deposit)]))
    CALL same_bits('integrated resumed erosion',RESHAPE(erosion,[SIZE(full_erosion)]), &
         RESHAPE(full_erosion,[SIZE(full_erosion)]))
    CALL same_bits('integrated resumed erodible',RESHAPE(erodible,[SIZE(full_erodible)]), &
         RESHAPE(full_erodible,[SIZE(full_erodible)]))
    CALL source_registry%finalize;CALL equation_partition%finalize
    rheology_flag=.FALSE.;rheology_model='none'
    thermal_cooling_flag=.FALSE.;lava_deposition_flag=.FALSE.
    lava_deposition_rate=0.0_wp
    thermal_cadence='STAGE';thermal_group=0
    WRITE(*,*) 'PASS D18 lava Newton/STEP deposition/STEP cooling/restart integration'
  END SUBROUTINE test_integrated_lava_deposition_cooling_restart

  !> \brief Exercise prescribed liquid loss with zero transported solid fraction.
  SUBROUTINE test_liquid_loss_without_solids
    REAL(wp) :: before(nq,nx,ny), bed_before(nx+1,ny+1)
    REAL(wp) :: expected_mass,expected_energy,flow_temperature,carrier_loss
    REAL(wp), PARAMETER :: delta=0.1_wp
    INTEGER :: j,k,l
    CALL reset_fixture
    IF (.NOT.ALLOCATED(loss_rate)) ERROR STOP 'D18 liquid loss fixture unallocated'
    loss_rate=0.02_wp
    DO k=1,ny
       DO j=1,nx
          CALL set_cell(j,k,0.5_wp,[0.0_wp,0.0_wp],350.0_wp)
       END DO
    END DO
    state%q(:,5,4)=0.0_wp
    before=state%q
    bed_before=B_vertex
    CALL configure_source_registry
    CALL apply_mass_exchange_step(state%q,state%qp,delta,domain,source_registry)
    DO l=1,domain%solve_cells
       j=domain%j_cent(l); k=domain%k_cent(l)
       IF (before(1,j,k)<=0.0_wp) THEN
          IF (ANY(state%q(:,j,k)/=before(:,j,k))) &
               ERROR STOP 'D18 liquid loss altered dry cell'
          CYCLE
       END IF
       carrier_loss=rho_l*loss_rate*delta
       expected_mass=before(1,j,k)-carrier_loss
       flow_temperature=350.0_wp
       expected_energy=before(4,j,k)-carrier_loss*sp_heat_l*flow_temperature
       CALL assert_small('D18 zero-solid liquid mass',ABS(state%q(1,j,k)-expected_mass),expected_mass)
       CALL assert_small('D18 zero-solid liquid energy',ABS(state%q(4,j,k)-expected_energy),expected_energy)
       CALL assert_small('D18 zero-solid liquid momentum x', &
            ABS(state%q(2,j,k)-before(2,j,k)*expected_mass/before(1,j,k)),before(2,j,k))
       CALL assert_small('D18 zero-solid liquid momentum y', &
            ABS(state%q(3,j,k)-before(3,j,k)*expected_mass/before(1,j,k)),before(3,j,k))
       IF (ANY(state%q(5:6,j,k)/=0.0_wp)) ERROR STOP 'D18 liquid loss created solids'
    END DO
    CALL same_bits('D18 zero-solid unchanged bed',RESHAPE(B_vertex,[SIZE(B_vertex)]), &
         RESHAPE(bed_before,[SIZE(bed_before)]))
    IF (ANY(deposit/=0.0_wp) .OR. ANY(erosion/=0.0_wp)) &
         ERROR STOP 'D18 liquid loss changed solid inventories'
    CALL source_registry%finalize
    CALL equation_partition%finalize
    loss_rate=0.0_wp
    WRITE(*,*) 'PASS D18 prescribed carrier loss without solid fraction'
  END SUBROUTINE test_liquid_loss_without_solids

  !> \brief Ensure an empty registered post-RK schedule suppresses legacy-active sources.
  SUBROUTINE test_empty_step_registry
    TYPE(source_group_type) :: empty_step
    TYPE(source_group_type) :: no_stages(0)
    REAL(wp) :: q_before(nq,nx,ny), deposit_before(nx,ny,ns)
    REAL(wp) :: erosion_before(nx,ny,ns), inventory_before(ns,nx,ny)
    REAL(wp) :: bed_before(nx+1,ny+1)
    CHARACTER(LEN=256) :: error
    INTEGER :: i

    CALL reset_fixture
    erosion_coeff=100.0_wp
    empty_step%id=1
    empty_step%cadence=SOURCE_CADENCE_STEP
    empty_step%method=SOURCE_METHOD_EXPLICIT

    empty_step%source_names=[CHARACTER(LEN=64)::]
    empty_step%equation_ids=[(i,i=1,nq)]
    CALL source_registry%initialize(nq,0,no_stages,[empty_step],error)
    IF (LEN_TRIM(error)>0) ERROR STOP 'failed to initialize empty step registry'
    q_before=state%q; deposit_before=deposit; erosion_before=erosion
    inventory_before=erodible; bed_before=B_vertex
    CALL apply_mass_exchange_step(state%q,state%qp,1.0_wp,domain,source_registry)
    CALL same_bits('empty step q',RESHAPE(state%q,[SIZE(q_before)]),RESHAPE(q_before,[SIZE(q_before)]))
    CALL same_bits('empty step deposit',RESHAPE(deposit,[SIZE(deposit_before)]), &
         RESHAPE(deposit_before,[SIZE(deposit_before)]))
    CALL same_bits('empty step erosion',RESHAPE(erosion,[SIZE(erosion_before)]), &
         RESHAPE(erosion_before,[SIZE(erosion_before)]))
    CALL same_bits('empty step inventory',RESHAPE(erodible,[SIZE(inventory_before)]), &
         RESHAPE(inventory_before,[SIZE(inventory_before)]))
    CALL same_bits('empty step bed',RESHAPE(B_vertex,[SIZE(bed_before)]),RESHAPE(bed_before,[SIZE(bed_before)]))
    CALL source_registry%finalize()
  END SUBROUTINE test_empty_step_registry

  ! D18: explicit group members are evaluated against the same state.
  ! Exercise the shared STAGE adapter and the production STEP executor.
  SUBROUTINE test_explicit_multi_source_groups
    TYPE(source_registry_type) :: test_registry
    TYPE(source_group_type) :: joint, reversed, thermal_only, curvature_only
    REAL(wp) :: rate_joint(nq),rate_reverse(nq),rate_thermal(nq)
    REAL(wp) :: rate_curvature(nq),rate_reference(nq),before(nq,nx,ny)
    REAL(wp) :: qp_cell(nq+2),p_dyn
    REAL(wp), PARAMETER :: delta=0.01_wp
    CHARACTER(LEN=256) :: error
    INTEGER :: j,k,l,mode
    CALL reset_fixture
    thermal_cooling_flag=.TRUE.
    curvature_term_flag=.TRUE.
    radiative_term_coeff=5.103E-8_wp
    convective_term_coeff=9.0_wp
    T_env=300.0_wp
    joint%id=1
    joint%method=SOURCE_METHOD_EXPLICIT
    joint%equation_ids=[2,3,4]
    joint%source_names=[CHARACTER(LEN=64):: &
         'thermal_cooling','geometry_curvature']
    reversed=joint
    reversed%source_names=[CHARACTER(LEN=64):: &
         'geometry_curvature','thermal_cooling']
    thermal_only=joint
    thermal_only%source_names=[CHARACTER(LEN=64)::'thermal_cooling']
    thermal_only%equation_ids=[4]
    curvature_only=joint
    curvature_only%source_names=[CHARACTER(LEN=64)::'geometry_curvature']
    curvature_only%equation_ids=[2,3]
    DO mode=1,2
       IF (mode==1) THEN
          joint%cadence=SOURCE_CADENCE_STAGE
       ELSE
          joint%cadence=SOURCE_CADENCE_STEP
       END IF
       reversed%cadence=joint%cadence
       thermal_only%cadence=joint%cadence
       curvature_only%cadence=joint%cadence
       IF (mode==1) THEN
          CALL test_registry%initialize(nq,0,[joint], &
               [source_group_type::],error)
       ELSE
          CALL test_registry%initialize(nq,0,[source_group_type::], &
               [joint],error)
       END IF
       IF (LEN_TRIM(error)>0) ERROR STOP TRIM(error)
       before=state%q
       DO l=1,domain%solve_cells
          j=domain%j_cent(l);k=domain%k_cent(l)
          IF (before(1,j,k)<=0.0_wp) CYCLE
          CALL qc_to_qp(before(:,j,k),qp_cell,p_dyn)
          CALL eval_cell_source_explicit(j,k,qp_cell,0.0_wp, &
               joint,rate_joint)
          CALL eval_cell_source_explicit(j,k,qp_cell,0.0_wp, &
               reversed,rate_reverse)
          CALL eval_cell_source_explicit(j,k,qp_cell,0.0_wp, &
               thermal_only,rate_thermal)
          CALL eval_cell_source_explicit(j,k,qp_cell,0.0_wp, &
               curvature_only,rate_curvature)
          rate_reference=rate_thermal+rate_curvature
          IF (MAXVAL(ABS(rate_joint-rate_reference))> &
               tol*MAX(1.0_wp,MAXVAL(ABS(rate_reference)))) &
               ERROR STOP 'D18: explicit group did not sum initial-state rates'
          IF (ANY(rate_joint/=rate_reverse)) &
               ERROR STOP 'D18: explicit group depends on member order'
       END DO
       IF (mode==2) THEN
          DO l=1,domain%solve_cells
             j=domain%j_cent(l);k=domain%k_cent(l)
             IF (before(1,j,k)<=0.0_wp) CYCLE
             CALL qc_to_qp(before(:,j,k),state%qp(:,j,k),p_dyn)
          END DO
          CALL apply_mass_exchange_step(state%q,state%qp,delta, &
               domain,test_registry)
          DO l=1,domain%solve_cells
             j=domain%j_cent(l);k=domain%k_cent(l)
             IF (before(1,j,k)<=0.0_wp) CYCLE
             CALL qc_to_qp(before(:,j,k),qp_cell,p_dyn)
             CALL eval_cell_source_explicit(j,k,qp_cell,0.0_wp, &
                  joint,rate_reference)
             IF (MAXVAL(ABS(state%q(:,j,k)-before(:,j,k)- &
                  delta*rate_reference))> &
                  tol*MAX(1.0_wp,MAXVAL(ABS(before(:,j,k))))) &
                  ERROR STOP 'D18: explicit STEP executor mismatch'
          END DO
       END IF
       CALL test_registry%finalize
    END DO
    thermal_cooling_flag=.FALSE.
    curvature_term_flag=.FALSE.
    WRITE(*,*) 'PASS D18 explicit multi-source STAGE/STEP common-state rates'
  END SUBROUTINE test_explicit_multi_source_groups

  !> \\brief Check sequential STEP groups with simultaneous cooling and erosion.
  !> \\details Compare the combined production executor with separately
  !>          executed registry groups and retain bed/inventory comparisons.
  SUBROUTINE test_combined_step_sources
    TYPE(source_registry_type) :: first_group, second_group
    TYPE(source_group_type), ALLOCATABLE :: groups(:)
    REAL(wp) :: q_full(nq,nx,ny), bed_full(nx+1,ny+1)
    REAL(wp) :: dep_full(nx,ny,ns), ers_full(nx,ny,ns), inv_full(ns,nx,ny)
    REAL(wp) :: initial_bed(nx+1,ny+1)
    CHARACTER(LEN=256) :: error
    INTEGER :: i,mass_idx,cooling_idx

    CALL reset_fixture
    erosion_coeff=0.1_wp
    thermal_cooling_flag=.TRUE.
    thermal_method='ANALYTIC'; thermal_cadence='STEP'; thermal_group=1
    radiative_term_coeff=5.103E-8_wp
    convective_term_coeff=9.0_wp
    T_env=300.0_wp
    CALL configure_source_registry
    groups=source_registry%step_groups()
    mass_idx=0; cooling_idx=0
    DO i=1,SIZE(groups)
       IF (ANY(groups(i)%source_names=='mass_exchange_erosion')) mass_idx=i
       IF (ANY(groups(i)%source_names=='thermal_cooling')) cooling_idx=i
    END DO
    IF (mass_idx==0 .OR. cooling_idx==0) &
         ERROR STOP 'D18 combined STEP missing mass exchange or cooling'
    IF (mass_idx>=cooling_idx) ERROR STOP 'D18 combined STEP group ordering'
    CALL first_group%initialize(nq,0,[source_group_type::],[groups(mass_idx)],error)
    IF (LEN_TRIM(error)>0) ERROR STOP TRIM(error)
    CALL second_group%initialize(nq,0,[source_group_type::],[groups(cooling_idx)],error)
    IF (LEN_TRIM(error)>0) ERROR STOP TRIM(error)
    initial_bed=B_vertex
    CALL apply_mass_exchange_step(state%q,state%qp,0.1_wp,domain,source_registry)
    q_full=state%q; bed_full=B_vertex
    dep_full=deposit; ers_full=erosion; inv_full=erodible
    IF (ALL(bed_full==initial_bed)) ERROR STOP 'D18 combined STEP exchange inactive'

    CALL reset_fixture
    erosion_coeff=0.1_wp
    CALL apply_mass_exchange_step(state%q,state%qp,0.1_wp,domain,first_group)
    CALL apply_mass_exchange_step(state%q,state%qp,0.1_wp,domain,second_group)
    CALL same_bits('D18 combined STEP Q',RESHAPE(state%q,[SIZE(q_full)]),RESHAPE(q_full,[SIZE(q_full)]))
    CALL same_bits('D18 combined STEP bed',RESHAPE(B_vertex,[SIZE(bed_full)]),RESHAPE(bed_full,[SIZE(bed_full)]))
    CALL same_bits('D18 combined STEP deposit',RESHAPE(deposit,[SIZE(dep_full)]),RESHAPE(dep_full,[SIZE(dep_full)]))
    CALL same_bits('D18 combined STEP erosion',RESHAPE(erosion,[SIZE(ers_full)]),RESHAPE(ers_full,[SIZE(ers_full)]))
    CALL same_bits('D18 combined STEP inventory',RESHAPE(erodible,[SIZE(inv_full)]),RESHAPE(inv_full,[SIZE(inv_full)]))
    CALL first_group%finalize
    CALL second_group%finalize
    CALL source_registry%finalize
    CALL equation_partition%finalize
    thermal_cooling_flag=.FALSE.
    thermal_method='ANALYTIC';thermal_cadence='STAGE';thermal_group=0
    WRITE(*,*) 'PASS D18 combined STEP mass exchange and analytic cooling'
  END SUBROUTINE test_combined_step_sources

  ! D18: official production restart round trip with analytical STEP cooling.
  SUBROUTINE test_official_cooling_restart
    REAL(wp) :: checkpoint_state(nq,nx,ny),continuous(nq,nx,ny)
    REAL(wp) :: checkpoint_bed(nx+1,ny+1), continuous_bed(nx+1,ny+1)
    REAL(wp) :: checkpoint_deposit(nx,ny,ns),continuous_deposit(nx,ny,ns)
    REAL(wp) :: checkpoint_erosion(nx,ny,ns),continuous_erosion(nx,ny,ns)
    REAL(wp) :: checkpoint_inventory(ns,nx,ny),continuous_inventory(ns,nx,ny)
    REAL(wp), PARAMETER :: delta=0.05_wp
    CHARACTER(LEN=256) :: error
    INTEGER :: l,j,k
    REAL(wp) :: pressure
    CALL reset_fixture
    erosion_coeff=0.1_wp
    thermal_cooling_flag=.TRUE.
    thermal_method='ANALYTIC'
    thermal_cadence='STEP'
    thermal_group=1
    radiative_term_coeff=5.103E-8_wp
    convective_term_coeff=9.0_wp
    T_env=300.0_wp
    CALL initialize_configured_sources(source_registry,equation_partition,error)
    IF (LEN_TRIM(error)>0) ERROR STOP TRIM(error)
    DO l=1,domain%solve_cells
       j=domain%j_cent(l);k=domain%k_cent(l)
       IF (state%q(1,j,k)<=0.0_wp) CYCLE
       CALL qc_to_qp(state%q(:,j,k),state%qp(:,j,k),pressure)
    END DO
    CALL apply_mass_exchange_step(state%q,state%qp,delta,domain,source_registry)
    checkpoint_state=state%q
    checkpoint_bed=B_vertex;checkpoint_deposit=deposit
    checkpoint_erosion=erosion;checkpoint_inventory=erodible
    IF (ALL(checkpoint_erosion==0.0_wp)) ERROR STOP 'D18 restart erosion inactive'
    runtime%t=delta;runtime%dt=delta
    CALL write_restart_file('cooling_checkpoint.bin',runtime,stochastic,state,domain)
    CALL apply_mass_exchange_step(state%q,state%qp,delta,domain,source_registry)
    continuous=state%q
    continuous_bed=B_vertex;continuous_deposit=deposit
    continuous_erosion=erosion;continuous_inventory=erodible
    state%q=-999.0_wp
    state%qp=-999.0_wp
    CALL read_restart_file('cooling_checkpoint.bin',runtime,stochastic,state,domain)
    CALL same_bits('official cooling restart state', &
         RESHAPE(state%q,[SIZE(checkpoint_state)]), &
         RESHAPE(checkpoint_state,[SIZE(checkpoint_state)]))
    CALL same_bits('official restart B_vertex',RESHAPE(B_vertex,[SIZE(checkpoint_bed)]), &
         RESHAPE(checkpoint_bed,[SIZE(checkpoint_bed)]))
    CALL same_bits('official restart deposit',RESHAPE(deposit,[SIZE(checkpoint_deposit)]), &
         RESHAPE(checkpoint_deposit,[SIZE(checkpoint_deposit)]))
    CALL same_bits('official restart erosion',RESHAPE(erosion,[SIZE(checkpoint_erosion)]), &
         RESHAPE(checkpoint_erosion,[SIZE(checkpoint_erosion)]))
    CALL same_bits('official restart erodible',RESHAPE(erodible,[SIZE(checkpoint_inventory)]), &
         RESHAPE(checkpoint_inventory,[SIZE(checkpoint_inventory)]))
    CALL source_registry%finalize
    CALL equation_partition%finalize
    CALL initialize_configured_sources(source_registry,equation_partition,error)
    IF (LEN_TRIM(error)>0) ERROR STOP TRIM(error)
    state%qp=0.0_wp
    state%qp(4,:,:)=T_ambient
    DO l=1,domain%solve_cells
       j=domain%j_cent(l);k=domain%k_cent(l)
       IF (state%q(1,j,k)<=0.0_wp) CYCLE
       CALL qc_to_qp(state%q(:,j,k),state%qp(:,j,k),pressure)
    END DO
    CALL apply_mass_exchange_step(state%q,state%qp,delta,domain,source_registry)
    CALL same_bits('official cooling restart continuation', &
         RESHAPE(state%q,[SIZE(continuous)]), &
         RESHAPE(continuous,[SIZE(continuous)]))
    CALL same_bits('official restart resumed bed',RESHAPE(B_vertex,[SIZE(continuous_bed)]), &
         RESHAPE(continuous_bed,[SIZE(continuous_bed)]))
    CALL same_bits('official restart resumed deposit',RESHAPE(deposit,[SIZE(continuous_deposit)]), &
         RESHAPE(continuous_deposit,[SIZE(continuous_deposit)]))
    CALL same_bits('official restart resumed erosion',RESHAPE(erosion,[SIZE(continuous_erosion)]), &
         RESHAPE(continuous_erosion,[SIZE(continuous_erosion)]))
    CALL same_bits('official restart resumed erodible',RESHAPE(erodible,[SIZE(continuous_inventory)]), &
         RESHAPE(continuous_inventory,[SIZE(continuous_inventory)]))
    CALL source_registry%finalize
    CALL equation_partition%finalize
    erosion_coeff=0.0_wp
    thermal_cooling_flag=.FALSE.
    thermal_method='ANALYTIC'
    thermal_cadence='STAGE'
    thermal_group=0
    WRITE(*,*) 'PASS D18 production restart with analytic cooling STEP'
  END SUBROUTINE test_official_cooling_restart

  !> D18: production post-RK dispatcher applies exact cooling once over dt.
  SUBROUTINE test_analytic_cooling_step_executor
    REAL(wp) :: prior(nq,nx,ny), prior_qp(nq+2,nx,ny)
    REAL(wp) :: temperature,cap,expected_temperature,expected_energy,p_dyn
    REAL(wp) :: before_stage(nq),after_stage(nq),im_rate(nq),residual(nq)
    INTEGER :: iterations,linear_info
    LOGICAL :: converged,line_search_failed
    REAL(wp), PARAMETER :: delta=0.15_wp
    INTEGER :: j,k,l
    CHARACTER(LEN=256) :: error

    CALL reset_fixture
    ! Retain the fixture's liquid carrier; zero loss_rate makes the
    ! mass-exchange transaction a no-op in this cooling-only experiment.
    thermal_cooling_flag=.TRUE.
    rheology_flag=.TRUE.;rheology_model='lava'
    rheology_method='NEWTON';rheology_cadence='STAGE';rheology_group=1
    nu_ref=0.02_wp;visc_par=0.003_wp;tau0=0.0_wp;T_ref=900.0_wp
    radiative_term_coeff=5.103E-8_wp
    convective_term_coeff=9.0_wp
    T_env=300.0_wp
    thermal_method='ANALYTIC'
    thermal_cadence='STEP'
    thermal_group=1
    CALL initialize_configured_sources(source_registry,equation_partition,error)
    IF (LEN_TRIM(error)>0) ERROR STOP TRIM(error)
    IF (equation_partition%n_IM_groups/=1) &
         ERROR STOP 'D18: expected only lava Newton in IMEX stage groups'
    CALL initialize_nonlinear_solver
    ! Apply the actual stage coordinator to a wet cell, then feed its output
    ! into the production post-RK executor on the full domain fixture.
    before_stage=state%q(:,2,2)
    residual=0.0_wp
    CALL solve_imex_stage(equation_partition,source_registry%stage_groups(), &
         before_stage,before_stage,delta,0.75_wp,residual, &
         B_prime_x_geom(2,2),B_prime_y_geom(2,2), &
         B_second_xx_geom(2,2),B_second_xy_geom(2,2), &
         B_second_yy_geom(2,2),grav_coeff(2,2),0.0_wp, &
         .FALSE.,.FALSE.,after_stage,im_rate,iterations,converged, &
         linear_info,line_search_failed)
    IF (.NOT.converged .OR. linear_info/=0 .OR. line_search_failed) &
         ERROR STOP 'D18: lava NEWTON/STAGE stage did not converge'
    IF (MAXVAL(ABS(after_stage(2:3)-before_stage(2:3)))<=0.0_wp) &
         ERROR STOP 'D18: lava NEWTON/STAGE did not alter momentum'
    IF (after_stage(4)/=before_stage(4)) &
         ERROR STOP 'D18: STEP cooling leaked into Newton stage'
    state%q(:,2,2)=after_stage
    CALL finalize_nonlinear_solver
    DO l=1,domain%solve_cells
      j=domain%j_cent(l); k=domain%k_cent(l)
      IF (state%q(1,j,k)<=0.0_wp) CYCLE
      CALL qc_to_qp(state%q(:,j,k),state%qp(:,j,k),p_dyn)
    END DO
    prior=state%q;prior_qp=state%qp
    CALL apply_mass_exchange_step(state%q,state%qp,delta,domain,source_registry)
    DO l=1,domain%solve_cells
      j=domain%j_cent(l);k=domain%k_cent(l)
      IF (prior(1,j,k)<=0.0_wp) THEN
        IF (ANY(state%q(:,j,k)/=prior(:,j,k))) &
             ERROR STOP 'D18: analytic STEP modified dry conservative state'
        CYCLE
      END IF
      temperature=prior_qp(4,j,k)
      cap=prior(4,j,k)/temperature
      CALL solve_thermal_cooling_analytic(temperature,delta,cap, &
           radiative_term_coeff,convective_term_coeff,T_env,expected_temperature)
      expected_energy=cap*expected_temperature
      IF (ABS(state%q(4,j,k)-expected_energy)> &
           1.0E-12_wp*MAX(1.0_wp,ABS(expected_energy))) &
           ERROR STOP 'D18: STEP cooling does not match exact full-step map'
      IF (ANY(state%q([1,2,3,5,6],j,k)/=prior([1,2,3,5,6],j,k))) &
           ERROR STOP 'D18: STEP cooling changed unrelated conservative equations'
      IF (ABS(state%qp(4,j,k)-expected_temperature)> &
           1.0E-12_wp*MAX(1.0_wp,ABS(expected_temperature))) &
           ERROR STOP 'D18: STEP cooling did not synchronize primitive state'
    END DO
    IF (ANY(state%q(:,nx,ny)/=prior(:,nx,ny))) &
         ERROR STOP 'D18: STEP cooling modified inactive cell'
    ! Exact cooling is a flow map: two full steps must equal one
    ! analytic application over twice the physical interval.
    CALL apply_mass_exchange_step(state%q,state%qp,delta,domain,source_registry)
    DO l=1,domain%solve_cells
      j=domain%j_cent(l);k=domain%k_cent(l)
      IF (prior(1,j,k)<=0.0_wp) CYCLE
      temperature=prior_qp(4,j,k)
      cap=prior(4,j,k)/temperature
      CALL solve_thermal_cooling_analytic(temperature,2.0_wp*delta,cap, &
           radiative_term_coeff,convective_term_coeff,T_env,expected_temperature)
      expected_energy=cap*expected_temperature
      IF (ABS(state%q(4,j,k)-expected_energy)> &
           2.0E-12_wp*MAX(1.0_wp,ABS(expected_energy))) &
           ERROR STOP 'D18: cooling STEP violates exact semigroup property'
      IF (ANY(state%q([1,2,3,5,6],j,k)/=prior([1,2,3,5,6],j,k))) &
           ERROR STOP 'D18: repeated cooling altered unrelated equations'
    END DO
    CALL source_registry%finalize
    CALL equation_partition%finalize
    thermal_cooling_flag=.FALSE.;rheology_flag=.FALSE.;rheology_model='none'
    rheology_method='NEWTON';rheology_cadence='STAGE';rheology_group=0
    thermal_method='ANALYTIC';thermal_cadence='STAGE';thermal_group=0
    WRITE(*,*) 'PASS D18 lava Newton stage followed by exact cooling STEP executor'
  END SUBROUTINE test_analytic_cooling_step_executor

  !> \brief Reset conservative inventories and an analytic initial bed before each independent case.
  SUBROUTINE reset_fixture
    INTEGER :: j,k
    REAL(wp) :: h, phi(ns), temperature, x,y
    erosion_coeff=0.0_wp; settling_flag=.FALSE.
    IF (ALLOCATED(loss_rate)) loss_rate=0.0_wp
    alphastot_min=0.0_wp; erodible_porosity=0.25_wp
    coeff_porosity=erodible_porosity/(1.0_wp-erodible_porosity)
    entrainment_flag=.FALSE.
    erodible_deposit_flag=.FALSE.; bottom_radial_source_flag=.FALSE.
    bottom_fissural_source_flag=.FALSE.
    T_erodible=317.0_wp; verbose_level=2
    deposit=0.0_wp; erosion=0.0_wp; erodible=0.0_wp;lava_deposit=0.0_wp
    B_zone=0; cell_source_fractions=0.0_wp
    DO k=1,ny+1
       y=REAL(k-1,wp)*dy
       DO j=1,nx+1
          x=REAL(j-1,wp)*dx
          B_vertex(j,k)=2.0_wp+0.03_wp*x-0.02_wp*y+0.002_wp*x*x+0.001_wp*x*y
       END DO
    END DO
    CALL refresh_topography_geometry
    DO k=1,ny
       DO j=1,nx
          h=0.5_wp+0.02_wp*REAL(j+k,wp)
          phi=[0.10_wp,0.04_wp]
          temperature=300.0_wp+2.0_wp*REAL(j,wp)+3.0_wp*REAL(k,wp)
          CALL set_cell(j,k,h,phi,temperature)
          erodible(:,j,k)=[0.002_wp,0.001_wp]*(1.0_wp+0.01_wp*REAL(j+k,wp))
       END DO
    END DO
    state%q(:,5,4)=0.0_wp; state%q(:,7,6)=0.0_wp
  END SUBROUTINE reset_fixture

  !> \brief Construct conserved masses directly from known phase volumes, avoiding inverse test logic.
  !> \param[in] j,k Cell indices.
  !> \param[in] h Thickness [m].
  !> \param[in] phi Solid volume fractions.
  !> \param[in] temperature Flow temperature [K].
  SUBROUTINE set_cell(j,k,h,phi,temperature)
    INTEGER, INTENT(IN) :: j,k
    REAL(wp), INTENT(IN) :: h, phi(ns), temperature
    REAL(wp) :: carrier_mass, carrier_density, carrier_heat
    carrier_density=rho_l; carrier_heat=sp_heat_l
    IF (gas_flag) THEN
       carrier_density=pres/(sp_gas_const_a*temperature)
       carrier_heat=sp_heat_a
    END IF
    carrier_mass=carrier_density*h*(1.0_wp-SUM(phi))
    state%q(5:6,j,k)=rho_s*h*phi
    state%q(1,j,k)=carrier_mass+SUM(state%q(5:6,j,k))
    state%q(2,j,k)=state%q(1,j,k)*(1.0_wp+0.01_wp*REAL(j,wp))
    state%q(3,j,k)=state%q(1,j,k)*(-0.3_wp+0.01_wp*REAL(k,wp))
    state%q(4,j,k)=temperature*(carrier_mass*carrier_heat+DOT_PRODUCT(state%q(5:6,j,k),sp_heat_s))
    IF (pore_pressure_flag) state%q(idx_pore,j,k)=state%q(1,j,k)*1000.0_wp
  END SUBROUTINE set_cell

  !> \brief Select a saturation/guard case and apply the real production transaction.
  !> \param[in] name Case name, also used by the diagnostic comparator.
  SUBROUTINE run_case(name)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(wp) :: dt
    INTEGER :: j,k
    CALL reset_fixture
    dt=1.0_wp
    SELECT CASE(name)
    CASE('deposition')
       settling_flag=.TRUE.; dt=1.0E6_wp
    CASE('erosion')
       erosion_coeff=100.0_wp
    CASE('combined','masked','masked_fissural')
       erosion_coeff=100.0_wp; settling_flag=.TRUE.
       erodible_deposit_flag=.TRUE.; dt=1.0E6_wp
       IF (name=='masked' .OR. name=='masked_fissural') THEN
          bottom_radial_source_flag=name=='masked'
          bottom_fissural_source_flag=name=='masked_fissural'
          cell_source_fractions(1,1)=1.0_wp
          cell_source_fractions(2,1)=0.5_wp
          cell_source_fractions(3,1)=0.25_wp
       END IF
    CASE('cutoff')
       erosion_coeff=100.0_wp; alphastot_min=0.3_wp
    CASE('loss_only')
       loss_rate=1000.0_wp
    CASE('loss_reserve')
       loss_rate=1000.0_wp
       ! Activate the exchange evaluator through erosion as well, with an
       ! empty substrate: the carrier reserve must never become a gain term.
       erosion_coeff=1.0_wp; erodible=0.0_wp
       erodible_porosity=0.7_wp
       coeff_porosity=erodible_porosity/(1.0_wp-erodible_porosity)
       DO k=1,ny
          DO j=1,nx
             IF (state%q(1,j,k)>0.0_wp) CALL set_cell(j,k,0.8_wp,[0.38_wp,0.17_wp],350.0_wp)
          END DO
       END DO
    CASE('thermal')
       erosion_coeff=100.0_wp; erodible_deposit_flag=.TRUE.
    END SELECT
    CALL configure_source_registry
    CALL check_update(name,dt)
    IF (name=='thermal') CALL check_thermal_entrainment
  END SUBROUTINE run_case

  !> \brief Verify that production thermal assembly includes entrained ambient-air energy.
  SUBROUTINE check_thermal_entrainment
    REAL(wp) :: qp_local(n_vars+2), dynamic_pressure
    REAL(wp) :: richardson, mixture_density, carrier_density, reduced_gravity
    REAL(wp) :: carrier_heat_capacity, mixture_heat_capacity
    REAL(wp) :: erosion_rate(ns), deposition_rate(ns)
    REAL(wp) :: carrier_erosion_rate, carrier_loss_rate
    REAL(wp) :: equation_rate(n_eqns), topo_rate, air_rate, expected
    REAL(wp) :: speed_squared

    CALL reset_fixture
    erosion_coeff=100.0_wp; erodible_deposit_flag=.TRUE.
    entrainment_flag=.TRUE.
    CALL qc_to_qp(state%q(:,1,1), qp_local, dynamic_pressure)
    CALL mixt_var(qp_local, richardson, mixture_density, carrier_density,   &
         reduced_gravity, carrier_heat_capacity, mixture_heat_capacity)
    IF (MIN(mixture_density,carrier_density,reduced_gravity,                 &
         carrier_heat_capacity,mixture_heat_capacity)<=0.0_wp)               &
         ERROR STOP 'invalid thermal entrainment fixture properties'
    speed_squared=qp_local(idx_u)**2+qp_local(idx_v)**2
    air_rate=eval_air_entrainment_rate(richardson, speed_squared)
    IF (air_rate<=0.0_wp) ERROR STOP 'thermal entrainment fixture is inactive'

    CALL eval_mass_exchange_terms(qp_local,0,0.0_wp,0.0_wp,erodible(:,1,1), &
         1.0_wp,erosion_rate,deposition_rate,carrier_erosion_rate,          &
         carrier_loss_rate,equation_rate,topo_rate)
    expected=-qp_local(4)*(SUM(rho_s*sp_heat_s*deposition_rate) +           &
         carrier_density*carrier_heat_capacity*carrier_loss_rate) +        &
         qp_local(4)*(SUM(rho_s*sp_heat_s*erosion_rate) +                   &
         rho_c_sub*carrier_heat_capacity*carrier_erosion_rate) +            &
         T_ambient*sp_heat_a*rho_a_amb*air_rate
    CALL assert_small('thermal entrainment energy source',                  &
         ABS(equation_rate(4)-expected),MAX(1.0_wp,ABS(expected)))
  END SUBROUTINE check_thermal_entrainment

  !> \brief Exercise the existing signed pore-pressure loss law with an independent gas fixture.
  SUBROUTINE run_gas_cases
    INTEGER :: j,k,c
    REAL(wp) :: dt,phi(ns)
    CHARACTER(LEN=16), PARAMETER :: names(6)=[CHARACTER(LEN=16) :: &
         'gas_packing','gas_reserve','gas_exhausted','gas_inflow','gas_combined','gas_masked']
    CALL state%finalize
    CALL layout%finalize
    nq=7; n_vars=nq; n_eqns=nq; n_pore_vars=1
    idx_pore=7; idx_poreEqn=7; idx_u=8; idx_v=9
    gas_flag=.TRUE.; liquid_flag=.FALSE.; pore_pressure_flag=.TRUE.; gas_loss_flag=.TRUE.
    dynamic_permeability_flag=.FALSE.; f_inhibit_mode='OFF'
    kin_visc_c=1.0E-5_wp; hydraulic_permeability=1.0E-6_wp
    rho_c_sub=1.2_wp; pi_g=ACOS(-1.0_wp)
    ! Gas-only execution must not require the optional liquid-loss allocation.
    DEALLOCATE(loss_rate)
    CALL layout%initialize(1,nq,nq,ns,0,0,1,.FALSE.)
    CALL state%initialize(layout)
    OPEN(NEWUNIT=snapshot_unit,FILE='gas_snapshot.bin',ACCESS='STREAM', &
         FORM='UNFORMATTED',STATUS='REPLACE')
    DO c=1,SIZE(names)
       CALL reset_fixture
       dt=1.0_wp
       SELECT CASE(TRIM(names(c)))
       CASE('gas_reserve','gas_exhausted')
          phi=[0.25_wp,0.15_wp]; erodible_porosity=0.55_wp
          IF (names(c)=='gas_exhausted') THEN
             phi=[0.38_wp,0.17_wp]; erodible_porosity=0.7_wp
          END IF
          coeff_porosity=erodible_porosity/(1.0_wp-erodible_porosity)
          DO k=1,ny
             DO j=1,nx
                IF (state%q(1,j,k)>0.0_wp) CALL set_cell(j,k,0.8_wp,phi,350.0_wp)
             END DO
          END DO
       CASE('gas_inflow')
          state%q(idx_pore,:,:)=-0.01_wp*state%q(1,:,:)
       CASE('gas_combined','gas_masked')
          erosion_coeff=100.0_wp; settling_flag=.TRUE.; erodible_deposit_flag=.TRUE.
          dt=1.0E6_wp
          IF (names(c)=='gas_masked') THEN
             bottom_radial_source_flag=.TRUE.
             cell_source_fractions(1,1)=1.0_wp
             cell_source_fractions(2,1)=0.5_wp
             cell_source_fractions(3,1)=0.25_wp
          END IF
       END SELECT
       CALL configure_source_registry
       CALL check_update(TRIM(names(c)),dt)
    END DO
    CLOSE(snapshot_unit)
  END SUBROUTINE run_gas_cases

  !> \brief Compare one flat-cell evaluation with four identical vertex evaluations, including guards.
  !> \param[in] name Case label.
  !> \param[in] dt Transaction duration [s].
  SUBROUTINE check_flat_evaluations(name,dt)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(wp), INTENT(IN) :: dt
    REAL(wp) :: qp0(nq+2), qp_saved(nq+2), inventory(ns), saved_inventory(ns), p_dyn
    REAL(wp) :: one(2*ns+2+nq+1), four(2*ns+2+nq+1,4), config(8), saved_config(8)
    INTEGER :: vertex, guard
    DO guard=1,2
       CALL qc_to_qp(state%q(:,1,1),qp0,p_dyn)
       IF (guard==2) qp0(1)=0.0_wp
       inventory=erodible(:,1,1); qp_saved=qp0; saved_inventory=inventory
       saved_config=[T_erodible,kin_visc_c,hydraulic_permeability,rho_c_sub, &
            erosion_coeff,erodible_porosity,coeff_porosity,maximum_solid_packing]
       CALL evaluate_local(qp0,inventory,dt,one)
       !$OMP PARALLEL DO DEFAULT(NONE) SHARED(qp0,inventory,dt,four) PRIVATE(vertex)
       DO vertex=1,4
          CALL evaluate_local(qp0,inventory,dt,four(:,vertex))
       END DO
       !$OMP END PARALLEL DO
       DO vertex=1,4
          CALL same_bits('flat identical vertex',four(:,vertex),one)
       END DO
       CALL assert_small('flat quarter-weighted sum',MAXVAL(ABS(0.25_wp*SUM(four,DIM=2)-one)), &
            MAXVAL(ABS(one)))
       IF (guard==2) CALL assert_small('dry flat outputs',MAXVAL(ABS(one)),1.0_wp)
       config=[T_erodible,kin_visc_c,hydraulic_permeability,rho_c_sub, &
            erosion_coeff,erodible_porosity,coeff_porosity,maximum_solid_packing]
       CALL same_bits('read-only flat configuration',config,saved_config)
       CALL same_bits('read-only physical input',qp0,qp_saved)
       CALL same_bits('read-only erodible input',inventory,saved_inventory)
       WRITE(extension_unit) one,four
    END DO
    WRITE(*,*) 'N8_FLAT_EQUIVALENCE ',name
  END SUBROUTINE check_flat_evaluations

  !> \brief Pack all production outputs evaluated with zero slope and identical cell inputs.
  !> \param[in] qp0 Physical base state, including velocity and optional pore pressure.
  !> \param[in] inventory Available solid substrate [m].
  !> \param[in] dt Transaction duration [s].
  !> \param[out] packed Erosion, deposition, carrier rates, equation rates and bed rate.
  SUBROUTINE evaluate_local(qp0,inventory,dt,packed)
    REAL(wp), INTENT(IN) :: qp0(n_vars+2),inventory(ns),dt
    REAL(wp), INTENT(OUT) :: packed(2*ns+2+n_eqns+1)
    CALL eval_mass_exchange_terms(qp0,0,0.0_wp,0.0_wp,inventory,dt, &
         packed(1:ns),packed(ns+1:2*ns),packed(2*ns+1),packed(2*ns+2), &
         packed(2*ns+3:2*ns+2+n_eqns),packed(2*ns+3+n_eqns))
  END SUBROUTINE evaluate_local

  !> \brief Independently limit two-solid liquid or signed pore-pressure gas exchange increments.
  !> \param[in] q0 Conservative base state.
  !> \param[in] inventory Available solid volume per area [m].
  !> \param[in] bx,by Filtered slopes used in the existing erosion law.
  !> \param[in] dt Step duration [s].
  !> \param[out] dep,ers Limited solid volume increments [m].
  !> \param[out] loss Limited carrier volume loss [m].
  !> \param[out] h,u,v,temperature Physical base state recovered analytically.
  SUBROUTINE reference_exchange(q0,inventory,bx,by,dt,dep,ers,loss,h,u,v,temperature)
    REAL(wp), INTENT(IN) :: q0(nq),inventory(ns),bx,by,dt
    REAL(wp), INTENT(OUT) :: dep(ns),ers(ns),loss,h,u,v,temperature
    REAL(wp) :: solid_volume(ns), phi(ns), carrier_mass, speed, settling, reserve, carrier_density, carrier_heat
    INTEGER :: i
    dep=0.0_wp; ers=0.0_wp; loss=0.0_wp
    h=0.0_wp; u=0.0_wp; v=0.0_wp; temperature=T_ambient
    IF (q0(1)<=0.0_wp) RETURN
    solid_volume=q0(5:6)/rho_s
    carrier_mass=q0(1)-SUM(q0(5:6))
    carrier_density=rho_l; carrier_heat=sp_heat_l
    IF (gas_flag) THEN
       carrier_heat=sp_heat_a
       temperature=q0(4)/(carrier_mass*carrier_heat+DOT_PRODUCT(q0(5:6),sp_heat_s))
       carrier_density=pres/(sp_gas_const_a*temperature)
    END IF
    h=SUM(solid_volume)+carrier_mass/carrier_density
    phi=solid_volume/h
    u=q0(2)/q0(1); v=q0(3)/q0(1)
    temperature=q0(4)/(carrier_mass*carrier_heat+DOT_PRODUCT(q0(5:6),sp_heat_s))
    IF (SUM(phi)<=alphastot_min) RETURN
    speed=SQRT(u*u+v*v+(u*bx+v*by)**2)
    ers=MAX(0.0_wp,MIN(inventory,dt*erosion_coeff*speed*h*(1.0_wp-SUM(phi)) &
         *(1.0_wp-erodible_porosity)*erodible_fract))
    IF (settling_flag) THEN
       DO i=1,ns
          settling=settling_velocity(diam_s(i),rho_s(i),carrier_density,1.0_wp/kin_visc_c)
          dep(i)=MIN(solid_volume(i),dt*phi(i)*settling                    &
               *(1.0_wp-MIN(1.0_wp,SUM(phi)/0.6_wp))**4.65_wp)
       END DO
    END IF
    ! Independent volume-balance oracle for the coupled deposition and loss.
    BLOCK
       REAL(wp) :: ds,es,ec,dc,requested,inventory,margin,k,factor
       ds=SUM(dep)
       es=SUM(ers)
       ec=coeff_porosity*es
       dc=coeff_porosity*ds

       IF (gas_flag) THEN

          requested=MAX(0.0_wp,dt*hydraulic_permeability/(kin_visc_c*carrier_density) &
               /MAX(1.0E-5_wp,h)*0.5_wp*ACOS(-1.0_wp)*q0(idx_pore)/q0(1))

       ELSE

          requested=MAX(0.0_wp,dt*loss_rate)
       END IF

       inventory=h*MAX(0.0_wp,1.0_wp-SUM(phi))
       factor=1.0_wp

       IF (dc+requested>0.0_wp) factor=MIN(factor,inventory/(dc+requested))

       margin=maximum_solid_packing*(h+es+ec)-h*SUM(phi)-es
       k=maximum_solid_packing*(ds+dc+requested)-ds

       IF (k>0.0_wp) factor=MIN(factor,MAX(0.0_wp,margin)/k)

       factor=MAX(0.0_wp,MIN(1.0_wp,factor))
       dep=dep*factor
       loss=(dc+requested)*factor
    END BLOCK
  END SUBROUTINE reference_exchange

  !> \brief Assert the complete production transaction and projection after all cell-local limits.
  !> \param[in] name Case label for failures and diagnostic output.
  !> \param[in] dt Duration of the production update [s].
  SUBROUTINE check_update(name,dt)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(wp), INTENT(IN) :: dt
    REAL(wp) :: q0(nq,nx,ny), dep0(nx,ny,ns), ers0(nx,ny,ns), inv0(ns,nx,ny)
    REAL(wp) :: bed0(nx+1,ny+1), q_expected(nq,nx,ny), inv_expected(ns,nx,ny)
    REAL(wp) :: dep_expected(nx,ny,ns), ers_expected(nx,ny,ns), cell_delta(nx,ny)
    REAL(wp) :: vertex_delta(nx+1,ny+1), geom_delta(nx,ny), mismatch(nx,ny)
    REAL(wp) :: dep(ns),ers(ns),loss,h,u,v,temperature,substrate_temperature,mask
    REAL(wp) :: outgoing_mass, entering_mass, weight, nodal_volume, diagnostics(5), carrier_density, carrier_heat
    INTEGER :: j,k,l,j0,j1,k0,k1,i
    CALL check_flat_evaluations(name,dt)
    q0=state%q; bed0=B_vertex; dep0=deposit; ers0=erosion; inv0=erodible
    q_expected=q0; dep_expected=dep0; ers_expected=ers0; inv_expected=inv0
    cell_delta=0.0_wp
    DO l=1,domain%solve_cells
       j=domain%j_cent(l); k=domain%k_cent(l)
       CALL reference_exchange(q0(:,j,k),inv0(:,j,k),B_prime_x_geom(j,k), &
            B_prime_y_geom(j,k),dt,dep,ers,loss,h,u,v,temperature)
       IF (j==1 .AND. k==1) THEN
          SELECT CASE(name)
          CASE('deposition','combined','masked','masked_fissural')
             CALL assert_small(name//' saturated deposition',MAXVAL(ABS(dep-q0(5:6,j,k)/rho_s)),1.0_wp)
          END SELECT
          SELECT CASE(name)
          CASE('erosion','combined','masked','masked_fissural','thermal')
             CALL assert_small(name//' saturated erosion',MAXVAL(ABS(ers-inv0(:,j,k))),1.0_wp)
          END SELECT
          IF (name=='loss_only' .AND. loss<=tol) ERROR STOP 'carrier loss not exercised'
          IF (name=='loss_reserve' .AND. loss/=0.0_wp) ERROR STOP 'carrier reserve not exhausted'
          SELECT CASE(name)
          CASE('gas_packing')
             IF (loss<=tol) ERROR STOP 'gas packing loss not exercised'
             IF (SUM(q0(5:6,j,k)/rho_s)>maximum_solid_packing*(h-loss)+tol) &
                  ERROR STOP 'gas packing limit violated'
          CASE('gas_reserve')
             IF (loss<=tol) ERROR STOP 'gas available reserve not exercised'
             IF (loss>h*(1.0_wp-0.4_wp)+tol) ERROR STOP 'gas inventory overdrawn'
          CASE('gas_exhausted')
             IF (loss/=0.0_wp) ERROR STOP 'gas exhausted reserve became a gain'
          CASE('gas_inflow')
             IF (loss<0.0_wp) ERROR STOP 'negative gas withdrawal'
          CASE('gas_combined','gas_masked')
             CALL assert_small('gas saturated deposition',MAXVAL(ABS(dep-q0(5:6,j,k)/rho_s)),1.0_wp)
             CALL assert_small('gas saturated erosion',MAXVAL(ABS(ers-inv0(:,j,k))),1.0_wp)
             IF (loss<=tol) ERROR STOP 'combined gas loss not exercised'
          END SELECT
       END IF
       mask=1.0_wp
       IF (bottom_radial_source_flag .OR. bottom_fissural_source_flag) &
            mask=1.0_wp-cell_source_fractions(j,k)
       dep=dep*mask; ers=ers*mask; loss=loss*mask
       dep_expected(j,k,:)=dep0(j,k,:)+dep
       ers_expected(j,k,:)=ers0(j,k,:)+ers
       inv_expected(:,j,k)=inv0(:,j,k)-ers
       IF (erodible_deposit_flag) inv_expected(:,j,k)=inv_expected(:,j,k)+dep
       carrier_density=rho_l; carrier_heat=sp_heat_l
       IF (gas_flag) THEN
          carrier_density=pres/(sp_gas_const_a*temperature); carrier_heat=sp_heat_a
       END IF
       outgoing_mass=DOT_PRODUCT(rho_s,dep)+carrier_density*loss
       entering_mass=DOT_PRODUCT(rho_s,ers)+rho_c_sub*coeff_porosity*SUM(ers)
       q_expected(1,j,k)=q0(1,j,k)+entering_mass-outgoing_mass
       q_expected(2:3,j,k)=q0(2:3,j,k)-outgoing_mass*[u,v]
       substrate_temperature=T_erodible
       IF (erodible_deposit_flag) substrate_temperature=temperature
       q_expected(4,j,k)=q0(4,j,k)                                      &
            -temperature*(DOT_PRODUCT(rho_s*sp_heat_s,dep)+carrier_density*carrier_heat*loss) &
            +substrate_temperature*(DOT_PRODUCT(rho_s*sp_heat_s,ers)     &
            +rho_c_sub*carrier_heat*coeff_porosity*SUM(ers))
       q_expected(5:6,j,k)=q0(5:6,j,k)+rho_s*(ers-dep)
       IF (pore_pressure_flag .AND. q0(1,j,k)>0.0_wp) &
            q_expected(idx_pore,j,k)=q0(idx_pore,j,k)+(entering_mass-outgoing_mass)*q0(idx_pore,j,k)/q0(1,j,k)
       cell_delta(j,k)=SUM(dep-ers)/(1.0_wp-erodible_porosity)
    END DO
    DO k=1,ny+1
       DO j=1,nx+1
          j0=MAX(1,j-1); j1=MIN(nx,j); k0=MAX(1,k-1); k1=MIN(ny,k)
          vertex_delta(j,k)=SUM(cell_delta(j0:j1,k0:k1))/REAL((j1-j0+1)*(k1-k0+1),wp)
       END DO
    END DO
    CALL apply_mass_exchange_step(state%q,state%qp,dt,domain,source_registry)
    CALL assert_small(name//' conservative transaction',MAXVAL(ABS(state%q-q_expected)),MAXVAL(ABS(q0)))
    IF (gas_flag) THEN
       ! Compare each equation separately: thermal-energy magnitudes must not
       ! hide a momentum or carrier error. Dry/inactive states are included.
       DO i=1,nq
          CALL assert_small(name//' individual equation',MAXVAL(ABS(state%q(i,:,:)-q_expected(i,:,:))), &
               MAXVAL(ABS(q0(i,:,:))))
       END DO
       IF (MINVAL(state%q(1,:,:)-SUM(state%q(5:6,:,:),DIM=1))<-tol) ERROR STOP 'negative gas carrier'
    END IF
    CALL assert_small(name//' deposit',MAXVAL(ABS(deposit-dep_expected)),1.0_wp)
    CALL assert_small(name//' erosion',MAXVAL(ABS(erosion-ers_expected)),1.0_wp)
    CALL assert_small(name//' erodible',MAXVAL(ABS(erodible-inv_expected)),1.0_wp)
    CALL assert_small(name//' read-only substrate temperature',ABS(T_erodible-317.0_wp),317.0_wp)
    IF (MINVAL(erodible)<-tol .OR. MINVAL(state%q(5:6,:,:))<-tol) ERROR STOP 'negative inventory'
    CALL assert_small(name//' limited nodal proposal',MAXVAL(ABS(B_vertex-bed0-vertex_delta)),MAXVAL(ABS(bed0)))
    ! Sum the solid flow budget and the actual accumulated inventory exchange.
    DO j=1,ns
       CALL assert_small(name//' solid budget',                          &
            ABS(SUM(state%q(4+j,:,:)-q0(4+j,:,:))+rho_s(j)*              &
            SUM(deposit(:,:,j)-dep0(:,:,j)-erosion(:,:,j)+ers0(:,:,j))), &
            SUM(ABS(q0(4+j,:,:))))
    END DO
    nodal_volume=0.0_wp
    DO k=1,ny+1
       DO j=1,nx+1
          weight=dx*dy
          IF (j==1 .OR. j==nx+1) weight=0.5_wp*weight
          IF (k==1 .OR. k==ny+1) weight=0.5_wp*weight
          nodal_volume=nodal_volume+weight*(B_vertex(j,k)-bed0(j,k))
       END DO
    END DO
    CALL assert_small(name//' Q1 weighted volume',ABS(nodal_volume-dx*dy*SUM(cell_delta)), &
         MAX(1.0_wp,dx*dy*SUM(ABS(cell_delta))))
    DO k=1,ny
       DO j=1,nx
          geom_delta(j,k)=0.25_wp*SUM(vertex_delta(j:j+1,k:k+1))
       END DO
    END DO
    mismatch=geom_delta-cell_delta
    IF (SUM(ABS(cell_delta))>tol) THEN
       IF (ABS(mismatch(5,4))<=tol .OR. ABS(mismatch(nx,ny))<=tol) &
            ERROR STOP 'missing dry/inactive geometric mismatch'
    END IF
    diagnostics=[dx*dy*SUM(cell_delta),dx*dy*SUM(geom_delta),dx*dy*SUM(mismatch), &
         dx*dy*SUM(ABS(mismatch)),MAXVAL(ABS(mismatch))]/dt
    WRITE(*,'(A,1X,A,5(1X,ES24.16))') 'N8_EXPECTED_DIAGNOSTICS',name,diagnostics
    CALL check_geometry
    WRITE(snapshot_unit) state%q,state%qp,B_vertex,B_cent,B_face_x,B_face_y, &
         B_prime_x_geom,B_prime_y_geom,B_second_xx_geom,B_second_xy_geom,   &
         B_second_yy_geom,grav_coeff,grav_coeff_stag_x,grav_coeff_stag_y,   &
         deposit,erosion,erodible
    WRITE(*,*) 'PASS: ',name
  END SUBROUTINE check_update

  !> \brief Independently evaluate Q1 identities and the retained five-point LS derivative filter.
  SUBROUTINE check_geometry
    REAL(wp), PARAMETER :: c1(5)=[-2.0_wp,-1.0_wp,0.0_wp,1.0_wp,2.0_wp]
    REAL(wp), PARAMETER :: c2(5)=[2.0_wp,-1.0_wp,-2.0_wp,-1.0_wp,2.0_wp]
    REAL(wp) :: centers(nx,ny), bx,by,bxx,byy,bxy,G
    INTEGER :: j,k,jc,kc,a,b
    DO k=1,ny
       DO j=1,nx
          centers(j,k)=0.25_wp*SUM(B_vertex(j:j+1,k:k+1))
       END DO
    END DO
    CALL assert_small('fresh Q1 centers',MAXVAL(ABS(B_cent-centers)),MAXVAL(ABS(centers)))
    CALL assert_small('fresh Q1 x faces',MAXVAL(ABS(B_face_x-             &
         0.5_wp*(B_vertex(:,1:ny)+B_vertex(:,2:ny+1)))),MAXVAL(ABS(B_vertex)))
    CALL assert_small('fresh Q1 y faces',MAXVAL(ABS(B_face_y-             &
         0.5_wp*(B_vertex(1:nx,:)+B_vertex(2:nx+1,:)))),MAXVAL(ABS(B_vertex)))
    DO k=1,ny
       kc=MAX(3,MIN(ny-2,k))
       DO j=1,nx
          jc=MAX(3,MIN(nx-2,j))
          bx=DOT_PRODUCT(c1,centers(jc-2:jc+2,kc))/(10.0_wp*dx)
          by=DOT_PRODUCT(c1,centers(jc,kc-2:kc+2))/(10.0_wp*dy)
          bxx=DOT_PRODUCT(c2,centers(jc-2:jc+2,kc))/(7.0_wp*dx*dx)
          byy=DOT_PRODUCT(c2,centers(jc,kc-2:kc+2))/(7.0_wp*dy*dy)
          bxy=0.0_wp
          DO b=1,5
             DO a=1,5
                bxy=bxy+c1(a)*c1(b)*centers(jc+a-3,kc+b-3)
             END DO
          END DO
          bxy=bxy/(100.0_wp*dx*dy)
          G=1.0_wp/(1.0_wp+bx*bx+by*by)
          CALL assert_small('fresh slope x',ABS(B_prime_x_geom(j,k)-bx),1.0_wp)
          CALL assert_small('fresh slope y',ABS(B_prime_y_geom(j,k)-by),1.0_wp)
          CALL assert_small('fresh curvature xx',ABS(B_second_xx_geom(j,k)-bxx),1.0_wp)
          CALL assert_small('fresh curvature yy',ABS(B_second_yy_geom(j,k)-byy),1.0_wp)
          CALL assert_small('fresh curvature xy',ABS(B_second_xy_geom(j,k)-bxy),1.0_wp)
          CALL assert_small('fresh center G',ABS(grav_coeff(j,k)-G),1.0_wp)
       END DO
    END DO
    CALL assert_small('fresh interior x-face G',MAXVAL(ABS(grav_coeff_stag_x(2:nx,:)- &
         0.5_wp*(grav_coeff(1:nx-1,:)+grav_coeff(2:nx,:)))),1.0_wp)
    CALL assert_small('fresh interior y-face G',MAXVAL(ABS(grav_coeff_stag_y(:,2:ny)- &
         0.5_wp*(grav_coeff(:,1:ny-1)+grav_coeff(:,2:ny)))),1.0_wp)
    CALL assert_small('fresh boundary x-face G',MAXVAL(ABS(grav_coeff_stag_x(1,:)-grav_coeff(1,:))),1.0_wp)
    CALL assert_small('fresh boundary y-face G',MAXVAL(ABS(grav_coeff_stag_y(:,1)-grav_coeff(:,1))),1.0_wp)
    CALL assert_small('fresh east x-face G',MAXVAL(ABS(grav_coeff_stag_x(nx+1,:)-grav_coeff(nx,:))),1.0_wp)
    CALL assert_small('fresh north y-face G',MAXVAL(ABS(grav_coeff_stag_y(:,ny+1)-grav_coeff(:,ny))),1.0_wp)
  END SUBROUTINE check_geometry

  !> \brief Restore a real checkpoint after corrupting caches, then compare a nonzero second update.
  SUBROUTINE run_restart
    REAL(wp) :: q1(nq,nx,ny),bed1(nx+1,ny+1),dep1(nx,ny,ns),ers1(nx,ny,ns),inv1(ns,nx,ny)
    REAL(wp) :: q2(nq,nx,ny),bed2(nx+1,ny+1),dep2(nx,ny,ns),ers2(nx,ny,ns),inv2(ns,nx,ny)
    CALL reset_fixture
    settling_flag=.TRUE.; erosion_coeff=0.001_wp; erodible_deposit_flag=.TRUE.
    CALL configure_source_registry
    CALL check_update('restart_first',0.2_wp)
    q1=state%q; bed1=B_vertex; dep1=deposit; ers1=erosion; inv1=erodible
    IF (SUM(deposit)+SUM(erosion)<=0.0_wp) ERROR STOP 'restart bed did not evolve'
    runtime%t=0.2_wp; runtime%dt=0.2_wp
    CALL write_restart_file('checkpoint.bin',runtime,stochastic,state,domain)
    CALL check_update('restart_continuous',0.2_wp)
    q2=state%q; bed2=B_vertex; dep2=deposit; ers2=erosion; inv2=erodible
    IF (ALL(bed2==bed1)) ERROR STOP 'restart continuation is a no-op'
    state%q=-999.0_wp; state%qp=-999.0_wp; B_vertex=-999.0_wp
    deposit=-999.0_wp; erosion=-999.0_wp; erodible=-999.0_wp
    B_cent=-999.0_wp; B_face_x=-999.0_wp; B_face_y=-999.0_wp
    B_prime_x_geom=-999.0_wp; B_prime_y_geom=-999.0_wp
    B_second_xx_geom=-999.0_wp; B_second_xy_geom=-999.0_wp; B_second_yy_geom=-999.0_wp
    grav_coeff=-999.0_wp; grav_coeff_stag_x=-999.0_wp; grav_coeff_stag_y=-999.0_wp
    CALL read_restart_file('checkpoint.bin',runtime,stochastic,state,domain)
    CALL same_bits('restart q',RESHAPE(state%q,[SIZE(q1)]),RESHAPE(q1,[SIZE(q1)]))
    CALL same_bits('restart B_vertex',RESHAPE(B_vertex,[SIZE(bed1)]),RESHAPE(bed1,[SIZE(bed1)]))
    CALL same_bits('restart deposit',RESHAPE(deposit,[SIZE(dep1)]),RESHAPE(dep1,[SIZE(dep1)]))
    CALL same_bits('restart erosion',RESHAPE(erosion,[SIZE(ers1)]),RESHAPE(ers1,[SIZE(ers1)]))
    CALL same_bits('restart erodible',RESHAPE(erodible,[SIZE(inv1)]),RESHAPE(inv1,[SIZE(inv1)]))
    CALL check_geometry
    ! qp is intentionally not checkpointed; the main program normally rebuilds
    ! it from q. Reset its inactive cache before the continuation fingerprint.
    state%qp=0.0_wp; state%qp(4,:,:)=T_ambient
    CALL check_update('restart_resumed',0.2_wp)
    CALL same_bits('continued q',RESHAPE(state%q,[SIZE(q2)]),RESHAPE(q2,[SIZE(q2)]))
    CALL same_bits('continued bed',RESHAPE(B_vertex,[SIZE(bed2)]),RESHAPE(bed2,[SIZE(bed2)]))
    CALL same_bits('continued deposit',RESHAPE(deposit,[SIZE(dep2)]),RESHAPE(dep2,[SIZE(dep2)]))
    CALL same_bits('continued erosion',RESHAPE(erosion,[SIZE(ers2)]),RESHAPE(ers2,[SIZE(ers2)]))
    CALL same_bits('continued inventory',RESHAPE(erodible,[SIZE(inv2)]),RESHAPE(inv2,[SIZE(inv2)]))
  END SUBROUTINE run_restart

  !> \brief Reject nonfinite values and scaled errors beyond the specified roundoff tolerance.
  !> \param[in] label Assertion name.
  !> \param[in] error Absolute error.
  !> \param[in] scale Reference magnitude for the tolerance.
  SUBROUTINE assert_small(label,error,scale)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: error,scale
    IF (.NOT.ieee_is_finite(error) .OR. error>tol*MAX(1.0_wp,scale)) THEN
       WRITE(*,*) 'FAIL: ',label,error,tol*MAX(1.0_wp,scale)
       ERROR STOP 1
    END IF
  END SUBROUTINE assert_small

  !> \brief Compare serialized REAL fields including their exact floating-point bit patterns.
  !> \param[in] label Assertion name.
  !> \param[in] actual,reference Flattened REAL fields to compare.
  SUBROUTINE same_bits(label,actual,reference)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: actual(:),reference(:)
    IF (ANY(TRANSFER(actual,[0_int64],SIZE(actual)) /=                   &
         TRANSFER(reference,[0_int64],SIZE(reference)))) THEN
       WRITE(*,*) 'FAIL: ',label
       ERROR STOP 1
    END IF
  END SUBROUTINE same_bits

END PROGRAM test_mass_exchange
