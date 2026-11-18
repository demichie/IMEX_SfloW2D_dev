!> \brief Exercise the production spatial operator, active worksets and CFL API without time stepping.
PROGRAM test_spatial_operator
  USE thermal_parameters_2d, ONLY : thermal_l1
  USE rheology_parameters_2d, ONLY : rheology_l1
  USE parameters_2d
  USE substrate_parameters_2d, ONLY : substrate,substrate_config_type,erodible_fract
  USE mass_exchange_parameters_2d, ONLY : mass_exchange_l1,mass_exchange_config_type
  USE pore_pressure_parameters_2d, ONLY : pore_pressure_l1
  USE constitutive_parameters_2d
  USE solid_transport_parameters_2d
  USE geometry_2d
  USE domain_2d, ONLY : domain_type
  USE spatial_operator_2d, ONLY : spatial_operator_type
  USE reconstruction_2d, ONLY : reconstruction_workspace_type
  USE state_conversions, ONLY : primitive_to_conservative
  USE state_fields, ONLY : state_type
  USE init_2d, ONLY : init_empty
  USE hyperbolic_2d, ONLY : hyperbolic_workspace_type, add_interlayer_terms, &
      eval_total_depth_face_speed_bounds
  USE pccu_2d, ONLY : eval_pccu_advective_flux, eval_central_upwind_flux, eval_oriented_pccu_pair, &
      eval_interlayer_pressure_path, eval_interlayer_pressure_correction, eval_hydrostatic_path, &
      eval_interlayer_oriented_face
  USE source_terms_2d, ONLY : eval_source_explicit
  USE equation_terms_2d, ONLY : eval_hydrostatic_coefficient, eval_local_speeds_direction
  USE state_composition, ONLY : conservative_mass_admissible
  USE mass_injection_coefficients_2d, ONLY : injection_material_set_type, prepare_injection_materials
  USE source_catalog_2d, ONLY : SOURCE_METHOD_EXPLICIT, SOURCE_ID_GEOMETRY_CURVATURE, &
      SOURCE_ID_EXCHANGE_DEPOSITION,SOURCE_ID_EXCHANGE_EROSION
  USE source_registry_2d, ONLY : source_group_type
  USE time_integration_2d, ONLY : time_integration_workspace_type
  USE source_registry_2d, ONLY : source_registry_type, &
      SOURCE_SCOPE_L1,SOURCE_SCOPE_L2,SOURCE_SCOPE_INTERLAYER
  USE source_catalog_2d, ONLY : SOURCE_METHOD_ANALYTIC
  USE layer_physics_2d, ONLY : layer_physics_type
  USE thermal_parameters_2d, ONLY : thermal_l2
  USE rheology_parameters_2d, ONLY : rheology_l2
  USE mass_exchange_parameters_2d, ONLY : mass_exchange_l2
  USE pore_pressure_parameters_2d, ONLY : pore_pressure_l2,pore_pressure_cache_l1,pore_pressure_cache_l2
  USE interlayer_mass_exchange_parameters_2d, ONLY : interlayer_exchange,interlayer_mass_exchange_config_type
  USE interlayer_mass_exchange_2d, ONLY : apply_interlayer_lie_group
  USE state_types, ONLY : physical_state_real,initialize_physical_state_real
  USE interlayer_momentum_exchange_2d, ONLY : interlayer_momentum_exchange, &
      interlayer_momentum_exchange_config_type,apply_interlayer_momentum_exchange_analytic
  USE stochastic_module, ONLY : stochastic_workspace_type,stochastic_l1,stochastic_l2
  USE omp_lib
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
  USE, INTRINSIC :: iso_fortran_env, ONLY : int64
  IMPLICIT NONE

  TYPE(domain_type) :: domain
    TYPE(domain_type) :: slot2_domain !< Independent L2 active-cell and face worksets.
  TYPE(spatial_operator_type) :: spatial
    TYPE(spatial_operator_type) :: slot2_spatial !< Independent inactive-slot HP-PCCU scratch.
  TYPE(reconstruction_workspace_type) :: reconstruction
  TYPE(hyperbolic_workspace_type) :: hyperbolic
  TYPE(source_group_type) :: explicit_group
  TYPE(injection_material_set_type) :: materials
  REAL(wp), ALLOCATABLE :: q_faceL(:,:,:), q_faceR(:,:,:), q_faceB(:,:,:), q_faceT(:,:,:)
  REAL(wp), ALLOCATABLE :: q(:,:,:), q0(:,:,:), qp(:,:,:), term(:,:,:), direct(:,:,:), cell_sources(:,:,:)
    REAL(wp), ALLOCATABLE :: slot2_qp(:,:,:),slot2_term(:,:,:) !< Independent slot buffers.
    REAL(wp), ALLOCATABLE :: slot_q(:,:,:,:) !< Test-only two-slot conservative states.
    LOGICAL, ALLOCATABLE :: l1_solve_mask_before(:,:) !< Check L2 refresh does not modify L1.
  REAL(wp) :: dt, repeated_dt, saved_cfl, capped_dt, no_fissures(0),slot2_dt,coupled_dt
  !> Independent rate derived from the public backend face payload [1/s].
  REAL(wp) :: lambda_x, lambda_y, expected_dt, combined_dt, saved_max_dt, rate
  CHARACTER(LEN=32) :: argument
  INTEGER :: nx,ny,slope,curvature,limiter_id,input_unit,output_unit,j,k,l,threads,actual_threads

  CALL get_command_argument(1,argument)
  READ(argument,*) threads
  actual_threads=0
  !$OMP PARALLEL
  !$OMP SINGLE
  actual_threads=omp_get_num_threads()
  !$OMP END SINGLE
  !$OMP END PARALLEL
  IF (actual_threads/=threads) ERROR STOP 'incorrect actual OpenMP team'
  WRITE(*,*) 'Actual OpenMP team:',actual_threads

  OPEN(NEWUNIT=input_unit,FILE='fixture.inp',STATUS='OLD',ACTION='READ')
  READ(input_unit,*) nx,ny,slope,curvature,limiter_id
  READ(input_unit,*) dx,dy
  CALL initialize_fixture
  CALL check_uniform_minmod_kernel
  CALL check_optional_pccu_flux_slots
  CALL check_interlayer_pressure_path
  CALL check_interlayer_metric_no_double_count
  CALL check_interlayer_variable_gamma_overlap
  CALL check_interlayer_oriented_face_pair
  CALL check_interlayer_coupled_face
  CALL check_degenerate_interlayer_face
  CALL check_interlayer_cell_residual
  explicit_group%method=SOURCE_METHOD_EXPLICIT
  explicit_group%source_names=[CHARACTER(LEN=64)::]
  IF (curvature_term_flag) THEN
       explicit_group%source_names=[CHARACTER(LEN=64)::'geometry_curvature']
       explicit_group%source_mask=IBSET(0,SOURCE_ID_GEOMETRY_CURVATURE-1)
  END IF
  READ(input_unit,*) B_vertex
  READ(input_unit,*) q
  CLOSE(input_unit)
  q0=q
  CALL refresh_topography_geometry
  CALL check_interlayer_face_residual
  CALL check_interlayer_cell_assembly
  CALL check_two_layer_union_workset
  CALL check_dry_dirichlet_inlet
  CALL domain%check_solve(q,0.0_wp,.FALSE.)
  CALL check_coupled_cfl_policy
  IF (domain%solve_cells<=0 .OR. domain%reconstruction_cells<domain%solve_cells) &
       ERROR STOP 'invalid production worksets'

  ! Only the read halo is deliberately poisoned. The timestep API must refresh
  ! wet states and dry halo entries; unrelated inactive caches are not inputs.
  qp=0.0_wp; qp(4,:,:)=T_ambient
  DO k=1,ny
     DO j=1,nx
        IF (domain%reconstruction_mask(j,k)) qp(:,j,k)=-999.0_wp
     END DO
  END DO
  CALL spatial%compute_timestep(q,qp,0.0_wp,dt,domain)
  CALL check_total_depth_cfl_equivalence
  IF (.NOT.ieee_is_finite(dt) .OR. dt<=0.0_wp) ERROR STOP 'invalid production timestep'
  IF (ANY(qp(1,:,:) < 0.0_wp)) ERROR STOP 'physical cache was not refreshed'
  CALL same_bits('timestep input is read-only',RESHAPE(q,[SIZE(q)]),RESHAPE(q0,[SIZE(q0)]))
  CALL spatial%evaluate(qp,term,0.0_wp,domain)
  CALL hyperbolic%evaluate_terms(reconstruction,qp,direct,0.0_wp,          &
       domain%solve_cells,domain%j_cent,domain%k_cent,                    &
       domain%solve_interfaces_x,domain%j_stag_x,domain%k_stag_x,          &
       domain%solve_interfaces_y,domain%j_stag_y,domain%k_stag_y)
  ! The API defines the returned term only on solve_cells. Initialize the
  ! remaining diagnostic payload explicitly, rather than reading undefined
  ! inactive entries from an INTENT(OUT) array.
  DO k=1,ny
     DO j=1,nx
        IF (.NOT.domain%solve_mask(j,k)) THEN
           term(:,j,k)=0.0_wp; direct(:,j,k)=0.0_wp
        END IF
     END DO
  END DO
  CALL same_bits('public operator delegates complete spatial term',     &
       RESHAPE(term,[SIZE(term)]),RESHAPE(direct,[SIZE(direct)]))

    ! The two independently owned HP-PCCU workspaces must reproduce the
    ! same single-layer flux and CFL for identical states. Neither one
    ! borrows the other's face traces, reconstruction or speed arrays.
    ALLOCATE(slot2_qp(n_vars,nx,ny),slot2_term(n_eqns,nx,ny))
    ALLOCATE(slot_q(n_vars,2,nx,ny))
    slot_q(:,1,:,:)=q
    slot_q(:,2,:,:)=q
    l1_solve_mask_before=domain%solve_mask
    CALL slot2_domain%initialize
    CALL slot2_domain%check_solve(slot_q(:,2,:,:),0.0_wp,.FALSE.)
    IF (slot2_domain%solve_cells/=domain%solve_cells) ERROR STOP 'identical layer worksets differ'
    IF (ANY(slot2_domain%solve_mask .NEQV. domain%solve_mask)) ERROR STOP 'L2 mask differs for same state'
    slot2_qp=qp
    CALL slot2_spatial%initialize(2)
    CALL slot2_spatial%compute_timestep(slot_q(:,2,:,:),slot2_qp,0.0_wp,slot2_dt,slot2_domain)
    CALL same_bits('slot-local CFL matches L1',[slot2_dt],[dt])

    ! Exercise the model-facing optional L2 geometry path without enabling
    ! coupled dynamics. L1 was prepared by the production CFL pass above.
    CALL slot2_spatial%compute_timestep(slot_q(:,2,:,:),slot2_qp,0.0_wp,slot2_dt,slot2_domain, &
        lower_operator=spatial,h_lower=qp(1,:,:))
    IF (.NOT.ieee_is_finite(slot2_dt) .OR. slot2_dt<=0.0_wp) &
        ERROR STOP 'Invalid L2 reconstructed isolated CFL'

    ! Complete coupled spatial residual and CFL on the SAME prepared workset.
    ! The public L1 and L2 entry points are unchanged.
    CALL spatial%evaluate_coupled(slot2_spatial,qp,slot2_qp,direct,slot2_term, &
        0.0_wp,domain,coupled_dt)
    IF (.NOT.ieee_is_finite(coupled_dt) .OR. coupled_dt<=0.0_wp) &
        ERROR STOP 'Invalid total-depth coupled CFL'
    CALL check_coupled_spatial_trial
    DO l=1,domain%solve_cells
        j=domain%j_cent(l)
        k=domain%k_cent(l)
        IF (.NOT.ALL(ieee_is_finite(direct(:,j,k))) .OR. &
            .NOT.ALL(ieee_is_finite(slot2_term(:,j,k)))) &
            ERROR STOP 'Nonfinite coupled HP-PCCU spatial residual'
    END DO

    ! Restore independent-layer traces for the existing bitwise regression.
    CALL slot2_spatial%compute_timestep(slot_q(:,2,:,:),slot2_qp,0.0_wp,slot2_dt,slot2_domain)
    CALL same_bits('slot-local CFL survives optional L2 geometry',[slot2_dt],[dt])
    CALL slot2_spatial%evaluate(slot2_qp,slot2_term,0.0_wp,slot2_domain)

    DO k=1,ny
        DO j=1,nx
            IF (.NOT.slot2_domain%solve_mask(j,k)) slot2_term(:,j,k)=0.0_wp
        END DO
    END DO

    CALL same_bits('slot-local HP-PCCU preserves L1 result', &
        RESHAPE(slot2_term,[SIZE(slot2_term)]),RESHAPE(term,[SIZE(term)]))
    ! Only slot 2 changes wet/dry state; the already built L1 lists remain intact.
    slot_q(:,2,2:nx-1,2:ny-1)=0.0_wp
    CALL slot2_domain%check_solve(slot_q(:,2,:,:),0.0_wp,.FALSE.)
    IF (ANY(domain%solve_mask .NEQV. l1_solve_mask_before)) &
        ERROR STOP 'L2 wet/dry domain rebuild modified L1 solve mask'
    CALL same_bits('L2 wet/dry domain rebuild did not change L1 state', &
        RESHAPE(slot_q(:,1,:,:),[SIZE(q)]),RESHAPE(q,[SIZE(q)]))
    IF (slot2_domain%reconstruction_cells<slot2_domain%solve_cells) &
        ERROR STOP 'Invalid independent L2 reconstruction workset'

    ! Exercise the existing isolated single-layer kernel on a different dry/wet
    ! slot. This is an operator test, not physical two-layer coupling.
    slot2_qp=0.0_wp
    slot2_qp(4,:,:)=T_ambient
    CALL slot2_spatial%compute_timestep(slot_q(:,2,:,:),slot2_qp,0.0_wp,slot2_dt,slot2_domain)
    IF (.NOT.ieee_is_finite(slot2_dt) .OR. slot2_dt<=0.0_wp) &
        ERROR STOP 'Invalid independent L2 dry/wet CFL'
    CALL slot2_spatial%evaluate(slot2_qp,slot2_term,0.0_wp,slot2_domain)
    DO k=1,ny
        DO j=1,nx
            IF (.NOT.slot2_domain%solve_mask(j,k)) CYCLE
            IF (.NOT.ALL(ieee_is_finite(slot2_term(:,j,k)))) &
                ERROR STOP 'Non-finite independent L2 dry/wet spatial term'
        END DO
    END DO

    slot2_domain%solve_mask_time(1,1)=12.5_wp
    IF (domain%solve_mask_time(1,1)==slot2_domain%solve_mask_time(1,1)) &
        ERROR STOP 'L2 domain workset aliases the L1 domain'
    CALL check_isolated_gas_l2
    CALL slot2_spatial%finalize
    CALL slot2_domain%finalize
    DEALLOCATE(slot2_qp,slot2_term,slot_q,l1_solve_mask_before)

  IF (.NOT.ALL(ieee_is_finite(term))) ERROR STOP 'nonfinite spatial term'
  CALL prepare_injection_materials(materials)
  cell_sources=0.0_wp
  DO k=1,ny
     DO j=1,nx
        IF (.NOT.domain%solve_mask(j,k)) CYCLE
        CALL eval_source_explicit(B_prime_x_geom(j,k),B_prime_y_geom(j,k),    &
             B_second_xx_geom(j,k),B_second_xy_geom(j,k),              &
             B_second_yy_geom(j,k),grav_coeff(j,k),qp(:,j,k),           &
             cell_sources(:,j,k),0.0_wp,0.0_wp,no_fissures,            &
             0.0_wp,0.0_wp,0.0_wp,dx*dy,explicit_group,0.0_wp,materials=materials,thermal=thermal_l1, &
             rheology=rheology_l1, &
             pore_enabled=pore_pressure_l1%pore_pressure_flag,stochastic_enabled=rheology_l1%stochastic_rheology_flag)
     END DO
  END DO
  IF (.NOT.ALL(ieee_is_finite(cell_sources))) ERROR STOP 'nonfinite local source'
  CALL spatial%compute_timestep(q,qp,0.0_wp,repeated_dt,domain)
  CALL same_bits('CFL refresh/evaluation is repeatable',[dt],[repeated_dt])

  ! Test the existing max_dt cap and prescribed-step sentinel independently
  ! of the reference's explicit SSPRK2 timestep combination policy.
  max_dt=0.5_wp*dt
  CALL spatial%compute_timestep(q,qp,0.0_wp,capped_dt,domain)
  CALL same_bits('max_dt cap',[capped_dt],[max_dt])
  saved_cfl=cfl; cfl=-1.0_wp
  CALL spatial%compute_timestep(q,qp,0.0_wp,capped_dt,domain)
  CALL same_bits('prescribed timestep sentinel',[capped_dt],[max_dt])
  cfl=saved_cfl

  ! Both policies consume the same face speeds. Check their algebra separately
  ! from the pinned reference, including collapsed axes, caps and no propagation.
  saved_max_dt=max_dt
  max_dt=1.0E6_wp
  lambda_x=0.0_wp; lambda_y=0.0_wp
  DO l=1,domain%solve_cells
     j=domain%j_cent(l); k=domain%k_cent(l)
     lambda_x=MAX(lambda_x,hyperbolic%a_interface_xPos(j,k), &
          -hyperbolic%a_interface_xNeg(j,k),hyperbolic%a_interface_xPos(j+1,k), &
          -hyperbolic%a_interface_xNeg(j+1,k))
     lambda_y=MAX(lambda_y,hyperbolic%a_interface_yPos(j,k), &
          -hyperbolic%a_interface_yNeg(j,k),hyperbolic%a_interface_yPos(j,k+1), &
          -hyperbolic%a_interface_yNeg(j,k+1))
  END DO
  rate=0.0_wp
  IF (nx>1) rate=lambda_x/dx
  IF (ny>1) rate=rate+lambda_y/dy
  expected_dt=max_dt
  IF (rate>0.0_wp) expected_dt=MIN(max_dt,cfl/rate)
  IF (nx==1 .OR. ny==1) expected_dt=dt
  cfl_policy='combined'
  CALL spatial%compute_timestep(q,qp,0.0_wp,combined_dt,domain)
  CALL same_bits('combined CFL algebra',[combined_dt],[expected_dt])
  IF (nx>1 .AND. ny>1) THEN
     IF (combined_dt>dt) ERROR STOP 'combined step exceeds directional minimum'
  ELSE
     CALL same_bits('one-dimensional policies agree',[combined_dt],[dt])
  END IF
  max_dt=0.5_wp*combined_dt
  CALL spatial%compute_timestep(q,qp,0.0_wp,capped_dt,domain)
  CALL same_bits('combined max_dt cap',[capped_dt],[max_dt])
  cfl=-1.0_wp
  qp=-999.0_wp
  CALL spatial%compute_timestep(q,qp,0.0_wp,capped_dt,domain)
  CALL same_bits('combined prescribed timestep sentinel',[capped_dt],[max_dt])
  IF (ANY(qp/=-999.0_wp)) ERROR STOP 'sentinel unexpectedly refreshed face states'
  cfl=saved_cfl; cfl_policy='directional_min'; max_dt=1.0E6_wp
  CALL spatial%compute_timestep(q,qp,0.0_wp,repeated_dt,domain)
  CALL same_bits('policy roundtrip preserves default',[dt],[repeated_dt])
  ! A wholly dry workset has no propagation in either direction. Empty lists
  ! must return the cap, not divide by zero, independent of policy.
  q=0.0_wp
  CALL domain%check_solve(q,0.0_wp,.FALSE.)
  DO j=1,2
     IF (j==1) cfl_policy='directional_min'
     IF (j==2) cfl_policy='combined'
     CALL spatial%compute_timestep(q,qp,0.0_wp,capped_dt,domain)
     CALL same_bits('empty propagation returns max_dt',[capped_dt],[max_dt])
  END DO
  q=q0; cfl_policy='directional_min'; max_dt=saved_max_dt
  CALL domain%check_solve(q,0.0_wp,.FALSE.)
  CALL spatial%compute_timestep(q,qp,0.0_wp,repeated_dt,domain)

  OPEN(NEWUNIT=output_unit,FILE='cfl_policies.bin',ACCESS='STREAM', &
       FORM='UNFORMATTED',STATUS='REPLACE')
  WRITE(output_unit) dt,combined_dt
  CLOSE(output_unit)

  CALL normalize_inactive_face_payloads
  ! Preserve the historic test payload layout, without obsolete production arrays.
  ALLOCATE(q_faceL(n_vars,nx+1,ny),q_faceR(n_vars,nx+1,ny))
  ALLOCATE(q_faceB(n_vars,nx,ny+1),q_faceT(n_vars,nx,ny+1))
  q_faceL=0.0_wp; q_faceR=0.0_wp
  q_faceB=0.0_wp; q_faceT=0.0_wp
  DO k=1,ny
    DO j=1,nx+1
      IF (.NOT.domain%solve_mask_x(j,k)) CYCLE
      CALL primitive_to_conservative(reconstruction%qp_adv_interfaceL(:,j,k),q_faceL(:,j,k))
      CALL primitive_to_conservative(reconstruction%qp_adv_interfaceR(:,j,k),q_faceR(:,j,k))
    END DO
  END DO
  DO k=1,ny+1
    DO j=1,nx
      IF (.NOT.domain%solve_mask_y(j,k)) CYCLE
      CALL primitive_to_conservative(reconstruction%qp_adv_interfaceB(:,j,k),q_faceB(:,j,k))
      CALL primitive_to_conservative(reconstruction%qp_adv_interfaceT(:,j,k),q_faceT(:,j,k))
    END DO
  END DO
  OPEN(NEWUNIT=output_unit,FILE='result.bin',ACCESS='STREAM',           &
       FORM='UNFORMATTED',STATUS='REPLACE')
  ! Version-two test-only stream order, decoded by compare_reference.py.
  ! No production output or private spatial-operator interface is changed.
  WRITE(output_unit) dt
  WRITE(output_unit) q,qp,term,cell_sources
  WRITE(output_unit) MERGE(1.0_wp,0.0_wp,domain%solve_mask)
  WRITE(output_unit) B_vertex,B_cent,B_face_x,B_face_y
  WRITE(output_unit) B_prime_x_geom,B_prime_y_geom,B_second_xx_geom, &
       B_second_xy_geom,B_second_yy_geom,grav_coeff,                &
       grav_coeff_stag_x,grav_coeff_stag_y
  WRITE(output_unit) hyperbolic%a_interface_xNeg,hyperbolic%a_interface_xPos, &
       hyperbolic%a_interface_yNeg,hyperbolic%a_interface_yPos
  WRITE(output_unit) q_faceL,q_faceR,q_faceB,q_faceT
  WRITE(output_unit) reconstruction%qp_adv_interfaceL,reconstruction%qp_adv_interfaceR, &
       reconstruction%qp_adv_interfaceB,reconstruction%qp_adv_interfaceT
  WRITE(output_unit) reconstruction%eta_interfaceL,reconstruction%eta_interfaceR, &
       reconstruction%eta_interfaceB,reconstruction%eta_interfaceT
  WRITE(output_unit) reconstruction%h_froude_interfaceL, &
       reconstruction%h_froude_interfaceR, reconstruction%h_froude_interfaceB, &
       reconstruction%h_froude_interfaceT
  CLOSE(output_unit)
  CALL check_stratified_rest_geometry
  CALL check_closed_coupled_spatial_budget
  ! Keep this state-mutating BC fixture separate from the reference payload.
  ! It runs after serialization so it cannot perturb the measured spatial fields.
  CALL check_layer_local_dirichlet_inlet
  WRITE(*,*) 'Production worksets:',domain%solve_cells,domain%reconstruction_cells
  CALL spatial%finalize
  CALL hyperbolic%finalize
  CALL reconstruction%finalize
  CALL domain%finalize
  CALL check_optional_layer_spatial_slots
  WRITE(*,*) 'PASS: production spatial API, active caches, repeatable CFL and caps'

CONTAINS

    !> @brief Validate an Euler trial on the existing two-layer spatial fixture.
    !> @details The fixture has Neumann boundaries, so its global residual sum
    !> need not vanish. A separate compact-support case checks conservation.
    SUBROUTINE check_coupled_spatial_trial
        REAL(wp) :: step,updated(n_eqns)
        INTEGER :: l,j,k

        step=MIN(coupled_dt*0.01_wp,1.0E-5_wp)
        IF (step<=0.0_wp) ERROR STOP 'Coupled spatial trial has invalid timestep'
        DO l=1,domain%solve_cells
            j=domain%j_cent(l)
            k=domain%k_cent(l)
            updated=slot_q(1:n_eqns,1,j,k)-step*direct(:,j,k)
            IF (.NOT.conservative_mass_admissible(updated)) &
                ERROR STOP 'Coupled L1 spatial candidate is inadmissible'
            updated=slot_q(1:n_eqns,2,j,k)-step*slot2_term(:,j,k)
            IF (.NOT.conservative_mass_admissible(updated)) &
                ERROR STOP 'Coupled L2 spatial candidate is inadmissible'
        END DO
        WRITE(*,*) 'PASS: coupled Euler trial admissibility for both layers'
    END SUBROUTINE check_coupled_spatial_trial

    !> @brief Execute the production spatial-Lie-commit path on both layers.
    !> @details A zero-drag reference and an identical Chezy step must have
    !>          the same mass, solids and thermal energy per layer; their
    !>          momentum differences must cancel pairwise. This test bypasses
    !>          only the model-layout production gate, not the timestep code.
    SUBROUTINE check_integrated_coupled_lie_timestep(paired,q1,q2,p1,p2,dt_bound)
        TYPE(domain_type), INTENT(IN) :: paired
        REAL(wp), INTENT(IN) :: q1(n_vars,nx,ny),q2(n_vars,nx,ny)
        REAL(wp), INTENT(IN) :: p1(n_vars,nx,ny),p2(n_vars,nx,ny),dt_bound
        TYPE(spatial_operator_type) :: operators(2)
        TYPE(time_integration_workspace_type) :: workspace
        TYPE(source_registry_type) :: registry
        TYPE(source_group_type) :: groups(4),external_groups(6)
        TYPE(layer_physics_type) :: physics(2)
        TYPE(state_type) :: coupled_state
        TYPE(stochastic_workspace_type) :: stochastic
        TYPE(interlayer_momentum_exchange_config_type) :: old_drag
        TYPE(interlayer_mass_exchange_config_type) :: old_mass
        TYPE(mass_exchange_config_type) :: old_external_l1,old_external_l2
        TYPE(substrate_config_type) :: old_substrate
        TYPE(physical_state_real) :: reference_physical(4)
        REAL(wp) :: base_q(n_vars,2,nx,ny),base_qp(n_vars,2,nx,ny)
        REAL(wp) :: expected_q(n_vars,2),expected_qp(n_vars,2)
        REAL(wp) :: before_q(n_vars,2,nx,ny)
        REAL(wp) :: before_qp(n_vars,2,nx,ny)
        REAL(wp) :: pt1(n_vars,nx,ny),pt2(n_vars,nx,ny)
        REAL(wp) :: bound1,bound2,dt_step,momentum_difference,drag_difference,scale
        REAL(wp) :: mass_difference,component_difference,kin_before,vonK_before,primitive(n_vars)
        REAL(wp) :: curvature_rate(n_eqns),curvature_effect,solid_budget,expected_deposition
        INTEGER :: original_layers,j,k,l,eq,scenario,integrated_unit,layer,group_mode,source_status
        LOGICAL :: created_zone,created_arc,created_diam,created_erodible_fraction,accepted
        LOGICAL :: created_source_fraction,created_fissure_fraction,created_arc_normals
        REAL(wp), ALLOCATABLE :: original_diam(:),original_erodible_fraction(:)
        CHARACTER(LEN=256) :: error

        original_layers=n_layers
        n_layers=2
        old_drag=interlayer_momentum_exchange
        old_mass=interlayer_exchange
        old_external_l1=mass_exchange_l1
        old_external_l2=mass_exchange_l2
        old_substrate=substrate
        DO l=1,4
            CALL initialize_physical_state_real(reference_physical(l))
        END DO
        ! Both layers are configured from the canonical physics objects.
        CALL physics(1)%initialize(.FALSE.,rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1,error)
        IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid L1 physics in integrated Lie fixture'
        CALL physics(2)%initialize(.FALSE.,rheology_l2,thermal_l2,pore_pressure_l2,mass_exchange_l2,error)
        IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid L2 physics in integrated Lie fixture'

        groups(1)%id=0
        groups(1)%scope=SOURCE_SCOPE_L1
        groups(1)%method=SOURCE_METHOD_EXPLICIT
        groups(1)%source_names=[CHARACTER(LEN=64)::]
        groups(1)%equation_ids=[(eq,eq=1,n_eqns)]
        groups(2)=groups(1)
        groups(2)%scope=SOURCE_SCOPE_L2
        groups(3)%id=4
        groups(3)%scope=SOURCE_SCOPE_INTERLAYER
        groups(3)%method=SOURCE_METHOD_EXPLICIT
        groups(3)%source_names=[CHARACTER(LEN=64)::'interlayer_mass_exchange']
        groups(3)%equation_ids=[(eq,eq=1,n_eqns)]
        groups(4)%id=5
        groups(4)%scope=SOURCE_SCOPE_INTERLAYER
        groups(4)%method=SOURCE_METHOD_ANALYTIC
        groups(4)%source_names=[CHARACTER(LEN=64)::'interlayer_drag']
        groups(4)%equation_ids=[INTEGER::]
        CALL registry%initialize(n_eqns,0,groups,error)
        IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid integrated Lie group registry'

        ALLOCATE(coupled_state%q(n_vars,2,nx,ny),coupled_state%qp(n_vars,2,nx,ny))
        ALLOCATE(coupled_state%deposit(nx,ny,n_solid),coupled_state%erosion(nx,ny,n_solid))
        ALLOCATE(coupled_state%erodible(n_solid,nx,ny))
        ALLOCATE(coupled_state%deposit_tot(nx,ny),coupled_state%erosion_tot(nx,ny))
        coupled_state%deposit=0.0_wp
        coupled_state%erosion=0.0_wp
        coupled_state%erodible=0.0_wp
        coupled_state%deposit_tot=0.0_wp
        coupled_state%erosion_tot=0.0_wp

        created_zone=.NOT.ALLOCATED(B_zone)
        IF (created_zone) ALLOCATE(B_zone(nx,ny))
        B_zone=0
        created_arc=.NOT.ALLOCATED(cell_arc_perim)
        IF (created_arc) ALLOCATE(cell_arc_perim(nx,ny))
        cell_arc_perim=0.0_wp
        ! The full source executor consumes these geometry inputs even when
        ! the configured external process has no injection or fissure source.
        created_source_fraction=.NOT.ALLOCATED(cell_source_fractions)
        IF (created_source_fraction) ALLOCATE(cell_source_fractions(nx,ny))
        cell_source_fractions=0.0_wp
        created_fissure_fraction=.NOT.ALLOCATED(cell_fissure_fractions)
        IF (created_fissure_fraction) ALLOCATE(cell_fissure_fractions(nx,ny,n_fissures))
        cell_fissure_fractions=0.0_wp
        created_arc_normals=.NOT.ALLOCATED(cell_arc_n_x)
        IF (created_arc_normals) THEN
            ALLOCATE(cell_arc_n_x(nx,ny),cell_arc_n_y(nx,ny))
        END IF
        cell_arc_n_x=0.0_wp
        cell_arc_n_y=0.0_wp

        ! The existing single-layer OU field is merely a zero-valued input;
        ! no stochastic source is active in either layer in this fixture.
        stochastic_l1%stochastic_transport_flag=.FALSE.
        stochastic_l2%stochastic_transport_flag=.FALSE.
        stochastic%config=>stochastic_l1
        stochastic%rheology=>rheology_l1
        ALLOCATE(stochastic%effective_Z(nx,ny))
        stochastic%effective_Z=0.0_wp

        CALL operators(1)%initialize(1)
        CALL operators(2)%initialize(2)
        CALL operators(1)%prepare_hp_workset(paired)
        CALL operators(2)%prepare_hp_workset(paired)
        CALL workspace%initialize(registry)
        ! The drag-only reference must remain a no-op for mass transfer.
        interlayer_exchange%settling_flag=.FALSE.
        interlayer_exchange%elutriation_flag=.FALSE.
        interlayer_momentum_exchange%model='CHEZY'
        interlayer_momentum_exchange%method='ANALYTIC'
        interlayer_momentum_exchange%group=5

        OPEN(NEWUNIT=integrated_unit,FILE='integrated_lie.bin',ACCESS='STREAM', &
            FORM='UNFORMATTED',STATUS='REPLACE')
        ! The original union workset safely covers each subset of wet cells.
        ! Test both wet, L1 only, L2 only, and both dry without rebuilding
        ! the source schedule or replacing the production timestep.
        DO scenario=1,4
            before_q(:,1,:,:)=q1
            before_q(:,2,:,:)=q2
            before_qp(:,1,:,:)=p1
            before_qp(:,2,:,:)=p2
            IF (scenario==2 .OR. scenario==4) THEN
                before_q(:,2,:,:)=0.0_wp
                before_qp(:,2,:,:)=0.0_wp
                before_qp(idx_T,2,:,:)=T_ambient
            END IF
            IF (scenario==3 .OR. scenario==4) THEN
                before_q(:,1,:,:)=0.0_wp
                before_qp(:,1,:,:)=0.0_wp
                before_qp(idx_T,1,:,:)=T_ambient
            END IF
            pt1=before_qp(:,1,:,:)
            pt2=before_qp(:,2,:,:)
            CALL operators(1)%compute_timestep(before_q(:,1,:,:),pt1,0.0_wp,bound1,paired)
            CALL operators(2)%compute_timestep(before_q(:,2,:,:),pt2,0.0_wp,bound2,paired, &
                lower_operator=operators(1),h_lower=pt1(idx_h,:,:))
            dt_step=MIN(1.0E-5_wp,0.01_wp*MIN(dt_bound,bound1,bound2))
            IF (dt_step<=0.0_wp) ERROR STOP 'Invalid integrated two-layer timestep'
            before_qp(:,1,:,:)=pt1
            before_qp(:,2,:,:)=pt2

            interlayer_momentum_exchange%chezy_drag_coeff=0.0_wp
            coupled_state%q=before_q
            coupled_state%qp=before_qp
            CALL workspace%advance(coupled_state,0.0_wp,dt_step,stochastic,paired,operators, &
                pre_spatial_primitives_ready=.TRUE.,layer_physics=physics, &
                pore_caches=[pore_pressure_cache_l1,pore_pressure_cache_l2])
            base_q=coupled_state%q

            interlayer_momentum_exchange%chezy_drag_coeff=0.05_wp
            coupled_state%q=before_q
            coupled_state%qp=before_qp
            CALL workspace%advance(coupled_state,0.0_wp,dt_step,stochastic,paired,operators, &
                pre_spatial_primitives_ready=.TRUE.,layer_physics=physics, &
                pore_caches=[pore_pressure_cache_l1,pore_pressure_cache_l2])

            drag_difference=0.0_wp
            DO l=1,paired%solve_cells
                j=paired%j_cent(l)
                k=paired%k_cent(l)
                DO eq=1,n_eqns
                    IF (eq==idx_uEqn .OR. eq==idx_vEqn) CYCLE
                    IF (coupled_state%q(eq,1,j,k)/=base_q(eq,1,j,k) .OR. &
                        coupled_state%q(eq,2,j,k)/=base_q(eq,2,j,k)) &
                        ERROR STOP 'Interlayer drag changed mass, composition or thermal energy'
                END DO
                DO eq=idx_uEqn,idx_vEqn
                    momentum_difference=SUM(coupled_state%q(eq,:,j,k)-base_q(eq,:,j,k))
                    scale=MAX(1.0_wp,SUM(ABS(base_q(eq,:,j,k))))
                    IF (ABS(momentum_difference)>1.0E-10_wp*scale) &
                        ERROR STOP 'Integrated interlayer drag did not conserve paired momentum'
                    drag_difference=drag_difference+SUM(ABS(coupled_state%q(eq,:,j,k)-base_q(eq,:,j,k)))
                END DO
                IF (.NOT.ALL(ieee_is_finite(coupled_state%qp(:,:,j,k)))) &
                    ERROR STOP 'Non-finite primitive state after integrated two-layer Lie timestep'
            END DO
            IF (scenario==1 .AND. drag_difference<=0.0_wp) &
                ERROR STOP 'Integrated wet/wet interlayer drag had no momentum effect'
            IF (scenario/=1 .AND. drag_difference/=0.0_wp) &
                ERROR STOP 'Integrated drag acted without two wet layers'
            ! Compare the complete accepted timestep payload across 1/4 threads.
            WRITE(integrated_unit) coupled_state%q,coupled_state%qp
        END DO

        ! Test simultaneous L1/L2 group-zero geometry sources on the same
        ! two-layer timestep. The no-source step is the post-spatial recipient;
        ! each curvature rate is evaluated on the pre-spatial primitive state.
        IF (curvature_term_flag) THEN
            before_q(:,1,:,:)=q1
            before_q(:,2,:,:)=q2
            pt1=p1
            pt2=p2
            CALL operators(1)%compute_timestep(before_q(:,1,:,:),pt1,0.0_wp,bound1,paired)
            CALL operators(2)%compute_timestep(before_q(:,2,:,:),pt2,0.0_wp,bound2,paired, &
                lower_operator=operators(1),h_lower=pt1(idx_h,:,:))
            dt_step=MIN(1.0E-5_wp,0.01_wp*MIN(dt_bound,bound1,bound2))
            IF (dt_step<=0.0_wp) ERROR STOP 'Invalid two-layer curvature timestep'
            before_qp(:,1,:,:)=pt1
            before_qp(:,2,:,:)=pt2
            interlayer_momentum_exchange%chezy_drag_coeff=0.0_wp
            interlayer_exchange%settling_flag=.FALSE.
            interlayer_exchange%elutriation_flag=.FALSE.
            coupled_state%q=before_q
            coupled_state%qp=before_qp
            CALL workspace%advance(coupled_state,0.0_wp,dt_step,stochastic,paired,operators, &
                pre_spatial_primitives_ready=.TRUE.,layer_physics=physics, &
                pore_caches=[pore_pressure_cache_l1,pore_pressure_cache_l2])
            base_q=coupled_state%q

            CALL workspace%finalize()
            CALL registry%finalize()
            groups(1)%source_names=[CHARACTER(LEN=64)::'geometry_curvature']
            groups(2)%source_names=[CHARACTER(LEN=64)::'geometry_curvature']
            groups(1)%source_mask=IBSET(0,SOURCE_ID_GEOMETRY_CURVATURE-1)
            groups(2)%source_mask=groups(1)%source_mask
            CALL registry%initialize(n_eqns,0,groups,error)
            IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid paired curvature Lie registry'
            CALL workspace%initialize(registry)
            coupled_state%q=before_q
            coupled_state%qp=before_qp
            CALL workspace%advance(coupled_state,0.0_wp,dt_step,stochastic,paired,operators, &
                pre_spatial_primitives_ready=.TRUE.,layer_physics=physics, &
                pore_caches=[pore_pressure_cache_l1,pore_pressure_cache_l2])

            curvature_effect=0.0_wp
            DO l=1,paired%solve_cells
                j=paired%j_cent(l)
                k=paired%k_cent(l)
                expected_q=base_q(:,:,j,k)
                DO layer=1,2
                    CALL eval_source_explicit(B_prime_x_geom(j,k),B_prime_y_geom(j,k), &
                        B_second_xx_geom(j,k),B_second_xy_geom(j,k),B_second_yy_geom(j,k), &
                        grav_coeff(j,k),before_qp(:,layer,j,k),curvature_rate,0.0_wp,0.0_wp, &
                        no_fissures,0.0_wp,0.0_wp,0.0_wp,dx*dy,groups(layer),0.0_wp, &
                        materials=materials,thermal=physics(layer)%thermal_configuration(), &
                        rheology=physics(layer)%rheology_configuration(),pore_enabled=.FALSE., &
                        stochastic_enabled=.FALSE.)
                    expected_q(1:n_eqns,layer)=expected_q(1:n_eqns,layer)+dt_step*curvature_rate
                END DO
                scale=MAX(1.0_wp,MAXVAL(ABS(expected_q)))
                IF (MAXVAL(ABS(expected_q-coupled_state%q(:,:,j,k)))>1.0E-10_wp*scale) &
                    ERROR STOP 'Paired group-zero curvature did not use pre-spatial rates and post-spatial recipients'
                curvature_effect=curvature_effect+SUM(ABS(coupled_state%q(:,:,j,k)-base_q(:,:,j,k)))
                IF (.NOT.ALL(ieee_is_finite(coupled_state%qp(:,:,j,k)))) &
                    ERROR STOP 'Nonfinite paired curvature primitive state'
            END DO
            IF (.NOT.ieee_is_finite(curvature_effect)) ERROR STOP 'Invalid paired curvature effect'
            WRITE(integrated_unit) coupled_state%q,coupled_state%qp

            CALL workspace%finalize()
            CALL registry%finalize()
            groups(1)%source_names=[CHARACTER(LEN=64)::]
            groups(2)%source_names=[CHARACTER(LEN=64)::]
            groups(1)%source_mask=0
            groups(2)%source_mask=0
            CALL registry%initialize(n_eqns,0,groups,error)
            IF (LEN_TRIM(error)>0) ERROR STOP 'Could not restore paired Lie registry'
            CALL workspace%initialize(registry)
        END IF

        ! Exercise the exact production Lie mass source on the same compact
        ! two-layer workset. Keep external erosion and deposition inactive:
        ! differences from the reference are purely interlayer transfers.
        created_diam=.NOT.ALLOCATED(diam_s)
        IF (created_diam) THEN
            ALLOCATE(diam_s(n_solid))
        ELSE
            original_diam=diam_s
        END IF
        diam_s=2.0E-5_wp
        kin_before=kin_visc_c
        vonK_before=vonK
        kin_visc_c=1.0E-6_wp
        vonK=0.41_wp
        interlayer_exchange%method='EXPLICIT'
        interlayer_exchange%group=4
        interlayer_exchange%elutriation_solid_volume_fraction=0.20_wp
        interlayer_exchange%elutriation_carrier_mode='LIMIT_SOLIDS'
        interlayer_exchange%elutriation_law='LINEAR_USTAR'
        interlayer_exchange%rouse_critical=100.0_wp
        interlayer_exchange%elutriation_efficiency=0.1_wp
        interlayer_momentum_exchange%chezy_drag_coeff=0.0_wp

        ! First run a later Lie mass group, then the same physical fixture
        ! with group zero. Check each accepted cell against the existing
        ! paired physical transaction, not a separately coded rate formula.
        DO group_mode=1,2
            IF (group_mode==2) THEN
                CALL workspace%finalize()
                CALL registry%finalize()
                groups(3)%id=0
                interlayer_exchange%group=0
                CALL registry%initialize(n_eqns,0,groups,error)
                IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid interlayer group-zero registry'
                CALL workspace%initialize(registry)
            END IF
            DO scenario=1,4
                before_q=0.0_wp
                before_qp=0.0_wp
                before_qp(idx_T,:,:,:)=T_ambient
                DO k=5,ny-4
                    DO j=5,nx-4
                        DO layer=1,2
                            IF (scenario==4) CYCLE
                            IF (scenario==2 .AND. layer==2) CYCLE
                            IF (scenario==3 .AND. layer==1) CYCLE
                            primitive=0.0_wp
                            primitive(idx_T)=T_ambient
                            primitive(idx_h)=0.3_wp
                            primitive(idx_u)=2.0_wp
                            IF (layer==1) THEN
                                primitive(idx_h)=0.5_wp
                                primitive(idx_u)=20.0_wp
                            END IF
                            primitive(idx_solid_first)=0.10_wp
                            CALL primitive_to_conservative(primitive,before_q(:,layer,j,k),accepted=accepted)
                            IF (.NOT.accepted) ERROR STOP 'Invalid solid-bearing interlayer cold-start fixture'
                            before_qp(:,layer,j,k)=primitive
                        END DO
                    END DO
                END DO

                pt1=before_qp(:,1,:,:)
                pt2=before_qp(:,2,:,:)
                CALL operators(1)%compute_timestep(before_q(:,1,:,:),pt1,0.0_wp,bound1,paired)
                CALL operators(2)%compute_timestep(before_q(:,2,:,:),pt2,0.0_wp,bound2,paired, &
                    lower_operator=operators(1),h_lower=pt1(idx_h,:,:))
                dt_step=MIN(1.0E-5_wp,0.01_wp*MIN(dt_bound,bound1,bound2))
                IF (dt_step<=0.0_wp) ERROR STOP 'Invalid mass-transfer integrated timestep'
                before_qp(:,1,:,:)=pt1
                before_qp(:,2,:,:)=pt2

                interlayer_exchange%settling_flag=.FALSE.
                interlayer_exchange%elutriation_flag=.FALSE.
                interlayer_momentum_exchange%chezy_drag_coeff=0.0_wp
                coupled_state%q=before_q
                coupled_state%qp=before_qp
                CALL workspace%advance(coupled_state,0.0_wp,dt_step,stochastic,paired,operators, &
                    pre_spatial_primitives_ready=.TRUE.,layer_physics=physics, &
                    pore_caches=[pore_pressure_cache_l1,pore_pressure_cache_l2])
                base_q=coupled_state%q
                base_qp=coupled_state%qp

                interlayer_exchange%settling_flag=.TRUE.
                interlayer_exchange%elutriation_flag=.TRUE.
                ! In group-zero mode, the wet/wet scenario also verifies that
                ! analytic drag follows the paired mass transfer in Lie order.
                interlayer_momentum_exchange%chezy_drag_coeff=0.0_wp
                IF (group_mode==2 .AND. scenario==1) interlayer_momentum_exchange%chezy_drag_coeff=0.05_wp
                coupled_state%q=before_q
                coupled_state%qp=before_qp
                CALL workspace%advance(coupled_state,0.0_wp,dt_step,stochastic,paired,operators, &
                    pre_spatial_primitives_ready=.TRUE.,layer_physics=physics, &
                    pore_caches=[pore_pressure_cache_l1,pore_pressure_cache_l2])

                mass_difference=0.0_wp
                DO l=1,paired%solve_cells
                    j=paired%j_cent(l)
                    k=paired%k_cent(l)
                    expected_q=base_q(:,:,j,k)
                    expected_qp=base_qp(:,:,j,k)
                    CALL apply_interlayer_lie_group(interlayer_exchange,rheology_l2,mass_exchange_l2, &
                        before_qp(:,:,j,k),group_mode==2,dt_step,expected_q,expected_qp, &
                        reference_physical,source_status)
                    IF (source_status/=0) ERROR STOP 'Rejected reference interlayer Lie transaction'
                    IF (group_mode==2 .AND. scenario==1) THEN
                        CALL apply_interlayer_momentum_exchange_analytic(interlayer_momentum_exchange,dt_step, &
                            expected_q,expected_qp,accepted)
                        IF (.NOT.accepted) ERROR STOP 'Rejected reference ordered interlayer drag'
                    END IF
                    IF (ANY(expected_q/=coupled_state%q(:,:,j,k)) .OR. &
                        ANY(expected_qp/=coupled_state%qp(:,:,j,k))) &
                        ERROR STOP 'Integrated Lie group order or group-zero recipient/snapshot mismatch'
                    DO eq=1,n_eqns
                        component_difference=SUM(coupled_state%q(eq,:,j,k)-base_q(eq,:,j,k))
                        scale=MAX(1.0_wp,SUM(ABS(base_q(eq,:,j,k))))
                        IF (ABS(component_difference)>1.0E-10_wp*scale) &
                            ERROR STOP 'Integrated interlayer mass exchange violates paired conservative budget'
                    END DO
                    DO layer=1,2
                        IF (.NOT.conservative_mass_admissible(coupled_state%q(:,layer,j,k))) &
                            ERROR STOP 'Integrated interlayer mass exchange produced inadmissible state'
                        IF (.NOT.ALL(ieee_is_finite(coupled_state%qp(:,layer,j,k)))) &
                            ERROR STOP 'Integrated mass exchange produced non-finite primitives'
                    END DO
                    mass_difference=mass_difference+SUM(ABS(coupled_state%q(idx_solidEqn_first,:,j,k) &
                        -base_q(idx_solidEqn_first,:,j,k)))
                END DO
                IF (scenario/=4 .AND. mass_difference<=0.0_wp) &
                    ERROR STOP 'Integrated interlayer settling or elutriation did not transfer solids'
                IF (scenario==4 .AND. mass_difference/=0.0_wp) &
                    ERROR STOP 'Dry column generated interlayer mass exchange'
                IF (scenario==2 .AND. .NOT.ANY(coupled_state%q(idx_totMassEqn,2,:,:)> &
                    base_q(idx_totMassEqn,2,:,:))) &
                    ERROR STOP 'Integrated elutriation failed to cold-start L2'
                IF (scenario==3 .AND. .NOT.ANY(coupled_state%q(idx_totMassEqn,1,:,:)> &
                    base_q(idx_totMassEqn,1,:,:))) &
                    ERROR STOP 'Integrated settling failed to cold-start L1'
                WRITE(integrated_unit) coupled_state%q,coupled_state%qp
            END DO
        END DO
        ! Integrate external deposition and internal exchange in the same
        ! cell-outer timestep. The shared terrain inventory receives deposits
        ! only from the exposed layer; the paired transfer conserves solids.
        CALL workspace%finalize()
        CALL registry%finalize()
        groups(3)%id=4
        interlayer_exchange%group=4
        external_groups(1:2)=groups(1:2)
        external_groups(3)=groups(1)
        external_groups(3)%id=4
        external_groups(3)%source_names=[CHARACTER(LEN=64)::'mass_exchange_deposition']
        external_groups(3)%source_mask=IBSET(0,SOURCE_ID_EXCHANGE_DEPOSITION-1)
        external_groups(4)=external_groups(3)
        external_groups(4)%scope=SOURCE_SCOPE_L2
        external_groups(5)=groups(3)
        external_groups(6)=groups(4)
        CALL registry%initialize(n_eqns,0,external_groups,error)
        IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid paired external/interlayer Lie registry'
        CALL workspace%initialize(registry)
        ! External source groups refresh aggregate bed thickness even with
        ! topographic evolution disabled. Use a physically valid fixture bed.
        substrate%erodible_porosity=0.25_wp
        substrate%T_ground=T_ambient
        mass_exchange_l1%deposition%settling_flag=.TRUE.
        mass_exchange_l2%deposition%settling_flag=.TRUE.
        CALL physics(1)%finalize()
        CALL physics(2)%finalize()
        CALL physics(1)%initialize(.FALSE.,rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1,error)
        IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid L1 deposition Lie physics'
        CALL physics(2)%initialize(.FALSE.,rheology_l2,thermal_l2,pore_pressure_l2,mass_exchange_l2,error)
        IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid L2 deposition Lie physics'
        interlayer_exchange%settling_flag=.TRUE.
        interlayer_exchange%elutriation_flag=.TRUE.
        interlayer_momentum_exchange%chezy_drag_coeff=0.0_wp

        DO scenario=1,4
            before_q=0.0_wp
            before_qp=0.0_wp
            before_qp(idx_T,:,:,:)=T_ambient
            DO k=5,ny-4
                DO j=5,nx-4
                    DO layer=1,2
                        IF (scenario==4) CYCLE
                        IF (scenario==2 .AND. layer==2) CYCLE
                        IF (scenario==3 .AND. layer==1) CYCLE
                        primitive=0.0_wp
                        primitive(idx_T)=T_ambient
                        primitive(idx_h)=0.5_wp
                        primitive(idx_u)=2.0_wp
                        primitive(idx_solid_first)=0.10_wp
                        CALL primitive_to_conservative(primitive,before_q(:,layer,j,k),accepted=accepted)
                        IF (.NOT.accepted) ERROR STOP 'Invalid external exchange fixture state'
                        before_qp(:,layer,j,k)=primitive
                    END DO
                END DO
            END DO
            pt1=before_qp(:,1,:,:)
            pt2=before_qp(:,2,:,:)
            CALL operators(1)%compute_timestep(before_q(:,1,:,:),pt1,0.0_wp,bound1,paired)
            CALL operators(2)%compute_timestep(before_q(:,2,:,:),pt2,0.0_wp,bound2,paired, &
                lower_operator=operators(1),h_lower=pt1(idx_h,:,:))
            dt_step=MIN(1.0E-5_wp,0.01_wp*MIN(dt_bound,bound1,bound2))
            IF (dt_step<=0.0_wp) ERROR STOP 'Invalid external exchange timestep'
            before_qp(:,1,:,:)=pt1
            before_qp(:,2,:,:)=pt2

            ! Rebuild the same cell workset and source schedule with the
            ! deposition closures temporarily disabled, as a paired reference.
            mass_exchange_l1%deposition%settling_flag=.FALSE.
            mass_exchange_l2%deposition%settling_flag=.FALSE.
            CALL physics(1)%finalize()
            CALL physics(2)%finalize()
            CALL physics(1)%initialize(.FALSE.,rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1,error)
            IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid L1 reference exchange physics'
            CALL physics(2)%initialize(.FALSE.,rheology_l2,thermal_l2,pore_pressure_l2,mass_exchange_l2,error)
            IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid L2 reference exchange physics'
            coupled_state%q=before_q
            coupled_state%qp=before_qp
            coupled_state%deposit=0.0_wp
            CALL workspace%advance(coupled_state,0.0_wp,dt_step,stochastic,paired,operators, &
                pre_spatial_primitives_ready=.TRUE.,layer_physics=physics, &
                pore_caches=[pore_pressure_cache_l1,pore_pressure_cache_l2])
            base_q=coupled_state%q

            mass_exchange_l1%deposition%settling_flag=.TRUE.
            mass_exchange_l2%deposition%settling_flag=.TRUE.
            CALL physics(1)%finalize()
            CALL physics(2)%finalize()
            CALL physics(1)%initialize(.FALSE.,rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1,error)
            IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid active L1 exchange physics'
            CALL physics(2)%initialize(.FALSE.,rheology_l2,thermal_l2,pore_pressure_l2,mass_exchange_l2,error)
            IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid active L2 exchange physics'
            coupled_state%q=before_q
            coupled_state%qp=before_qp
            coupled_state%deposit=0.0_wp
            CALL workspace%advance(coupled_state,0.0_wp,dt_step,stochastic,paired,operators, &
                pre_spatial_primitives_ready=.TRUE.,layer_physics=physics, &
                pore_caches=[pore_pressure_cache_l1,pore_pressure_cache_l2])

            expected_deposition=SUM(coupled_state%deposit(:,:,1))
            IF (scenario==4) THEN
                IF (expected_deposition/=0.0_wp) &
                    ERROR STOP 'Dry paired column deposited solids'
            ELSE
                IF (expected_deposition<=0.0_wp) &
                    ERROR STOP 'Exposed layer did not deposit into shared substrate'
            END IF
            DO l=1,paired%solve_cells
                j=paired%j_cent(l)
                k=paired%k_cent(l)
                solid_budget=SUM(coupled_state%q(idx_solidEqn_first,:,j,k) &
                    -base_q(idx_solidEqn_first,:,j,k)) &
                    +rho_s(1)*coupled_state%deposit(j,k,1)
                scale=MAX(1.0_wp,SUM(ABS(base_q(idx_solidEqn_first,:,j,k))))
                IF (ABS(solid_budget)>1.0E-9_wp*scale) &
                    ERROR STOP 'External deposition and interlayer exchange do not conserve paired solids'
                DO layer=1,2
                    IF (.NOT.conservative_mass_admissible(coupled_state%q(:,layer,j,k))) &
                        ERROR STOP 'Paired external exchange produced inadmissible state'
                END DO
            END DO
            WRITE(integrated_unit) coupled_state%q,coupled_state%qp,coupled_state%deposit
        END DO
        ! A second external exchange fixture covers erosion from the one
        ! shared substrate, with both external groups and the paired mass
        ! exchange present in the exact same production Lie sequence.
        CALL workspace%finalize()
        CALL registry%finalize()
        external_groups(3)%source_names=[CHARACTER(LEN=64)::'mass_exchange_erosion']
        external_groups(4)%source_names=[CHARACTER(LEN=64)::'mass_exchange_erosion']
        external_groups(3)%source_mask=IBSET(0,SOURCE_ID_EXCHANGE_EROSION-1)
        external_groups(4)%source_mask=external_groups(3)%source_mask
        CALL registry%initialize(n_eqns,0,external_groups,error)
        IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid paired erosion/interlayer Lie registry'
        CALL workspace%initialize(registry)
        mass_exchange_l1%deposition%settling_flag=.FALSE.
        mass_exchange_l2%deposition%settling_flag=.FALSE.
        mass_exchange_l1%erosion%coefficient=0.3_wp
        mass_exchange_l2%erosion%coefficient=0.3_wp
        substrate%T_ground=T_ambient
        substrate%erodible_porosity=0.25_wp
        created_erodible_fraction=.NOT.ALLOCATED(erodible_fract)
        IF (created_erodible_fraction) THEN
            ALLOCATE(erodible_fract(n_solid))
        ELSE
            original_erodible_fraction=erodible_fract
        END IF
        erodible_fract=1.0_wp
        CALL physics(1)%finalize()
        CALL physics(2)%finalize()
        CALL physics(1)%initialize(.FALSE.,rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1,error)
        IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid L1 erosion Lie physics'
        CALL physics(2)%initialize(.FALSE.,rheology_l2,thermal_l2,pore_pressure_l2,mass_exchange_l2,error)
        IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid L2 erosion Lie physics'

        DO scenario=1,4
            before_q=0.0_wp
            before_qp=0.0_wp
            before_qp(idx_T,:,:,:)=T_ambient
            DO k=5,ny-4
                DO j=5,nx-4
                    DO layer=1,2
                        IF (scenario==4) CYCLE
                        IF (scenario==2 .AND. layer==2) CYCLE
                        IF (scenario==3 .AND. layer==1) CYCLE
                        primitive=0.0_wp
                        primitive(idx_T)=T_ambient
                        primitive(idx_h)=0.5_wp
                        primitive(idx_u)=2.0_wp
                        primitive(idx_solid_first)=0.10_wp
                        CALL primitive_to_conservative(primitive,before_q(:,layer,j,k),accepted=accepted)
                        IF (.NOT.accepted) ERROR STOP 'Invalid paired erosion fixture state'
                        before_qp(:,layer,j,k)=primitive
                    END DO
                END DO
            END DO
            pt1=before_qp(:,1,:,:)
            pt2=before_qp(:,2,:,:)
            CALL operators(1)%compute_timestep(before_q(:,1,:,:),pt1,0.0_wp,bound1,paired)
            CALL operators(2)%compute_timestep(before_q(:,2,:,:),pt2,0.0_wp,bound2,paired, &
                lower_operator=operators(1),h_lower=pt1(idx_h,:,:))
            dt_step=MIN(1.0E-5_wp,0.01_wp*MIN(dt_bound,bound1,bound2))
            IF (dt_step<=0.0_wp) ERROR STOP 'Invalid paired erosion timestep'
            before_qp(:,1,:,:)=pt1
            before_qp(:,2,:,:)=pt2

            mass_exchange_l1%erosion%coefficient=-1.0_wp
            mass_exchange_l2%erosion%coefficient=-1.0_wp
            CALL physics(1)%finalize()
            CALL physics(2)%finalize()
            CALL physics(1)%initialize(.FALSE.,rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1,error)
            IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid L1 no-erosion reference'
            CALL physics(2)%initialize(.FALSE.,rheology_l2,thermal_l2,pore_pressure_l2,mass_exchange_l2,error)
            IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid L2 no-erosion reference'
            coupled_state%q=before_q
            coupled_state%qp=before_qp
            coupled_state%erodible=0.2_wp
            coupled_state%erosion=0.0_wp
            CALL workspace%advance(coupled_state,0.0_wp,dt_step,stochastic,paired,operators, &
                pre_spatial_primitives_ready=.TRUE.,layer_physics=physics, &
                pore_caches=[pore_pressure_cache_l1,pore_pressure_cache_l2])
            base_q=coupled_state%q

            mass_exchange_l1%erosion%coefficient=0.3_wp
            mass_exchange_l2%erosion%coefficient=0.3_wp
            CALL physics(1)%finalize()
            CALL physics(2)%finalize()
            CALL physics(1)%initialize(.FALSE.,rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1,error)
            IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid active L1 erosion physics'
            CALL physics(2)%initialize(.FALSE.,rheology_l2,thermal_l2,pore_pressure_l2,mass_exchange_l2,error)
            IF (LEN_TRIM(error)>0) ERROR STOP 'Invalid active L2 erosion physics'
            coupled_state%q=before_q
            coupled_state%qp=before_qp
            coupled_state%erodible=0.2_wp
            coupled_state%erosion=0.0_wp
            CALL workspace%advance(coupled_state,0.0_wp,dt_step,stochastic,paired,operators, &
                pre_spatial_primitives_ready=.TRUE.,layer_physics=physics, &
                pore_caches=[pore_pressure_cache_l1,pore_pressure_cache_l2])
            expected_deposition=SUM(coupled_state%erosion(:,:,1))
            IF (scenario==4) THEN
                IF (expected_deposition/=0.0_wp) ERROR STOP 'Dry paired column eroded the substrate'
            ELSE
                IF (expected_deposition<=0.0_wp) ERROR STOP 'Exposed layer did not erode shared substrate'
            END IF
            DO l=1,paired%solve_cells
                j=paired%j_cent(l)
                k=paired%k_cent(l)
                solid_budget=SUM(coupled_state%q(idx_solidEqn_first,:,j,k) &
                    -base_q(idx_solidEqn_first,:,j,k)) &
                    +rho_s(1)*(coupled_state%erodible(1,j,k)-0.2_wp)
                scale=MAX(1.0_wp,SUM(ABS(base_q(idx_solidEqn_first,:,j,k))))
                IF (ABS(solid_budget)>1.0E-9_wp*scale) &
                    ERROR STOP 'External erosion and interlayer exchange do not conserve paired solids'
                IF (coupled_state%erodible(1,j,k)<0.0_wp) &
                    ERROR STOP 'Paired Lie erosion exhausted substrate inventory'
                DO layer=1,2
                    IF (.NOT.conservative_mass_admissible(coupled_state%q(:,layer,j,k))) &
                        ERROR STOP 'Paired external erosion produced inadmissible fluid state'
                END DO
            END DO
            WRITE(integrated_unit) coupled_state%q,coupled_state%qp,coupled_state%erodible
        END DO
        substrate=old_substrate
        IF (created_erodible_fraction) THEN
            DEALLOCATE(erodible_fract)
        ELSE
            erodible_fract=original_erodible_fraction
        END IF
        mass_exchange_l1=old_external_l1
        mass_exchange_l2=old_external_l2
        IF (created_diam) THEN
            DEALLOCATE(diam_s)
        ELSE
            diam_s=original_diam
        END IF
        kin_visc_c=kin_before
        vonK=vonK_before
        CLOSE(integrated_unit)

        CALL workspace%finalize()
        CALL operators(1)%finalize()
        CALL operators(2)%finalize()
        CALL registry%finalize()
        CALL physics(1)%finalize()
        CALL physics(2)%finalize()
        IF (created_zone) DEALLOCATE(B_zone)
        IF (created_arc) DEALLOCATE(cell_arc_perim)
        IF (created_source_fraction) DEALLOCATE(cell_source_fractions)
        IF (created_fissure_fraction) DEALLOCATE(cell_fissure_fractions)
        IF (created_arc_normals) DEALLOCATE(cell_arc_n_x,cell_arc_n_y)
        interlayer_momentum_exchange=old_drag
        interlayer_exchange=old_mass
        n_layers=original_layers
        WRITE(*,*) 'PASS: coupled spatial-Lie-commit, interlayer group 0/later, ordered drag, wet/dry'
    END SUBROUTINE check_integrated_coupled_lie_timestep

    !> @brief Check exact global conservative budgets for compact-support layers.
    !> @details Both layers are dry in a collar, so no boundary flux is present.
    !> Only the mass, thermal and constituent equations have no pressure paths.
    SUBROUTINE check_closed_coupled_spatial_budget
        TYPE(domain_type) :: paired
        TYPE(spatial_operator_type) :: lower,upper
        REAL(wp) :: q1(n_vars,nx,ny),q2(n_vars,nx,ny)
        REAL(wp) :: p1(n_vars,nx,ny),p2(n_vars,nx,ny)
        REAL(wp) :: r1(n_eqns,nx,ny),r2(n_eqns,nx,ny)
        REAL(wp) :: primitive(n_vars),dt1,dt2,dt_pair
        REAL(wp) :: balance,scale,tol
        INTEGER :: j,k,l,layer,i,eq
        LOGICAL :: accepted

        IF (nx<12 .OR. ny<12) RETURN
        q1=0.0_wp
        q2=0.0_wp
        p1=0.0_wp
        p2=0.0_wp
        p1(idx_T,:,:)=T_ambient
        p2(idx_T,:,:)=T_ambient
        primitive=0.0_wp
        primitive(idx_T)=T_ambient
        DO k=5,ny-4
            DO j=5,nx-4
                primitive(idx_h)=0.5_wp
                primitive(idx_u)=0.2_wp
                CALL primitive_to_conservative(primitive,q1(:,j,k),accepted=accepted)
                IF (.NOT.accepted) ERROR STOP 'Invalid compact lower-layer fixture'
                primitive(idx_h)=0.3_wp
                primitive(idx_u)=-0.1_wp
                CALL primitive_to_conservative(primitive,q2(:,j,k),accepted=accepted)
                IF (.NOT.accepted) ERROR STOP 'Invalid compact upper-layer fixture'
            END DO
        END DO
        CALL paired%initialize
        CALL lower%initialize(1)
        CALL upper%initialize(2)
        CALL paired%check_solve(q1,0.0_wp,.FALSE.,q_other=q2)
        CALL lower%prepare_hp_workset(paired)
        CALL upper%prepare_hp_workset(paired)
        CALL lower%compute_timestep(q1,p1,0.0_wp,dt1,paired)
        CALL upper%compute_timestep(q2,p2,0.0_wp,dt2,paired, &
            lower_operator=lower,h_lower=p1(idx_h,:,:))
        CALL lower%evaluate_coupled(upper,p1,p2,r1,r2,0.0_wp,paired,dt_pair)

        DO layer=1,2
            DO i=1,3
                SELECT CASE(i)
                CASE(1)
                    eq=idx_totMassEqn
                CASE(2)
                    eq=idx_engyEqn
                CASE(3)
                    eq=idx_solidEqn_first
                END SELECT
                balance=0.0_wp
                scale=0.0_wp
                DO l=1,paired%solve_cells
                    j=paired%j_cent(l)
                    k=paired%k_cent(l)
                    IF (layer==1) THEN
                        balance=balance+r1(eq,j,k)
                        scale=scale+ABS(r1(eq,j,k))
                    ELSE
                        balance=balance+r2(eq,j,k)
                        scale=scale+ABS(r2(eq,j,k))
                    END IF
                END DO
                tol=4096.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,scale)
                IF (ABS(balance)>tol) &
                    ERROR STOP 'Compact two-layer spatial conservative budget mismatch'
            END DO
        END DO
        CALL check_integrated_coupled_lie_timestep(paired,q1,q2,p1,p2,dt_pair)
        CALL lower%finalize
        CALL upper%finalize
        CALL paired%finalize
        WRITE(*,*) 'PASS: compact two-layer spatial mass, thermal and constituent budgets'
    END SUBROUTINE check_closed_coupled_spatial_budget

    !> @brief Quantify stratified-rest residuals on flat, affine and curved beds.
    !> @details Two constant free surfaces eta1 and eta2 give stationary layers.
    !> Flat and affine geometry must balance within roundoff; curved geometry
    !> is measured without assuming exact momentum well-balancing.
    !> Original shared geometry is restored before other tests run.
    SUBROUTINE check_stratified_rest_geometry
        TYPE(domain_type) :: paired
        TYPE(spatial_operator_type) :: lower,upper
        REAL(wp) :: saved_b(nx+1,ny+1)
        REAL(wp) :: q1(n_vars,nx,ny),q2(n_vars,nx,ny)
        REAL(wp) :: p1(n_vars,nx,ny),p2(n_vars,nx,ny)
        REAL(wp) :: r1(n_eqns,nx,ny),r2(n_eqns,nx,ny)
        REAL(wp) :: primitive(n_vars),eta1,dt1,dt2,dt_pair
        REAL(wp) :: max_m1,max_m2,max_other,scale,tolerance
        INTEGER :: terrain,j,k,l,i
        LOGICAL :: accepted

        IF (nx<8 .OR. ny<8) RETURN
        saved_b=B_vertex
        CALL paired%initialize
        CALL lower%initialize(1)
        CALL upper%initialize(2)
        primitive=0.0_wp
        primitive(idx_T)=T_ambient

        DO terrain=1,3
            DO k=1,ny+1
                DO j=1,nx+1
                    SELECT CASE(terrain)
                    CASE(1)
                        B_vertex(j,k)=1.0_wp
                    CASE(2)
                        B_vertex(j,k)=1.0_wp+0.01_wp*x_stag(j)+0.02_wp*y_stag(k)
                    CASE(3)
                        B_vertex(j,k)=1.0_wp+0.002_wp*x_stag(j)**2+0.003_wp*y_stag(k)**2
                    END SELECT
                END DO
            END DO
            CALL refresh_topography_geometry

            eta1=MAXVAL(B_cent)+2.0_wp
            DO k=1,ny
                DO j=1,nx
                    primitive(idx_h)=eta1-B_cent(j,k)
                    CALL primitive_to_conservative(primitive,q1(:,j,k),accepted=accepted)
                    IF (.NOT.accepted) ERROR STOP 'Invalid lower stratified-rest state'
                    primitive(idx_h)=0.7_wp
                    CALL primitive_to_conservative(primitive,q2(:,j,k),accepted=accepted)
                    IF (.NOT.accepted) ERROR STOP 'Invalid upper stratified-rest state'
                END DO
            END DO

            CALL paired%check_solve(q1,0.0_wp,.FALSE.,q_other=q2)
            CALL lower%prepare_hp_workset(paired)
            CALL upper%prepare_hp_workset(paired)
            p1=0.0_wp
            p2=0.0_wp
            CALL lower%compute_timestep(q1,p1,0.0_wp,dt1,paired)
            CALL upper%compute_timestep(q2,p2,0.0_wp,dt2,paired, &
                lower_operator=lower,h_lower=p1(idx_h,:,:))
            CALL lower%evaluate_coupled(upper,p1,p2,r1,r2,0.0_wp,paired,dt_pair)
            IF (.NOT.ieee_is_finite(dt_pair) .OR. dt_pair<=0.0_wp) &
                ERROR STOP 'Invalid stratified-rest coupled timestep'

            max_m1=0.0_wp
            max_m2=0.0_wp
            max_other=0.0_wp
            DO l=1,paired%solve_cells
                j=paired%j_cent(l)
                k=paired%k_cent(l)
                IF (.NOT.ALL(ieee_is_finite(r1(:,j,k))) .OR. &
                    .NOT.ALL(ieee_is_finite(r2(:,j,k)))) &
                    ERROR STOP 'Nonfinite stratified-rest spatial residual'
                ! Exclude the imposed boundary collar from equilibrium estimates.
                IF (j<=3 .OR. j>=nx-2 .OR. k<=3 .OR. k>=ny-2) CYCLE
                max_m1=MAX(max_m1,MAXVAL(ABS(r1(idx_uEqn:idx_vEqn,j,k))))
                max_m2=MAX(max_m2,MAXVAL(ABS(r2(idx_uEqn:idx_vEqn,j,k))))
                DO i=1,n_eqns
                    IF (i==idx_uEqn .OR. i==idx_vEqn) CYCLE
                    max_other=MAX(max_other,ABS(r1(i,j,k)),ABS(r2(i,j,k)))
                END DO
            END DO

            scale=grav*MAXVAL(eta1-B_cent+0.7_wp)*MAX(one_by_dx,one_by_dy)
            tolerance=1.0E-9_wp*MAX(1.0_wp,scale)
            IF (max_other>tolerance) &
                ERROR STOP 'Stratified rest created scalar transport'
            IF (terrain<=2 .AND. (max_m1>tolerance .OR. max_m2>tolerance)) &
                ERROR STOP 'Flat/affine stratified equilibrium failed'
            WRITE(*,'(A,I0,A,3(ES13.5,1X))') 'STRATIFIED_REST terrain=',terrain, &
                ' momentum_L1_L2_other=',max_m1,max_m2,max_other
        END DO

        B_vertex=saved_b
        CALL refresh_topography_geometry
        CALL upper%finalize
        CALL lower%finalize
        CALL paired%finalize
        WRITE(*,*) 'PASS: stratified rest on flat, affine, curved Q1 beds'
    END SUBROUTINE check_stratified_rest_geometry

    !> @brief Verify an interlayer face on two active cells for all four layer states.
    !> @details Compare directly with the oriented PCCU pressure kernel, check
    !> that only momentum changes and that the sparse workset is respected.
    !> This checks spatial assembly, not the eventual coupled time integrator.
    SUBROUTINE check_interlayer_face_residual
        TYPE(reconstruction_workspace_type) :: l1,l2
        TYPE(hyperbolic_workspace_type) :: hp2
        REAL(wp) :: r1(n_eqns,nx,ny),r2(n_eqns,nx,ny)
        REAL(wp) :: fl1(n_eqns),fr1(n_eqns),fl2(n_eqns),fr2(n_eqns)
        REAL(wp) :: h1a,h1b,h2a,h2b,am,ap,tol,scale
        LOGICAL :: active(nx,ny)
        INTEGER :: case_number,i

        IF (nx<3 .OR. ny<3) RETURN

        CALL l1%initialize(layer=1)
        CALL l2%initialize(layer=2)
        CALL hp2%initialize(clamp_buoyant_pressure=.TRUE.)
        active=.FALSE.
        active(1,2)=.TRUE.
        active(2,2)=.TRUE.

        DO case_number=1,4
            h1a=1.2_wp
            h1b=1.4_wp
            h2a=0.8_wp
            h2b=0.6_wp
            am=-2.0_wp
            ap=3.0_wp

            SELECT CASE(case_number)
            CASE(1)
                ! Both layers present with nonuniform upper pressure.
                CONTINUE
            CASE(2)
                ! Only L1 present, including its shoreline pressure path.
                h2a=0.0_wp
                h2b=0.0_wp
            CASE(3)
                ! Only L2 present, with no lower-layer pressure.
                h1a=0.0_wp
                h1b=0.0_wp
            CASE(4)
                ! Both layers dry: the characteristic interval is degenerate.
                h1a=0.0_wp
                h1b=0.0_wp
                h2a=0.0_wp
                h2b=0.0_wp
                am=0.0_wp
                ap=0.0_wp
            END SELECT

            r1=0.0_wp
            r2=0.0_wp
            l1%h_selected_interfaceL(2,2)=h1a
            l1%h_selected_interfaceR(2,2)=h1b
            l2%h_selected_interfaceL(2,2)=h2a
            l2%h_selected_interfaceR(2,2)=h2b
            l2%eta_interfaceL(2,2)=h1a+h2a
            l2%eta_interfaceR(2,2)=h1b+h2b
            hp2%gamma_interface_xL(2,2)=4.0_wp
            hp2%gamma_interface_xR(2,2)=6.0_wp
            hp2%a_interface_xNeg(2,2)=am
            hp2%a_interface_xPos(2,2)=ap

            CALL eval_interlayer_oriented_face(1, &
                h1a,h2a,4.0_wp,h1a+h2a,h1b,h2b,6.0_wp,h1b+h2b, &
                grav_coeff_stag_x(2,2),am,ap,fl1,fr1,fl2,fr2)
            CALL add_interlayer_terms(l1,l2,hp2,r1,r2,active,0,[INTEGER::],[INTEGER::],qp, &
                1,[2],[2],0,[INTEGER::],[INTEGER::], &
                hp2%a_interface_xNeg,hp2%a_interface_xPos, &
                hp2%a_interface_yNeg,hp2%a_interface_yPos)

            scale=MAX(1.0_wp,MAXVAL(ABS(fl1)),MAXVAL(ABS(fr1)), &
                MAXVAL(ABS(fl2)),MAXVAL(ABS(fr2)))
            tol=128.0_wp*EPSILON(1.0_wp)*scale*one_by_dx
            IF (ANY(ABS(r1(:,1,2)-fl1*one_by_dx)>tol) .OR. &
                ANY(ABS(r2(:,1,2)-fl2*one_by_dx)>tol) .OR. &
                ANY(ABS(r1(:,2,2)+fr1*one_by_dx)>tol) .OR. &
                ANY(ABS(r2(:,2,2)+fr2*one_by_dx)>tol)) &
                ERROR STOP 'Interlayer oriented face not assembled in adjacent cells'

            ! The interlayer operator exchanges pressure momentum only.
            ! Neither layer can gain mass, energy or scalar constituents here.
            DO i=1,n_eqns
                IF (i==2) CYCLE
                IF (r1(i,1,2)/=0.0_wp .OR. r1(i,2,2)/=0.0_wp .OR. &
                    r2(i,1,2)/=0.0_wp .OR. r2(i,2,2)/=0.0_wp) &
                    ERROR STOP 'Interlayer face changed a non-momentum equation'
            END DO

            IF (case_number>1) THEN
                IF (ANY(r1/=0.0_wp) .OR. ANY(r2/=0.0_wp)) &
                    ERROR STOP 'Absent interlayer contact generated pressure'
            END IF

            r1(:,1:2,2)=0.0_wp
            r2(:,1:2,2)=0.0_wp
            IF (ANY(r1/=0.0_wp) .OR. ANY(r2/=0.0_wp)) &
                ERROR STOP 'Interlayer face modified an inactive cell'
        END DO

        CALL hp2%finalize
        CALL l1%finalize
        CALL l2%finalize
        WRITE(*,*) 'PASS: four layer-presence cases, face momentum and sparse workset'
    END SUBROUTINE check_interlayer_face_residual

    !> @brief Check within-cell pressure assembly without any active coupling face.
    !> @details Uses two wet endpoints and variable Gamma2 and metric. It
    !> compares the L1/L2 residuals with the same nonconservative path kernel.
    SUBROUTINE check_interlayer_cell_assembly
        TYPE(reconstruction_workspace_type) :: lower,upper
        TYPE(hyperbolic_workspace_type) :: upper_hp
        REAL(wp) :: r1(n_eqns,nx,ny),r2(n_eqns,nx,ny),qp_upper(n_vars,nx,ny)
        REAL(wp) :: path1(n_eqns),path2(n_eqns),tol,scale
        LOGICAL :: active(nx,ny)
        INTEGER, PARAMETER :: j=2,k=2

        IF (nx<3 .OR. ny<3) RETURN

        CALL lower%initialize(layer=1)
        CALL upper%initialize(layer=2)
        CALL upper_hp%initialize(clamp_buoyant_pressure=.TRUE.)
        r1=0.0_wp
        r2=0.0_wp
        qp_upper=0.0_wp
        qp_upper(1,j,k)=1.0_wp
        active=.FALSE.
        active(j,k)=.TRUE.

        ! No face jumps: isolate the within-cell x path only.
        lower%eta_cellW(j,k)=B_face_x(j,k)+1.2_wp
        lower%eta_cellE(j,k)=B_face_x(j+1,k)+1.7_wp
        upper%eta_cellW(j,k)=lower%eta_cellW(j,k)+0.8_wp
        upper%eta_cellE(j,k)=lower%eta_cellE(j,k)+1.1_wp
        upper%qp_adv_interfaceR(:,j,k)=0.0_wp
        upper%qp_adv_interfaceL(:,j+1,k)=0.0_wp
        upper%h_selected_interfaceR(j,k)=0.8_wp
        upper%h_selected_interfaceL(j+1,k)=1.1_wp
        upper_hp%gamma_interface_xR(j,k)=4.0_wp
        upper_hp%gamma_interface_xL(j+1,k)=6.0_wp
        upper%qp_adv_interfaceR(1,j,k)=0.8_wp
        upper%qp_adv_interfaceL(1,j+1,k)=1.1_wp

        CALL eval_interlayer_pressure_correction(1, &
            1.2_wp,0.8_wp,4.0_wp,upper%eta_cellW(j,k),grav_coeff_stag_x(j,k), &
            1.7_wp,1.1_wp,6.0_wp,upper%eta_cellE(j,k),grav_coeff_stag_x(j+1,k),path1,path2)

        ! Keep y dimensionality active but use a dry, constant y cell path.
        lower%eta_cellS(j,k)=B_face_y(j,k)+1.2_wp
        lower%eta_cellN(j,k)=B_face_y(j,k+1)+1.2_wp
        upper%eta_cellS(j,k)=lower%eta_cellS(j,k)
        upper%eta_cellN(j,k)=lower%eta_cellN(j,k)
        upper%qp_adv_interfaceT(:,j,k)=0.0_wp
        upper%qp_adv_interfaceB(:,j,k+1)=0.0_wp
        upper%h_selected_interfaceT(j,k)=0.0_wp
        upper%h_selected_interfaceB(j,k+1)=0.0_wp
        upper_hp%gamma_interface_yT(j,k)=0.0_wp
        upper_hp%gamma_interface_yB(j,k+1)=0.0_wp

        CALL add_interlayer_terms(lower,upper,upper_hp,r1,r2,active,1,[j],[k],qp_upper, &
            0,[INTEGER::],[INTEGER::],0,[INTEGER::],[INTEGER::], &
            upper_hp%a_interface_xNeg,upper_hp%a_interface_xPos, &
            upper_hp%a_interface_yNeg,upper_hp%a_interface_yPos)

        scale=MAX(1.0_wp,MAXVAL(ABS(path1)),MAXVAL(ABS(path2)))
        tol=128.0_wp*EPSILON(1.0_wp)*scale*one_by_dx
        IF (ANY(ABS(r1(:,j,k)+path1*one_by_dx)>tol) .OR. &
            ANY(ABS(r2(:,j,k)+path2*one_by_dx)>tol)) &
            ERROR STOP 'Within-cell interlayer path has wrong residual sign'
        r1(:,j,k)=0.0_wp
        r2(:,j,k)=0.0_wp
        IF (ANY(r1/=0.0_wp) .OR. ANY(r2/=0.0_wp)) &
            ERROR STOP 'Within-cell interlayer path touched inactive cells'

        CALL upper_hp%finalize
        CALL upper%finalize
        CALL lower%finalize
        WRITE(*,*) 'PASS: within-cell interlayer pressure residual and active mask'
    END SUBROUTINE check_interlayer_cell_assembly

    !> @brief Verify identical-material total-depth bounds and CFL equivalence.
    !> @details Reuses a physically converted wet primitive state, including
    !> its temperature and composition. Tests both directions and dry limits.
    SUBROUTINE check_total_depth_cfl_equivalence
        REAL(wp) :: total(n_vars),lower(n_vars),upper(n_vars),dry(n_vars)
        REAL(wp) :: alo,ahi,reflo,refhi
        REAL(wp) :: xn(comp_interfaces_x,comp_cells_y),xp(comp_interfaces_x,comp_cells_y)
        REAL(wp) :: yn(comp_cells_x,comp_interfaces_y),yp(comp_cells_x,comp_interfaces_y)
        REAL(wp) :: step,expected,cfl_limit,rate,xspeed,yspeed
        REAL(wp) :: xlo,xhi,ylo,yhi,face_lo,face_hi
        TYPE(reconstruction_workspace_type) :: rec1,rec2
        TYPE(hyperbolic_workspace_type) :: speed_cache
        CHARACTER(LEN=LEN(cfl_policy)) :: old_policy
        INTEGER :: j,k,jwet,kwet,dir,case_number,normal

        jwet=0
        kwet=0
        DO k=1,ny
            DO j=1,nx
                IF (jwet==0 .AND. qp(1,j,k)>10.0_wp*dry_thickness_tolerance) THEN
                    jwet=j
                    kwet=k
                END IF
            END DO
        END DO
        IF (jwet==0) RETURN

        total=qp(:,jwet,kwet)
        dry=0.0_wp
        lower=total
        upper=total
        lower(1)=0.5_wp*total(1)
        upper(1)=total(1)-lower(1)

        DO dir=1,2
            normal=idx_u
            IF (dir==2) normal=idx_v

            DO case_number=1,4
                SELECT CASE(case_number)
                CASE(1)
                    CONTINUE
                CASE(2)
                    lower=total
                    upper=dry
                CASE(3)
                    lower=dry
                    upper=total
                CASE(4)
                    lower=dry
                    upper=dry
                END SELECT

                CALL eval_total_depth_face_speed_bounds(lower,upper,lower,upper, &
                    1.0_wp,normal,alo,ahi)
                IF (case_number==4) THEN
                    IF (alo/=0.0_wp .OR. ahi/=0.0_wp) &
                        ERROR STOP 'Dry two-layer face generated wave speeds'
                ELSE
                    CALL eval_local_speeds_direction(total,1.0_wp,normal,reflo,refhi)
                    CALL same_bits('total-depth negative wave bound', &
                        [alo],[MIN(0.0_wp,reflo)])
                    CALL same_bits('total-depth positive wave bound', &
                        [ahi],[MAX(0.0_wp,refhi)])
                END IF
            END DO

            lower=total
            upper=total
            lower(1)=0.5_wp*total(1)
            upper(1)=total(1)-lower(1)
        END DO

        CALL eval_total_depth_face_speed_bounds(lower,upper,lower,upper, &
            1.0_wp,idx_u,xlo,xhi)
        CALL eval_total_depth_face_speed_bounds(lower,upper,lower,upper, &
            1.0_wp,idx_v,ylo,yhi)
        ! Verify the face pass reuses the canonical workspace, not a shadow array.
        CALL rec1%initialize(layer=1)
        CALL rec2%initialize(layer=2)
        CALL speed_cache%initialize
        IF (nx>1) THEN
            rec1%qp_adv_interfaceL(:,2,1)=lower
            rec1%qp_adv_interfaceR(:,2,1)=lower
            rec2%qp_adv_interfaceL(:,2,1)=upper
            rec2%qp_adv_interfaceR(:,2,1)=upper
            rec1%h_selected_interfaceL(2,1)=lower(1)
            rec1%h_selected_interfaceR(2,1)=lower(1)
            rec2%h_selected_interfaceL(2,1)=upper(1)
            rec2%h_selected_interfaceR(2,1)=upper(1)
            CALL speed_cache%evaluate_total_depth_speeds(rec1,rec2,1,[2],[1], &
                0,[INTEGER::],[INTEGER::])
            CALL eval_total_depth_face_speed_bounds(lower,upper,lower,upper, &
                grav_coeff_stag_x(2,1),idx_u,face_lo,face_hi)
            CALL same_bits('shared x-face total-depth bounds', &
                [speed_cache%a_interface_xNeg(2,1),speed_cache%a_interface_xPos(2,1)], &
                [face_lo,face_hi])
        END IF
        IF (ny>1) THEN
            rec1%qp_adv_interfaceB(:,1,2)=lower
            rec1%qp_adv_interfaceT(:,1,2)=lower
            rec2%qp_adv_interfaceB(:,1,2)=upper
            rec2%qp_adv_interfaceT(:,1,2)=upper
            rec1%h_selected_interfaceB(1,2)=lower(1)
            rec1%h_selected_interfaceT(1,2)=lower(1)
            rec2%h_selected_interfaceB(1,2)=upper(1)
            rec2%h_selected_interfaceT(1,2)=upper(1)
            CALL speed_cache%evaluate_total_depth_speeds(rec1,rec2,0,[INTEGER::],[INTEGER::], &
                1,[1],[2])
            CALL eval_total_depth_face_speed_bounds(lower,upper,lower,upper, &
                grav_coeff_stag_y(1,2),idx_v,face_lo,face_hi)
            CALL same_bits('shared y-face total-depth bounds', &
                [speed_cache%a_interface_yNeg(1,2),speed_cache%a_interface_yPos(1,2)], &
                [face_lo,face_hi])
        END IF
        CALL speed_cache%finalize
        CALL rec2%finalize
        CALL rec1%finalize

        xn=xlo
        xp=xhi
        yn=ylo
        yp=yhi
        xspeed=MAX(-xlo,xhi)
        yspeed=MAX(-ylo,yhi)
        old_policy=cfl_policy
        cfl_limit=MIN(cfl,0.25_wp)

        DO case_number=1,2
            expected=max_dt
            cfl_policy='directional_min'
            IF (case_number==2) cfl_policy='combined'
            IF (case_number==1) THEN
                IF (nx>1 .AND. xspeed>0.0_wp) expected=MIN(expected,cfl_limit*dx/xspeed)
                IF (ny>1 .AND. yspeed>0.0_wp) expected=MIN(expected,cfl_limit*dy/yspeed)
            ELSEIF (ny==1) THEN
                IF (nx>1 .AND. xspeed>0.0_wp) expected=MIN(expected,cfl_limit*dx/xspeed)
            ELSEIF (nx==1) THEN
                IF (yspeed>0.0_wp) expected=MIN(expected,cfl_limit*dy/yspeed)
            ELSE
                rate=xspeed/dx+yspeed/dy
                IF (rate>0.0_wp) expected=MIN(expected,cfl_limit/rate)
            END IF

            CALL spatial%compute_coupled_timestep(domain,xn,xp,yn,yp,step)
            CALL same_bits('homogeneous coupled CFL equals total-depth CFL',[step],[expected])
        END DO
        cfl_policy=old_policy
        WRITE(*,*) 'PASS: total-depth wave bounds and homogeneous CFL equivalence'
    END SUBROUTINE check_total_depth_cfl_equivalence

    !> @brief Check shared coupled-face CFL reduction against both policy formulas.
    !> @details This tests face-rate reduction only, not coupled eigenvalues.
    !> The same suite runs 2D, x-only, y-only and single-cell fixtures.
    SUBROUTINE check_coupled_cfl_policy
        REAL(wp) :: axneg(comp_interfaces_x,comp_cells_y)
        REAL(wp) :: axpos(comp_interfaces_x,comp_cells_y)
        REAL(wp) :: ayneg(comp_cells_x,comp_interfaces_y)
        REAL(wp) :: aypos(comp_cells_x,comp_interfaces_y)
        REAL(wp) :: proposed,reference,lambda_x,lambda_y,rate
        CHARACTER(LEN=LEN(cfl_policy)) :: saved_policy
        REAL(wp) :: saved_cap
        INTEGER :: policy_case

        saved_policy=cfl_policy
        saved_cap=max_dt
        axneg=-2.0_wp
        axpos=3.0_wp
        ayneg=-4.0_wp
        aypos=5.0_wp
        lambda_x=0.0_wp
        lambda_y=0.0_wp
        IF (nx>1) lambda_x=3.0_wp
        IF (ny>1) lambda_y=5.0_wp

        DO policy_case=1,2
            IF (policy_case==1) THEN
                cfl_policy='directional_min'
                reference=max_dt
                IF (lambda_x>0.0_wp) reference=MIN(reference,MIN(cfl,0.25_wp)*dx/lambda_x)
                IF (lambda_y>0.0_wp) reference=MIN(reference,MIN(cfl,0.25_wp)*dy/lambda_y)
            ELSE
                cfl_policy='combined'
                reference=max_dt
                IF (ny==1) THEN
                    IF (nx>1) reference=MIN(reference,MIN(cfl,0.25_wp)*dx/lambda_x)
                ELSEIF (nx==1) THEN
                    reference=MIN(reference,MIN(cfl,0.25_wp)*dy/lambda_y)
                ELSE
                    rate=lambda_x/dx+lambda_y/dy
                    IF (rate>0.0_wp) reference=MIN(reference,MIN(cfl,0.25_wp)/rate)
                END IF
            END IF

            CALL spatial%compute_coupled_timestep(domain,axneg,axpos,ayneg,aypos,proposed)
            CALL same_bits('coupled-speed CFL policy',[proposed],[reference])
        END DO

        cfl_policy=saved_policy
        max_dt=saved_cap
        WRITE(*,*) 'PASS: externally supplied coupled speed bounds use canonical CFL'
    END SUBROUTINE check_coupled_cfl_policy

    !> @brief Validate the shared workset as the union of independent L1/L2 masks.
    !> @details Checks active cells, faces and HP reconstruction halo without
    !> allocating any persistent two-layer domain or enabling coupled physics.
    SUBROUTINE check_two_layer_union_workset
        TYPE(domain_type) :: lower,upper,coupled
        REAL(wp) :: lower_q(n_vars,nx,ny),upper_q(n_vars,nx,ny)
        REAL(wp) :: lower_q_before(n_vars,nx,ny),upper_q_before(n_vars,nx,ny)

        ! One-dimensional and single-cell fixtures have no two interior cells.
        IF (nx<3 .OR. ny<3) RETURN

        lower_q=0.0_wp
        upper_q=0.0_wp
        lower_q(idx_totMassEqn,2,2)=1.0_wp
        upper_q(idx_totMassEqn,nx-1,ny-1)=2.0_wp
        lower_q_before=lower_q
        upper_q_before=upper_q

        CALL lower%initialize
        CALL upper%initialize
        CALL coupled%initialize
        CALL lower%check_solve(lower_q,0.0_wp,.FALSE.)
        CALL upper%check_solve(upper_q,0.0_wp,.FALSE.)
        CALL coupled%check_solve(lower_q,0.0_wp,.FALSE.,q_other=upper_q)

        IF (ANY(coupled%solve_mask .NEQV. (lower%solve_mask .OR. upper%solve_mask))) &
            ERROR STOP 'Coupled active cells differ from layer union'
        IF (ANY(coupled%solve_mask_x .NEQV. (lower%solve_mask_x .OR. upper%solve_mask_x))) &
            ERROR STOP 'Coupled x faces differ from layer union'
        IF (ANY(coupled%solve_mask_y .NEQV. (lower%solve_mask_y .OR. upper%solve_mask_y))) &
            ERROR STOP 'Coupled y faces differ from layer union'
        IF (ANY(coupled%hp_eta_mask .NEQV. (lower%hp_eta_mask .OR. upper%hp_eta_mask))) &
            ERROR STOP 'Coupled HP cells differ from layer union'
        IF (ANY(coupled%reconstruction_mask .NEQV. &
            (lower%reconstruction_mask .OR. upper%reconstruction_mask))) &
            ERROR STOP 'Coupled reconstruction halo differs from layer union'
        IF (coupled%solve_cells/=COUNT(coupled%solve_mask) .OR. &
            coupled%hp_eta_cells/=COUNT(coupled%hp_eta_mask) .OR. &
            coupled%solve_interfaces_x/=COUNT(coupled%solve_mask_x) .OR. &
            coupled%solve_interfaces_y/=COUNT(coupled%solve_mask_y)) &
            ERROR STOP 'Coupled workset contains duplicate or missing entries'
        IF (ANY(lower_q/=lower_q_before) .OR. ANY(upper_q/=upper_q_before)) &
            ERROR STOP 'Coupled workset modified layer state'

        CALL lower%finalize
        CALL upper%finalize
        CALL coupled%finalize
        WRITE(*,*) 'PASS: two-layer union workset, HP halo, faces and read-only states'
    END SUBROUTINE check_two_layer_union_workset

    !> Verify the common interlayer pressure quadrature and its conservative sum.
    !> This exercises the path algebra only; the full coupled operator remains gated.
    SUBROUTINE check_interlayer_pressure_path
        REAL(wp) :: path_l1(5),path_l2(5),reverse_l1(5),reverse_l2(5)
        REAL(wp) :: h_left,h_right,load_left,load_right,metric_left,metric_right
        REAL(wp) :: expected,scale,tolerance
        INTEGER :: direction,normal_momentum,case_number

        DO direction=1,2
            normal_momentum=direction+1
            DO case_number=1,4
                SELECT CASE(case_number)
                CASE(1)
                    ! Vary every field, including G across a cell-internal path.
                    h_left=2.0_wp; h_right=3.0_wp
                    load_left=5.0_wp; load_right=7.0_wp
                    metric_left=0.7_wp; metric_right=0.9_wp
                CASE(2)
                    ! Shared metric on an interface, load gradient only.
                    h_left=2.0_wp; h_right=2.0_wp
                    load_left=5.0_wp; load_right=7.0_wp
                    metric_left=0.8_wp; metric_right=0.8_wp
                CASE(3)
                    ! Upper layer disappears at the right endpoint.
                    h_left=2.0_wp; h_right=1.0_wp
                    load_left=5.0_wp; load_right=0.0_wp
                    metric_left=0.6_wp; metric_right=0.7_wp
                CASE(4)
                    ! Zero interlayer load: no pressure force in either layer.
                    h_left=2.0_wp; h_right=3.0_wp
                    load_left=0.0_wp; load_right=0.0_wp
                    metric_left=0.7_wp; metric_right=0.9_wp
                END SELECT

                CALL eval_interlayer_pressure_path(direction,h_left,load_left,metric_left, &
                    h_right,load_right,metric_right,path_l1,path_l2)
                expected=metric_left*h_left*load_left-metric_right*h_right*load_right
                scale=MAX(1.0_wp,ABS(metric_left*h_left*load_left),ABS(metric_right*h_right*load_right))
                tolerance=128.0_wp*EPSILON(1.0_wp)*scale
                IF (ABS(path_l1(normal_momentum)+path_l2(normal_momentum)-expected)>tolerance) &
                    ERROR STOP 'Interlayer pressure path violates product balance'
                IF (ANY(path_l1(1:normal_momentum-1)/=0.0_wp) .OR. &
                    ANY(path_l2(1:normal_momentum-1)/=0.0_wp) .OR. &
                    ANY(path_l1(normal_momentum+1:)/=0.0_wp) .OR. &
                    ANY(path_l2(normal_momentum+1:)/=0.0_wp)) &
                    ERROR STOP 'Interlayer path introduced a non-momentum source'

                CALL eval_interlayer_pressure_path(direction,h_right,load_right,metric_right, &
                    h_left,load_left,metric_left,reverse_l1,reverse_l2)
                IF (ABS(path_l1(normal_momentum)+reverse_l1(normal_momentum))>tolerance .OR. &
                    ABS(path_l2(normal_momentum)+reverse_l2(normal_momentum))>tolerance) &
                    ERROR STOP 'Interlayer pressure path is not antisymmetric'
                IF (case_number==4 .AND. (ANY(path_l1/=0.0_wp) .OR. ANY(path_l2/=0.0_wp))) &
                    ERROR STOP 'Zero upper-layer load created interlayer pressure'
            END DO
        END DO
        WRITE(*,*) 'PASS: interlayer pressure path product balance, dry limit and reversal'
    END SUBROUTINE check_interlayer_pressure_path

    !> @brief Verify that using eta2=B+h1+h2 does not count -G*L*dh1 twice.
    !> @details With constant Gamma2, both the standard hydrostatic path and
    !> the interlayer path interpolate the same L=Gamma2*h2. Their difference
    !> must equal the embedded lower-surface contribution; only -h1*L*dG
    !> remains to be added to the L2 momentum equation.
    SUBROUTINE check_interlayer_metric_no_double_count
        REAL(wp) :: full_l1(5),full_l2(5),metric_l2(5),reverse_metric(5)
        REAL(wp) :: hydro_total(5),hydro_without_h1(5)
        REAL(wp) :: h1_l,h1_r,h2_l,h2_r,load_l,load_r,g_l,g_r,gamma
        REAL(wp) :: eta_l,eta_r,expected,tolerance,scale,embedded
        INTEGER :: direction,normal_momentum,case_number

        gamma=4.25_wp
        DO direction=1,2
            normal_momentum=direction+1
            DO case_number=1,4
                h1_l=1.4_wp
                h1_r=2.0_wp
                h2_l=0.7_wp
                h2_r=1.1_wp
                g_l=0.72_wp
                g_r=0.91_wp
                SELECT CASE(case_number)
                CASE(1)
                    ! Both thicknesses and terrain metric vary.
                    CONTINUE
                CASE(2)
                    ! A shared face metric has no explicit L2 correction.
                    g_r=g_l
                CASE(3)
                    ! Dry upper endpoint and unequal lower-layer thickness.
                    h2_r=0.0_wp
                CASE(4)
                    ! Entire upper layer dry.
                    h2_l=0.0_wp
                    h2_r=0.0_wp
                END SELECT
                load_l=gamma*h2_l
                load_r=gamma*h2_r
                CALL eval_interlayer_pressure_path(direction,h1_l,load_l,g_l, &
                    h1_r,load_r,g_r,full_l1,full_l2,metric_l2)
                ! Integral of -h1(s)*L(s)*dG with linear h1,L,G.
                expected=-(g_r-g_l)*(h1_l*load_l &
                    +0.5_wp*(h1_l*(load_r-load_l)+load_l*(h1_r-h1_l)) &
                    +(h1_r-h1_l)*(load_r-load_l)/3.0_wp)
                scale=MAX(1.0_wp,ABS(g_l*h1_l*load_l),ABS(g_r*h1_r*load_r))
                tolerance=128.0_wp*EPSILON(1.0_wp)*scale
                IF (ABS(metric_l2(normal_momentum)-expected)>tolerance) &
                    ERROR STOP 'Metric-only L2 path has incorrect quadrature'
                IF (ANY(metric_l2(:normal_momentum-1)/=0.0_wp) .OR. &
                    ANY(metric_l2(normal_momentum+1:)/=0.0_wp)) &
                    ERROR STOP 'Metric-only L2 path contaminated other components'
                IF (case_number==2 .AND. ANY(metric_l2/=0.0_wp)) &
                    ERROR STOP 'Shared metric created duplicate L2 face pressure'

                ! B is constant in this identity; total eta includes h1.
                eta_l=h1_l+h2_l
                eta_r=h1_r+h2_r
                CALL eval_hydrostatic_path(direction,h2_l,gamma,eta_l,g_l, &
                    h2_r,gamma,eta_r,g_r,hydro_total)
                CALL eval_hydrostatic_path(direction,h2_l,gamma,eta_l-h1_l,g_l, &
                    h2_r,gamma,eta_r-h1_r,g_r,hydro_without_h1)
                embedded=hydro_total(normal_momentum)-hydro_without_h1(normal_momentum)
                IF (ABS(embedded-(full_l2(normal_momentum)-metric_l2(normal_momentum)))>tolerance) &
                    ERROR STOP 'L2 hydrostatic path already contains inconsistent interlayer pressure'

                CALL eval_interlayer_pressure_path(direction,h1_r,load_r,g_r, &
                    h1_l,load_l,g_l,full_l1,full_l2,reverse_metric)
                IF (ABS(metric_l2(normal_momentum)+reverse_metric(normal_momentum))>tolerance) &
                    ERROR STOP 'Metric-only L2 pressure is not path-antisymmetric'
            END DO
        END DO
        WRITE(*,*) 'PASS: L2 embedded h1 pressure and metric-only interlayer correction'
    END SUBROUTINE check_interlayer_metric_no_double_count

    !> @brief Verify the exact discrete overlap for variable upper pressure coefficient.
    !> @details The interlayer path is straight in L=Gamma*h2, whereas the
    !> existing hydrostatic path is straight in Gamma and h2 separately.
    SUBROUTINE check_interlayer_variable_gamma_overlap
        REAL(wp) :: path_l1(5),correction_l2(5),full_l1(5),full_l2(5)
        REAL(wp) :: hydro_total(5),hydro_without_h1(5)
        REAL(wp) :: h1_l,h1_r,h2_l,h2_r,gamma_l,gamma_r,eta_l,eta_r,g_l,g_r
        REAL(wp) :: balance,scale,tolerance
        INTEGER :: direction,normal_momentum,case_number

        h1_l=1.3_wp
        h1_r=1.9_wp
        h2_l=0.6_wp
        gamma_l=4.0_wp
        gamma_r=7.0_wp
        g_l=0.7_wp
        g_r=0.9_wp
        DO direction=1,2
            normal_momentum=direction+1
            DO case_number=1,3
                h2_r=1.2_wp
                IF (case_number==2) h2_r=0.0_wp
                IF (case_number==3) THEN
                    h2_r=0.0_wp
                    h2_l=0.0_wp
                END IF
                eta_l=0.4_wp+h1_l+h2_l
                eta_r=0.4_wp+h1_r+h2_r
                CALL eval_interlayer_pressure_correction(direction, &
                    h1_l,h2_l,gamma_l,eta_l,g_l, &
                    h1_r,h2_r,gamma_r,eta_r,g_r,path_l1,correction_l2)
                CALL eval_interlayer_pressure_path(direction, &
                    h1_l,gamma_l*h2_l,g_l,h1_r,gamma_r*h2_r,g_r,full_l1,full_l2)
                CALL eval_hydrostatic_path(direction,h2_l,gamma_l,eta_l,g_l, &
                    h2_r,gamma_r,eta_r,g_r,hydro_total)
                CALL eval_hydrostatic_path(direction,h2_l,gamma_l,eta_l-h1_l,g_l, &
                    h2_r,gamma_r,eta_r-h1_r,g_r,hydro_without_h1)
                scale=MAX(1.0_wp,ABS(g_l*h1_l*gamma_l*h2_l),ABS(g_r*h1_r*gamma_r*h2_r))
                tolerance=128.0_wp*EPSILON(1.0_wp)*scale
                IF (ANY(ABS(path_l1-full_l1)>tolerance)) &
                    ERROR STOP 'Upper pressure correction modified L1 path'
                IF (ANY(ABS(correction_l2+hydro_total-hydro_without_h1-full_l2)>tolerance)) &
                    ERROR STOP 'Variable Gamma L2 pressure was counted twice'
                balance=path_l1(normal_momentum)+full_l2(normal_momentum) &
                    +g_r*h1_r*gamma_r*h2_r-g_l*h1_l*gamma_l*h2_l
                IF (ABS(balance)>tolerance) &
                    ERROR STOP 'Variable Gamma interlayer pressure product is unbalanced'
            END DO
            h2_l=0.6_wp
        END DO
        WRITE(*,*) 'PASS: variable Gamma L2 pressure overlap and dry upper layer'
    END SUBROUTINE check_interlayer_variable_gamma_overlap

    !> @brief Check interlayer path splitting into the two oriented PCCU face values.
    !> @details The numerical flux is set to zero to isolate the nonconservative
    !> contribution. The production oriented split must reconstruct the same
    !> left-to-right jump for each layer, including variable G and a dry L2 side.
    SUBROUTINE check_interlayer_oriented_face_pair
        REAL(wp) :: path_l1(5),path_l2(5)
        REAL(wp) :: left_l1(5),right_l1(5),left_l2(5),right_l2(5)
        REAL(wp) :: zero_flux(5),expected,tolerance,scale
        REAL(wp) :: a_minus,a_plus,h1_left,h1_right,load_left,load_right
        REAL(wp) :: metric_left,metric_right
        INTEGER :: direction,case_number,speed_case,normal_momentum

        zero_flux=0.0_wp
        DO direction=1,2
            normal_momentum=direction+1
            DO case_number=1,3
                SELECT CASE(case_number)
                CASE(1)
                    ! Shared metric at a face, unequal upper pressure.
                    h1_left=1.8_wp; h1_right=1.2_wp
                    load_left=7.0_wp; load_right=3.0_wp
                    metric_left=0.8_wp; metric_right=0.8_wp
                CASE(2)
                    ! Variable metric along a cell-internal quadrature path.
                    h1_left=1.8_wp; h1_right=1.2_wp
                    load_left=7.0_wp; load_right=3.0_wp
                    metric_left=0.6_wp; metric_right=0.9_wp
                CASE(3)
                    ! Wet L1 with a vanishing upper layer.
                    h1_left=1.8_wp; h1_right=1.2_wp
                    load_left=0.0_wp; load_right=0.0_wp
                    metric_left=0.6_wp; metric_right=0.9_wp
                END SELECT
                CALL eval_interlayer_pressure_path(direction,h1_left,load_left,metric_left, &
                    h1_right,load_right,metric_right,path_l1,path_l2)
                expected=metric_left*h1_left*load_left-metric_right*h1_right*load_right
                scale=MAX(1.0_wp,ABS(metric_left*h1_left*load_left),ABS(metric_right*h1_right*load_right))
                tolerance=128.0_wp*EPSILON(1.0_wp)*scale

                DO speed_case=1,3
                    SELECT CASE(speed_case)
                    CASE(1)
                        a_minus=-2.0_wp; a_plus=3.0_wp
                    CASE(2)
                        a_minus=0.0_wp; a_plus=3.0_wp
                    CASE(3)
                        a_minus=-2.0_wp; a_plus=0.0_wp
                    END SELECT
                    CALL eval_oriented_pccu_pair(zero_flux,path_l1,a_minus,a_plus,left_l1,right_l1)
                    CALL eval_oriented_pccu_pair(zero_flux,path_l2,a_minus,a_plus,left_l2,right_l2)
                    IF (ABS(right_l1(normal_momentum)-left_l1(normal_momentum) &
                        -path_l1(normal_momentum))>tolerance) &
                        ERROR STOP 'Oriented interlayer L1 face jump differs from path'
                    IF (ABS(right_l2(normal_momentum)-left_l2(normal_momentum) &
                        -path_l2(normal_momentum))>tolerance) &
                        ERROR STOP 'Oriented interlayer L2 face jump differs from path'
                    IF (ABS(right_l1(normal_momentum)+right_l2(normal_momentum) &
                        -left_l1(normal_momentum)-left_l2(normal_momentum)-expected)>tolerance) &
                        ERROR STOP 'Oriented interlayer face pair violates pressure product balance'
                    IF (ANY(left_l1(1:normal_momentum-1)/=0.0_wp) .OR. &
                        ANY(right_l1(1:normal_momentum-1)/=0.0_wp) .OR. &
                        ANY(left_l2(1:normal_momentum-1)/=0.0_wp) .OR. &
                        ANY(right_l2(1:normal_momentum-1)/=0.0_wp)) &
                        ERROR STOP 'Oriented interlayer face pair created a non-momentum flux'
                    IF (case_number==3) THEN
                        IF (ANY(left_l1/=0.0_wp) .OR. ANY(right_l1/=0.0_wp) .OR. &
                            ANY(left_l2/=0.0_wp) .OR. ANY(right_l2/=0.0_wp)) &
                            ERROR STOP 'Dry L2 created an oriented interlayer face force'
                    END IF
                END DO
            END DO
        END DO
        WRITE(*,*) 'PASS: oriented interlayer face pair, one-sided speeds and dry upper layer'
    END SUBROUTINE check_interlayer_oriented_face_pair

    !> @brief Check a stationary, pressure-free face with zero coupled speeds.
    !> @details No speed floor is introduced and all four oriented values
    !> must vanish. The nondegenerate-pressure case remains guarded.
    SUBROUTINE check_degenerate_interlayer_face
        REAL(wp) :: l1left(n_eqns),l1right(n_eqns)
        REAL(wp) :: l2left(n_eqns),l2right(n_eqns)

        CALL eval_interlayer_oriented_face(1, &
            0.0_wp,0.0_wp,0.0_wp,0.0_wp, &
            0.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp, &
            0.0_wp,0.0_wp,l1left,l1right,l2left,l2right)

        IF (ANY(l1left/=0.0_wp) .OR. ANY(l1right/=0.0_wp) .OR. &
            ANY(l2left/=0.0_wp) .OR. ANY(l2right/=0.0_wp)) &
            ERROR STOP 'Stationary pressure-free interlayer face is not zero'

        WRITE(*,*) 'PASS: stationary pressure-free interlayer face'
    END SUBROUTINE check_degenerate_interlayer_face

    !> @brief Validate oriented coupled pressure corrections with supplied common speeds.
    SUBROUTINE check_interlayer_coupled_face
        REAL(wp) :: path_l1(5),correction_l2(5),left_l1(5),right_l1(5)
        REAL(wp) :: left_l2(5),right_l2(5),metric,eta_l,eta_r
        REAL(wp) :: tol,scale,a_minus,a_plus
        INTEGER :: direction,case_number,normal

        metric=0.78_wp
        eta_l=2.3_wp
        eta_r=2.7_wp
        DO direction=1,2
            normal=direction+1
            DO case_number=1,3
                SELECT CASE(case_number)
                CASE(1)
                    a_minus=-2.0_wp
                    a_plus=3.0_wp
                CASE(2)
                    a_minus=0.0_wp
                    a_plus=3.0_wp
                CASE(3)
                    a_minus=-2.0_wp
                    a_plus=0.0_wp
                END SELECT
                CALL eval_interlayer_oriented_face(direction, &
                    1.2_wp,1.1_wp,4.0_wp,eta_l,1.8_wp,0.9_wp,6.0_wp,eta_r,metric, &
                    a_minus,a_plus,left_l1,right_l1,left_l2,right_l2)
                CALL eval_interlayer_pressure_correction(direction, &
                    1.2_wp,1.1_wp,4.0_wp,eta_l,metric, &
                    1.8_wp,0.9_wp,6.0_wp,eta_r,metric,path_l1,correction_l2)
                scale=MAX(1.0_wp,ABS(path_l1(normal)),ABS(correction_l2(normal)))
                tol=128.0_wp*EPSILON(1.0_wp)*scale
                IF (ANY(ABS(right_l1-left_l1-path_l1)>tol) .OR. &
                    ANY(ABS(right_l2-left_l2-correction_l2)>tol)) &
                    ERROR STOP 'Interlayer face correction differs from path jump'
                IF (case_number==2) THEN
                    IF (ANY(left_l1/=0.0_wp) .OR. ANY(left_l2/=0.0_wp)) &
                        ERROR STOP 'One-sided interlayer face has wrong negative-side split'
                ELSEIF (case_number==3) THEN
                    IF (ANY(right_l1/=0.0_wp) .OR. ANY(right_l2/=0.0_wp)) &
                        ERROR STOP 'One-sided interlayer face has wrong positive-side split'
                END IF
                IF (ANY(left_l1(:normal-1)/=0.0_wp) .OR. ANY(right_l1(normal+1:)/=0.0_wp) .OR. &
                    ANY(left_l2(:normal-1)/=0.0_wp) .OR. ANY(right_l2(normal+1:)/=0.0_wp)) &
                    ERROR STOP 'Interlayer face pressure contaminated non-momentum components'
            END DO
        END DO
        WRITE(*,*) 'PASS: shared-speed interlayer oriented face with embedded L2 pressure'
    END SUBROUTINE check_interlayer_coupled_face

    !> @brief Check the HP-PCCU cell assembly sign for interlayer pressure.
    !> @details Both physical faces have identical traces on their two sides,
    !> so the extra face paths vanish. The cell pressure contribution is the
    !> complete interlayer quadrature with the L2 hydrostatic overlap removed.
    !> Checks all four local layer-presence cases and a variable Gamma2 path.
    SUBROUTINE check_interlayer_cell_residual
        REAL(wp) :: h1_w,h1_e,h2_w,h2_e,gamma_w,gamma_e
        REAL(wp) :: eta_w,eta_e,g_w,g_e,path_l1(5),path_l2(5)
        REAL(wp) :: full_l1(5),full_l2(5),hydro_total(5),hydro_without_h1(5)
        REAL(wp) :: residual_l1(5),residual_l2(5),reference(5),tol,scale
        INTEGER :: direction,case_number,normal

        DO direction=1,2
            normal=direction+1
            DO case_number=1,4
                h1_w=1.6_wp
                h1_e=1.9_wp
                h2_w=0.5_wp
                h2_e=0.8_wp
                gamma_w=4.0_wp
                gamma_e=6.0_wp
                g_w=0.72_wp
                g_e=0.89_wp
                SELECT CASE(case_number)
                CASE(1)
                    ! Both layers present with varying material and metric.
                    CONTINUE
                CASE(2)
                    ! Only lower layer present: no interlayer pressure.
                    h2_w=0.0_wp
                    h2_e=0.0_wp
                CASE(3)
                    ! Only upper layer present: h1=0 on both sides.
                    h1_w=0.0_wp
                    h1_e=0.0_wp
                CASE(4)
                    ! Both layers dry.
                    h1_w=0.0_wp
                    h1_e=0.0_wp
                    h2_w=0.0_wp
                    h2_e=0.0_wp
                END SELECT

                eta_w=0.3_wp+h1_w+h2_w
                eta_e=0.3_wp+h1_e+h2_e
                CALL eval_interlayer_pressure_correction(direction, &
                    h1_w,h2_w,gamma_w,eta_w,g_w, &
                    h1_e,h2_e,gamma_e,eta_e,g_e,path_l1,path_l2)
                CALL eval_interlayer_pressure_path(direction, &
                    h1_w,gamma_w*h2_w,g_w,h1_e,gamma_e*h2_e,g_e,full_l1,full_l2)
                CALL eval_hydrostatic_path(direction,h2_w,gamma_w,eta_w,g_w, &
                    h2_e,gamma_e,eta_e,g_e,hydro_total)
                CALL eval_hydrostatic_path(direction,h2_w,gamma_w,eta_w-h1_w,g_w, &
                    h2_e,gamma_e,eta_e-h1_e,g_e,hydro_without_h1)

                ! In the production divergence the within-cell path is subtracted.
                ! Equal physical states at both faces give zero face path jumps.
                residual_l1=-path_l1
                residual_l2=-path_l2
                reference=-full_l2+hydro_total-hydro_without_h1
                scale=MAX(1.0_wp,ABS(g_w*h1_w*gamma_w*h2_w),ABS(g_e*h1_e*gamma_e*h2_e))
                tol=128.0_wp*EPSILON(1.0_wp)*scale
                IF (ANY(ABS(residual_l1+full_l1)>tol)) &
                    ERROR STOP 'Lower-layer cell pressure has wrong sign'
                IF (ANY(ABS(residual_l2-reference)>tol)) &
                    ERROR STOP 'Upper-layer cell pressure was counted twice'
                IF (case_number>1) THEN
                    IF (ANY(residual_l1/=0.0_wp) .OR. ANY(residual_l2/=0.0_wp)) &
                        ERROR STOP 'Absent interlayer contact created cell pressure'
                END IF
                IF (ANY(residual_l1(:normal-1)/=0.0_wp) .OR. &
                    ANY(residual_l1(normal+1:)/=0.0_wp) .OR. &
                    ANY(residual_l2(:normal-1)/=0.0_wp) .OR. &
                    ANY(residual_l2(normal+1:)/=0.0_wp)) &
                    ERROR STOP 'Interlayer cell pressure contaminated transported variables'
            END DO
        END DO
        WRITE(*,*) 'PASS: coupled cell path residual signs and four layer-presence cases'
    END SUBROUTINE check_interlayer_cell_residual

    !> @brief Verify inactive optional components remain zero in the production PCCU algebra.
    !> @details Both layer vectors use the same seven component positions. L1 carries only
    !>          the pore component (7); L2 carries only the stochastic component (6).
    !>          Test central-upwind, factored advection and oriented face values with
    !>          two-sided and one-sided wave bounds. This is a local face-algebra test,
    !>          not a coupled reconstruction or boundary-condition integration test.
    SUBROUTINE check_optional_pccu_flux_slots
        REAL(wp) :: q_left(7),q_right(7),f_left(7),f_right(7)
        REAL(wp) :: flux(7),left_value(7),right_value(7),path(7)
        REAL(wp) :: a_minus,a_plus,u_left,u_right
        INTEGER :: layer,case_number,active_index,inactive_index

        u_left=0.5_wp
        u_right=0.8_wp
        path=0.0_wp

        DO layer=1,2
            IF (layer==1) THEN
                active_index=7
                inactive_index=6
            ELSE
                active_index=6
                inactive_index=7
            END IF

            q_left=0.0_wp
            q_right=0.0_wp
            q_left(1)=2.0_wp
            q_right(1)=1.5_wp
            q_left(active_index)=0.4_wp
            q_right(active_index)=0.6_wp
            f_left=u_left*q_left
            f_right=u_right*q_right

            DO case_number=1,3
                SELECT CASE(case_number)
                CASE(1)
                    a_minus=-0.4_wp
                    a_plus=1.3_wp
                CASE(2)
                    a_minus=0.0_wp
                    a_plus=1.3_wp
                CASE(3)
                    a_minus=-0.4_wp
                    a_plus=0.0_wp
                END SELECT

                CALL eval_pccu_advective_flux(a_minus,a_plus,u_left,u_right, &
                    q_left,q_right,f_left,f_right,flux)
                IF (ABS(flux(inactive_index))>0.0_wp) &
                    ERROR STOP 'PCCU advection generated an inactive optional component'
                IF (ABS(flux(active_index))<=0.0_wp) &
                    ERROR STOP 'PCCU advection lost an active optional component'

                CALL eval_central_upwind_flux(a_minus,a_plus,f_left,f_right, &
                    q_left,q_right,flux)
                IF (ABS(flux(inactive_index))>0.0_wp) &
                    ERROR STOP 'PCCU central-upwind flux generated an inactive optional component'
                IF (ABS(flux(active_index))<=0.0_wp) &
                    ERROR STOP 'PCCU central-upwind flux lost an active optional component'

                CALL eval_oriented_pccu_pair(flux,path,a_minus,a_plus,left_value,right_value)
                IF (ABS(left_value(inactive_index))>0.0_wp .OR. &
                    ABS(right_value(inactive_index))>0.0_wp) &
                    ERROR STOP 'PCCU oriented face generated an inactive optional component'
                IF (ABS(left_value(active_index))<=0.0_wp .OR. &
                    ABS(right_value(active_index))<=0.0_wp) &
                    ERROR STOP 'PCCU oriented face lost an active optional component'
            END DO
        END DO

        WRITE(*,*) 'PASS: PCCU inactive optional face slots remain exactly zero'
    END SUBROUTINE check_optional_pccu_flux_slots


    !> @brief Verify isolated L2 pressure, advective transport and conservative Euler update.
    !> @details A buoyant plume at rest must have zero residual on shared terrain.
    !>          A dense moving gas-solid layer must preserve total mass, solid mass
    !>          and energy through a CFL-limited spatial step with a quiet collar.
    !>          L1 is dry and unchanged; this does not activate coupled two-layer physics.
    SUBROUTINE check_isolated_gas_l2
        REAL(wp) :: primitive(n_vars),reduced_gravity,gamma,isolated_dt
        REAL(wp) :: initial_mass,initial_solids,initial_energy
        REAL(wp) :: final_mass,final_solids,final_energy,transport_norm
        LOGICAL :: conversion_ok
        INTEGER :: j,k,l

        gas_flag=.TRUE.
        liquid_flag=.FALSE.
        slot_q(:,1,:,:)=0.0_wp
        slot2_domain%solve_mask_time(1,1)=0.0_wp

        primitive=0.0_wp
        primitive(idx_h)=0.7_wp
        primitive(idx_T)=600.0_wp
        primitive(idx_solid_first)=0.01_wp
        CALL eval_hydrostatic_coefficient(primitive,reduced_gravity,gamma)
        IF (reduced_gravity>=0.0_wp .OR. gamma>=0.0_wp) &
            ERROR STOP 'Expected a positively buoyant dilute gas layer'

        ! Vary depth at rest: without the L2 pressure clamp even a flat-bed
        ! case would create a nonzero signed-gravity hydrostatic residual.
        DO k=1,ny
            DO j=1,nx
                primitive(idx_h)=0.5_wp+0.01_wp*REAL(j+k,wp)
                CALL primitive_to_conservative(primitive,slot_q(:,2,j,k),accepted=conversion_ok)
                IF (.NOT.conversion_ok) ERROR STOP 'Rejected buoyant L2 gas composition'
            END DO
        END DO

        CALL slot2_domain%check_solve(slot_q(:,2,:,:),0.0_wp,.FALSE.)
        slot2_qp=0.0_wp
        CALL slot2_spatial%compute_timestep(slot_q(:,2,:,:),slot2_qp,0.0_wp,isolated_dt,slot2_domain)
        IF (.NOT.ieee_is_finite(isolated_dt) .OR. isolated_dt<=0.0_wp) &
            ERROR STOP 'Invalid buoyant L2 timestep'
        CALL slot2_spatial%evaluate(slot2_qp,slot2_term,0.0_wp,slot2_domain)

        DO l=1,slot2_domain%solve_cells
            j=slot2_domain%j_cent(l)
            k=slot2_domain%k_cent(l)
            IF (ANY(slot2_term(:,j,k)/=0.0_wp)) &
                ERROR STOP 'Buoyant L2 created hydrostatic or transport forces at rest'
        END DO

        ! Reuse the same L2 operator and geometry for a dense, moving ash cloud.
        primitive(idx_T)=330.0_wp
        primitive(idx_solid_first)=0.30_wp
        CALL eval_hydrostatic_coefficient(primitive,reduced_gravity,gamma)
        IF (reduced_gravity<=0.0_wp .OR. gamma<=0.0_wp) &
            ERROR STOP 'Expected positive reduced gravity for dense gas-solid L2'

        ! A constant boundary collar has no normal transport: all component
        ! budgets should therefore be conserved by the isolated L2 spatial step.
        DO k=1,ny
            DO j=1,nx
                primitive(idx_h)=0.7_wp
                primitive(idx_u)=0.0_wp
                primitive(idx_v)=0.0_wp
                IF (j>3 .AND. j<nx-2 .AND. k>3 .AND. k<ny-2) THEN
                    primitive(idx_h)=0.7_wp+0.02_wp*SIN(REAL(j+k,wp))
                    primitive(idx_u)=0.10_wp*SIN(REAL(j,wp))
                    primitive(idx_v)=0.08_wp*COS(REAL(k,wp))
                END IF
                CALL primitive_to_conservative(primitive,slot_q(:,2,j,k),accepted=conversion_ok)
                IF (.NOT.conversion_ok) ERROR STOP 'Rejected dense L2 gas composition'
            END DO
        END DO

        initial_mass=SUM(slot_q(idx_totMassEqn,2,:,:))
        initial_solids=SUM(slot_q(idx_solidEqn_first,2,:,:))
        initial_energy=SUM(slot_q(idx_engyEqn,2,:,:))

        CALL slot2_domain%check_solve(slot_q(:,2,:,:),0.0_wp,.FALSE.)
        slot2_qp=0.0_wp
        CALL slot2_spatial%compute_timestep(slot_q(:,2,:,:),slot2_qp,0.0_wp,isolated_dt,slot2_domain)
        IF (.NOT.ieee_is_finite(isolated_dt) .OR. isolated_dt<=0.0_wp) &
            ERROR STOP 'Invalid dense L2 CFL'
        CALL slot2_spatial%evaluate(slot2_qp,slot2_term,0.0_wp,slot2_domain)

        transport_norm=0.0_wp
        DO l=1,slot2_domain%solve_cells
            j=slot2_domain%j_cent(l)
            k=slot2_domain%k_cent(l)
            IF (.NOT.ALL(ieee_is_finite(slot2_term(:,j,k)))) &
                ERROR STOP 'Non-finite moving L2 gas transport'
            transport_norm=MAX(transport_norm,ABS(slot2_term(idx_totMassEqn,j,k)))
            slot_q(1:n_eqns,2,j,k)=slot_q(1:n_eqns,2,j,k)-isolated_dt*slot2_term(:,j,k)
            IF (.NOT.conservative_mass_admissible(slot_q(1:n_eqns,2,j,k))) &
                ERROR STOP 'Isolated L2 spatial update lost admissible mass fractions'
            IF (slot_q(idx_engyEqn,2,j,k)<=0.0_wp) &
                ERROR STOP 'Isolated L2 spatial update lost positive thermal energy'
        END DO

        IF (nx>8 .AND. ny>8 .AND. transport_norm==0.0_wp) &
            ERROR STOP 'Moving dense L2 has no conservative transport'
        final_mass=SUM(slot_q(idx_totMassEqn,2,:,:))
        final_solids=SUM(slot_q(idx_solidEqn_first,2,:,:))
        final_energy=SUM(slot_q(idx_engyEqn,2,:,:))
        IF (ABS(final_mass-initial_mass)>4096.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,initial_mass)) &
            ERROR STOP 'Isolated L2 total mass not conserved'
        IF (ABS(final_solids-initial_solids)>4096.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,initial_solids)) &
            ERROR STOP 'Isolated L2 solid mass not conserved'
        IF (ABS(final_energy-initial_energy)>4096.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,initial_energy)) &
            ERROR STOP 'Isolated L2 energy not conserved'

        IF (ANY(slot_q(:,1,:,:)/=0.0_wp)) ERROR STOP 'Isolated L2 modified dry L1'
        CALL same_bits('L2 gas evaluation does not modify L1',RESHAPE(q,[SIZE(q)]),RESHAPE(q0,[SIZE(q0)]))
        gas_flag=.FALSE.
        liquid_flag=.TRUE.
        WRITE(*,*) 'PASS: isolated buoyant/dense L2 gas pressure, conservation and transport'
    END SUBROUTINE check_isolated_gas_l2


  !> \brief Configure a one-solid liquid model with zero solid fraction and no physical source except curvature.
  !> \details Initializes shared read-only configuration and owns all fixture arrays;
  !> called serially before OpenMP evaluation. The finalizer releases backend storage.
  SUBROUTINE initialize_fixture
    INTEGER :: i
    comp_cells_x=nx; comp_cells_y=ny; comp_interfaces_x=nx+1; comp_interfaces_y=ny+1
    cell_size=dx; dx2=0.5_wp*dx; dy2=0.5_wp*dy
    one_by_dx=1.0_wp/dx; one_by_dy=1.0_wp/dy
    n_vars=5; n_eqns=5; n_solid=1; n_add_gas=0; n_stoch_vars=0; n_pore_vars=0
    kt_trace_reconstruction = 'froude_blend'
    CALL configure_canonical_state_indices()
    liquid_flag=.TRUE.; gas_flag=.FALSE.; rheology_flag=.FALSE.
    slope_correction_flag=slope/=0; curvature_term_flag=curvature/=0
    rheology_l1%stochastic_rheology_flag=.FALSE.; pore_pressure_l1%pore_pressure_flag=.FALSE.
    pore_pressure_l1%gas_loss_flag=.FALSE.
    mass_exchange_l1%vaporization%liquid_vaporization_flag=.FALSE.
    mass_exchange_l1%entrainment%flag=.FALSE.
    radial_source_flag=.FALSE.; lateral_source_flag=.FALSE.
    bottom_radial_source_flag=.FALSE.; bottom_fissural_source_flag=.FALSE.
    sutherland_flag=.FALSE.; substrate%topo_change_flag=.FALSE.; verbose_level=-1
    n_fissures=0; n_intervals=0
    rho_l=1000.0_wp; inv_rho_l=1.0_wp/rho_l; rho_a_amb=1.2_wp; grav=9.81_wp
    sp_heat_l=4180.0_wp; sp_heat_a=998.0_wp; sp_gas_const_a=287.05_wp
    T_ambient=300.0_wp; pres=101300.0_wp; inv_pres=1.0_wp/pres
    eps_sing=1.0E-8_wp; eps_sing4=eps_sing**4
    limiter=limiter_id; theta=1.3_wp; reconstr_coeff=1.0_wp
    cfl=0.24_wp; max_dt=1.0E6_wp; maximum_solid_packing=0.6_wp
    ALLOCATE(rho_s(1),inv_rho_s(1),sp_heat_s(1),sp_heat_g(0),sp_gas_const_g(0))
    rho_s=2500.0_wp; inv_rho_s=1.0_wp/rho_s; sp_heat_s=1100.0_wp
    ALLOCATE(bcW(n_vars,2),bcE(n_vars,2),bcS(n_vars,2),bcN(n_vars,2))
    bcW%flag=1; bcW%value=0.0_wp; bcE%flag=1; bcE%value=0.0_wp
    bcS%flag=1; bcS%value=0.0_wp; bcN%flag=1; bcN%value=0.0_wp
    ALLOCATE(q(n_vars,nx,ny),q0(n_vars,nx,ny),qp(n_vars,nx,ny))
    ALLOCATE(term(n_eqns,nx,ny),direct(n_eqns,nx,ny),cell_sources(n_eqns,nx,ny))
    ALLOCATE(B_vertex(nx+1,ny+1),B_cent(nx,ny),B_face_x(nx+1,ny),B_face_y(nx,ny+1))
    ALLOCATE(B_prime_x_geom(nx,ny),B_prime_y_geom(nx,ny))
    ALLOCATE(B_second_xx_geom(nx,ny),B_second_xy_geom(nx,ny),B_second_yy_geom(nx,ny))
    ALLOCATE(grav_coeff(nx,ny),grav_coeff_stag_x(nx+1,ny),grav_coeff_stag_y(nx,ny+1))
    ALLOCATE(source_cell(nx,ny),sourceW(nx,ny),sourceE(nx,ny),sourceS(nx,ny),sourceN(nx,ny))
    source_cell=0; sourceW=.FALSE.; sourceE=.FALSE.; sourceS=.FALSE.; sourceN=.FALSE.
    ALLOCATE(x_comp(nx),y_comp(ny),x_stag(nx+1),y_stag(ny+1))
    DO i=1,nx+1
       x_stag(i)=REAL(i-1,wp)*dx
       IF (i<=nx) x_comp(i)=(REAL(i,wp)-0.5_wp)*dx
    END DO
    DO i=1,ny+1
       y_stag(i)=REAL(i-1,wp)*dy
       IF (i<=ny) y_comp(i)=(REAL(i,wp)-0.5_wp)*dy
    END DO
    CALL domain%initialize
    CALL spatial%initialize
    CALL reconstruction%initialize
    CALL hyperbolic%initialize
  END SUBROUTINE initialize_fixture

  !> Verify that a prescribed wet inlet remains wet in HP reconstruction,
  !> even if the entire initial flow is dry (as in EXAMPLE_BW).
  !> The same BC depth must be selected by all supported trace modes.
  SUBROUTINE check_dry_dirichlet_inlet
    TYPE(bc) :: saved_depth_bc,saved_discharge_bc
    CHARACTER(LEN=32) :: saved_mode
    CHARACTER(LEN=16), PARAMETER :: modes(4)=[CHARACTER(LEN=16):: &
         'h_x','eta_x','minmod','froude_blend']
    REAL(wp) :: scale
    INTEGER :: mode

    IF (nx<=1) RETURN
    saved_depth_bc=bcW(idx_h,1)
    saved_discharge_bc=bcW(idx_u,1)
    saved_mode=kt_trace_reconstruction

    bcW(idx_h,1)%flag=0
    bcW(idx_h,1)%value=4.0_wp
    bcW(idx_u,1)%flag=0
    bcW(idx_u,1)%value=10.0_wp
    q=0.0_wp
    qp=0.0_wp
    qp(idx_T,:,:)=T_ambient
    CALL domain%check_solve(q,0.0_wp,.FALSE.)
    IF (domain%solve_cells==0) ERROR STOP 'dry inlet not included in active domain'
    CALL reconstruction%prepare_hp_workset(domain%hp_eta_cells,domain%hp_eta_j,domain%hp_eta_k)

    DO mode=1,SIZE(modes)
       kt_trace_reconstruction=TRIM(modes(mode))
       CALL reconstruction%reconstruct(qp,0.0_wp,domain%solve_cells,domain%j_cent,domain%k_cent,.TRUE.)
       CALL reconstruction%reconstruct_kt_eta(qp)

       scale=MAX(1.0_wp,ABS(B_face_x(1,1)))
       IF (ABS(reconstruction%h_selected_interfaceR(1,1)-4.0_wp)> &
            64.0_wp*EPSILON(1.0_wp)) ERROR STOP 'wet Dirichlet HP depth lost'
       IF (ABS(reconstruction%qp_adv_interfaceR(idx_h,1,1)-4.0_wp)> &
            64.0_wp*EPSILON(1.0_wp)) ERROR STOP 'Dirichlet primitive depth lost'
       IF (ABS(reconstruction%qp_adv_interfaceR(idx_u,1,1)-2.5_wp)> &
            64.0_wp*EPSILON(1.0_wp)) ERROR STOP 'inlet discharge not preserved'
       IF (ABS(reconstruction%eta_interfaceR(1,1)-B_face_x(1,1)-4.0_wp)> &
            64.0_wp*EPSILON(1.0_wp)*scale) ERROR STOP 'boundary eta inconsistent with bed'
       IF (ABS(reconstruction%h_selected_interfaceL(1,1)-4.0_wp)> &
            64.0_wp*EPSILON(1.0_wp)) ERROR STOP 'boundary ghost depth not copied'
    END DO

    ! A dry Dirichlet boundary with zero discharge remains dry.
    bcW(idx_h,1)%value=0.0_wp
    bcW(idx_u,1)%value=0.0_wp
    kt_trace_reconstruction='froude_blend'
    CALL reconstruction%reconstruct(qp,0.0_wp,domain%solve_cells,domain%j_cent,domain%k_cent,.TRUE.)
    CALL reconstruction%reconstruct_kt_eta(qp)
    IF (reconstruction%h_selected_interfaceR(1,1)/=0.0_wp) &
         ERROR STOP 'zero prescribed boundary depth did not stay dry'
    IF (reconstruction%qp_adv_interfaceR(idx_u,1,1)/=0.0_wp) &
         ERROR STOP 'zero prescribed discharge did not stay zero'

    bcW(idx_h,1)=saved_depth_bc
    bcW(idx_u,1)=saved_discharge_bc
    kt_trace_reconstruction=saved_mode
    q=q0
    CALL reconstruction%finalize
    CALL reconstruction%initialize
    WRITE(*,*) 'PASS: wet Dirichlet inflow into a dry domain survives HP selection'
  END SUBROUTINE check_dry_dirichlet_inlet

  !> Verify that the same reconstruction kernel reads distinct L1/L2 boundary data.
  !> Both layers start dry; no coupled simulation or interlayer physics is enabled.
  SUBROUTINE check_layer_local_dirichlet_inlet
    TYPE(reconstruction_workspace_type) :: reconstruction_l2
    TYPE(bc) :: saved_l1_depth,saved_l1_discharge,saved_l2_depth,saved_l2_discharge
    CHARACTER(LEN=32) :: saved_mode
    CHARACTER(LEN=16), PARAMETER :: modes(4)=[CHARACTER(LEN=16):: &
         'h_x','eta_x','minmod','froude_blend']
    REAL(wp) :: l1_depth,l1_velocity,l1_eta,scale
    INTEGER :: mode

    IF (nx<=1) RETURN

    saved_l1_depth=bcW(idx_h,1)
    saved_l1_discharge=bcW(idx_u,1)
    saved_l2_depth=bcW(idx_h,2)
    saved_l2_discharge=bcW(idx_u,2)
    saved_mode=kt_trace_reconstruction

    bcW(idx_h,1)%flag=0
    bcW(idx_h,1)%value=4.0_wp
    bcW(idx_u,1)%flag=0
    bcW(idx_u,1)%value=10.0_wp
    bcW(idx_h,2)%flag=0
    bcW(idx_h,2)%value=2.0_wp
    bcW(idx_u,2)%flag=0
    bcW(idx_u,2)%value=6.0_wp

    q=0.0_wp
    qp=0.0_wp
    qp(idx_T,:,:)=T_ambient
    CALL domain%check_solve(q,0.0_wp,.FALSE.)
    IF (domain%solve_cells==0) ERROR STOP 'L1 dry inlet has no active workset'
    CALL reconstruction%prepare_hp_workset(domain%hp_eta_cells,domain%hp_eta_j,domain%hp_eta_k)
    CALL reconstruction_l2%initialize(2)
    CALL reconstruction_l2%prepare_hp_workset(domain%hp_eta_cells,domain%hp_eta_j,domain%hp_eta_k)

    DO mode=1,SIZE(modes)
       kt_trace_reconstruction=TRIM(modes(mode))
       CALL reconstruction%reconstruct(qp,0.0_wp,domain%solve_cells,domain%j_cent,domain%k_cent,.TRUE.)
       CALL reconstruction%reconstruct_kt_eta(qp)
       l1_depth=reconstruction%h_selected_interfaceR(1,1)
       l1_velocity=reconstruction%qp_adv_interfaceR(idx_u,1,1)
       l1_eta=reconstruction%eta_interfaceR(1,1)
       IF (ABS(l1_depth-4.0_wp)>64.0_wp*EPSILON(1.0_wp)) &
            ERROR STOP 'Layer-local test: L1 boundary depth mismatch'
       IF (ABS(l1_velocity-2.5_wp)>64.0_wp*EPSILON(1.0_wp)) &
            ERROR STOP 'Layer-local test: L1 boundary velocity mismatch'

       CALL reconstruction_l2%reconstruct(qp,0.0_wp,domain%solve_cells,domain%j_cent,domain%k_cent,.TRUE.)
       CALL reconstruction_l2%reconstruct_kt_eta(qp)
       IF (ABS(reconstruction_l2%h_selected_interfaceR(1,1)-2.0_wp)> &
            64.0_wp*EPSILON(1.0_wp)) ERROR STOP 'L2 boundary depth read from L1'
       IF (ABS(reconstruction_l2%qp_adv_interfaceR(idx_u,1,1)-3.0_wp)> &
            64.0_wp*EPSILON(1.0_wp)) ERROR STOP 'L2 boundary discharge read from L1'
       scale=MAX(1.0_wp,ABS(B_face_x(1,1)))
       IF (ABS(reconstruction_l2%eta_interfaceR(1,1)-B_face_x(1,1)-2.0_wp)> &
            64.0_wp*EPSILON(1.0_wp)*scale) ERROR STOP 'L2 boundary eta inconsistent with bed'

       IF (reconstruction%h_selected_interfaceR(1,1)/=l1_depth .OR. &
           reconstruction%qp_adv_interfaceR(idx_u,1,1)/=l1_velocity .OR. &
           reconstruction%eta_interfaceR(1,1)/=l1_eta) &
           ERROR STOP 'L2 reconstruction changed L1 boundary traces'
    END DO

    CALL reconstruction_l2%finalize
    bcW(idx_h,1)=saved_l1_depth
    bcW(idx_u,1)=saved_l1_discharge
    bcW(idx_h,2)=saved_l2_depth
    bcW(idx_u,2)=saved_l2_discharge
    kt_trace_reconstruction=saved_mode
    q=q0
    CALL reconstruction%finalize
    CALL reconstruction%initialize
    WRITE(*,*) 'PASS: layer-local Dirichlet inlet and independent HP traces'
  END SUBROUTINE check_layer_local_dirichlet_inlet

  !> Verify disabled optional slots in the complete reconstructed HP-PCCU operator.
  !> Each isolated layer uses the same seven-component layout: L1 has only pore
  !> pressure and L2 has only stochastic transport. The production N_LAYERS=2
  !> execution gate remains unchanged.
  SUBROUTINE check_optional_layer_spatial_slots
    TYPE(domain_type) :: test_domain
    TYPE(reconstruction_workspace_type) :: test_reconstruction
    TYPE(reconstruction_workspace_type) :: lower_reconstruction
    TYPE(hyperbolic_workspace_type) :: test_hyperbolic
    TYPE(state_type) :: empty_state
    TYPE(source_group_type) :: geometric_group
    REAL(wp), ALLOCATABLE :: primitive(:,:,:),conservative(:,:,:),residual(:,:,:),source_rate(:)
    REAL(wp), ALLOCATABLE :: upper_center(:,:),upper_w(:,:),upper_e(:,:),upper_s(:,:),upper_n(:,:)
    INTEGER :: layer,active_index,inactive_index,ix,iy
    LOGICAL :: accepted,active_flux_found
    REAL(wp) :: amplitude,saved_e,adjacent_w

    ! Some reference fixtures collapse one axis; only use an actual 2-D case.
    IF (nx<3 .OR. ny<3) RETURN

    n_stoch_vars=1
    n_pore_vars=1
    n_vars=4+n_solid+n_add_gas+n_stoch_vars+n_pore_vars
    n_eqns=n_vars
    CALL configure_canonical_state_indices()
    IF (idx_stoch==0 .OR. idx_pore==0) ERROR STOP 'Optional layout indices missing'

    DEALLOCATE(bcW,bcE,bcS,bcN)
    ALLOCATE(bcW(n_vars,2),bcE(n_vars,2),bcS(n_vars,2),bcN(n_vars,2))
    bcW%flag=1; bcW%value=0.0_wp
    bcE%flag=1; bcE%value=0.0_wp
    bcS%flag=1; bcS%value=0.0_wp
    bcN%flag=1; bcN%value=0.0_wp

    ! Both allocated layers must start with precisely zero conservative transport.
    ! Bypass the production model-layout gate: this is a structural test only.
    ALLOCATE(empty_state%q(n_vars,2,1,1))
    empty_state%q=1.0_wp
    CALL init_empty(empty_state)
    IF (ANY(empty_state%q/=0.0_wp)) ERROR STOP 'Dry initial state retained an optional layer'
    DEALLOCATE(empty_state%q)

    geometric_group%method=SOURCE_METHOD_EXPLICIT
    geometric_group%source_mask=IBSET(0,SOURCE_ID_GEOMETRY_CURVATURE-1)
    ALLOCATE(primitive(n_vars,nx,ny),conservative(n_vars,nx,ny),residual(n_eqns,nx,ny))
    ALLOCATE(source_rate(n_eqns))
    DO layer=1,2
        active_index=idx_pore
        inactive_index=idx_stoch
        IF (layer==2) THEN
            active_index=idx_stoch
            inactive_index=idx_pore
        END IF

        primitive=0.0_wp
        primitive(idx_h,:,:)=0.8_wp
        primitive(idx_u,:,:)=0.3_wp
        primitive(idx_v,:,:)=0.1_wp
        primitive(idx_T,:,:)=T_ambient
        DO iy=1,ny
            DO ix=1,nx
                amplitude=0.1_wp+0.01_wp*REAL(ix,wp)
                primitive(active_index,ix,iy)=amplitude
                CALL primitive_to_conservative(primitive(:,ix,iy),conservative(:,ix,iy),accepted)
                IF (.NOT.accepted) ERROR STOP 'Optional layer state conversion rejected'
            END DO
        END DO
        IF (ANY(conservative(inactive_index,:,:)/=0.0_wp)) &
            ERROR STOP 'Inactive optional conservative state is nonzero'

        ! Active curvature accelerates momentum without creating optional mass.
        ! The same cell-local source kernel is reused with each layer selection.
        CALL eval_source_explicit(0.2_wp,0.1_wp,0.05_wp,0.0_wp,0.02_wp,1.0_wp, &
            primitive(:,2,2),source_rate,0.0_wp,0.0_wp,no_fissures, &
            0.0_wp,0.0_wp,0.0_wp,dx*dy,geometric_group,0.0_wp, &
            materials=materials,thermal=thermal_l1,rheology=rheology_l1, &
            pore_enabled=layer==1,stochastic_enabled=layer==2)
        IF (source_rate(idx_stoch)/=0.0_wp .OR. source_rate(idx_pore)/=0.0_wp) &
            ERROR STOP 'Explicit geometry source created optional transport'
        IF (source_rate(idx_uEqn)==0.0_wp .AND. source_rate(idx_vEqn)==0.0_wp) &
            ERROR STOP 'Optional-layer source check did not activate curvature'

        CALL test_domain%initialize
        CALL test_domain%check_solve(conservative,0.0_wp,.TRUE.)
        CALL test_reconstruction%initialize(layer)
        CALL test_reconstruction%prepare_hp_workset( &
            test_domain%hp_eta_cells,test_domain%hp_eta_j,test_domain%hp_eta_k)
        IF (layer==2) THEN
            ! Reuse the existing HP reconstruction with the completed L1 traces.
            ! The ordinary L1 reference run above never supplies these arguments.
            CALL lower_reconstruction%initialize(1)
            CALL lower_reconstruction%prepare_hp_workset( &
                test_domain%hp_eta_cells,test_domain%hp_eta_j,test_domain%hp_eta_k)
            CALL lower_reconstruction%reconstruct(primitive,0.0_wp,test_domain%solve_cells, &
                test_domain%j_cent,test_domain%k_cent,.TRUE.)
            CALL lower_reconstruction%reconstruct_kt_eta(primitive)
            CALL test_reconstruction%reconstruct(primitive,0.0_wp,test_domain%solve_cells, &
                test_domain%j_cent,test_domain%k_cent,.TRUE., &
                lower_recon=lower_reconstruction,h_lower=primitive(idx_h,:,:))
            CALL test_reconstruction%reconstruct_kt_eta(primitive,lower_recon=lower_reconstruction)
            ! Exercise the complete upper-layer backend with the lower HP geometry.
            ! No two-layer momentum coupling is enabled by this isolated test.
            CALL test_hyperbolic%initialize(clamp_buoyant_pressure=.TRUE.)
            CALL test_hyperbolic%evaluate_terms(test_reconstruction,primitive,residual,0.0_wp, &
                test_domain%solve_cells,test_domain%j_cent,test_domain%k_cent, &
                test_domain%solve_interfaces_x,test_domain%j_stag_x,test_domain%k_stag_x, &
                test_domain%solve_interfaces_y,test_domain%j_stag_y,test_domain%k_stag_y, &
                lower_recon=lower_reconstruction,h_lower=primitive(idx_h,:,:))
            DO iy=1,ny
                DO ix=1,nx
                    IF (.NOT.test_domain%solve_mask(ix,iy)) CYCLE
                    IF (.NOT.ALL(ieee_is_finite(residual(:,ix,iy)))) &
                        ERROR STOP 'Nonfinite upper-layer hydrostatic path with L1 geometry'
                END DO
            END DO
            CALL test_hyperbolic%finalize
            IF (ANY(test_reconstruction%eta_cellW(2:3,2:3) < &
                lower_reconstruction%eta_cellW(2:3,2:3))) &
                ERROR STOP 'Upper surface fell below lower-layer reconstructed interface'
            IF (ANY(test_reconstruction%h_selected_interfaceR(2:3,2:3)<0.0_wp)) &
                ERROR STOP 'Upper layer reconstructed negative depth'
            CALL lower_reconstruction%finalize
        END IF
        CALL test_hyperbolic%initialize(clamp_buoyant_pressure=layer==2)
        CALL test_hyperbolic%evaluate_terms(test_reconstruction,primitive,residual,0.0_wp, &
            test_domain%solve_cells,test_domain%j_cent,test_domain%k_cent, &
            test_domain%solve_interfaces_x,test_domain%j_stag_x,test_domain%k_stag_x, &
            test_domain%solve_interfaces_y,test_domain%j_stag_y,test_domain%k_stag_y)

        IF (layer==1) THEN
            ALLOCATE(upper_center(nx,ny),upper_w(nx,ny),upper_e(nx,ny),upper_s(nx,ny),upper_n(nx,ny))
            CALL test_reconstruction%prepare_upper_layer_bed(primitive(idx_h,:,:), &
                upper_center,upper_w,upper_e,upper_s,upper_n)
            IF (ANY(upper_center/=B_cent+primitive(idx_h,:,:))) &
                ERROR STOP 'Upper-layer center base is not B+h1'

            DO iy=1,ny
                DO ix=1,nx
                    IF (.NOT.test_domain%hp_eta_mask(ix,iy)) CYCLE
                    IF (upper_w(ix,iy)/=test_reconstruction%eta_cellW(ix,iy) .OR. &
                        upper_e(ix,iy)/=test_reconstruction%eta_cellE(ix,iy) .OR. &
                        upper_s(ix,iy)/=test_reconstruction%eta_cellS(ix,iy) .OR. &
                        upper_n(ix,iy)/=test_reconstruction%eta_cellN(ix,iy)) &
                        ERROR STOP 'Upper-layer face base differs from L1 traces'
                END DO
            END DO

            ! Adjacent cell-side bases must not be silently averaged.
            saved_e=test_reconstruction%eta_cellE(2,2)
            adjacent_w=test_reconstruction%eta_cellW(3,2)
            test_reconstruction%eta_cellE(2,2)=saved_e+0.125_wp
            CALL test_reconstruction%prepare_upper_layer_bed(primitive(idx_h,:,:), &
                upper_center,upper_w,upper_e,upper_s,upper_n)
            IF (upper_e(2,2)/=saved_e+0.125_wp .OR. upper_w(3,2)/=adjacent_w) &
                ERROR STOP 'Upper-layer base averaged opposite L1 traces'
            test_reconstruction%eta_cellE(2,2)=saved_e
            DEALLOCATE(upper_center,upper_w,upper_e,upper_s,upper_n)
        END IF

        ! Inspect only active faces: unvisited scratch is outside the operator contract.
        active_flux_found=.FALSE.
        DO iy=1,ny
            DO ix=1,nx+1
                IF (.NOT.test_domain%solve_mask_x(ix,iy)) CYCLE
                IF (test_reconstruction%qp_adv_interfaceL(inactive_index,ix,iy)/=0.0_wp .OR. &
                    test_reconstruction%qp_adv_interfaceR(inactive_index,ix,iy)/=0.0_wp) &
                    ERROR STOP 'Inactive optional x trace is nonzero'
                IF (test_hyperbolic%G_interface_xL(inactive_index,ix,iy)/=0.0_wp .OR. &
                    test_hyperbolic%G_interface_xR(inactive_index,ix,iy)/=0.0_wp) &
                    ERROR STOP 'Inactive optional x flux is nonzero'
                IF (ABS(test_hyperbolic%G_interface_xL(active_index,ix,iy))>0.0_wp) active_flux_found=.TRUE.
            END DO
        END DO
        DO iy=1,ny+1
            DO ix=1,nx
                IF (.NOT.test_domain%solve_mask_y(ix,iy)) CYCLE
                IF (test_reconstruction%qp_adv_interfaceB(inactive_index,ix,iy)/=0.0_wp .OR. &
                    test_reconstruction%qp_adv_interfaceT(inactive_index,ix,iy)/=0.0_wp) &
                    ERROR STOP 'Inactive optional y trace is nonzero'
                IF (test_hyperbolic%G_interface_yB(inactive_index,ix,iy)/=0.0_wp .OR. &
                    test_hyperbolic%G_interface_yT(inactive_index,ix,iy)/=0.0_wp) &
                    ERROR STOP 'Inactive optional y flux is nonzero'
            END DO
        END DO
        IF (.NOT.active_flux_found) ERROR STOP 'Active optional component not transported'
        DO iy=1,ny
            DO ix=1,nx
                IF (.NOT.test_domain%solve_mask(ix,iy)) CYCLE
                IF (residual(inactive_index,ix,iy)/=0.0_wp) &
                    ERROR STOP 'Inactive optional spatial residual is nonzero'
                IF (test_hyperbolic%P_cell_x(inactive_index,ix,iy)/=0.0_wp .OR. &
                    test_hyperbolic%P_cell_y(inactive_index,ix,iy)/=0.0_wp) &
                    ERROR STOP 'Inactive optional spatial path is nonzero'
            END DO
        END DO

        CALL test_hyperbolic%finalize
        CALL test_reconstruction%finalize
        CALL test_domain%finalize
    END DO
    DEALLOCATE(primitive,conservative,residual,source_rate)
    WRITE(*,*) 'PASS: inactive optional L1/L2 initialization, sources, traces and fluxes'
  END SUBROUTINE check_optional_layer_spatial_slots

  !> \brief Verify the vector generalized-minmod kernel against the scalar limiter.
  !> \details Uses mixed-sign and near-zero component stencils with a
  !> local coordinate triplet of spacing dx, independent of domain dimensions.
  !> This guards the numerical
  !> formula independently from the spatial reference comparison.
  SUBROUTINE check_uniform_minmod_kernel
    REAL(wp) :: west(5),center(5),east(5),a(5),b(5),c(5),vector_slope(5)
    REAL(wp) :: scalar_slope,coords(3),scale
    INTEGER :: i
    west=[0.1_wp,-0.5_wp,1.0E-45_wp,3.0_wp,4.0_wp]
    center=[0.5_wp,0.2_wp,0.0_wp,3.0_wp,2.0_wp]
    east=[0.9_wp,-0.1_wp,-1.0E-45_wp,1.0_wp,3.0_wp]
    ! The CFL-axis tests also construct grids with nx=1 or ny=1. This
    ! mathematical unit test must not index nonexistent mesh neighbours.
    coords=[-dx,0.0_wp,dx]
    a=(east-center)/(coords(3)-coords(2))
    b=(center-west)/(coords(2)-coords(1))
    c=(east-west)/(coords(3)-coords(1))
    vector_slope=minmod(c,theta*minmod(a,b))
    DO i=1,5
       CALL limit([west(i),center(i),east(i)],coords,3,scalar_slope)
       scale=MAX(1.0_wp,ABS(scalar_slope))
       IF (ABS(vector_slope(i)-scalar_slope)>32.0_wp*EPSILON(1.0_wp)*scale) THEN
          ERROR STOP 'generalized minmod array kernel differs from scalar limiter'
       END IF
    END DO
  END SUBROUTINE check_uniform_minmod_kernel

  !> \brief Zero diagnostic entries outside the production active-face contract before serialization.
  SUBROUTINE normalize_inactive_face_payloads
    DO k=1,ny
       DO j=1,nx+1
          IF (domain%solve_mask_x(j,k)) CYCLE
          reconstruction%qp_adv_interfaceL(:,j,k)=0.0_wp; reconstruction%qp_adv_interfaceR(:,j,k)=0.0_wp
          reconstruction%eta_interfaceL(j,k)=0.0_wp; reconstruction%eta_interfaceR(j,k)=0.0_wp
       END DO
    END DO
    DO k=1,ny+1
       DO j=1,nx
          IF (domain%solve_mask_y(j,k)) CYCLE
          reconstruction%qp_adv_interfaceB(:,j,k)=0.0_wp; reconstruction%qp_adv_interfaceT(:,j,k)=0.0_wp
          reconstruction%eta_interfaceB(j,k)=0.0_wp; reconstruction%eta_interfaceT(j,k)=0.0_wp
       END DO
    END DO
  END SUBROUTINE normalize_inactive_face_payloads

  !> \brief Compare exact bit patterns to test immutable inputs and repeated public API evaluations.
  !> \param[in] label Assertion name.
  !> \param[in] actual,reference Flattened REAL fields to compare.
  SUBROUTINE same_bits(label,actual,reference)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: actual(:),reference(:)
    IF (ANY(TRANSFER(actual,[0_int64],SIZE(actual)) /= TRANSFER(reference,[0_int64],SIZE(reference)))) THEN
       WRITE(*,*) 'FAIL: ',label
       ERROR STOP 1
    END IF
  END SUBROUTINE same_bits
END PROGRAM test_spatial_operator
