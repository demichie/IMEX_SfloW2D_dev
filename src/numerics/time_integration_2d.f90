!> @brief Coordinate first-order Lie integration of spatial fluxes and physical sources.
!> @details The spatial operator advances first without rate limiting. Ordered
!>          source groups subsequently modify a transaction candidate; only a
!>          fully accepted timestep publishes fluid, substrate and geometry.
MODULE time_integration_2d
    USE parameters_2d, ONLY : wp,n_eqns,n_vars,idx_totMassEqn,idx_solidEqn_first,idx_solidEqn_last,                  &
            idx_add_gas_first,idx_add_gas_last,sea_mass_removal_flag,n_layers,idx_h
    USE parameters_2d, ONLY : cfl,max_dt,lateral_source_flag
    USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
    USE geometry_2d, ONLY : comp_cells_x,comp_cells_y,                                                                 &
            begin_geometry_transaction,commit_geometry_transaction,                                                    &
            B_zone
    USE state_conversions, ONLY : conservative_to_primitive
    USE state_fields, ONLY : state_type,initialize_dry_physical_cell
    USE state_composition, ONLY : conservative_mass_admissible
    USE diagnostics_2d, ONLY : fatal_error
    USE thermal_parameters_2d, ONLY : thermal_config_type
    USE rheology_parameters_2d, ONLY : rheology_config_type
    USE pore_pressure_parameters_2d, ONLY : pore_pressure_config_type, pore_pressure_cache_type
    USE mass_exchange_parameters_2d, ONLY : mass_exchange_config_type
    USE layer_physics_2d, ONLY : layer_physics_type
    USE source_registry_2d, ONLY : source_registry_type
    USE stochastic_module, ONLY : stochastic_workspace_type
    USE spatial_operator_2d, ONLY : spatial_operator_type
    USE domain_2d, ONLY : domain_type
    USE lie_source_executor_2d, ONLY : execute_lie_source_groups,initialize_lie_source_executor, &
        finalize_lie_source_executor,report_lie_source_profile,lie_source_executor_type
    IMPLICIT NONE
    PRIVATE
    ABSTRACT INTERFACE
        !> @brief stochastic stage observer operation.
        SUBROUTINE stochastic_stage_observer()
        END SUBROUTINE stochastic_stage_observer
    END INTERFACE
    PUBLIC :: time_integration_workspace_type
    !> @brief Scratch and transactional state buffers for one Lie integration.
    TYPE :: time_integration_workspace_type
        !> Spatial residual, dimensioned by conservative equation and cell.
        REAL(wp), ALLOCATABLE :: spatial_rate(:,:,:,:)
        !> Optional pre-spatial primitive halo when CFL traces cannot be reused.
        REAL(wp), ALLOCATABLE :: qp_before(:,:,:,:)
        !> Candidate fields: only worklist entries are initialized or published.
        REAL(wp), ALLOCATABLE :: q_trial(:,:,:,:), qp_trial(:,:,:,:)
        !> Persistent status worklist for primitive and spatial checks.
        INTEGER, ALLOCATABLE :: cell_status(:)
        TYPE(lie_source_executor_type) :: lie_executor
        LOGICAL :: profile_enabled=.FALSE.,front_mass_diagnostics=.FALSE.
        INTEGER :: profile_calls=0
        !> Profile times for preparing the pre-step residual and applying the Lie step [s].
        REAL(wp) :: profile_spatial=0.0_wp,profile_prepare=0.0_wp,profile_sources=0.0_wp
        REAL(wp) :: profile_source_executor=0.0_wp,profile_publish=0.0_wp
    CONTAINS
        PROCEDURE :: initialize => initialize_time_integration
        PROCEDURE :: finalize => finalize_time_integration
        PROCEDURE :: advance => advance_lie
        PROCEDURE :: prepare_spatial => prepare_spatial_for_timestep
    END TYPE
CONTAINS
    !> @brief Allocate integration workspaces from the configured domain sizes.
    !> @param[in,out] this Integration workspace to allocate.
    !> @param[in] source_registry Validated Lie source groups for the executor.
    SUBROUTINE initialize_time_integration(this,source_registry)
        CLASS(time_integration_workspace_type), INTENT(INOUT) :: this
        TYPE(source_registry_type), INTENT(IN) :: source_registry
        BLOCK
            CHARACTER(LEN=8) :: value
            INTEGER :: status
            value=''
            CALL GET_ENVIRONMENT_VARIABLE('IMEX_PROFILE_LIE',value,STATUS=status)
            this%profile_enabled=status==0 .AND. TRIM(value)=='1'
        END BLOCK
        BLOCK
            CHARACTER(LEN=8) :: value
            INTEGER :: status
            value=''
            CALL GET_ENVIRONMENT_VARIABLE('IMEX_FRONT_MASS_DIAGNOSTICS',value,STATUS=status)
            this%front_mass_diagnostics=status==0 .AND. TRIM(value)=='1'
        END BLOCK
        ALLOCATE(this%spatial_rate(n_eqns,n_layers,comp_cells_x,comp_cells_y))
        ALLOCATE(this%q_trial(n_vars,n_layers,comp_cells_x,comp_cells_y))
        ALLOCATE(this%qp_trial(n_vars,n_layers,comp_cells_x,comp_cells_y))
        ALLOCATE(this%cell_status(comp_cells_x*comp_cells_y))
        CALL initialize_lie_source_executor(this%lie_executor,source_registry)
    END SUBROUTINE initialize_time_integration

    !> @brief Prepare the unique spatial residual and both pre-step time bounds.
    !> @details The conservative state is not advanced. The residual is independent
    !> of dt and is consumed by the subsequent Lie step without reevaluation.
    !> Isolated L1/L2 bounds are reduced before the paired bounds overwrite them;
    !> their original CFL policy is applied only after residual assembly.
    !> @param[in,out] this Integration workspace owning the unique spatial residual.
    !> @param[in,out] state Immutable conservative field and refreshed primitive cache.
    !> @param[in] t Beginning of the timestep [s].
    !> @param[in] domain Shared active cell/face and HP reconstruction workset.
    !> @param[in,out] spatial_operators Canonical L1/L2 PCCU workspaces.
    !> @param[out] dt_cfl Existing CFL policy bound, capped at MAX_DT [s].
    !> @param[out] dt_pos Positive mass-only forward-Euler bound, capped at MAX_DT [s].
    !> @param[out] max_speed_x,max_speed_y Isolated L1 directional speed maxima [m/s].
    !> @param[in] boundary_source_enabled Optional original CFL boundary source switch.
    !> @param[in] diagnose_max_x_face Optional sampled CFL-face diagnostic switch.
    !> @param[out] max_speed_x max speed x argument of prepare_spatial_for_timestep.
    !> @param[out] max_speed_y max speed y argument of prepare_spatial_for_timestep.
    SUBROUTINE prepare_spatial_for_timestep(this,state,t,domain,spatial_operators, &
        dt_cfl,dt_pos,max_speed_x,max_speed_y,boundary_source_enabled,diagnose_max_x_face)
        CLASS(time_integration_workspace_type), INTENT(INOUT) :: this
        CLASS(state_type), INTENT(INOUT) :: state
        REAL(wp), INTENT(IN) :: t
        CLASS(domain_type), INTENT(IN) :: domain
        CLASS(spatial_operator_type), INTENT(INOUT) :: spatial_operators(:)
        REAL(wp), INTENT(OUT) :: dt_cfl,dt_pos,max_speed_x,max_speed_y
        LOGICAL, INTENT(IN), OPTIONAL :: boundary_source_enabled,diagnose_max_x_face
        REAL(wp) :: l1_x,l1_y,l2_x,l2_y,paired_x,paired_y
        REAL(wp) :: dt_l2,dt_paired,mass,rate,local_dt
        INTEGER :: l,j,k,layer,tick,tock,tick_rate
        LOGICAL :: invalid_mass_rate

        IF (SIZE(spatial_operators)/=n_layers) &
            CALL fatal_error('Spatial preparation requires exactly one operator per layer')
        IF (.NOT.ieee_is_finite(max_dt) .OR. max_dt<=0.0_wp) &
            CALL fatal_error('Invalid MAX_DT in positivity preparation')
        IF (this%profile_enabled) CALL SYSTEM_CLOCK(tick,tick_rate)

        ! Refresh the primitive halo and the original isolated face speeds.
        CALL spatial_operators(1)%prepare_traces(state%q(:,1,:,:),state%qp(:,1,:,:),t,domain, &
            boundary_source_enabled=boundary_source_enabled)
        IF (n_layers==2) THEN
            CALL spatial_operators(2)%prepare_traces(state%q(:,2,:,:),state%qp(:,2,:,:),t,domain, &
                boundary_source_enabled=boundary_source_enabled, &
                lower_operator=spatial_operators(1),h_lower=state%qp(idx_h,1,:,:))
        END IF

        l1_x=0.0_wp
        l1_y=0.0_wp
        l2_x=0.0_wp
        l2_y=0.0_wp
        ! Preserve isolated CFL maxima before coupling overwrites face speeds.
        ! A single-layer lateral source also needs its original CFL traces:
        ! Lie reconstructs these with the source-side boundary defaults.
        IF (cfl/=-1.0_wp .AND. (n_layers==2 .OR. lateral_source_flag)) THEN
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l,j,k) REDUCTION(MAX:l1_x,l1_y,l2_x,l2_y)
            DO l=1,domain%solve_cells
                j=domain%j_cent(l)
                k=domain%k_cent(l)
                CALL spatial_operators(1)%accumulate_cfl_cell(j,k,l1_x,l1_y)
                IF (n_layers==2) CALL spatial_operators(2)%accumulate_cfl_cell(j,k,l2_x,l2_y)
            END DO
            !$OMP END PARALLEL DO
            CALL spatial_operators(1)%report_cfl_face(state%q(:,1,:,:),state%qp(:,1,:,:), &
                domain,l1_x,diagnose_max_x_face)
        END IF

        ! One authoritative residual buffer for both the timestep limits and Lie.
        IF (n_layers==2) THEN
            CALL spatial_operators(1)%evaluate_coupled(spatial_operators(2), &
                state%qp(:,1,:,:),state%qp(:,2,:,:), &
                this%spatial_rate(:,1,:,:),this%spatial_rate(:,2,:,:),t,domain)
        ELSE
            CALL spatial_operators(1)%evaluate(state%qp(:,1,:,:),this%spatial_rate(:,1,:,:),t,domain, &
                reuse_cfl_traces=.NOT.lateral_source_flag)
        END IF

        paired_x=0.0_wp
        paired_y=0.0_wp
        dt_pos=max_dt
        invalid_mass_rate=.FALSE.
        ! The same active-cell pass reduces final characteristic maxima and
        ! the forward-Euler conservative mixture-mass limit for both layers.
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l,j,k,layer,mass,rate,local_dt) &
        !$OMP& REDUCTION(MAX:paired_x,paired_y) REDUCTION(MIN:dt_pos) REDUCTION(.OR.:invalid_mass_rate)
        DO l=1,domain%solve_cells
            j=domain%j_cent(l)
            k=domain%k_cent(l)
            IF (cfl/=-1.0_wp .AND. (n_layers==2 .OR. .NOT.lateral_source_flag)) &
                CALL spatial_operators(1)%accumulate_cfl_cell(j,k,paired_x,paired_y,paired_bounds=n_layers==2)

            DO layer=1,n_layers
                mass=state%q(idx_totMassEqn,layer,j,k)
                rate=this%spatial_rate(idx_totMassEqn,layer,j,k)
                IF (.NOT.ieee_is_finite(mass) .OR. .NOT.ieee_is_finite(rate) .OR. mass<0.0_wp) THEN
                    invalid_mass_rate=.TRUE.
                    CYCLE
                END IF
                IF (rate<=0.0_wp) CYCLE
                IF (mass==0.0_wp) THEN
                    dt_pos=0.0_wp
                    CYCLE
                END IF
                ! Avoid overflow from mass/max_dt when MAX_DT is very small;
                ! such a cell cannot constrain a timestep already <= MAX_DT.
                IF (max_dt<1.0_wp) THEN
                    IF (mass>HUGE(1.0_wp)*max_dt) CYCLE
                END IF
                IF (rate<mass/max_dt) CYCLE
                local_dt=mass/rate
                ! A ratio rounded up can make m-dt*rate negative at equality.
                ! Use the immediately smaller representable positive timestep.
                IF (local_dt>0.0_wp) local_dt=NEAREST(local_dt,-1.0_wp)
                dt_pos=MIN(dt_pos,local_dt)
            END DO
        END DO
        !$OMP END PARALLEL DO

        IF (invalid_mass_rate) &
            CALL fatal_error('Negative/nonfinite mass or mass residual in timestep preparation')
        IF (.NOT.ieee_is_finite(dt_pos) .OR. dt_pos<=0.0_wp) THEN
            DO l=1,domain%solve_cells
                j=domain%j_cent(l)
                k=domain%k_cent(l)
                DO layer=1,n_layers
                    mass=state%q(idx_totMassEqn,layer,j,k)
                    rate=this%spatial_rate(idx_totMassEqn,layer,j,k)
                    IF (rate<=0.0_wp) CYCLE
                    IF (mass>0.0_wp .AND. mass/rate>0.0_wp) CYCLE
                    WRITE(*,'(A,3(I0,1X),A,2(ES24.16,1X))') &
                        'DT_POS_INVALID layer,j,k=',layer,j,k,' mass,rate=',mass,rate
                    EXIT
                END DO
                IF (layer<=n_layers) EXIT
            END DO
            CALL fatal_error('No positive mass-preserving timestep exists')
        END IF

        IF (n_layers==1 .AND. .NOT.lateral_source_flag) THEN
            l1_x=paired_x
            l1_y=paired_y
            IF (cfl/=-1.0_wp) CALL spatial_operators(1)%report_cfl_face( &
                state%q(:,1,:,:),state%qp(:,1,:,:),domain,l1_x,diagnose_max_x_face)
        END IF
        max_speed_x=l1_x
        max_speed_y=l1_y
        CALL spatial_operators(1)%limit_from_speeds(l1_x,l1_y,dt_cfl)
        IF (n_layers==2) THEN
            CALL spatial_operators(2)%limit_from_speeds(l2_x,l2_y,dt_l2)
            CALL spatial_operators(1)%limit_from_speeds(paired_x,paired_y,dt_paired)
            dt_cfl=MIN(dt_cfl,dt_l2,dt_paired)
        END IF

        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            this%profile_prepare=this%profile_prepare+REAL(tock-tick,wp)/REAL(tick_rate,wp)
        END IF
    END SUBROUTINE prepare_spatial_for_timestep

    !> @brief Advance one timestep with spatial update followed by Lie source groups.
    !> @details Reads group-zero physical rates from the pre-spatial snapshot,
    !>          checks spatial admissibility without limiting the spatial rate,
    !>          and commits state and geometry only on success.
    !> @param[in,out] this Owns pre-step and trial workspace arrays.
    !> @param[in,out] state Fluid and substrate state, committed after acceptance.
    !> @param[in] t Start time of timestep [s].
    !> @param[in] dt Positive timestep duration [s].
    !> @param[in,out] stochastic Independent OU workspace for each active layer.
    !> @param[in] pre_spatial_primitives_ready Reuse CFL-prepared primitive traces when true.
    !> @param[in] spatial_residual_ready Reuse the unique residual prepared before selecting dt.
    !> @param[in] stochastic_prepared_hook Optional observer after OU preparation.
    !> @param[in,out] domain Shared worklists and wetness mask updated on accepted publication.
    !> @param[in,out] spatial_operators Per-layer HP-PCCU workspaces; both slots execute for L1/L2.
    !> @param[in] layer_physics Canonical per-layer physical descriptors.
    !> @param[in] pore_caches Derived Pascal coefficients for active layers.
    SUBROUTINE advance_lie(this,state,t,dt,stochastic,domain,spatial_operators,pre_spatial_primitives_ready, &
        spatial_residual_ready,stochastic_prepared_hook,layer_physics,pore_caches)
        CLASS(time_integration_workspace_type), TARGET, INTENT(INOUT) :: this
        CLASS(state_type), TARGET, INTENT(INOUT) :: state
        REAL(wp), INTENT(IN) :: t,dt
        TYPE(stochastic_workspace_type), INTENT(INOUT) :: stochastic(:)
        PROCEDURE(stochastic_stage_observer), OPTIONAL :: stochastic_prepared_hook
        CLASS(domain_type), INTENT(INOUT) :: domain
        CLASS(spatial_operator_type), INTENT(INOUT) :: spatial_operators(:)
        TYPE(layer_physics_type), INTENT(IN) :: layer_physics(:)
        TYPE(pore_pressure_cache_type), INTENT(IN) :: pore_caches(:)
        TYPE(thermal_config_type) :: thermal(n_layers)
        TYPE(rheology_config_type) :: rheology(n_layers)
        TYPE(pore_pressure_config_type) :: pore_pressure(n_layers)
        TYPE(mass_exchange_config_type) :: mass_config(n_layers)
        ! These are read-only timestep-local snapshots, not persistent state:
        ! they supply the existing per-layer array arguments to the OpenMP Lie executor.
        LOGICAL :: transport_stochastic(n_layers),pore_enabled(n_layers),stochastic_enabled(n_layers)
        LOGICAL, INTENT(IN), OPTIONAL :: pre_spatial_primitives_ready,spatial_residual_ready
        LOGICAL :: reuse_cfl_traces,residual_prepared
        ! The full 4D buffer is contiguous; fixing the layer index produces a
        ! noncontiguous 3D section when N_LAYERS=2. Pass those sections directly
        ! so explicit-shape callees receive correctly packed values.
        REAL(wp), POINTER, CONTIGUOUS :: qp_pre_all(:,:,:,:)
        REAL(wp) :: coupled_dt
        INTEGER :: l,j,k,first_failure,eq,failed_layer,layer_status,layer,front_j,front_k
        REAL(wp) :: face_left,face_right
        INTEGER :: tick,tock,tick_rate
        LOGICAL :: source_step_failed,conversion_ok
        CHARACTER(LEN=256) :: source_failure
        LOGICAL :: needs_geometry_update
        IF (SIZE(layer_physics)<n_layers .OR. SIZE(pore_caches)<n_layers) &
            CALL fatal_error('TL07: missing physical descriptor or Darcy cache for layer')
        DO layer=1,n_layers
            IF (.NOT.layer_physics(layer)%is_initialized()) &
                CALL fatal_error('TL07: uninitialized layer physical configuration')
            thermal(layer)=layer_physics(layer)%thermal_configuration()
            rheology(layer)=layer_physics(layer)%rheology_configuration()
            pore_pressure(layer)=layer_physics(layer)%pore_pressure_configuration()
            mass_config(layer)=layer_physics(layer)%mass_exchange_configuration()
            pore_enabled(layer)=pore_pressure(layer)%pore_pressure_flag
            stochastic_enabled(layer)=rheology(layer)%stochastic_rheology_flag
        END DO
        IF (SIZE(stochastic)<n_layers) CALL fatal_error('Missing stochastic workspace for layer')
        DO layer=1,n_layers
            IF (.NOT.ASSOCIATED(stochastic(layer)%config)) &
                CALL fatal_error('Uninitialized stochastic layer configuration')
            transport_stochastic(layer)=stochastic(layer)%config%stochastic_transport_flag
        END DO
        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tick,tick_rate)
            this%profile_calls=this%profile_calls+1
        END IF
        ! The conservative state remains immutable until commit. Trial buffers
        ! are initialized only on active cells, avoiding four full-grid copies.
        ! Production CFL already converted this unchanged pre-spatial state.
        ! Direct callers retain the original preflight unless they explicitly opt in.
        reuse_cfl_traces=.FALSE.
        IF (PRESENT(pre_spatial_primitives_ready)) reuse_cfl_traces=pre_spatial_primitives_ready
        residual_prepared=.FALSE.
        IF (PRESENT(spatial_residual_ready)) residual_prepared=spatial_residual_ready
        IF (residual_prepared .AND. .NOT.reuse_cfl_traces) &
            CALL fatal_error('Prepared spatial residual requires prepared primitive state')
        IF (n_layers==2 .AND. .NOT.reuse_cfl_traces) &
            CALL fatal_error('TL07 coupled Lie requires prepared L1/L2 CFL traces on a shared workset')
        IF (reuse_cfl_traces) THEN
            qp_pre_all=>state%qp
        ELSE
            ! This full halo buffer exists only for direct/non-CFL callers.
            ! Preserve their original reconstruction-stencil semantics.
            IF (.NOT.ALLOCATED(this%qp_before))                                                                     &
                ALLOCATE(this%qp_before(n_vars,n_layers,comp_cells_x,comp_cells_y))
            this%qp_before(:,1,:,:)=state%qp(:,1,:,:)
            qp_pre_all=>this%qp_before
        END IF
        IF (.NOT.reuse_cfl_traces) THEN
            ! Each cell initializes its own status in the conversion sweep.
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l,j,k,conversion_ok)
        DO l=1,domain%reconstruction_cells
            this%cell_status(l)=0
            j=domain%j_reconstruction(l)
            k=domain%k_reconstruction(l)

            IF (state%q(idx_totMassEqn,1,j,k)>0.0_wp) THEN
                CALL conservative_to_primitive(state%q(:,1,j,k),this%qp_before(:,1,j,k),                              &
                    accepted=conversion_ok)
                IF (.NOT.conversion_ok) this%cell_status(l)=1
            ELSE
                CALL initialize_dry_physical_cell(this%qp_before(:,1,j,k))

            END IF

        END DO
        !$OMP END PARALLEL DO
        first_failure=MINLOC(this%cell_status(1:domain%reconstruction_cells),1,                                      &
            MASK=this%cell_status(1:domain%reconstruction_cells)/=0)
        IF (first_failure/=0)                                                                                          &
            CALL fatal_error('Lie pre-spatial state conversion rejected')
        END IF

        ! Evaluate the shared spatial operator once, using the prepared OR workset.
        ! The production L1 path retains its original evaluation and floating-point order.
        IF (.NOT.residual_prepared) THEN
            IF (n_layers==2) THEN
                IF (SIZE(spatial_operators)/=2) CALL fatal_error('TL07 coupled Lie needs two spatial operators')
                CALL spatial_operators(1)%evaluate_coupled(spatial_operators(2), &
                    qp_pre_all(:,1,:,:),qp_pre_all(:,2,:,:), &
                    this%spatial_rate(:,1,:,:),this%spatial_rate(:,2,:,:),t,domain,coupled_dt)
                IF (dt>coupled_dt) CALL fatal_error('TL07 timestep exceeds preliminary coupled CFL bound')
            ELSE
                CALL spatial_operators(1)%evaluate(qp_pre_all(:,1,:,:),this%spatial_rate(:,1,:,:),t,domain, &
                    reuse_cfl_traces=reuse_cfl_traces)
            END IF
        END IF

        ! Select the first wet-to-dry x-face of L1 from the published state.
        ! The right cell sees the mass influx from G_interface_xR.
        front_j=0
        front_k=0
        IF (this%front_mass_diagnostics .AND. comp_cells_x>1) THEN
            find_front: DO k=1,comp_cells_y
                DO j=2,comp_cells_x
                    IF (state%q(idx_totMassEqn,1,j-1,k)<=0.0_wp) CYCLE
                    IF (state%q(idx_totMassEqn,1,j,k)>0.0_wp) CYCLE
                    IF (.NOT.domain%solve_mask(j,k)) CYCLE
                    front_j=j
                    front_k=k
                    EXIT find_front
                END DO
            END DO find_front
            IF (front_j>0) THEN
                WRITE(*,'(A,L1,A,2(ES24.16,1X))') 'FRONT_TRACE face_active=', &
                    domain%solve_mask_x(front_j,front_k), ' center_mass_L,R=', &
                    state%q(idx_totMassEqn,1,front_j-1,front_k), &
                    state%q(idx_totMassEqn,1,front_j,front_k)
                CALL spatial_operators(1)%front_mass_flux(front_j,front_k,face_left,face_right)
                WRITE(*,'(A,2(I0,1X),A,2(ES24.16,1X),A,2(ES24.16,1X))') &
                    'FRONT_MASS face(j,k)=',front_j,front_k,' t,dt=',t,dt, &
                    ' flux_to_left,right=',face_left,face_right
                WRITE(*,'(A,3(ES24.16,1X))') 'FRONT_MASS before,spatial_rate,expected_after=', &
                    state%q(idx_totMassEqn,1,front_j,front_k), &
                    this%spatial_rate(idx_totMassEqn,1,front_j,front_k), &
                    state%q(idx_totMassEqn,1,front_j,front_k)- &
                    dt*this%spatial_rate(idx_totMassEqn,1,front_j,front_k)
            ELSE
                WRITE(*,'(A,ES24.16)') 'FRONT_MASS no wet-to-dry x-face at t=',t
            END IF
        END IF

        ! Both layer candidates are checked before any authoritative fluid publication.
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l,j,k,layer,layer_status)
        DO l=1,domain%solve_cells
            j=domain%j_cent(l)
            k=domain%k_cent(l)
            this%cell_status(l)=0
            DO layer=1,n_layers
                CALL update_spatial_cell(state%q(:,layer,j,k),qp_pre_all(:,layer,j,k), &
                    this%spatial_rate(:,layer,j,k),dt,sea_mass_removal_flag .AND. B_zone(j,k)==1, &
                    this%q_trial(:,layer,j,k),this%qp_trial(:,layer,j,k),layer_status)
                IF (this%cell_status(l)==0 .AND. layer_status/=0) &
                    this%cell_status(l)=layer_status+2*(layer-1)
            END DO
        END DO
        !$OMP END PARALLEL DO
        first_failure=MINLOC(this%cell_status(1:domain%solve_cells),1,                                                &
            MASK=this%cell_status(1:domain%solve_cells)/=0)

        IF (first_failure/=0) THEN
            ! Diagnose the first failing cell outside the OpenMP region.
            ! Do not change, limit or repair the rejected conservative candidate.
            j=domain%j_cent(first_failure)
            k=domain%k_cent(first_failure)
            failed_layer=1
            layer_status=this%cell_status(first_failure)
            IF (layer_status>2) THEN
                failed_layer=2
                layer_status=layer_status-2
            END IF
            WRITE(*,'(A,2(I0,1X),A,ES23.15,A,ES23.15)') &
                'LIE_SPATIAL_FAILURE cell(j,k)=',j,k,' t=',t,' dt=',dt
            WRITE(*,'(A,I0)') 'LIE_SPATIAL_FAILURE status=',layer_status
            IF (n_layers==2) WRITE(*,'(A,I0)') 'LIE_SPATIAL_FAILURE layer=',failed_layer
            WRITE(*,'(A,ES23.15,A,ES23.15)') &
                'LIE_SPATIAL_FAILURE mass_before=',state%q(idx_totMassEqn,failed_layer,j,k), &
                ' mass_after=',this%q_trial(idx_totMassEqn,failed_layer,j,k)
            WRITE(*,'(A,ES23.15,A,ES23.15)') &
                'LIE_SPATIAL_FAILURE components_before=', &
                SUM(state%q(idx_solidEqn_first:idx_solidEqn_last,failed_layer,j,k))+ &
                SUM(state%q(idx_add_gas_first:idx_add_gas_last,failed_layer,j,k)), &
                ' components_after=', &
                SUM(this%q_trial(idx_solidEqn_first:idx_solidEqn_last,failed_layer,j,k))+ &
                SUM(this%q_trial(idx_add_gas_first:idx_add_gas_last,failed_layer,j,k))
            DO eq=1,n_eqns
                WRITE(*,'(A,I0,A,3(1X,ES23.15))') 'LIE_SPATIAL_FAILURE equation=',eq, &
                    ' before/rate/after=',state%q(eq,failed_layer,j,k), &
                    this%spatial_rate(eq,failed_layer,j,k),this%q_trial(eq,failed_layer,j,k)
                CALL spatial_operators(failed_layer)%diagnose_cell(j,k,eq)
            END DO
            IF (n_layers==2 .AND. layer_status==1) THEN
                IF (failed_layer==2) THEN
                    CALL spatial_operators(2)%diagnose_positivity(j,k,qp_pre_all(:,2,j,k), &
                        state%q(idx_totMassEqn,2,j,k),dt, &
                        lower=spatial_operators(1),lower_h_cell=qp_pre_all(idx_h,1,j,k))
                ELSE
                    CALL spatial_operators(1)%diagnose_positivity(j,k,qp_pre_all(:,1,j,k), &
                        state%q(idx_totMassEqn,1,j,k),dt)
                END IF
            END IF
            IF (layer_status==1)                                                                                      &
                CALL fatal_error('Lie spatial update: inadmissible conservative state')
            CALL fatal_error('Lie spatial state conversion rejected')

        END IF

        IF (front_j>0) THEN
            WRITE(*,'(A,ES24.16)') 'FRONT_MASS after_spatial L1=', &
                this%q_trial(idx_totMassEqn,1,front_j,front_k)
        END IF

        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            this%profile_spatial=this%profile_spatial+REAL(tock-tick,wp)/REAL(tick_rate,wp)
            tick=tock
        END IF
        ! Both one- and two-layer spatial candidates enter the same ordered
        ! Lie executor. There is no alternative or layer-specific commit path.

        ! OU is now the post-spatial fractional stage.  q_trial/qp_trial
        ! are valid only on the solve worklist, so stochastic%prepare_timestep
        ! uses the published state for inactive cells. Keep both published
        ! fluid arrays untouched until the final Lie transaction commits.
        DO layer=1,n_layers
            IF (.NOT.stochastic_enabled(layer) .AND. .NOT.transport_stochastic(layer)) CYCLE
            CALL stochastic(layer)%prepare_timestep(state,dt,transport_stochastic(layer), &
                this%q_trial(:,layer,:,:),this%qp_trial(:,layer,:,:),domain)
        END DO
        IF (PRESENT(stochastic_prepared_hook)) CALL stochastic_prepared_hook()
        ! Group 0 evaluates q^n but limits and applies to post-spatial q*.
        ! The geometry transaction supports deferred topography updates;
        ! it is not a recoverable timestep or retry mechanism.
        needs_geometry_update=this%lie_executor%has_bed_group
        IF (needs_geometry_update) CALL begin_geometry_transaction(sparse=.TRUE.)
        IF (n_layers==2) THEN
            CALL execute_lie_source_groups(this%lie_executor,state,this%q_trial,this%qp_trial,dt,domain, &
                stochastic(1)%effective_Z,t,state%q,qp_pre_all,source_step_failed,source_failure, &
                thermal,rheology,pore_pressure,pore_caches,mass_config, &
                transport_stochastic,pore_enabled,stochastic_enabled, &
                stochastic_values_l2=stochastic(2)%effective_Z)
        ELSE
            CALL execute_lie_source_groups(this%lie_executor,state,this%q_trial,this%qp_trial,dt,domain, &
                stochastic(1)%effective_Z,t,state%q,qp_pre_all,source_step_failed,source_failure, &
                thermal,rheology,pore_pressure,pore_caches,mass_config, &
                transport_stochastic,pore_enabled,stochastic_enabled)
        END IF
        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            this%profile_source_executor=this%profile_source_executor+REAL(tock-tick,wp)/REAL(tick_rate,wp)
        END IF

        ! A source failure is fatal: no restoration or timestep retry.
        IF (source_step_failed) CALL fatal_error('Lie timestep failed: '//TRIM(source_failure))

        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            this%profile_sources=this%profile_sources+REAL(tock-tick,wp)/REAL(tick_rate,wp)
        END IF
        IF (front_j>0) THEN
            WRITE(*,'(A,ES24.16)') 'FRONT_MASS after_sources L1=', &
                this%q_trial(idx_totMassEqn,1,front_j,front_k)
        END IF

        ! Publish the accepted fluid fields only after all Lie source groups return.
        ! The integrator owns evolving inventories; geometry owns only bed caches.
        IF (this%profile_enabled) CALL SYSTEM_CLOCK(tick,tick_rate)
        ! Publish both layer slots from one accepted transaction, using
        ! only the shared solve worklist after every source group succeeded.
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l,j,k)
        DO l=1,domain%solve_cells
            j=domain%j_cent(l)
            k=domain%k_cent(l)
            state%q(:,:,j,k)=this%q_trial(:,:,j,k)
            state%qp(:,:,j,k)=this%qp_trial(:,:,j,k)
            ! The fixed outer boundary is not reset by check_solve.
            IF (j>1 .AND. j<comp_cells_x .AND. k>1 .AND. k<comp_cells_y) &
                domain%solve_mask(j,k)=ANY(this%q_trial(idx_totMassEqn,1:n_layers,j,k)>0.0_wp)
        END DO
        !$OMP END PARALLEL DO
        domain%committed_mask_ready=.TRUE.
        IF (front_j>0) THEN
            WRITE(*,'(A,ES24.16)') 'FRONT_MASS published L1=', &
                state%q(idx_totMassEqn,1,front_j,front_k)
        END IF
        ! The source executor refreshed substrate aggregates after its cell sweep.
        ! Recomputing them here would traverse the whole domain a second time.
        IF (needs_geometry_update) CALL commit_geometry_transaction
        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            this%profile_publish=this%profile_publish+REAL(tock-tick,wp)/REAL(tick_rate,wp)
        END IF
    END SUBROUTINE advance_lie

    !> @brief Assemble and validate one layer's post-spatial cell candidate.
    !> @details Identical kernel for either layer; no committed state is changed.
    !> An invalid conservative candidate is rejected before sea removal, as in
    !> the original Lie stage. Status is 0 (accepted), 1 (inadmissible mass) or
    !> 2 (failed primitive conversion).
    !> @param[in] q_before q before argument of update_spatial_cell.
    !> @param[in] qp_before qp before argument of update_spatial_cell.
    !> @param[in] rate rate argument of update_spatial_cell.
    !> @param[in] dt Integration time increment.
    !> @param[in] sea_cell sea cell argument of update_spatial_cell.
    !> @param[out] q_candidate q candidate argument of update_spatial_cell.
    !> @param[out] qp_candidate qp candidate argument of update_spatial_cell.
    !> @param[out] status status argument of update_spatial_cell.
    SUBROUTINE update_spatial_cell(q_before,qp_before,rate,dt,sea_cell,q_candidate,qp_candidate,status)
        REAL(wp), INTENT(IN) :: q_before(n_vars),qp_before(n_vars),rate(n_eqns),dt
        LOGICAL, INTENT(IN) :: sea_cell
        REAL(wp), INTENT(OUT) :: q_candidate(n_vars),qp_candidate(n_vars)
        INTEGER, INTENT(OUT) :: status
        LOGICAL :: conversion_ok

        status=0
        q_candidate(1:n_eqns)=q_before(1:n_eqns)-dt*rate
        IF (n_vars>n_eqns) q_candidate(n_eqns+1:n_vars)=q_before(n_eqns+1:n_vars)
        IF (.NOT.conservative_mass_admissible(q_candidate(1:n_eqns))) THEN
            status=1
            RETURN
        END IF

        IF (sea_cell) THEN
            q_candidate=0.0_wp
            CALL initialize_dry_physical_cell(qp_candidate)
            RETURN
        END IF

        IF (q_candidate(idx_totMassEqn)>0.0_wp) THEN
            IF (ALL(q_candidate==q_before)) THEN
                qp_candidate=qp_before
            ELSE
                CALL conservative_to_primitive(q_candidate,qp_candidate,accepted=conversion_ok)
                IF (.NOT.conversion_ok) status=2
            END IF
        ELSE
            CALL initialize_dry_physical_cell(qp_candidate)
        END IF
    END SUBROUTINE update_spatial_cell

    !> @brief Release previously allocated integration workspaces.
    !> @param[in,out] this Integration workspace to release.
    SUBROUTINE finalize_time_integration(this)
        CLASS(time_integration_workspace_type), INTENT(INOUT) :: this
        IF (this%profile_enabled) CALL report_lie_source_profile()
        IF (this%profile_enabled) THEN
            WRITE(*,'(A,I0)') 'LIE_PROFILE steps=',this%profile_calls
            WRITE(*,'(A,ES15.7)') 'LIE_PROFILE spatial_prepare_seconds=',this%profile_prepare
            WRITE(*,'(A,ES15.7)') 'LIE_PROFILE spatial_update_seconds=',this%profile_spatial
            WRITE(*,'(A,ES15.7)') 'LIE_PROFILE source_groups_seconds=',this%profile_sources
            WRITE(*,'(A,ES15.7)') 'LIE_PROFILE source_executor_seconds=',this%profile_source_executor
            WRITE(*,'(A,ES15.7)') 'LIE_PROFILE publish_seconds=',this%profile_publish
        END IF
        CALL finalize_lie_source_executor(this%lie_executor)
        IF (ALLOCATED(this%spatial_rate)) DEALLOCATE(this%spatial_rate)
        IF (ALLOCATED(this%qp_before)) DEALLOCATE(this%qp_before)
        IF (ALLOCATED(this%q_trial)) DEALLOCATE(this%q_trial)
        IF (ALLOCATED(this%qp_trial)) DEALLOCATE(this%qp_trial)
        IF (ALLOCATED(this%cell_status)) DEALLOCATE(this%cell_status)
    END SUBROUTINE finalize_time_integration

END MODULE time_integration_2d
