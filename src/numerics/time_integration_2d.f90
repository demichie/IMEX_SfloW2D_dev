!> @brief Coordinate first-order Lie integration of spatial fluxes and physical sources.
!> @details The spatial operator advances first without rate limiting. Ordered
!>          source groups subsequently modify a transaction candidate; only a
!>          fully accepted timestep publishes fluid, substrate and geometry.
MODULE time_integration_2d
    USE parameters_2d, ONLY : wp,n_eqns,n_vars,idx_totMassEqn,idx_solidEqn_first,idx_solidEqn_last,                  &
            idx_add_gas_first,idx_add_gas_last,sea_mass_removal_flag,n_layers
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
        LOGICAL :: profile_enabled=.FALSE.
        INTEGER :: profile_calls=0
        REAL(wp) :: profile_spatial=0.0_wp,profile_sources=0.0_wp
        REAL(wp) :: profile_source_executor=0.0_wp,profile_publish=0.0_wp
    CONTAINS
        PROCEDURE :: initialize => initialize_time_integration
        PROCEDURE :: finalize => finalize_time_integration
        PROCEDURE :: advance => advance_lie
    END TYPE
CONTAINS
    !> @brief Allocate integration workspaces from the configured domain sizes.
    !> @param[in,out] this Integration workspace to allocate.
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
        ALLOCATE(this%spatial_rate(n_eqns,n_layers,comp_cells_x,comp_cells_y))
        ALLOCATE(this%q_trial(n_vars,n_layers,comp_cells_x,comp_cells_y))
        ALLOCATE(this%qp_trial(n_vars,n_layers,comp_cells_x,comp_cells_y))
        ALLOCATE(this%cell_status(comp_cells_x*comp_cells_y))
        CALL initialize_lie_source_executor(this%lie_executor,source_registry)
    END SUBROUTINE initialize_time_integration

    !> @brief Advance one timestep with spatial update followed by Lie source groups.
    !> @details Reads group-zero physical rates from the pre-spatial snapshot,
    !>          checks spatial admissibility without limiting the spatial rate,
    !>          and commits state and geometry only on success.
    !> @param[in,out] this Owns pre-step and trial workspace arrays.
    !> @param[in,out] state Fluid and substrate state, committed after acceptance.
    !> @param[in] t Start time of timestep [s].
    !> @param[in] dt Positive timestep duration [s].
    !> @param[in] stochastic_values Per-cell stochastic input array.
    !> @param[in] source_registry Ordered Lie source groups.
    !> @param[in] domain Active and reconstruction cell worklists.
    !> @param[in,out] spatial_operators Per-layer HP-PCCU workspaces; only slot 1 executes today.
    !> @param[in] layer_physics Canonical per-layer physical descriptors.
    !> @param[in] pore_caches Derived Pascal coefficients for active layers.
    SUBROUTINE advance_lie(this,state,t,dt,stochastic,domain,spatial_operators,pre_spatial_primitives_ready, &
        stochastic_prepared_hook,layer_physics,pore_caches)
        CLASS(time_integration_workspace_type), TARGET, INTENT(INOUT) :: this
        CLASS(state_type), TARGET, INTENT(INOUT) :: state
        REAL(wp), INTENT(IN) :: t,dt
        CLASS(stochastic_workspace_type), INTENT(INOUT) :: stochastic
        PROCEDURE(stochastic_stage_observer), OPTIONAL :: stochastic_prepared_hook
        CLASS(domain_type), INTENT(IN) :: domain
        CLASS(spatial_operator_type), INTENT(INOUT) :: spatial_operators(:)
        TYPE(layer_physics_type), INTENT(IN) :: layer_physics(:)
        TYPE(pore_pressure_cache_type), INTENT(IN) :: pore_caches(:)
        TYPE(thermal_config_type) :: thermal(n_layers)
        TYPE(rheology_config_type) :: rheology(n_layers)
        TYPE(pore_pressure_config_type) :: pore_pressure(n_layers)
        TYPE(mass_exchange_config_type) :: mass_config(n_layers)
        LOGICAL :: transport_stochastic(n_layers),pore_enabled(n_layers),stochastic_enabled(n_layers)
        LOGICAL, INTENT(IN), OPTIONAL :: pre_spatial_primitives_ready
        LOGICAL :: reuse_cfl_traces
        REAL(wp), POINTER, CONTIGUOUS :: qp_pre(:,:,:),qp_pre2(:,:,:),qp_pre_all(:,:,:,:)
        REAL(wp) :: coupled_dt
        INTEGER :: l,j,k,first_failure,eq,failed_layer,layer_status,layer
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
            transport_stochastic(layer)=.FALSE.
            pore_enabled(layer)=pore_pressure(layer)%pore_pressure_flag
            stochastic_enabled(layer)=rheology(layer)%stochastic_rheology_flag
        END DO
        transport_stochastic(1)=stochastic%config%stochastic_transport_flag
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
        qp_pre=>qp_pre_all(:,1,:,:)
        IF (n_layers==2) qp_pre2=>qp_pre_all(:,2,:,:)
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
        IF (n_layers==2) THEN
            IF (SIZE(spatial_operators)/=2) CALL fatal_error('TL07 coupled Lie needs two spatial operators')
            CALL spatial_operators(1)%evaluate_coupled(spatial_operators(2),qp_pre,qp_pre2, &
                this%spatial_rate(:,1,:,:),this%spatial_rate(:,2,:,:),t,domain,coupled_dt)
            IF (dt>coupled_dt) CALL fatal_error('TL07 timestep exceeds preliminary coupled CFL bound')
        ELSE
            CALL spatial_operators(1)%evaluate(qp_pre,this%spatial_rate(:,1,:,:),t,domain, &
                reuse_cfl_traces=reuse_cfl_traces)
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
            IF (layer_status==1)                                                                                      &
                CALL fatal_error('Lie spatial update: inadmissible conservative state')
            CALL fatal_error('Lie spatial state conversion rejected')

        END IF

        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            this%profile_spatial=this%profile_spatial+REAL(tock-tick,wp)/REAL(tick_rate,wp)
            tick=tock
        END IF
        ! TL07-B must extend the existing Lie groups to L1/L2/INTERLAYER before
        ! the coupled trial may reach the source executor or the atomic commit.
        IF (n_layers==2) CALL fatal_error('TL07 coupled spatial trial ready; coupled Lie sources not yet implemented')

        ! OU is now the post-spatial fractional stage.  q_trial/qp_trial
        ! are valid only on the solve worklist, so stochastic%prepare_timestep
        ! uses the published state for inactive cells. Keep both published
        ! fluid arrays untouched until the final Lie transaction commits.
        IF (rheology(1)%stochastic_rheology_flag .OR. stochastic%config%stochastic_transport_flag) THEN
            CALL stochastic%prepare_timestep(state,dt,stochastic%config%stochastic_transport_flag, &
                this%q_trial(:,1,:,:),this%qp_trial(:,1,:,:),domain)
        END IF
        IF (PRESENT(stochastic_prepared_hook)) CALL stochastic_prepared_hook()
        ! Group 0 evaluates q^n but limits and applies to post-spatial q*.
        ! The geometry transaction supports deferred topography updates;
        ! it is not a recoverable timestep or retry mechanism.
        needs_geometry_update=this%lie_executor%has_bed_group
        IF (needs_geometry_update) CALL begin_geometry_transaction(sparse=.TRUE.)
        CALL execute_lie_source_groups(this%lie_executor,state,this%q_trial,this%qp_trial,dt,domain, &
            stochastic%effective_Z,t,state%q,qp_pre_all,source_step_failed,source_failure, &
            thermal,rheology,pore_pressure,pore_caches,mass_config, &
            transport_stochastic,pore_enabled,stochastic_enabled)
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
        ! Publish the accepted fluid fields only after all Lie source groups return.
        ! The integrator owns evolving inventories; geometry owns only bed caches.
        IF (this%profile_enabled) CALL SYSTEM_CLOCK(tick,tick_rate)
        ! Publish all active layer slots from one accepted transaction, using
        ! only the shared solve worklist. The TL07 coupled source guard above
        ! prevents reaching this commit before both layers are source-complete.
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l,j,k)
        DO l=1,domain%solve_cells
            j=domain%j_cent(l)
            k=domain%k_cent(l)
            state%q(:,:,j,k)=this%q_trial(:,:,j,k)
            state%qp(:,:,j,k)=this%qp_trial(:,:,j,k)
        END DO
        !$OMP END PARALLEL DO
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
