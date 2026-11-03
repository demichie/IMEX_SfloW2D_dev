!> @brief Coordinate first-order Lie integration of spatial fluxes and physical sources.
!> @details The spatial operator advances first without rate limiting. Ordered
!>          source groups subsequently modify a transaction candidate; only a
!>          fully accepted timestep publishes fluid, substrate and geometry.
MODULE time_integration_2d
    USE parameters_2d, ONLY : wp,n_eqns,n_vars,idx_totMassEqn,idx_solidEqn_first,idx_solidEqn_last,                  &
            idx_add_gas_first,idx_add_gas_last,topo_change_flag,sea_mass_removal_flag
    USE geometry_2d, ONLY : comp_cells_x,comp_cells_y,                                                                 &
            begin_geometry_transaction,commit_geometry_transaction,                                                    &
            rollback_geometry_transaction,B_vertex,B_cent,B_face_x,B_face_y,                                           &
            B_prime_x_geom,B_prime_y_geom,grav_coeff,B_zone,B_nodata
    USE state_conversion_2d, ONLY : qc_to_qp
    USE state_2d, ONLY : state_type,substrate_snapshot_type,                                                           &
            initialize_dry_physical_cell,conservative_mass_admissible
    USE diagnostics_2d, ONLY : fatal_error
    USE source_registry_2d, ONLY : source_registry_type
    USE spatial_operator_2d, ONLY : spatial_operator_type
    USE domain_2d, ONLY : domain_type
    USE lie_source_executor_2d, ONLY : execute_lie_source_groups,initialize_lie_source_executor, &
        finalize_lie_source_executor,report_lie_source_profile,lie_source_executor_type
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: time_integration_workspace_type
    !> @brief Scratch and transactional state buffers for one Lie integration.
    TYPE :: time_integration_workspace_type
        !> Spatial residual, dimensioned by conservative equation and cell.
        REAL(wp), ALLOCATABLE :: spatial_rate(:,:,:)
        !> Optional pre-spatial primitive halo when CFL traces cannot be reused.
        REAL(wp), ALLOCATABLE :: qp_before(:,:,:)
        !> Candidate fields: only worklist entries are initialized or published.
        REAL(wp), ALLOCATABLE :: q_trial(:,:,:), qp_trial(:,:,:)
        !> Persistent status worklist for primitive and spatial checks.
        INTEGER, ALLOCATABLE :: cell_status(:)
        TYPE(lie_source_executor_type) :: lie_executor
        ! Retain rollback storage between steps to avoid repeated allocation/free.
        TYPE(substrate_snapshot_type) :: substrate_before
        LOGICAL :: profile_enabled=.FALSE.
        INTEGER :: profile_calls=0
        REAL(wp) :: profile_spatial=0.0_wp,profile_sources=0.0_wp
        REAL(wp) :: profile_state_snapshot=0.0_wp,profile_substrate_snapshot=0.0_wp
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
        ALLOCATE(this%spatial_rate(n_eqns,comp_cells_x,comp_cells_y))
        ALLOCATE(this%q_trial(n_vars,comp_cells_x,comp_cells_y))
        ALLOCATE(this%qp_trial(n_vars+2,comp_cells_x,comp_cells_y))
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
    !> @param[in,out] spatial_operator Spatial flux operator and scratch buffers.
    SUBROUTINE advance_lie(this,state,t,dt,stochastic_values,domain,spatial_operator,pre_spatial_primitives_ready)
        CLASS(time_integration_workspace_type), TARGET, INTENT(INOUT) :: this
        CLASS(state_type), TARGET, INTENT(INOUT) :: state
        REAL(wp), INTENT(IN) :: t,dt,stochastic_values(comp_cells_x,comp_cells_y)
        CLASS(domain_type), INTENT(IN) :: domain
        CLASS(spatial_operator_type), INTENT(INOUT) :: spatial_operator
        LOGICAL, INTENT(IN), OPTIONAL :: pre_spatial_primitives_ready
        LOGICAL :: reuse_cfl_traces
        REAL(wp), POINTER, CONTIGUOUS :: qp_pre(:,:,:)
        INTEGER :: l,j,k,first_failure,eq
        INTEGER :: tick,tock,tick_rate
        LOGICAL :: source_step_failed,conversion_ok
        CHARACTER(LEN=256) :: source_failure
        REAL(wp) :: pd
        LOGICAL :: needs_substrate_snapshot,needs_geometry_snapshot
        REAL(wp), ALLOCATABLE :: q_authoritative(:,:,:),qp_authoritative(:,:,:),bed_before(:,:),                      &
            bed_cent_before(:,:),bed_face_x_before(:,:),bed_face_y_before(:,:),                                      &
            slope_x_before(:,:),slope_y_before(:,:),gravity_before(:,:)
        LOGICAL :: check_injected_rollback,had_bed_change
        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tick,tick_rate)
            this%profile_calls=this%profile_calls+1
        END IF
        ! Retain full pre-step snapshots only for the D20 injected rollback test.
        check_injected_rollback=this%lie_executor%reject_active

        IF (check_injected_rollback) THEN
            q_authoritative=state%q
            qp_authoritative=state%qp
            bed_before=B_vertex
            bed_cent_before=B_cent
            bed_face_x_before=B_face_x
            bed_face_y_before=B_face_y
            slope_x_before=B_prime_x_geom
            slope_y_before=B_prime_y_geom
            gravity_before=grav_coeff

        END IF

        ! The conservative state remains immutable until commit. Trial buffers
        ! are initialized only on active cells, avoiding four full-grid copies.
        ! Production CFL already converted this unchanged pre-spatial state.
        ! Direct callers retain the original preflight unless they explicitly opt in.
        reuse_cfl_traces=.FALSE.
        IF (PRESENT(pre_spatial_primitives_ready)) reuse_cfl_traces=pre_spatial_primitives_ready
        IF (reuse_cfl_traces) THEN
            qp_pre=>state%qp
        ELSE
            ! This full halo buffer exists only for direct/non-CFL callers.
            ! Preserve their original reconstruction-stencil semantics.
            IF (.NOT.ALLOCATED(this%qp_before))                                                                     &
                ALLOCATE(this%qp_before(n_vars+2,comp_cells_x,comp_cells_y))
            this%qp_before=state%qp
            qp_pre=>this%qp_before
        END IF
        IF (.NOT.reuse_cfl_traces) THEN
            ! Each cell initializes its own status in the conversion sweep.
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l,j,k,pd,conversion_ok)
        DO l=1,domain%reconstruction_cells
            this%cell_status(l)=0
            j=domain%j_reconstruction(l)
            k=domain%k_reconstruction(l)

            IF (state%q(idx_totMassEqn,j,k)>0.0_wp) THEN
                CALL qc_to_qp(state%q(:,j,k),this%qp_before(:,j,k),pd,                                                 &
                    'Lie pre-spatial reconstruction state',conversion_ok)
                IF (.NOT.conversion_ok) this%cell_status(l)=1
            ELSE
                CALL initialize_dry_physical_cell(this%qp_before(:,j,k))

            END IF

        END DO
        !$OMP END PARALLEL DO
        first_failure=MINLOC(this%cell_status(1:domain%reconstruction_cells),1,                                      &
            MASK=this%cell_status(1:domain%reconstruction_cells)/=0)
        IF (first_failure/=0)                                                                                          &
            CALL fatal_error('Lie pre-spatial state conversion rejected')
        END IF

        ! The spatial step is conservative; only admissibility is checked.
        CALL spatial_operator%evaluate(qp_pre,this%spatial_rate,t,domain,                                            &
            reuse_cfl_traces=reuse_cfl_traces)
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l,j,k,pd,conversion_ok)
        DO l=1,domain%solve_cells
            this%cell_status(l)=0
            j=domain%j_cent(l)
            k=domain%k_cent(l)
            ! Fuse candidate initialization and conservative update in one sweep.
            this%q_trial(1:n_eqns,j,k)=state%q(1:n_eqns,j,k)-dt*this%spatial_rate(:,j,k)
            IF (n_vars>n_eqns) this%q_trial(n_eqns+1:n_vars,j,k)=state%q(n_eqns+1:n_vars,j,k)
            ! Spatial admissibility check. Do not limit or clip its conservative rate.
            IF (.NOT.conservative_mass_admissible(this%q_trial(1:n_eqns,j,k))) THEN
                this%cell_status(l)=1
                CYCLE

            END IF

            ! A sea cell is an absorbing sink after the spatial step, before ALL Lie sources.
            ! Keep this inside the existing OpenMP sweep and synchronize q/qp.
            IF (B_nodata(j,k) .OR. (sea_mass_removal_flag .AND. B_zone(j,k)==1)) THEN
                this%q_trial(:,j,k)=0.0_wp
                CALL initialize_dry_physical_cell(this%qp_trial(:,j,k))
                CYCLE
            END IF

            IF (this%q_trial(idx_totMassEqn,j,k)>0.0_wp) THEN
                CALL qc_to_qp(this%q_trial(:,j,k),this%qp_trial(:,j,k),pd,'Lie spatial update',conversion_ok)
                IF (.NOT.conversion_ok) this%cell_status(l)=2
            ELSE
                CALL initialize_dry_physical_cell(this%qp_trial(:,j,k))

            END IF

        END DO
        !$OMP END PARALLEL DO
        first_failure=MINLOC(this%cell_status(1:domain%solve_cells),1,                                                &
            MASK=this%cell_status(1:domain%solve_cells)/=0)

        IF (first_failure/=0) THEN
            ! Diagnose the first failing cell outside the OpenMP region.
            ! Do not change, limit or repair the rejected conservative candidate.
            j=domain%j_cent(first_failure)
            k=domain%k_cent(first_failure)
            WRITE(*,'(A,2(I0,1X),A,ES23.15,A,ES23.15)') &
                'LIE_SPATIAL_FAILURE cell(j,k)=',j,k,' t=',t,' dt=',dt
            WRITE(*,'(A,I0)') 'LIE_SPATIAL_FAILURE status=',this%cell_status(first_failure)
            WRITE(*,'(A,ES23.15,A,ES23.15)') &
                'LIE_SPATIAL_FAILURE mass_before=',state%q(idx_totMassEqn,j,k), &
                ' mass_after=',this%q_trial(idx_totMassEqn,j,k)
            WRITE(*,'(A,ES23.15,A,ES23.15)') &
                'LIE_SPATIAL_FAILURE components_before=', &
                SUM(state%q(idx_solidEqn_first:idx_solidEqn_last,j,k))+ &
                SUM(state%q(idx_add_gas_first:idx_add_gas_last,j,k)), &
                ' components_after=', &
                SUM(this%q_trial(idx_solidEqn_first:idx_solidEqn_last,j,k))+ &
                SUM(this%q_trial(idx_add_gas_first:idx_add_gas_last,j,k))
            DO eq=1,n_eqns
                WRITE(*,'(A,I0,A,3(1X,ES23.15))') 'LIE_SPATIAL_FAILURE equation=',eq, &
                    ' before/rate/after=',state%q(eq,j,k), &
                    this%spatial_rate(eq,j,k),this%q_trial(eq,j,k)
                CALL spatial_operator%diagnose_cell(j,k,eq)
            END DO
            IF (this%cell_status(first_failure)==1)                                                                  &
                CALL fatal_error('Lie spatial update: inadmissible conservative state')
            CALL fatal_error('Lie spatial state conversion rejected')

        END IF

        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            this%profile_spatial=this%profile_spatial+REAL(tock-tick,wp)/REAL(tick_rate,wp)
            tick=tock
        END IF
        ! GROUP-0 source rates use immutable q^n/qp^n, while applying to post-spatial q.
        ! Later groups evaluate from each predecessor's updated state.
        ! Geometry stays frozen for all source groups; retain authoritative
        ! fields until the full trial succeeds.
        ! Snapshot only fields a configured source can actually change.
        ! The metadata is immutable within a timestep; no per-cell scan is needed.
        needs_substrate_snapshot=this%lie_executor%has_exchange_group
        needs_geometry_snapshot=this%lie_executor%has_bed_group
        IF (needs_substrate_snapshot) CALL state%snapshot_substrate(this%substrate_before)
        IF (needs_geometry_snapshot) CALL begin_geometry_transaction
        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            this%profile_substrate_snapshot=this%profile_substrate_snapshot+REAL(tock-tick,wp)/REAL(tick_rate,wp)
            tick=tock
        END IF
        CALL execute_lie_source_groups(this%lie_executor,state,this%q_trial,this%qp_trial,dt,domain, &
            stochastic_values,t,state%q,qp_pre,source_step_failed,source_failure)
        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            this%profile_source_executor=this%profile_source_executor+REAL(tock-tick,wp)/REAL(tick_rate,wp)
        END IF

        IF (source_step_failed) THEN
            had_bed_change=.FALSE.
            IF (check_injected_rollback) had_bed_change=ANY(B_vertex/=bed_before)
            IF (needs_substrate_snapshot) CALL state%restore_substrate(this%substrate_before)
            IF (needs_geometry_snapshot) CALL rollback_geometry_transaction

            IF (check_injected_rollback) THEN
                IF (TRIM(source_failure)/='D20_INJECTED_LATE_LIE_FAILURE')                                             &
                    CALL fatal_error('D20 rollback: failpoint did not trigger')
                IF (.NOT.had_bed_change)                                                                               &
                    CALL fatal_error('D20 rollback: failpoint had no prior bed change')
                IF (.NOT.ALL(state%q==q_authoritative) .OR.                                                           &
                    .NOT.ALL(state%qp==qp_authoritative))                                                             &
                    CALL fatal_error('D20 rollback: published fluid state changed')
                IF (.NOT.ALL(state%deposit==this%substrate_before%deposit) .OR.                                             &
                    .NOT.ALL(state%erosion==this%substrate_before%erosion) .OR.                                             &
                    .NOT.ALL(state%erodible==this%substrate_before%erodible) .OR.                                           &
                    .NOT.ALL(state%lava_deposit==this%substrate_before%lava_deposit) .OR.                                   &
                    .NOT.ALL(state%deposit_tot==this%substrate_before%deposit_tot) .OR.                                     &
                    .NOT.ALL(state%erosion_tot==this%substrate_before%erosion_tot))                                         &
                    CALL fatal_error('D20 rollback: physical inventories not restored')
                IF (.NOT.ALL(B_vertex==bed_before) .OR.                                                                &
                    .NOT.ALL(B_cent==bed_cent_before) .OR.                                                             &
                    .NOT.ALL(B_face_x==bed_face_x_before) .OR.                                                         &
                    .NOT.ALL(B_face_y==bed_face_y_before) .OR.                                                         &
                    .NOT.ALL(B_prime_x_geom==slope_x_before) .OR.                                                      &
                    .NOT.ALL(B_prime_y_geom==slope_y_before) .OR.                                                      &
                    .NOT.ALL(grav_coeff==gravity_before))                                                              &
                    CALL fatal_error('D20 rollback: geometry not restored')
                WRITE(*,'(A)') 'D20_FULL_STEP_ROLLBACK_PASS'

            END IF

            CALL fatal_error('Lie timestep rolled back: '//TRIM(source_failure))

        END IF

        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            this%profile_sources=this%profile_sources+REAL(tock-tick,wp)/REAL(tick_rate,wp)
        END IF
        ! Publish the accepted fluid fields only after all Lie source groups return.
        ! The integrator owns evolving inventories; geometry owns only bed caches.
        IF (this%profile_enabled) CALL SYSTEM_CLOCK(tick,tick_rate)
        ! Publish only the active worklist; inactive cells remain unchanged.
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l,j,k)
        DO l=1,domain%solve_cells
            j=domain%j_cent(l)
            k=domain%k_cent(l)
            state%q(:,j,k)=this%q_trial(:,j,k)
            state%qp(:,j,k)=this%qp_trial(:,j,k)
        END DO
        !$OMP END PARALLEL DO
        ! The source executor refreshed substrate aggregates after its cell sweep.
        ! Recomputing them here would traverse the whole domain a second time.
        IF (needs_geometry_snapshot) CALL commit_geometry_transaction
        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            this%profile_publish=this%profile_publish+REAL(tock-tick,wp)/REAL(tick_rate,wp)
        END IF
    END SUBROUTINE advance_lie

    !> @brief Release previously allocated integration workspaces.
    !> @param[in,out] this Integration workspace to release.
    SUBROUTINE finalize_time_integration(this)
        CLASS(time_integration_workspace_type), INTENT(INOUT) :: this
        IF (this%profile_enabled) CALL report_lie_source_profile()
        IF (this%profile_enabled) THEN
            WRITE(*,'(A,I0)') 'LIE_PROFILE steps=',this%profile_calls
            WRITE(*,'(A,ES15.7)') 'LIE_PROFILE spatial_update_seconds=',this%profile_spatial
            WRITE(*,'(A,ES15.7)') 'LIE_PROFILE source_groups_seconds=',this%profile_sources
            WRITE(*,'(A,ES15.7)') 'LIE_PROFILE state_snapshot_seconds=',this%profile_state_snapshot
            WRITE(*,'(A,ES15.7)') 'LIE_PROFILE substrate_geometry_snapshot_seconds=',this%profile_substrate_snapshot
            WRITE(*,'(A,ES15.7)') 'LIE_PROFILE source_executor_seconds=',this%profile_source_executor
            WRITE(*,'(A,ES15.7)') 'LIE_PROFILE publish_seconds=',this%profile_publish
        END IF
        CALL finalize_lie_source_executor(this%lie_executor)
        IF (ALLOCATED(this%spatial_rate)) DEALLOCATE(this%spatial_rate)
        IF (ALLOCATED(this%qp_before)) DEALLOCATE(this%qp_before)
        IF (ALLOCATED(this%q_trial)) DEALLOCATE(this%q_trial)
        IF (ALLOCATED(this%qp_trial)) DEALLOCATE(this%qp_trial)
        IF (ALLOCATED(this%cell_status)) DEALLOCATE(this%cell_status)
        IF (ALLOCATED(this%substrate_before%deposit)) DEALLOCATE(this%substrate_before%deposit)
        IF (ALLOCATED(this%substrate_before%erosion)) DEALLOCATE(this%substrate_before%erosion)
        IF (ALLOCATED(this%substrate_before%erodible)) DEALLOCATE(this%substrate_before%erodible)
        IF (ALLOCATED(this%substrate_before%lava_deposit)) DEALLOCATE(this%substrate_before%lava_deposit)
        IF (ALLOCATED(this%substrate_before%deposit_tot)) DEALLOCATE(this%substrate_before%deposit_tot)
        IF (ALLOCATED(this%substrate_before%erosion_tot)) DEALLOCATE(this%substrate_before%erosion_tot)
    END SUBROUTINE finalize_time_integration

END MODULE time_integration_2d
