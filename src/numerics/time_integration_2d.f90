!> @brief Coordinate first-order Lie integration of spatial fluxes and physical sources.
!> @details The spatial operator advances first without rate limiting. Ordered
!>          source groups subsequently modify a transaction candidate; only a
!>          fully accepted timestep publishes fluid, substrate and geometry.
MODULE time_integration_2d
    USE parameters_2d, ONLY : wp,n_eqns,n_vars,idx_totMassEqn
    USE geometry_2d, ONLY : comp_cells_x,comp_cells_y, &
            begin_geometry_transaction,commit_geometry_transaction, &
            rollback_geometry_transaction
    USE state_conversion_2d, ONLY : qc_to_qp
    USE state_2d, ONLY : state_type,substrate_snapshot_type, &
            initialize_dry_physical_cell,conservative_mass_admissible
    USE diagnostics_2d, ONLY : fatal_error
    USE source_registry_2d, ONLY : source_registry_type
    USE spatial_operator_2d, ONLY : spatial_operator_type
    USE domain_2d, ONLY : domain_type
    USE lie_source_executor_2d, ONLY : execute_lie_source_groups
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: time_integration_workspace_type
    !> @brief Scratch and transactional state buffers for one Lie integration.
    TYPE :: time_integration_workspace_type
        !> Spatial residual, dimensioned by conservative equation and cell.
        REAL(wp), ALLOCATABLE :: spatial_rate(:,:,:)
        !> Immutable conservative and primitive pre-step snapshots.
        REAL(wp), ALLOCATABLE :: q_before(:,:,:), qp_before(:,:,:)
        !> Candidate cell states, published only after all groups succeed.
        REAL(wp), ALLOCATABLE :: q_trial(:,:,:), qp_trial(:,:,:)
    CONTAINS
        PROCEDURE :: initialize => initialize_time_integration
        PROCEDURE :: finalize => finalize_time_integration
        PROCEDURE :: advance => advance_lie
    END TYPE
CONTAINS
    !> @brief Allocate integration workspaces from the configured domain sizes.
    !> @param[in,out] this Integration workspace to allocate.
    SUBROUTINE initialize_time_integration(this)
        CLASS(time_integration_workspace_type), INTENT(INOUT) :: this
        ALLOCATE(this%spatial_rate(n_eqns,comp_cells_x,comp_cells_y))
        ALLOCATE(this%q_before(n_vars,comp_cells_x,comp_cells_y))
        ALLOCATE(this%qp_before(n_vars+2,comp_cells_x,comp_cells_y))
        ALLOCATE(this%q_trial(n_vars,comp_cells_x,comp_cells_y))
        ALLOCATE(this%qp_trial(n_vars+2,comp_cells_x,comp_cells_y))
    END SUBROUTINE initialize_time_integration

    !> @brief Release previously allocated integration workspaces.
    !> @param[in,out] this Integration workspace to release.
    SUBROUTINE finalize_time_integration(this)
        CLASS(time_integration_workspace_type), INTENT(INOUT) :: this
        IF (ALLOCATED(this%spatial_rate)) DEALLOCATE(this%spatial_rate)
        IF (ALLOCATED(this%q_before)) DEALLOCATE(this%q_before)
        IF (ALLOCATED(this%qp_before)) DEALLOCATE(this%qp_before)
        IF (ALLOCATED(this%q_trial)) DEALLOCATE(this%q_trial)
        IF (ALLOCATED(this%qp_trial)) DEALLOCATE(this%qp_trial)
    END SUBROUTINE finalize_time_integration

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
    SUBROUTINE advance_lie(this,state,t,dt,stochastic_values,source_registry,domain,spatial_operator)
        CLASS(time_integration_workspace_type), INTENT(INOUT) :: this
        CLASS(state_type), INTENT(INOUT) :: state
        REAL(wp), INTENT(IN) :: t,dt,stochastic_values(comp_cells_x,comp_cells_y)
        TYPE(source_registry_type), INTENT(IN) :: source_registry
        CLASS(domain_type), INTENT(IN) :: domain
        CLASS(spatial_operator_type), INTENT(INOUT) :: spatial_operator
        INTEGER :: l,j,k,first_failure
        INTEGER, ALLOCATABLE :: reconstruction_status(:),spatial_status(:)
        LOGICAL :: source_step_failed,conversion_ok
        CHARACTER(LEN=256) :: source_failure
        REAL(wp) :: pd
        TYPE(substrate_snapshot_type) :: substrate_before
        IF (.NOT.source_registry%is_initialized()) CALL fatal_error('Lie registry unavailable')
        ! Snapshot before evaluating any operator, including the primitive halo.
        this%q_before=state%q
        this%qp_before=state%qp
        this%q_trial=state%q
        this%qp_trial=state%qp
        ! Each cell writes a distinct primitive slice; collect failures by cell.
        ALLOCATE(reconstruction_status(domain%reconstruction_cells))
        reconstruction_status=0
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l,j,k,pd,conversion_ok)
        DO l=1,domain%reconstruction_cells
            j=domain%j_reconstruction(l)
            k=domain%k_reconstruction(l)
            IF (this%q_before(idx_totMassEqn,j,k)>0.0_wp) THEN
                CALL qc_to_qp(this%q_before(:,j,k),this%qp_before(:,j,k),pd, &
                    'Lie pre-spatial reconstruction state',conversion_ok)
                IF (.NOT.conversion_ok) reconstruction_status(l)=1
            ELSE
                CALL initialize_dry_physical_cell(this%qp_before(:,j,k))
            END IF
        END DO
        !$OMP END PARALLEL DO
        first_failure=MINLOC(reconstruction_status,1,MASK=reconstruction_status/=0)
        IF (first_failure/=0) &
            CALL fatal_error('Lie pre-spatial state conversion rejected')
        DEALLOCATE(reconstruction_status)

        ! The spatial step is conservative; only admissibility is checked.
        CALL spatial_operator%evaluate(this%qp_before,this%spatial_rate,t,domain)
        ALLOCATE(spatial_status(domain%solve_cells))
        spatial_status=0
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l,j,k,pd,conversion_ok)
        DO l=1,domain%solve_cells
            j=domain%j_cent(l)
            k=domain%k_cent(l)
            this%q_trial(1:n_eqns,j,k)=this%q_before(1:n_eqns,j,k)-dt*this%spatial_rate(:,j,k)
            ! Spatial admissibility check. Do not limit or clip its conservative rate.
            IF (.NOT.conservative_mass_admissible(this%q_trial(1:n_eqns,j,k))) THEN
                spatial_status(l)=1
                CYCLE
            END IF
            IF (this%q_trial(idx_totMassEqn,j,k)>0.0_wp) THEN
                CALL qc_to_qp(this%q_trial(:,j,k),this%qp_trial(:,j,k),pd,'Lie spatial update',conversion_ok)
                IF (.NOT.conversion_ok) spatial_status(l)=2
            ELSE
                CALL initialize_dry_physical_cell(this%qp_trial(:,j,k))
            END IF
        END DO
        !$OMP END PARALLEL DO
        first_failure=MINLOC(spatial_status,1,MASK=spatial_status/=0)
        IF (first_failure/=0) THEN
            IF (spatial_status(first_failure)==1) &
                CALL fatal_error('Lie spatial update: inadmissible conservative state')
            CALL fatal_error('Lie spatial state conversion rejected')
        END IF
        DEALLOCATE(spatial_status)

        ! GROUP-0 source rates use immutable q^n/qp^n, while applying to post-spatial q.
        ! Later groups evaluate from each predecessor's updated state.
        ! The groups may update live geometry for successor groups; retain the
        ! pre-step authoritative fields until the complete Lie sweep succeeds.
        ! Stage substrate and bed changes before running ordered groups.
        CALL state%snapshot_substrate(substrate_before)
        CALL begin_geometry_transaction
        CALL execute_lie_source_groups(state,this%q_trial,this%qp_trial,dt,domain, &
            source_registry,stochastic_values,t,this%q_before,this%qp_before, &
            source_step_failed,source_failure)
        IF (source_step_failed) THEN
            CALL state%restore_substrate(substrate_before)
            CALL rollback_geometry_transaction
            CALL fatal_error('Lie timestep rolled back: '//TRIM(source_failure))
        END IF

        ! Publish the accepted fluid fields only after all Lie source groups return.
        ! The integrator owns evolving inventories; geometry owns only bed caches.
        state%q=this%q_trial
        state%qp=this%qp_trial
        CALL state%refresh_substrate_totals()
        CALL commit_geometry_transaction
    END SUBROUTINE advance_lie

END MODULE time_integration_2d
