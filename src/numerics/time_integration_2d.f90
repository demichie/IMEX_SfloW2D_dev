! First-order Lie integration. GROUP 0: spatial update first; physical rates from q^n.
MODULE time_integration_2d
  USE parameters_2d, ONLY : wp,n_eqns,n_vars,idx_totMassEqn
  USE geometry_2d, ONLY : comp_cells_x,comp_cells_y, &
       begin_lie_geometry_transaction,commit_lie_geometry_transaction, &
       rollback_lie_geometry_transaction,deposit,erosion,erodible,lava_deposit
  USE state_conversion_2d, ONLY : qc_to_qp
  USE state_2d, ONLY : initialize_dry_physical_cell,conservative_mass_admissible
  USE diagnostics_2d, ONLY : fatal_error
  USE source_registry_2d, ONLY : source_registry_type
  USE spatial_operator_2d, ONLY : spatial_operator_type
  USE domain_2d, ONLY : domain_type
  USE lie_source_executor_2d, ONLY : execute_lie_source_groups
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: time_integration_workspace_type
  TYPE :: time_integration_workspace_type
    REAL(wp), ALLOCATABLE :: spatial_rate(:,:,:)
    REAL(wp), ALLOCATABLE :: q_before(:,:,:), qp_before(:,:,:)
    REAL(wp), ALLOCATABLE :: q_trial(:,:,:), qp_trial(:,:,:)
  CONTAINS
    PROCEDURE :: initialize => initialize_time_integration
    PROCEDURE :: finalize => finalize_time_integration
    PROCEDURE :: advance => advance_lie
  END TYPE
CONTAINS
  SUBROUTINE initialize_time_integration(this)
    CLASS(time_integration_workspace_type), INTENT(INOUT) :: this
    ALLOCATE(this%spatial_rate(n_eqns,comp_cells_x,comp_cells_y))
    ALLOCATE(this%q_before(n_vars,comp_cells_x,comp_cells_y))
    ALLOCATE(this%qp_before(n_vars+2,comp_cells_x,comp_cells_y))
    ALLOCATE(this%q_trial(n_vars,comp_cells_x,comp_cells_y))
    ALLOCATE(this%qp_trial(n_vars+2,comp_cells_x,comp_cells_y))
  END SUBROUTINE
  SUBROUTINE finalize_time_integration(this)
    CLASS(time_integration_workspace_type), INTENT(INOUT) :: this
    IF (ALLOCATED(this%spatial_rate)) DEALLOCATE(this%spatial_rate)
    IF (ALLOCATED(this%q_before)) DEALLOCATE(this%q_before)
    IF (ALLOCATED(this%qp_before)) DEALLOCATE(this%qp_before)
    IF (ALLOCATED(this%q_trial)) DEALLOCATE(this%q_trial)
    IF (ALLOCATED(this%qp_trial)) DEALLOCATE(this%qp_trial)
  END SUBROUTINE
  SUBROUTINE advance_lie(this,q,qp,t,dt,Z,source_registry,domain,spatial_operator)
    CLASS(time_integration_workspace_type), INTENT(INOUT) :: this
    REAL(wp), INTENT(INOUT) :: q(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(INOUT) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: t,dt,Z(comp_cells_x,comp_cells_y)
    TYPE(source_registry_type), INTENT(IN) :: source_registry
    CLASS(domain_type), INTENT(IN) :: domain
    CLASS(spatial_operator_type), INTENT(INOUT) :: spatial_operator
    INTEGER :: l,j,k
    LOGICAL :: source_step_failed,conversion_ok
    CHARACTER(LEN=256) :: source_failure
    REAL(wp) :: pd
    REAL(wp), ALLOCATABLE :: deposit_before(:,:,:),erosion_before(:,:,:)
    REAL(wp), ALLOCATABLE :: erodible_before(:,:,:),lava_before(:,:)
    IF (.NOT.source_registry%is_initialized()) CALL fatal_error('Lie registry unavailable')
    this%q_before=q
    this%qp_before=qp
    this%q_trial=q
    this%qp_trial=qp
    ! Reconstruct a defined primitive halo from the conserved pre-step state.
    DO l=1,domain%reconstruction_cells
      j=domain%j_reconstruction(l)
      k=domain%k_reconstruction(l)
      IF (this%q_before(1,j,k)>0.0_wp) THEN
        CALL qc_to_qp(this%q_before(:,j,k),this%qp_before(:,j,k),pd, &
             'Lie pre-spatial reconstruction state',conversion_ok)
        IF (.NOT.conversion_ok) &
             CALL fatal_error('Lie pre-spatial state conversion rejected')
      ELSE
        CALL initialize_dry_physical_cell(this%qp_before(:,j,k))
      END IF
    END DO
    CALL spatial_operator%evaluate(this%qp_before,this%spatial_rate,t,domain)
    DO l=1,domain%solve_cells
      j=domain%j_cent(l)
      k=domain%k_cent(l)
      this%q_trial(1:n_eqns,j,k)=this%q_before(1:n_eqns,j,k)-dt*this%spatial_rate(:,j,k)
      ! Spatial admissibility check. Do not limit or clip its conservative rate.
      IF (.NOT.conservative_mass_admissible(this%q_trial(1:n_eqns,j,k))) &
           CALL fatal_error('Lie spatial update: inadmissible conservative state')
      IF (this%q_trial(idx_totMassEqn,j,k)>0.0_wp) THEN
        CALL qc_to_qp(this%q_trial(:,j,k),this%qp_trial(:,j,k),pd,'Lie spatial update',conversion_ok)
        IF (.NOT.conversion_ok) &
             CALL fatal_error('Lie spatial state conversion rejected')
      ELSE
        CALL initialize_dry_physical_cell(this%qp_trial(:,j,k))
      END IF
    END DO
    ! GROUP-0 source rates use immutable q^n/qp^n, while applying to post-spatial q.
    ! Later groups evaluate from each predecessor's updated state.
    ! The groups may update live geometry for successor groups; retain the
    ! pre-step authoritative fields until the complete Lie sweep succeeds.
    deposit_before=deposit
    erosion_before=erosion
    erodible_before=erodible
    lava_before=lava_deposit
    CALL begin_lie_geometry_transaction
    CALL execute_lie_source_groups(this%q_trial,this%qp_trial,dt,domain, &
         source_registry,Z,t,this%q_before,this%qp_before,source_step_failed,source_failure)
    IF (source_step_failed) THEN
      deposit=deposit_before
      erosion=erosion_before
      erodible=erodible_before
      lava_deposit=lava_before
      CALL rollback_lie_geometry_transaction
      CALL fatal_error('Lie timestep rolled back: '//TRIM(source_failure))
    END IF
    ! Publish the accepted fluid fields only after all Lie source groups return.
    ! The integrator owns evolving inventories; geometry owns only bed caches.
    q=this%q_trial
    qp=this%qp_trial
    CALL commit_lie_geometry_transaction
  END SUBROUTINE
END MODULE time_integration_2d
