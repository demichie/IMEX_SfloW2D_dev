! First-order Lie integration. GROUP 0: spatial update first; physical rates from q^n.
MODULE time_integration_2d
  USE parameters_2d, ONLY : wp,n_eqns,n_vars,idx_solidEqn_first,idx_solidEqn_last, &
       idx_totMassEqn
  USE geometry_2d, ONLY : comp_cells_x,comp_cells_y
  USE state_conversion_2d, ONLY : qc_to_qp
  USE state_2d, ONLY : initialize_dry_physical_cell
  USE diagnostics_2d, ONLY : fatal_error
  USE source_registry_2d, ONLY : source_registry_type
  USE equation_metadata_2d, ONLY : equation_partition_type
  USE spatial_operator_2d, ONLY : spatial_operator_type
  USE domain_2d, ONLY : domain_type
  USE step_source_executor_2d, ONLY : execute_lie_source_groups
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: time_integration_workspace_type
  TYPE :: time_integration_workspace_type
    REAL(wp), ALLOCATABLE :: spatial_rate(:,:,:)
    REAL(wp), ALLOCATABLE :: q_before(:,:,:), qp_before(:,:,:)
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
  END SUBROUTINE
  SUBROUTINE finalize_time_integration(this)
    CLASS(time_integration_workspace_type), INTENT(INOUT) :: this
    IF (ALLOCATED(this%spatial_rate)) DEALLOCATE(this%spatial_rate)
    IF (ALLOCATED(this%q_before)) DEALLOCATE(this%q_before)
    IF (ALLOCATED(this%qp_before)) DEALLOCATE(this%qp_before)
  END SUBROUTINE
  SUBROUTINE advance_lie(this,q,qp,t,dt,Z,equation_partition,source_registry,domain,spatial_operator)
    CLASS(time_integration_workspace_type), INTENT(INOUT) :: this
    REAL(wp), INTENT(INOUT) :: q(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(INOUT) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: t,dt,Z(comp_cells_x,comp_cells_y)
    TYPE(equation_partition_type), INTENT(IN) :: equation_partition
    TYPE(source_registry_type), INTENT(IN) :: source_registry
    CLASS(domain_type), INTENT(IN) :: domain
    CLASS(spatial_operator_type), INTENT(INOUT) :: spatial_operator
    INTEGER :: l,j,k
    REAL(wp) :: pd
    IF (.NOT.source_registry%is_initialized()) CALL fatal_error('Lie registry unavailable')
    this%q_before=q
    this%qp_before=qp
    ! Reconstruct a defined primitive halo from the conserved pre-step state.
    DO l=1,domain%reconstruction_cells
      j=domain%j_reconstruction(l)
      k=domain%k_reconstruction(l)
      IF (this%q_before(1,j,k)>0.0_wp) THEN
        CALL qc_to_qp(this%q_before(:,j,k),this%qp_before(:,j,k),pd, &
             'Lie pre-spatial reconstruction state')
      ELSE
        CALL initialize_dry_physical_cell(this%qp_before(:,j,k))
      END IF
    END DO
    CALL spatial_operator%evaluate(this%qp_before,this%spatial_rate,t,domain)
    DO l=1,domain%solve_cells
      j=domain%j_cent(l)
      k=domain%k_cent(l)
      q(1:n_eqns,j,k)=this%q_before(1:n_eqns,j,k)-dt*this%spatial_rate(:,j,k)
      ! Spatial admissibility check. Do not limit or clip its conservative rate.
      IF (.NOT.ALL(q(1:n_eqns,j,k)==q(1:n_eqns,j,k))) &
           CALL fatal_error('Lie spatial update produced NaN')
      IF (q(idx_totMassEqn,j,k)<0.0_wp) &
           CALL fatal_error('Lie spatial update: negative mass')
      IF (ANY(q(idx_solidEqn_first:idx_solidEqn_last,j,k)<0.0_wp)) &
           CALL fatal_error('Lie spatial update: negative solid mass')
      IF (SUM(q(idx_solidEqn_first:idx_solidEqn_last,j,k))> &
           q(idx_totMassEqn,j,k)) &
           CALL fatal_error('Lie spatial update: solids exceed total mass')
      IF (q(idx_totMassEqn,j,k)>0.0_wp) THEN
        CALL qc_to_qp(q(:,j,k),qp(:,j,k),pd,'Lie spatial update')
      ELSE
        CALL initialize_dry_physical_cell(qp(:,j,k))
      END IF
    END DO
    ! GROUP-0 source rates use immutable q^n/qp^n, while applying to post-spatial q.
    ! Later groups evaluate from each predecessor's updated state.
    CALL execute_lie_source_groups(q,qp,dt,domain,source_registry,Z,t, &
         this%q_before,this%qp_before)
  END SUBROUTINE
END MODULE time_integration_2d
