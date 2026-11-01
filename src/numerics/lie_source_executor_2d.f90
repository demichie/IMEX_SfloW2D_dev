!********************************************************************************
!> \brief Execute ordered first-order Lie physical source groups.
!> \details Coordinates simultaneous cell-local rate evaluation, physical
!>          limitation and atomic fluid/substrate transfer. Rates come from source_terms_2d,
!>          while bed projection and geometry refresh remain in geometry_2d.
!********************************************************************************

MODULE lie_source_executor_2d
    USE parameters_2d, ONLY : idx_totMassEqn, idx_solidEqn_first, idx_solidEqn_last

  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
  USE diagnostics_2d, ONLY : fatal_error
  USE parameters_2d, ONLY : wp, n_eqns, n_vars, n_solid, mass_exchange_rates_type
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE domain_2d, ONLY : domain_type
  USE source_registry_2d, ONLY : source_registry_type, source_group_type
  USE state_2d, ONLY : state_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: execute_lie_source_groups

CONTAINS

  !> Evaluate and commit one complete Lie source group at a time.
  !> For GROUP 0, evaluate every physical rate at the original pre-spatial
  !> state, but cap and apply the combined source increment to the post-FV state.
  SUBROUTINE execute_lie_source_groups(state,q,qp,time_step,domain,source_registry, &
       stochastic_values,step_start_time,q_original,qp_original,step_failed,failure_message)
    USE geometry_2d, ONLY : B_prime_x_geom,B_prime_y_geom, &
         B_second_xx_geom,B_second_xy_geom,B_second_yy_geom,grav_coeff, &
         project_cell_field_to_vertices,assign_topography_vertices,B_vertex, &
         comp_interfaces_x,comp_interfaces_y
    USE source_terms_2d, ONLY : assemble_lie_explicit_transaction, &
         eval_source_analytic,cell_analytic_source_context, &
         analytic_term_context_type
    USE source_catalog_2d, ONLY : SOURCE_METHOD_ANALYTIC,SOURCE_METHOD_EXPLICIT, &
         SOURCE_METHOD_NEWTON
    USE nonlinear_solver_2d, ONLY : solve_backward_euler_group
    USE state_conversion_2d, ONLY : qc_to_qp
    USE state_2d, ONLY : initialize_dry_physical_cell,conservative_mass_admissible, &
         physical_cell_state_finite
    USE parameters_2d, ONLY : topo_change_flag,erodible_deposit_flag
    CLASS(state_type), INTENT(INOUT) :: state
    TYPE(source_group_type), ALLOCATABLE :: groups(:)
    TYPE(mass_exchange_rates_type) :: transfers
    TYPE(analytic_term_context_type) :: analytic_context
    REAL(wp), INTENT(INOUT) :: q(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(INOUT) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: time_step
    CLASS(domain_type), INTENT(IN) :: domain
    CLASS(source_registry_type), INTENT(IN) :: source_registry
    REAL(wp), INTENT(IN), OPTIONAL :: stochastic_values(:,:)
    REAL(wp), INTENT(IN), OPTIONAL :: step_start_time
    REAL(wp), INTENT(IN) :: q_original(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: qp_original(n_vars+2,comp_cells_x,comp_cells_y)
    LOGICAL, INTENT(OUT) :: step_failed
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: failure_message
    REAL(wp) :: qpentry(n_vars+2),q_candidate(n_vars),qp_candidate(n_vars+2)
    REAL(wp), ALLOCATABLE :: q_group(:,:,:),qp_group(:,:,:)
    REAL(wp) :: deposit_candidate(n_solid),erosion_candidate(n_solid)
    REAL(wp) :: erodible_candidate(n_solid),lava_candidate
    INTEGER, ALLOCATABLE :: group_j(:),group_k(:)
    REAL(wp), ALLOCATABLE :: bed_rate_cell(:,:),bed_rate_vertex(:,:),bed_candidate(:,:)
    REAL(wp), ALLOCATABLE :: staged_deposit(:,:),staged_erosion(:,:), &
         staged_erodible(:,:),staged_lava(:),staged_bed_rate(:)
    REAL(wp) :: total_rate(n_eqns)
    REAL(wp) :: analytic_rate(n_eqns),pressure,source_time,noise
    LOGICAL :: with_exchange,ok,failed,bed_accepted,conversion_ok
    INTEGER :: group_index,cell,j,k,status,iterations,linear_info

    CHARACTER(LEN=96) :: context
    CHARACTER(LEN=256) :: failure_reason

    step_failed=.FALSE.
    failure_reason=''
    IF (PRESENT(failure_message)) failure_message=''
    IF (.NOT.source_registry%is_initialized()) &
         CALL fatal_error('D19: uninitialized Lie source registry')
    groups=source_registry%groups()
    source_time=0.0_wp
    IF (PRESENT(step_start_time)) source_time=step_start_time
    IF (.NOT.ieee_is_finite(time_step)) CALL fatal_error('D19: non-finite Lie timestep')
    IF (time_step<=0.0_wp) CALL fatal_error('D19: nonpositive time interval')

    DO group_index=1,SIZE(groups)
       q_group=q
       qp_group=qp
       with_exchange=ANY(INDEX(groups(group_index)%source_names,'mass_exchange_')==1)
       IF (SIZE(groups(group_index)%source_names)==0) CYCLE
       IF (with_exchange) THEN
          IF (groups(group_index)%method/=SOURCE_METHOD_EXPLICIT) &
               CALL fatal_error('D19: mass exchange requires explicit group method')
          ALLOCATE(group_j(domain%solve_cells),group_k(domain%solve_cells))
          ALLOCATE(staged_deposit(n_solid,domain%solve_cells), &
               staged_erosion(n_solid,domain%solve_cells), &
               staged_erodible(n_solid,domain%solve_cells), &
               staged_lava(domain%solve_cells),staged_bed_rate(domain%solve_cells))
          IF (topo_change_flag) THEN
             ALLOCATE(bed_rate_cell(comp_cells_x,comp_cells_y))
             bed_rate_cell=0.0_wp
          END IF
       END IF

       cell_loop: DO cell=1,domain%solve_cells
          j=domain%j_cent(cell)
          k=domain%k_cent(cell)
          IF (groups(group_index)%id==0) THEN
             qpentry=qp_original(:,j,k)
          ELSE
             qpentry=qp_group(:,j,k)
          END IF
          context='Lie source group'
          noise=0.0_wp
          IF (PRESENT(stochastic_values)) noise=stochastic_values(j,k)

          ! Validate constitutive inputs before dispatch. Model/method errors
          ! remain configuration errors, not recoverable timestep rejections.
          IF (.NOT.physical_cell_state_finite(qpentry) .OR. &
               .NOT.ieee_is_finite(noise)) THEN
             CALL reject_cell('D19: non-finite source-group physical input')
             EXIT cell_loop
          END IF
          IF (.NOT.ieee_is_finite(B_prime_x_geom(j,k)) .OR. &
               .NOT.ieee_is_finite(B_prime_y_geom(j,k)) .OR. &
               .NOT.ieee_is_finite(B_second_xx_geom(j,k)) .OR. &
               .NOT.ieee_is_finite(B_second_xy_geom(j,k)) .OR. &
               .NOT.ieee_is_finite(B_second_yy_geom(j,k)) .OR. &
               .NOT.ieee_is_finite(grav_coeff(j,k))) THEN
             CALL reject_cell('D19: non-finite source-group geometry')
             EXIT cell_loop
          END IF
          IF (with_exchange) THEN
             ! All immobilized and erodible inventories are nonnegative.
             ! Validate the group input before computing any exchange rate.
             IF (.NOT.ALL(ieee_is_finite(state%erodible(1:n_solid,j,k))) .OR. &
                  .NOT.ALL(ieee_is_finite(state%deposit(j,k,1:n_solid))) .OR. &
                  .NOT.ALL(ieee_is_finite(state%erosion(j,k,1:n_solid))) .OR. &
                  .NOT.ieee_is_finite(state%lava_deposit(j,k)) .OR. &
                  ANY(state%erodible(1:n_solid,j,k)<0.0_wp) .OR. &
                  ANY(state%deposit(j,k,1:n_solid)<0.0_wp) .OR. &
                  ANY(state%erosion(j,k,1:n_solid)<0.0_wp) .OR. &
                  state%lava_deposit(j,k)<0.0_wp) THEN
                CALL reject_cell('D19: invalid substrate inventory before mass exchange')
                EXIT cell_loop
             END IF
          END IF

          SELECT CASE(groups(group_index)%method)
          CASE(SOURCE_METHOD_EXPLICIT)
             CALL assemble_lie_explicit_transaction(j,k,qpentry,q_group(1:n_eqns,j,k), &
                  time_step,source_time,groups(group_index),state%erodible(1:n_solid,j,k), &
                  total_rate,transfers,with_exchange,ok)
             IF (.NOT.ok) THEN
                CALL reject_cell('D19: physical source transaction rejected')
                EXIT cell_loop
             END IF
             ! Stage the entire local fluid/substrate transfer before committing it.
             q_candidate=q_group(:,j,k)
             q_candidate(1:n_eqns)=q_candidate(1:n_eqns)+time_step*total_rate
             IF (.NOT.conservative_mass_admissible(q_candidate(1:n_eqns))) THEN
                CALL reject_cell('D19: inadmissible conservative mass candidate')
                EXIT cell_loop
             END IF
             IF (with_exchange) THEN
                deposit_candidate=state%deposit(j,k,1:n_solid)+ &
                     time_step*transfers%solid_deposition_rate
                erosion_candidate=state%erosion(j,k,1:n_solid)+ &
                     time_step*transfers%solid_erosion_rate
                erodible_candidate=state%erodible(1:n_solid,j,k)- &
                     time_step*transfers%solid_erosion_rate
                IF (erodible_deposit_flag) erodible_candidate= &
                     erodible_candidate+time_step*transfers%solid_deposition_rate
                lava_candidate=state%lava_deposit(j,k)+time_step*transfers%lava_deposit_rate
                IF (.NOT.ALL(ieee_is_finite(deposit_candidate)) .OR. &
                     .NOT.ALL(ieee_is_finite(erosion_candidate)) .OR. &
                     .NOT.ALL(ieee_is_finite(erodible_candidate)) .OR. &
                     .NOT.ieee_is_finite(lava_candidate) .OR. &
                     ANY(erodible_candidate<0.0_wp) .OR. &
                     ANY(deposit_candidate<0.0_wp) .OR. &
                     ANY(erosion_candidate<0.0_wp) .OR. &
                     lava_candidate<0.0_wp) THEN
                   CALL reject_cell('D19: invalid candidate substrate inventory')
                   EXIT cell_loop
                END IF
             END IF
             IF (q_candidate(idx_totMassEqn)>0.0_wp) THEN
                CALL qc_to_qp(q_candidate,qp_candidate,pressure,'Lie explicit group candidate',conversion_ok)
                IF (.NOT.conversion_ok) THEN
                   CALL reject_cell('D19: explicit candidate conversion failed')
                   EXIT cell_loop
                END IF
             ELSE
                CALL initialize_dry_physical_cell(qp_candidate)
             END IF
             ! Commit only after the complete cell-local candidate is admissible.
             q_group(:,j,k)=q_candidate
             qp_group(:,j,k)=qp_candidate
             IF (with_exchange) THEN
                group_j(cell)=j
                group_k(cell)=k
                staged_deposit(:,cell)=deposit_candidate
                staged_erosion(:,cell)=erosion_candidate
                staged_erodible(:,cell)=erodible_candidate
                staged_lava(cell)=lava_candidate
                staged_bed_rate(cell)=transfers%dB_dt
             END IF
          CASE(SOURCE_METHOD_ANALYTIC)
             IF (q_group(idx_totMassEqn,j,k)<=0.0_wp) CYCLE
             analytic_context=cell_analytic_source_context(j,k,noise)
             CALL eval_source_analytic(groups(group_index),q_group(:,j,k), &
                  analytic_context,time_step,analytic_rate,status)
             IF (status/=0) THEN
                WRITE(failure_reason,'(A,I0)') 'D19: analytic group status ',status
                step_failed=.TRUE.
                EXIT cell_loop
             END IF
             IF (.NOT.ALL(ieee_is_finite(analytic_rate))) THEN
                CALL reject_cell('D19: analytic source returned non-finite rate')
                EXIT cell_loop
             END IF
             q(1:n_vars,j,k)=q(1:n_vars,j,k)+time_step*analytic_rate(1:n_vars)
          CASE(SOURCE_METHOD_NEWTON)
             IF (q_group(idx_totMassEqn,j,k)<=0.0_wp) CYCLE
             CALL solve_backward_euler_group(groups(group_index), &
                  q_group(:,j,k),time_step,B_prime_x_geom(j,k),B_prime_y_geom(j,k), &
                  noise,iterations,ok,linear_info,failed)
             IF (.NOT.ok .OR. failed .OR. linear_info/=0) THEN
                WRITE(failure_reason,'(A,I0,A,I0)') &
                     'D19: Newton failed, linear_info=',linear_info,', iterations=',iterations
                step_failed=.TRUE.
                EXIT cell_loop
             END IF
          CASE DEFAULT
             CALL reject_cell('D19: invalid group method')
             IF (step_failed) EXIT cell_loop
          END SELECT
          ! Only check admissibility: never clip the accepted source update.
          IF (.NOT.conservative_mass_admissible(q_group(1:n_eqns,j,k))) THEN
             CALL reject_cell('D19: inadmissible conservative source group state')
             EXIT cell_loop
          END IF
          IF (q_group(idx_totMassEqn,j,k)>0.0_wp) THEN
             CALL qc_to_qp(q_group(:,j,k),qp_group(:,j,k),pressure,TRIM(context),conversion_ok)
             IF (.NOT.conversion_ok) THEN
                CALL reject_cell('D19: source group conversion failed')
                EXIT cell_loop
             END IF
          ELSE
             CALL initialize_dry_physical_cell(qp_group(:,j,k))
          END IF
       END DO cell_loop
       IF (step_failed) THEN
          IF (PRESENT(failure_message)) failure_message=TRIM(failure_reason)
          RETURN
       END IF
       ! Publish inventory changes only after all cells in the group succeed.
       IF (with_exchange .AND. topo_change_flag) THEN
          ALLOCATE(bed_rate_vertex(comp_interfaces_x,comp_interfaces_y))
          CALL record_group_bed_rates
          bed_accepted=ALL(ieee_is_finite(bed_rate_cell))
          IF (bed_accepted) THEN
             CALL project_cell_field_to_vertices(bed_rate_cell,bed_rate_vertex)
             bed_accepted=ALL(ieee_is_finite(bed_rate_vertex))
          END IF
          IF (bed_accepted) THEN
             bed_candidate=B_vertex+time_step*bed_rate_vertex
             IF (.NOT.ALL(ieee_is_finite(bed_candidate))) bed_accepted=.FALSE.
          END IF
          IF (.NOT.bed_accepted) THEN
             CALL reject_cell('D19: invalid candidate projected bed')
             IF (PRESENT(failure_message)) failure_message=TRIM(failure_reason)
             RETURN
          END IF
       END IF
       IF (with_exchange) THEN
          ! Accept the geometric projection before publishing inventory changes.
          ! A rejected geometric refresh therefore leaves both inventories
          ! and the conserved/primitive group state unmodified.
          IF (topo_change_flag) THEN
             CALL assign_topography_vertices(bed_candidate,bed_accepted)
             IF (.NOT.bed_accepted) THEN
                CALL reject_cell('D19: invalid refreshed geometry')
                IF (PRESENT(failure_message)) failure_message=TRIM(failure_reason)
                RETURN
             END IF
             DEALLOCATE(bed_rate_cell,bed_rate_vertex,bed_candidate)
          END IF
          DO cell=1,domain%solve_cells
             j=group_j(cell)
             k=group_k(cell)
             state%deposit(j,k,1:n_solid)=staged_deposit(:,cell)
             state%erosion(j,k,1:n_solid)=staged_erosion(:,cell)
             state%erodible(1:n_solid,j,k)=staged_erodible(:,cell)
             state%lava_deposit(j,k)=staged_lava(cell)
          END DO
          DEALLOCATE(group_j,group_k,staged_deposit,staged_erosion, &
               staged_erodible,staged_lava,staged_bed_rate)
       END IF
       q=q_group
       qp=qp_group
    END DO
  CONTAINS
    SUBROUTINE record_group_bed_rates
      INTEGER :: n,j_cell,k_cell
      DO n=1,domain%solve_cells
         j_cell=group_j(n)
         k_cell=group_k(n)
         bed_rate_cell(j_cell,k_cell)=staged_bed_rate(n)
      END DO
    END SUBROUTINE record_group_bed_rates

    ! Return a recoverable cell/group rejection without leaving the executor
    ! before returning the source-group status.
    SUBROUTINE reject_cell(message)
      CHARACTER(*), INTENT(IN) :: message
      failure_reason=message
      step_failed=.TRUE.
    END SUBROUTINE reject_cell
  END SUBROUTINE execute_lie_source_groups

END MODULE lie_source_executor_2d
