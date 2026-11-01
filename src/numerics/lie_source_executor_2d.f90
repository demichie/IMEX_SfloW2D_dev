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

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: execute_lie_source_groups

CONTAINS

  !> Apply the same physical coverage to every rate of one exchange transaction.
  !> Fluid conservative source rates are scaled separately by the caller.
  SUBROUTINE scale_mass_exchange_transfer(transfer,coverage)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: transfer
    REAL(wp), INTENT(IN) :: coverage
    transfer%solid_deposition_rate=coverage*transfer%solid_deposition_rate
    transfer%solid_erosion_rate=coverage*transfer%solid_erosion_rate
    transfer%solid_mass_rate=coverage*transfer%solid_mass_rate
    transfer%add_gas_mass_rate=coverage*transfer%add_gas_mass_rate
    transfer%mixture_mass_rate=coverage*transfer%mixture_mass_rate
    transfer%energy_rate=coverage*transfer%energy_rate
    transfer%carrier_erosion_rate=coverage*transfer%carrier_erosion_rate
    transfer%carrier_deposition_rate=coverage*transfer%carrier_deposition_rate
    transfer%carrier_withdrawal_volume_rate=coverage*transfer%carrier_withdrawal_volume_rate
    transfer%lava_deposit_rate=coverage*transfer%lava_deposit_rate
    transfer%dB_dt=coverage*transfer%dB_dt
  END SUBROUTINE scale_mass_exchange_transfer

  !> Evaluate and commit one complete Lie source group at a time.
  !> For GROUP 0, evaluate every physical rate at the original pre-spatial
  !> state, but cap and apply the combined source increment to the post-FV state.
  SUBROUTINE execute_lie_source_groups(q,qp,time_step,domain,source_registry, &
       stochastic_values,step_start_time,q_original,qp_original,step_failed,failure_message)
    USE parameters_2d, ONLY : erodible_deposit_flag
    USE geometry_2d, ONLY : B_prime_x_geom,B_prime_y_geom,cell_source_fractions, &
         B_second_xx_geom,B_second_xy_geom,B_second_yy_geom,grav_coeff, &
         begin_topography_rate_update,record_topography_rate_cell, &
         apply_topography_rate_update,deposit,erosion,erodible,B_zone,lava_deposit
    USE source_terms_2d, ONLY : eval_source_explicit,eval_source_analytic, &
         eval_cell_source_explicit,cell_analytic_source_context, &
         analytic_term_context_type
    USE source_catalog_2d, ONLY : SOURCE_METHOD_ANALYTIC,SOURCE_METHOD_EXPLICIT, &
         SOURCE_METHOD_NEWTON
    USE nonlinear_solver_2d, ONLY : solve_backward_euler_group
    USE state_conversion_2d, ONLY : qc_to_qp
    USE state_2d, ONLY : initialize_dry_physical_cell,conservative_mass_admissible
    USE parameters_2d, ONLY : topo_change_flag,bottom_radial_source_flag, &
         bottom_fissural_source_flag
    TYPE(source_group_type), ALLOCATABLE :: groups(:)
    TYPE(source_group_type) :: other_sources
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
    REAL(wp) :: deposit_candidate(n_solid),erosion_candidate(n_solid)
    REAL(wp) :: erodible_candidate(n_solid),lava_candidate
    REAL(wp) :: physical_rate(n_eqns),exchange_rate(n_eqns),total_rate(n_eqns)
    REAL(wp) :: analytic_rate(n_eqns),factor,coverage,pressure,source_time,noise
    REAL(wp) :: remaining_mass,change_in_solid
    LOGICAL :: with_exchange,has_other,ok,failed,bed_accepted,conversion_ok
    INTEGER :: group_index,cell,j,k,i,solid_eq,status,iterations,linear_info

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
       other_sources=groups(group_index)
       with_exchange=ANY(INDEX(other_sources%source_names,'mass_exchange_')==1)
       IF (with_exchange) THEN
          IF (groups(group_index)%method/=SOURCE_METHOD_EXPLICIT) &
               CALL fatal_error('D19: mass exchange must use its supported explicit method')
          other_sources%source_names=PACK(other_sources%source_names, &
               INDEX(other_sources%source_names,'mass_exchange_')/=1)
       END IF
       has_other=SIZE(other_sources%source_names)>0
       IF (.NOT.has_other .AND. .NOT.with_exchange) CYCLE

       IF (with_exchange .AND. topo_change_flag) CALL begin_topography_rate_update

       cell_loop: DO cell=1,domain%solve_cells
          j=domain%j_cent(cell)
          k=domain%k_cent(cell)
          IF (groups(group_index)%id==0) THEN
             qpentry=qp_original(:,j,k)
          ELSE
             qpentry=qp(:,j,k)
          END IF
          context='Lie source group'
          noise=0.0_wp
          IF (PRESENT(stochastic_values)) noise=stochastic_values(j,k)

          ! Validate constitutive inputs before dispatch. Model/method errors
          ! remain configuration errors, not recoverable timestep rejections.
          IF (.NOT.ALL(ieee_is_finite(qpentry)) .OR. &
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
             IF (.NOT.ALL(ieee_is_finite(erodible(1:n_solid,j,k))) .OR. &
                  ANY(erodible(1:n_solid,j,k)<0.0_wp)) THEN
                CALL reject_cell('D19: invalid erodible inventory before mass exchange')
                EXIT cell_loop
             END IF
          END IF

          SELECT CASE(groups(group_index)%method)
          CASE(SOURCE_METHOD_EXPLICIT)
             physical_rate=0.0_wp
             exchange_rate=0.0_wp
             IF (has_other) CALL eval_cell_source_explicit(j,k,qpentry, &
                  source_time,other_sources,physical_rate)
             IF (with_exchange) THEN
                CALL eval_source_explicit(qpentry,B_zone(j,k), &
                     B_prime_x_geom(j,k),B_prime_y_geom(j,k), &
                     erodible(1:n_solid,j,k),time_step,transfers,exchange_rate)
                IF (bottom_radial_source_flag .OR. bottom_fissural_source_flag) THEN
                   coverage=1.0_wp-cell_source_fractions(j,k)
                   IF (.NOT.ieee_is_finite(coverage) .OR. &
                        coverage<0.0_wp .OR. coverage>1.0_wp) &
                        CALL reject_cell('D19: invalid mass-exchange source coverage')
                        IF (step_failed) EXIT cell_loop
                   exchange_rate=coverage*exchange_rate
                   CALL scale_mass_exchange_transfer(transfers,coverage)
                END IF
             END IF
             total_rate=physical_rate+exchange_rate
             IF (.NOT.ALL(ieee_is_finite(total_rate))) &
                  CALL reject_cell('D19: non-finite explicit source rate')
                  IF (step_failed) EXIT cell_loop
             IF (with_exchange) THEN
                IF (.NOT.ALL(ieee_is_finite(transfers%solid_deposition_rate)) .OR. &
                     .NOT.ALL(ieee_is_finite(transfers%solid_erosion_rate)) .OR. &
                     .NOT.ALL(ieee_is_finite(transfers%solid_mass_rate)) .OR. &
                     .NOT.ALL(ieee_is_finite(transfers%add_gas_mass_rate)) .OR. &
                     .NOT.ieee_is_finite(transfers%mixture_mass_rate) .OR. &
                     .NOT.ieee_is_finite(transfers%energy_rate) .OR. &
                     .NOT.ieee_is_finite(transfers%carrier_erosion_rate) .OR. &
                     .NOT.ieee_is_finite(transfers%carrier_deposition_rate) .OR. &
                     .NOT.ieee_is_finite(transfers%carrier_withdrawal_volume_rate) .OR. &
                     .NOT.ieee_is_finite(transfers%lava_deposit_rate) .OR. &
                     .NOT.ieee_is_finite(transfers%dB_dt)) &
                     CALL reject_cell('D19: non-finite coupled mass-exchange transfer')
                     IF (step_failed) EXIT cell_loop
             END IF
             ! One common limiter scales the coupled fluid + substrate transaction.
             ! It never touches the previously committed spatial FV rate.
             factor=1.0_wp
             IF (total_rate(idx_totMassEqn)<0.0_wp) factor=MIN(factor, &
                  MAX(0.0_wp,q(idx_totMassEqn,j,k)/ &
                  (-time_step*total_rate(idx_totMassEqn))))
             DO solid_eq=idx_solidEqn_first,idx_solidEqn_last
                IF (total_rate(solid_eq)<0.0_wp) factor=MIN(factor, &
                     MAX(0.0_wp,q(solid_eq,j,k)/(-time_step*total_rate(solid_eq))))
             END DO
             IF (with_exchange) THEN
                DO i=1,n_solid
                   IF (transfers%solid_erosion_rate(i)>0.0_wp) &
                        factor=MIN(factor,MAX(0.0_wp,erodible(i,j,k)/ &
                        (time_step*transfers%solid_erosion_rate(i))))
                END DO
             END IF
             remaining_mass=q(idx_totMassEqn,j,k)- &
                  SUM(q(idx_solidEqn_first:idx_solidEqn_last,j,k))
             change_in_solid=SUM(total_rate(idx_solidEqn_first:idx_solidEqn_last)) &
                  -total_rate(idx_totMassEqn)
             IF (change_in_solid>0.0_wp) factor=MIN(factor, &
                  MAX(0.0_wp,remaining_mass/(time_step*change_in_solid)))
             IF (.NOT.ieee_is_finite(factor) .OR. factor<0.0_wp .OR. factor>1.0_wp) &
                  CALL reject_cell('D19: invalid Lie source limiter factor')
                  IF (step_failed) EXIT cell_loop
             ! Stage the entire local fluid/substrate transfer before committing it.
             q_candidate=q(:,j,k)
             q_candidate(1:n_eqns)=q_candidate(1:n_eqns)+time_step*factor*total_rate
             IF (.NOT.ALL(ieee_is_finite(q_candidate(1:n_eqns)))) &
                  CALL reject_cell('D19: non-finite explicit group candidate')
                  IF (step_failed) EXIT cell_loop
             IF (.NOT.conservative_mass_admissible(q_candidate(1:n_eqns))) THEN
                CALL reject_cell('D19: inadmissible conservative mass candidate')
                EXIT cell_loop
             END IF
             IF (with_exchange) THEN
                deposit_candidate=deposit(j,k,1:n_solid)+ &
                     time_step*factor*transfers%solid_deposition_rate
                erosion_candidate=erosion(j,k,1:n_solid)+ &
                     time_step*factor*transfers%solid_erosion_rate
                erodible_candidate=erodible(1:n_solid,j,k)- &
                     time_step*factor*transfers%solid_erosion_rate
                IF (erodible_deposit_flag) erodible_candidate= &
                     erodible_candidate+time_step*factor*transfers%solid_deposition_rate
                lava_candidate=lava_deposit(j,k)+time_step*factor*transfers%lava_deposit_rate
                IF (.NOT.ALL(ieee_is_finite(deposit_candidate)) .OR. &
                     .NOT.ALL(ieee_is_finite(erosion_candidate)) .OR. &
                     .NOT.ALL(ieee_is_finite(erodible_candidate)) .OR. &
                     .NOT.ieee_is_finite(lava_candidate)) &
                     CALL reject_cell('D19: non-finite mass-exchange candidate inventory')
                     IF (step_failed) EXIT cell_loop
                IF (ANY(erodible_candidate<0.0_wp)) &
                     CALL reject_cell('D19: negative candidate erodible inventory')
                     IF (step_failed) EXIT cell_loop
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
             q(:,j,k)=q_candidate
             qp(:,j,k)=qp_candidate
             IF (with_exchange) THEN
                deposit(j,k,1:n_solid)=deposit_candidate
                erosion(j,k,1:n_solid)=erosion_candidate
                erodible(1:n_solid,j,k)=erodible_candidate
                lava_deposit(j,k)=lava_candidate
                IF (topo_change_flag) CALL record_topography_rate_cell(j,k, &
                     factor*transfers%dB_dt)
             END IF
          CASE(SOURCE_METHOD_ANALYTIC)
             IF (q(idx_totMassEqn,j,k)<=0.0_wp) CYCLE
             analytic_context=cell_analytic_source_context(j,k,noise)
             CALL eval_source_analytic(groups(group_index),q(:,j,k), &
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
             IF (q(idx_totMassEqn,j,k)<=0.0_wp) CYCLE
             CALL solve_backward_euler_group(groups(group_index), &
                  q(:,j,k),time_step,B_prime_x_geom(j,k),B_prime_y_geom(j,k), &
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
          IF (.NOT.conservative_mass_admissible(q(1:n_eqns,j,k))) THEN
             CALL reject_cell('D19: inadmissible conservative source group state')
             EXIT cell_loop
          END IF
          IF (q(idx_totMassEqn,j,k)>0.0_wp) THEN
             CALL qc_to_qp(q(:,j,k),qp(:,j,k),pressure,TRIM(context),conversion_ok)
             IF (.NOT.conversion_ok) THEN
                CALL reject_cell('D19: source group conversion failed')
                EXIT cell_loop
             END IF
          ELSE
             CALL initialize_dry_physical_cell(qp(:,j,k))
          END IF
       END DO cell_loop
       IF (step_failed) THEN
          IF (PRESENT(failure_message)) failure_message=TRIM(failure_reason)
          RETURN
       END IF
       IF (with_exchange .AND. topo_change_flag) THEN
          CALL apply_topography_rate_update(time_step,bed_accepted)
          IF (.NOT.bed_accepted) THEN
             CALL reject_cell('D19: projected topography update failed')
             IF (PRESENT(failure_message)) failure_message=TRIM(failure_reason)
             RETURN
          END IF
       END IF
    END DO
  CONTAINS
    ! Return a recoverable cell/group rejection without leaving the executor
    ! before returning the source-group status.
    SUBROUTINE reject_cell(message)
      CHARACTER(*), INTENT(IN) :: message
      failure_reason=message
      step_failed=.TRUE.
    END SUBROUTINE reject_cell
  END SUBROUTINE execute_lie_source_groups

END MODULE lie_source_executor_2d
