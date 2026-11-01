!********************************************************************************
!> \brief Execute validated STEP source groups after the RK update.
!> \details Coordinates the atomic mass-exchange transaction after IMEX
!>          integration. Physical rates are evaluated through source_terms_2d,
!>          while bed projection and geometry refresh remain in geometry_2d.
!********************************************************************************

MODULE step_source_executor_2d
    USE parameters_2d, ONLY : idx_totMassEqn, idx_solidEqn_first, idx_solidEqn_last

  USE diagnostics_2d, ONLY : fatal_error
  USE parameters_2d, ONLY : wp, n_eqns, n_vars, n_solid, mass_exchange_rates_type
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE domain_2d, ONLY : domain_type
  USE source_registry_2d, ONLY : source_registry_type, source_group_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: execute_lie_source_groups
  PUBLIC :: apply_mass_exchange_increment

CONTAINS

  !> \brief Apply and validate one atomic cell mass-exchange update.
  !> \details Applies conservative equation sources and material-inventory
  !>          increments for the completed time step, then enforces the
  !>          mass-exchange positivity checks and refreshes the primitive state.
  !> \param[in,out] q_cell Conservative cell state.
  !> \param[in,out] qp_cell Primitive cell state, refreshed from q_cell.
  !> \param[in] time_step Physical time increment [s].
  !> \param[in] equation_rate Conservative equation source rates.
  !> \param[in] deposition_rate Per-class deposition rates.
  !> \param[in] erosion_rate Per-class erosion rates.
  !> \param[in,out] deposit_inventory Accumulated deposited solid thickness.
  !> \param[in,out] erosion_inventory Accumulated eroded solid thickness.
  !> \param[in,out] erodible_inventory Remaining per-class substrate thickness.
  !> \param[in] redeposit_on_substrate Whether deposited solids replenish the substrate.
  !> \param[in] closure_context Cell context used in physical-state diagnostics.
  SUBROUTINE apply_mass_exchange_increment(q_cell, qp_cell, time_step,       &
       equation_rate, deposition_rate, erosion_rate, deposit_inventory,      &
       erosion_inventory, erodible_inventory, redeposit_on_substrate,        &
       closure_context)

    USE state_conversion_2d, ONLY : qc_to_qp, mixt_var
    USE state_2d, ONLY : initialize_dry_physical_cell

    REAL(wp), INTENT(INOUT) :: q_cell(n_vars), qp_cell(n_vars+2)
    REAL(wp), INTENT(IN) :: time_step, equation_rate(n_eqns)
    REAL(wp), INTENT(IN) :: deposition_rate(n_solid), erosion_rate(n_solid)
    REAL(wp), INTENT(INOUT) :: deposit_inventory(:)
    REAL(wp), INTENT(INOUT) :: erosion_inventory(:)
    REAL(wp), INTENT(INOUT) :: erodible_inventory(:)
    LOGICAL, INTENT(IN) :: redeposit_on_substrate
    CHARACTER(*), INTENT(IN) :: closure_context

    REAL(wp) :: p_dyn
    REAL(wp) :: r_Ri, r_rho_m, r_rho_c, r_red_grav
    REAL(wp) :: r_sp_heat_c, r_sp_heat_mix

    q_cell(1:n_eqns) = q_cell(1:n_eqns) + time_step * equation_rate
    q_cell(idx_solidEqn_first:idx_solidEqn_last) = MAX(0.0_wp, q_cell(idx_solidEqn_first:idx_solidEqn_last))

    deposit_inventory = deposit_inventory + time_step * deposition_rate
    erosion_inventory = erosion_inventory + time_step * erosion_rate
    erodible_inventory = erodible_inventory - time_step * erosion_rate

    IF (redeposit_on_substrate) THEN
       erodible_inventory = erodible_inventory + time_step * deposition_rate
    END IF

    IF (ANY(q_cell(idx_solidEqn_first:idx_solidEqn_last) .LT. 0.0_wp)) THEN
       CALL fatal_error('negative solid mass after mass exchange')
    END IF

    IF (q_cell(idx_totMassEqn) .LE. 0.0_wp) THEN
       IF (q_cell(idx_totMassEqn) .GT. -1.0E-10_wp) THEN
          q_cell(1:n_vars) = 0.0_wp
       ELSE
          CALL fatal_error('negative total mass after mass exchange')
       END IF
    END IF

    IF (SUM(q_cell(idx_solidEqn_first:idx_solidEqn_last)) .GT. q_cell(idx_totMassEqn)) THEN
       IF (q_cell(idx_totMassEqn) .LT. 1.0e-10_wp) THEN
          q_cell(idx_solidEqn_first:idx_solidEqn_last) = q_cell(idx_solidEqn_first:idx_solidEqn_last)                   &
               / SUM(q_cell(idx_solidEqn_first:idx_solidEqn_last)) * q_cell(idx_totMassEqn)
       ELSE
          CALL fatal_error('solid mass exceeds total mass after mass exchange')
       END IF
    END IF

    IF (q_cell(idx_totMassEqn) .GT. 0.0_wp) THEN
       CALL qc_to_qp(q_cell, qp_cell, p_dyn, TRIM(closure_context) // ' final state')
       CALL mixt_var(qp_cell, r_Ri, r_rho_m, r_rho_c, r_red_grav,       &
            r_sp_heat_c, r_sp_heat_mix)
    ELSE
       CALL initialize_dry_physical_cell(qp_cell)
       r_red_grav = 0.0_wp
    END IF

    IF (r_red_grav .LE. 0.0_wp) q_cell(1:n_vars) = 0.0_wp

  END SUBROUTINE apply_mass_exchange_increment

  !> Evaluate and commit one complete Lie source group at a time.
  !> For GROUP 0, evaluate every physical rate at the original pre-spatial
  !> state, but cap and apply the combined source increment to the post-FV state.
  SUBROUTINE execute_lie_source_groups(q,qp,time_step,domain,source_registry, &
       stochastic_values,step_start_time,q_original,qp_original)
    USE parameters_2d, ONLY : erodible_deposit_flag
    USE geometry_2d, ONLY : B_prime_x_geom,B_prime_y_geom,cell_source_fractions, &
         begin_topography_rate_update,record_topography_rate_cell, &
         apply_topography_rate_update,deposit,erosion,erodible,B_zone,lava_deposit
    USE source_terms_2d, ONLY : eval_source_explicit,eval_source_analytic, &
         eval_cell_source_explicit,cell_analytic_source_context, &
         analytic_term_context_type
    USE source_catalog_2d, ONLY : SOURCE_METHOD_ANALYTIC,SOURCE_METHOD_EXPLICIT, &
         SOURCE_METHOD_NEWTON
    USE equation_metadata_2d, ONLY : equation_partition_type,IM_GROUP_METHOD_NEWTON
    USE nonlinear_solver_2d, ONLY : solve_backward_euler_group
    USE state_conversion_2d, ONLY : qc_to_qp
    USE state_2d, ONLY : initialize_dry_physical_cell
    USE geometry_2d, ONLY : comp_cells_x,comp_cells_y
    USE parameters_2d, ONLY : topo_change_flag,bottom_radial_source_flag, &
         bottom_fissural_source_flag
    TYPE(source_group_type), ALLOCATABLE :: groups(:)
    TYPE(source_group_type) :: other_sources
    TYPE(mass_exchange_rates_type) :: transfers
    TYPE(analytic_term_context_type) :: analytic_context
    TYPE(equation_partition_type) :: local_partition
    REAL(wp), INTENT(INOUT) :: q(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(INOUT) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: time_step
    CLASS(domain_type), INTENT(IN) :: domain
    CLASS(source_registry_type), INTENT(IN) :: source_registry
    REAL(wp), INTENT(IN), OPTIONAL :: stochastic_values(:,:)
    REAL(wp), INTENT(IN), OPTIONAL :: step_start_time
    REAL(wp), INTENT(IN) :: q_original(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: qp_original(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp) :: qpentry(n_vars+2)
    REAL(wp) :: physical_rate(n_eqns),exchange_rate(n_eqns),total_rate(n_eqns)
    REAL(wp) :: analytic_rate(n_eqns),factor,coverage,pressure,source_time,noise
    REAL(wp) :: remaining_mass,change_in_solid
    LOGICAL :: with_exchange,has_other,implicit_mask(n_eqns),ok,failed
    INTEGER :: group_index,cell,j,k,i,solid_eq,status,iterations,linear_info
    INTEGER :: membership(n_eqns)
    CHARACTER(LEN=96) :: context

    IF (.NOT.source_registry%is_initialized()) &
         CALL fatal_error('D19: uninitialized Lie source registry')
    groups=source_registry%groups()
    source_time=0.0_wp
    IF (PRESENT(step_start_time)) source_time=step_start_time
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

       IF (groups(group_index)%method==SOURCE_METHOD_NEWTON) THEN
          implicit_mask=.FALSE.
          implicit_mask(groups(group_index)%equation_ids)=.TRUE.
          membership=0
          WHERE (implicit_mask) membership=1
          CALL local_partition%initialize(implicit_mask,membership,[IM_GROUP_METHOD_NEWTON])
       END IF
       IF (with_exchange .AND. topo_change_flag) CALL begin_topography_rate_update

       DO cell=1,domain%solve_cells
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
                   exchange_rate=coverage*exchange_rate
                   transfers%solid_deposition_rate=coverage*transfers%solid_deposition_rate
                   transfers%solid_erosion_rate=coverage*transfers%solid_erosion_rate
                   transfers%carrier_deposition_rate=coverage*transfers%carrier_deposition_rate
                   transfers%dB_dt=coverage*transfers%dB_dt
                   transfers%lava_deposit_rate=coverage*transfers%lava_deposit_rate
                END IF
             END IF
             total_rate=physical_rate+exchange_rate
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
             q(1:n_eqns,j,k)=q(1:n_eqns,j,k)+time_step*factor*total_rate
             IF (with_exchange) THEN
                deposit(j,k,1:n_solid)=deposit(j,k,1:n_solid)+ &
                     time_step*factor*transfers%solid_deposition_rate
                erosion(j,k,1:n_solid)=erosion(j,k,1:n_solid)+ &
                     time_step*factor*transfers%solid_erosion_rate
                erodible(1:n_solid,j,k)=erodible(1:n_solid,j,k)- &
                     time_step*factor*transfers%solid_erosion_rate
                IF (erodible_deposit_flag) erodible(1:n_solid,j,k)= &
                     erodible(1:n_solid,j,k)+time_step*factor*transfers%solid_deposition_rate
                lava_deposit(j,k)=lava_deposit(j,k)+time_step*factor*transfers%lava_deposit_rate
                IF (topo_change_flag) CALL record_topography_rate_cell(j,k, &
                     factor*transfers%dB_dt)
             END IF
          CASE(SOURCE_METHOD_ANALYTIC)
             IF (q(idx_totMassEqn,j,k)<=0.0_wp) CYCLE
             analytic_context=cell_analytic_source_context(j,k,noise)
             CALL eval_source_analytic(groups(group_index),q(:,j,k), &
                  analytic_context,time_step,analytic_rate,status)
             IF (status/=0) CALL fatal_error('D19: analytic group failed')
             q(1:n_vars,j,k)=q(1:n_vars,j,k)+time_step*analytic_rate(1:n_vars)
          CASE(SOURCE_METHOD_NEWTON)
             IF (q(idx_totMassEqn,j,k)<=0.0_wp) CYCLE
             CALL solve_backward_euler_group(local_partition,groups(group_index), &
                  q(:,j,k),time_step,B_prime_x_geom(j,k),B_prime_y_geom(j,k), &
                  noise,iterations,ok,linear_info,failed)
             IF (.NOT.ok .OR. failed .OR. linear_info/=0) &
                  CALL fatal_error('D19: Newton group failed')
          CASE DEFAULT
             CALL fatal_error('D19: invalid group method')
          END SELECT
          ! This is an admissibility check only, not an implicit post-solve clip.
          IF (.NOT.ALL(q(1:n_eqns,j,k)==q(1:n_eqns,j,k))) &
               CALL fatal_error('D19: source group produced NaN')
          IF (q(idx_totMassEqn,j,k)<0.0_wp) &
               CALL fatal_error('D19: negative mixture mass')
          IF (ANY(q(idx_solidEqn_first:idx_solidEqn_last,j,k)<0.0_wp)) &
               CALL fatal_error('D19: negative solid mass')
          IF (q(idx_totMassEqn,j,k)>0.0_wp) THEN
             CALL qc_to_qp(q(:,j,k),qp(:,j,k),pressure,TRIM(context))
          ELSE
             CALL initialize_dry_physical_cell(qp(:,j,k))
          END IF
       END DO
       IF (with_exchange .AND. topo_change_flag) CALL apply_topography_rate_update(time_step)
       IF (groups(group_index)%method==SOURCE_METHOD_NEWTON) CALL local_partition%finalize()
    END DO
  END SUBROUTINE execute_lie_source_groups

END MODULE step_source_executor_2d
