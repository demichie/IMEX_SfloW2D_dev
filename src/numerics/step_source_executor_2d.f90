!********************************************************************************
!> \brief Execute validated STEP source groups after the RK update.
!> \details Coordinates the atomic mass-exchange transaction after IMEX
!>          integration. Physical rates are evaluated through source_terms_2d,
!>          while bed projection and geometry refresh remain in geometry_2d.
!********************************************************************************

MODULE step_source_executor_2d

  USE diagnostics_2d, ONLY : fatal_error
  USE parameters_2d, ONLY : wp, n_eqns, n_vars, n_solid, mass_exchange_rates_type
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE domain_2d, ONLY : domain_type
  USE source_registry_2d, ONLY : source_registry_type, source_group_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: execute_step_source_groups
  PUBLIC :: apply_mass_exchange_increment

CONTAINS

  !> \brief Apply and validate a cell's post-IMEX mass-exchange update.
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
    q_cell(5:4+n_solid) = MAX(0.0_wp, q_cell(5:4+n_solid))

    deposit_inventory = deposit_inventory + time_step * deposition_rate
    erosion_inventory = erosion_inventory + time_step * erosion_rate
    erodible_inventory = erodible_inventory - time_step * erosion_rate

    IF (redeposit_on_substrate) THEN
       erodible_inventory = erodible_inventory + time_step * deposition_rate
    END IF

    IF (ANY(q_cell(5:4+n_solid) .LT. 0.0_wp)) THEN
       CALL fatal_error('negative solid mass after mass exchange')
    END IF

    IF (q_cell(1) .LE. 0.0_wp) THEN
       IF (q_cell(1) .GT. -1.0E-10_wp) THEN
          q_cell(1:n_vars) = 0.0_wp
       ELSE
          CALL fatal_error('negative total mass after mass exchange')
       END IF
    END IF

    IF (SUM(q_cell(5:4+n_solid)) .GT. q_cell(1)) THEN
       IF (q_cell(1) .LT. 1.0e-10_wp) THEN
          q_cell(5:4+n_solid) = q_cell(5:4+n_solid)                   &
               / SUM(q_cell(5:4+n_solid)) * q_cell(1)
       ELSE
          CALL fatal_error('solid mass exceeds total mass after mass exchange')
       END IF
    END IF

    IF (q_cell(1) .GT. 0.0_wp) THEN
       CALL qc_to_qp(q_cell, qp_cell, p_dyn, TRIM(closure_context) // ' final state')
       CALL mixt_var(qp_cell, r_Ri, r_rho_m, r_rho_c, r_red_grav,       &
            r_sp_heat_c, r_sp_heat_mix)
    ELSE
       CALL initialize_dry_physical_cell(qp_cell)
       r_red_grav = 0.0_wp
    END IF

    IF (r_red_grav .LE. 0.0_wp) q_cell(1:n_vars) = 0.0_wp

  END SUBROUTINE apply_mass_exchange_increment

  !> \brief Evaluate and apply the post-IMEX mass-exchange fractional step.
  !> \details Evaluates cell-local exchange sources, applies source-coverage
  !>          masks, updates the conservative state and material inventories,
  !>          then delegates any nodal bed update to geometry_2d. Call once
  !>          after the IMEX update and before advancing the simulation clock.
  !> \param[in,out] q Cell-centered conservative state.
  !> \param[in,out] qp Cell-centered primitive state, refreshed for solved cells.
  !> \param[in] time_step Physical time increment [s].
  !> \param[in] domain Domain descriptor for active solve-cell indices.
  !> \param[in] source_registry Optional validated schedule for post-RK sources.
  SUBROUTINE execute_step_source_groups(q, qp, time_step, domain, source_registry, stochastic_values, step_start_time)

    USE constitutive_parameters_2d, ONLY : erosion_coeff, settling_flag,    &
         entrainment_flag, loss_rate
    USE geometry_2d, ONLY : B_prime_x_geom, B_prime_y_geom, cell_source_fractions
    USE geometry_2d, ONLY : B_second_xx_geom, B_second_xy_geom, B_second_yy_geom, &
         grav_coeff, cell_fissure_fractions, cell_arc_perim, cell_arc_n_x, &
         cell_arc_n_y, dx, dy
    USE geometry_2d, ONLY : begin_topography_rate_update,                   &
         record_topography_rate_cell, apply_topography_rate_update
    USE geometry_2d, ONLY : deposit, erosion, erodible, B_zone, lava_deposit
    USE source_terms_2d, ONLY : eval_source_explicit, eval_source_analytic, &
         eval_cell_source_explicit, cell_analytic_source_context, &
         analytic_term_context_type
    USE source_catalog_2d, ONLY : SOURCE_METHOD_ANALYTIC, SOURCE_METHOD_EXPLICIT, &
         SOURCE_METHOD_NEWTON
    USE equation_metadata_2d, ONLY : equation_partition_type, IM_GROUP_METHOD_NEWTON
    USE nonlinear_solver_2d, ONLY : solve_rk_step
    USE state_conversion_2d, ONLY : qc_to_qp
    USE parameters_2d, ONLY : topo_change_flag, bottom_radial_source_flag,  &
         bottom_fissural_source_flag, erodible_deposit_flag,                &
         pore_pressure_flag, liquid_flag, lava_flag
    USE state_2d, ONLY : initialize_dry_physical_cell
    USE mass_exchange_lava_deposition_2d, ONLY : lava_deposition_flag

    REAL(wp), INTENT(INOUT) :: q(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(INOUT) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: time_step
    CLASS(domain_type), INTENT(IN) :: domain
    CLASS(source_registry_type), INTENT(IN), OPTIONAL :: source_registry
    REAL(wp), INTENT(IN), OPTIONAL :: stochastic_values(:,:)
    REAL(wp), INTENT(IN), OPTIONAL :: step_start_time

    TYPE(mass_exchange_rates_type) :: exchange_rates
    REAL(wp) :: eqns_term(n_eqns)
    REAL(wp) :: out_of_source_fraction
    REAL(wp) :: p_dyn
    CHARACTER(LEN=96) :: closure_context
    TYPE(source_group_type), ALLOCATABLE :: step_groups(:)
    INTEGER :: j, k, l
    INTEGER :: group_index, mass_group_count, analytic_status, cell_index
    TYPE(analytic_term_context_type) :: analytic_context
    REAL(wp) :: analytic_state(n_vars), old_state(n_vars), step_noise, &
         explicit_rate(n_eqns), source_time
    REAL(wp) :: empty_residual(n_eqns)
    LOGICAL :: implicit_mask(n_eqns), newton_converged, newton_search_failed
    INTEGER :: newton_iterations, newton_info
    TYPE(equation_partition_type) :: step_partition
    LOGICAL :: liquid_loss_active

    liquid_loss_active = .FALSE.
    IF (liquid_flag) THEN
       IF (ALLOCATED(loss_rate)) liquid_loss_active = loss_rate .GT. 0.0_wp
    END IF

    IF (PRESENT(source_registry)) THEN
       IF (.NOT.source_registry%is_initialized()) &
            CALL fatal_error('D17: post-RK source dispatch requires an initialized registry')
       step_groups=source_registry%step_groups()
       mass_group_count=0
       DO group_index=1,SIZE(step_groups)
          IF (ANY(INDEX(step_groups(group_index)%source_names,'mass_exchange_')==1)) &
               mass_group_count=mass_group_count+1
       END DO
       IF (mass_group_count>1) &
            CALL fatal_error('D17: atomic mass-exchange batch split across step groups')
       ! The absence of mass exchange does not suppress other STEP groups.
       empty_residual=0.0_wp
       source_time=0.0_wp
       IF (PRESENT(step_start_time)) source_time=step_start_time
       DO group_index=1,SIZE(step_groups)
          ! Honor GROUP ordering; exchange remains one atomic transaction.
          IF (SIZE(step_groups(group_index)%source_names)==0) CYCLE
          IF (ANY(INDEX(step_groups(group_index)%source_names,'mass_exchange_')==1)) THEN
             CALL execute_mass_exchange_transaction()
             CYCLE
          END IF
          IF (step_groups(group_index)%method==SOURCE_METHOD_NEWTON) THEN
             implicit_mask=.FALSE.
             implicit_mask(step_groups(group_index)%equation_ids)=.TRUE.
             BLOCK
               INTEGER :: group_membership(n_eqns)
               group_membership=0
               WHERE(implicit_mask) group_membership=1
               CALL step_partition%initialize(implicit_mask,group_membership, &
                    [IM_GROUP_METHOD_NEWTON])
             END BLOCK
          END IF
          DO cell_index=1,domain%solve_cells
             j=domain%j_cent(cell_index); k=domain%k_cent(cell_index)
             IF (q(1,j,k)<=0.0_wp) CYCLE
             IF (step_groups(group_index)%method==SOURCE_METHOD_EXPLICIT) THEN
                CALL eval_cell_source_explicit(j,k,qp(1:n_vars+2,j,k), &
                     source_time,step_groups(group_index),explicit_rate)
                q(1:n_eqns,j,k)=q(1:n_eqns,j,k)+time_step*explicit_rate
             ELSE IF (step_groups(group_index)%method==SOURCE_METHOD_ANALYTIC) THEN
                step_noise=0.0_wp
                IF (PRESENT(stochastic_values)) step_noise=stochastic_values(j,k)
                analytic_context=cell_analytic_source_context(j,k,step_noise)
                CALL eval_source_analytic(step_groups(group_index),q(1:n_vars,j,k), &
                     analytic_context,time_step,analytic_state,analytic_status)
                IF (analytic_status/=0) &
                     CALL fatal_error('D18: analytical STEP group solve failed')
                q(1:n_vars,j,k)=analytic_state
             ELSE
                old_state=q(1:n_vars,j,k)
                step_noise=0.0_wp
                IF (PRESENT(stochastic_values)) step_noise=stochastic_values(j,k)
                CALL solve_rk_step(step_partition,[step_groups(group_index)], &
                     q(1:n_vars,j,k),old_state,time_step,1.0_wp,empty_residual, &
                     B_prime_x_geom(j,k),B_prime_y_geom(j,k),step_noise, &
                     newton_iterations,newton_converged,newton_info,newton_search_failed)
                IF (.NOT.newton_converged .OR. newton_info/=0 .OR. newton_search_failed) &
                     CALL fatal_error('D18: Newton STEP group solve failed')
             END IF
             CALL qc_to_qp(q(1:n_vars,j,k),qp(1:n_vars+2,j,k),p_dyn, &
                  'analytic STEP source')
          END DO
       END DO
       CALL step_partition%finalize()
       RETURN
    END IF

    ! Preserve direct legacy calls without a registry.
    CALL execute_mass_exchange_transaction()

  CONTAINS

    !> \\brief Execute one atomic bed-and-inventory transaction.
    !> \\details Called in ascending STEP GROUP order when configured.
    SUBROUTINE execute_mass_exchange_transaction()
    IF ((erosion_coeff .EQ. 0.0_wp) .AND. (.NOT.settling_flag)            &
         .AND. (.NOT.pore_pressure_flag) .AND. (.NOT.entrainment_flag)     &
         .AND. (.NOT.liquid_loss_active) .AND. (.NOT.lava_deposition_flag)) RETURN

    IF (lava_deposition_flag .AND. .NOT.topo_change_flag) &
         CALL fatal_error('lava deposition requires TOPO_CHANGE_FLAG')
    IF (topo_change_flag) CALL begin_topography_rate_update

    !$OMP PARALLEL DO private(j,k,closure_context,exchange_rates,eqns_term, &
    !$OMP & out_of_source_fraction,p_dyn)
    DO l = 1, domain%solve_cells

       j = domain%j_cent(l)
       k = domain%k_cent(l)
       closure_context = 'mass exchange cell'

       IF (q(1,j,k) .GT. 0.0_wp) THEN
          CALL qc_to_qp(q(1:n_vars,j,k), qp(1:n_vars+2,j,k), p_dyn,         &
               TRIM(closure_context))
       ELSE
          CALL initialize_dry_physical_cell(qp(1:n_vars+2,j,k))
       END IF

       CALL eval_source_explicit(qp(1:n_vars+2,j,k), B_zone(j,k), &
            B_prime_x_geom(j,k), B_prime_y_geom(j,k), &
            erodible(1:n_solid,j,k), time_step, exchange_rates, eqns_term)

       IF (bottom_radial_source_flag .OR. bottom_fissural_source_flag) THEN
          out_of_source_fraction = 1.0_wp - cell_source_fractions(j,k)
          exchange_rates%solid_deposition_rate = exchange_rates%solid_deposition_rate*out_of_source_fraction
          exchange_rates%solid_erosion_rate = exchange_rates%solid_erosion_rate*out_of_source_fraction
          exchange_rates%carrier_deposition_rate = &
               exchange_rates%carrier_deposition_rate*out_of_source_fraction
          eqns_term = eqns_term * out_of_source_fraction
          exchange_rates%dB_dt = exchange_rates%dB_dt*out_of_source_fraction
          exchange_rates%lava_deposit_rate = exchange_rates%lava_deposit_rate*out_of_source_fraction
       END IF

       CALL apply_mass_exchange_increment(q(1:n_vars,j,k),                 &
            qp(1:n_vars+2,j,k), time_step, eqns_term, exchange_rates%solid_deposition_rate, &
            exchange_rates%solid_erosion_rate, deposit(j,k,1:n_solid), erosion(j,k,1:n_solid),  &
            erodible(1:n_solid,j,k), erodible_deposit_flag, closure_context)

       lava_deposit(j,k)=lava_deposit(j,k)+time_step*exchange_rates%lava_deposit_rate
       IF (topo_change_flag) CALL record_topography_rate_cell(j, k, exchange_rates%dB_dt)

    END DO
    !$OMP END PARALLEL DO

    IF (topo_change_flag) CALL apply_topography_rate_update(time_step)

    END SUBROUTINE execute_mass_exchange_transaction

  END SUBROUTINE execute_step_source_groups

END MODULE step_source_executor_2d
