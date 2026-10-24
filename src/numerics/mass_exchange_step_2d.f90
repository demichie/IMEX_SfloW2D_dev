!********************************************************************************
!> \brief Execute the registered post-RK mass-exchange step.
!> \details Coordinates the atomic mass-exchange transaction after IMEX
!>          integration. Physical rates are evaluated by equation_terms_2d,
!>          while bed projection and geometry refresh remain in geometry_2d.
!********************************************************************************

MODULE mass_exchange_step_2d

  USE diagnostics_2d, ONLY : fatal_error
  USE parameters_2d, ONLY : wp, n_eqns, n_vars, n_solid
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE domain_2d, ONLY : domain_type
  USE source_registry_2d, ONLY : source_registry_type, source_group_type
  USE source_catalog_2d, ONLY : SOURCE_CADENCE_STEP, &
       source_term_type, source_term_descriptor_by_name

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: apply_mass_exchange_step
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
    USE constitutive_parameters_2d, ONLY : T_ambient

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
       qp_cell = 0.0_wp
       qp_cell(4) = T_ambient
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
  SUBROUTINE apply_mass_exchange_step(q, qp, time_step, domain, source_registry)

    USE constitutive_parameters_2d, ONLY : erosion_coeff, settling_flag,    &
         entrainment_flag, loss_rate
    USE geometry_2d, ONLY : B_prime_x_geom, B_prime_y_geom, cell_source_fractions
    USE geometry_2d, ONLY : begin_topography_rate_update,                   &
         record_topography_rate_cell, apply_topography_rate_update
    USE geometry_2d, ONLY : deposit, erosion, erodible, B_zone
    USE equation_terms_2d, ONLY : eval_mass_exchange_terms
    USE state_conversion_2d, ONLY : qc_to_qp
    USE parameters_2d, ONLY : topo_change_flag, bottom_radial_source_flag,  &
         bottom_fissural_source_flag, erodible_deposit_flag,                &
         pore_pressure_flag, liquid_flag
    USE constitutive_parameters_2d, ONLY : T_ambient

    REAL(wp), INTENT(INOUT) :: q(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(INOUT) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: time_step
    CLASS(domain_type), INTENT(IN) :: domain
    CLASS(source_registry_type), INTENT(IN), OPTIONAL :: source_registry

    REAL(wp) :: erosion_term(n_solid)
    REAL(wp) :: deposition_term(n_solid)
    REAL(wp) :: continuous_phase_erosion_term
    REAL(wp) :: continuous_phase_loss_term
    REAL(wp) :: eqns_term(n_eqns)
    REAL(wp) :: topo_term
    REAL(wp) :: out_of_source_fraction
    REAL(wp) :: p_dyn
    CHARACTER(LEN=96) :: closure_context
    TYPE(source_group_type), ALLOCATABLE :: step_groups(:)
    TYPE(source_term_type) :: term
    INTEGER :: j, k, l
    INTEGER :: group_index, term_index, mass_group_count
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
          IF (step_groups(group_index)%cadence/=SOURCE_CADENCE_STEP) &
               CALL fatal_error('D17: invalid cadence in post-RK source schedule')
          IF (SIZE(step_groups(group_index)%source_names)>0) &
               mass_group_count=mass_group_count+1
          DO term_index=1,SIZE(step_groups(group_index)%source_names)
             term=source_term_descriptor_by_name(step_groups(group_index)%source_names(term_index))
             IF (term%target/='mass_exchange') &
                  CALL fatal_error('D17: unsupported non-mass post-RK source')
          END DO
       END DO
       IF (mass_group_count>1) &
            CALL fatal_error('D17: atomic mass-exchange batch split across step groups')
       IF (mass_group_count==0) RETURN
    END IF

    IF ((erosion_coeff .EQ. 0.0_wp) .AND. (.NOT.settling_flag)            &
         .AND. (.NOT.pore_pressure_flag) .AND. (.NOT.entrainment_flag)     &
         .AND. (.NOT.liquid_loss_active)) RETURN

    IF (topo_change_flag) CALL begin_topography_rate_update

    !$OMP PARALLEL DO private(j,k,closure_context,erosion_term,deposition_term, &
    !$OMP & eqns_term,                                                          &
    !$OMP & topo_term,                                                          &
    !$OMP & continuous_phase_erosion_term,continuous_phase_loss_term,           &
    !$OMP & out_of_source_fraction,p_dyn)
    DO l = 1, domain%solve_cells

       j = domain%j_cent(l)
       k = domain%k_cent(l)
       WRITE(closure_context,'(A,I0,A,I0)') 'mass exchange cell=',j,',',k

       IF (q(1,j,k) .GT. 0.0_wp) THEN
          CALL qc_to_qp(q(1:n_vars,j,k), qp(1:n_vars+2,j,k), p_dyn,         &
               TRIM(closure_context))
       ELSE
          qp(1:n_vars+2,j,k) = 0.0_wp
          qp(4,j,k) = T_ambient
       END IF

       CALL eval_mass_exchange_terms(qp(1:n_vars+2,j,k), B_zone(j,k),      &
            B_prime_x_geom(j,k), B_prime_y_geom(j,k),                      &
            erodible(1:n_solid,j,k), time_step, erosion_term,              &
            deposition_term, continuous_phase_erosion_term,                &
            continuous_phase_loss_term, eqns_term, topo_term)

       IF (bottom_radial_source_flag .OR. bottom_fissural_source_flag) THEN
          out_of_source_fraction = 1.0_wp - cell_source_fractions(j,k)
          deposition_term = deposition_term * out_of_source_fraction
          erosion_term = erosion_term * out_of_source_fraction
          eqns_term = eqns_term * out_of_source_fraction
          topo_term = topo_term * out_of_source_fraction
       END IF

       CALL apply_mass_exchange_increment(q(1:n_vars,j,k),                 &
            qp(1:n_vars+2,j,k), time_step, eqns_term, deposition_term,     &
            erosion_term, deposit(j,k,1:n_solid), erosion(j,k,1:n_solid),  &
            erodible(1:n_solid,j,k), erodible_deposit_flag, closure_context)

       IF (topo_change_flag) CALL record_topography_rate_cell(j, k, topo_term)

    END DO
    !$OMP END PARALLEL DO

    IF (topo_change_flag) CALL apply_topography_rate_update(time_step)

  END SUBROUTINE apply_mass_exchange_step

END MODULE mass_exchange_step_2d
