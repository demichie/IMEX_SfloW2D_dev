!********************************************************************************
!> \brief IMEX Runge-Kutta time integration
!>
!> This module owns the IMEX tableau and Runge-Kutta stage workspace,
!> orchestrates explicit/implicit stage advancement and applies post-IMEX
!> fractional-step mass-exchange increments to cell states and inventories.
!********************************************************************************

MODULE time_integration_2d

  USE diagnostics_2d, ONLY : debug_pause, fatal_error

  USE parameters_2d, ONLY : wp
  USE parameters_2d, ONLY : n_eqns, n_vars, n_RK, n_solid, n_add_gas
  USE parameters_2d, ONLY : verbose_level

  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE geometry_2d, ONLY : dx, dy
  USE geometry_2d, ONLY : B_cent
  USE geometry_2d, ONLY : B_prime_x_geom, B_prime_y_geom
  USE geometry_2d, ONLY : B_second_xx_geom, B_second_xy_geom, B_second_yy_geom
  USE geometry_2d, ONLY : grav_coeff
   USE geometry_2d, ONLY : cell_source_fractions, cell_fissure_fractions
  USE geometry_2d, ONLY : cell_arc_perim, cell_arc_n_x, cell_arc_n_y

  USE constitutive_parameters_2d, ONLY : T_ambient
  USE nonlinear_solver_2d, ONLY : solve_imex_stage

  USE equation_metadata_2d, ONLY : equation_partition_type

  USE source_registry_2d, ONLY : source_registry_type, source_group_type, &
       source_solver_type

  USE spatial_operator_2d, ONLY : spatial_operator_type

  USE domain_2d, ONLY : domain_type

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: validate_runge_kutta_stages
  PUBLIC :: apply_mass_exchange_increment
  PUBLIC :: apply_mass_exchange_step

  ABSTRACT INTERFACE
     !> \brief Observe an unmodified cell state at a specified IMEX audit point.
     !> \param[in] event 1: known stage state before closure; 2: solved stage; 3: raw final assembly.
     !> \param[in] stage Runge-Kutta stage index; zero for the final assembly.
     !> \param[in] j Cell x index in the active solve list.
     !> \param[in] k Cell y index in the active solve list.
     !> \param[in] q_cell Conservative cell state, before subsequent conversions or final clipping.
     !> \param[in] implicit_status 0: no local solve; 1: converged; -1: failed local solve.
     !> \note Called concurrently from OpenMP cell loops. Implementations must be
     !>       thread-safe and must not alter solver configuration or input states.
     SUBROUTINE imex_state_observer(event, stage, j, k, q_cell, implicit_status)
       IMPORT :: wp
       INTEGER, INTENT(IN) :: event, stage, j, k, implicit_status
       REAL(wp), INTENT(IN) :: q_cell(:)
     END SUBROUTINE imex_state_observer
  END INTERFACE

  !> \brief Persistent IMEX tableaux, current-stage states and stage source arrays.
  !> \details Only the current stage state is retained; previously evaluated
  !>          terms are stored along the Runge-Kutta stage index for assembly.
  TYPE, PUBLIC :: time_integration_workspace_type
     PRIVATE

     REAL(wp), ALLOCATABLE :: a_tilde_ij(:,:)
     REAL(wp), ALLOCATABLE :: a_dirk_ij(:,:)
     REAL(wp), ALLOCATABLE :: omega_tilde(:)
     REAL(wp), ALLOCATABLE :: omega(:)
     REAL(wp), ALLOCATABLE :: a_tilde(:)
     REAL(wp), ALLOCATABLE :: a_dirk(:)
     REAL(wp), ALLOCATABLE :: q_rk(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_rk(:,:,:)
     REAL(wp), ALLOCATABLE :: spatial_terms(:,:,:,:)
     REAL(wp), ALLOCATABLE :: IM_terms(:,:,:,:)
     REAL(wp), ALLOCATABLE :: EX_terms(:,:,:,:)
   CONTAINS
     PROCEDURE, PUBLIC :: initialize => initialize_time_integration
     PROCEDURE, PUBLIC :: finalize => finalize_time_integration
     PROCEDURE, PUBLIC :: advance => imex_RK_solver
  END TYPE time_integration_workspace_type

CONTAINS

  !> \brief Reject unsupported stage counts before any Runge-Kutta allocation.
  !> \param[in] stages Number of stages requested by input or a direct workspace caller.
  !> \note Only 2, 3 and 4 stages are supported; this is not the temporal order.
  SUBROUTINE validate_runge_kutta_stages(stages)
    INTEGER, INTENT(IN) :: stages

    IF (stages < 2 .OR. stages > 4) THEN
       WRITE(*,*) 'Invalid N_RK:', stages
       CALL fatal_error('N_RK must be 2, 3 or 4')
    END IF
  END SUBROUTINE validate_runge_kutta_stages

  !> \brief Apply and validate a cell's post-IMEX mass-exchange update.
  !> \details Applies conservative equation sources and material-inventory
  !>          increments for the completed time step, then enforces the
  !>          mass-exchange positivity checks and refreshes the primitive state.
  !>          The caller invokes this after the IMEX Runge-Kutta update and
  !>          before projecting any cell bed-rate proposals or advancing time.
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
  !******************************************************************************
  SUBROUTINE apply_mass_exchange_increment(q_cell, qp_cell, time_step,       &
       equation_rate, deposition_rate, erosion_rate, deposit_inventory,      &
       erosion_inventory, erodible_inventory, redeposit_on_substrate,        &
       closure_context)

    USE state_conversion_2d, ONLY : qc_to_qp, mixt_var

    IMPLICIT NONE

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
  !******************************************************************************
  SUBROUTINE apply_mass_exchange_step(q, qp, time_step, domain)

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

    IMPLICIT NONE

    REAL(wp), INTENT(INOUT) :: q(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(INOUT) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: time_step
    CLASS(domain_type), INTENT(IN) :: domain

    REAL(wp) :: erosion_term(n_solid)
    REAL(wp) :: deposition_term(n_solid)
    REAL(wp) :: continuous_phase_erosion_term
    REAL(wp) :: continuous_phase_loss_term
    REAL(wp) :: eqns_term(n_eqns)
    REAL(wp) :: topo_term
    REAL(wp) :: out_of_source_fraction
    REAL(wp) :: p_dyn
    CHARACTER(LEN=96) :: closure_context
    INTEGER :: j, k, l
    LOGICAL :: liquid_loss_active

    ! Carrier loss is a mass exchange even without erosion or settling. Use
    ! nested guards: Fortran does not require short-circuit evaluation, and
    ! gas-only configurations need not allocate the liquid loss parameter.
    liquid_loss_active = .FALSE.
    IF (liquid_flag) THEN
       IF (ALLOCATED(loss_rate)) liquid_loss_active = loss_rate .GT. 0.0_wp
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
          ! Apply the source-coverage mask before changing state, inventories,
          ! or the cell proposal passed to geometry_2d.
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

  !> \brief Build the selected IMEX tableau and allocate Runge-Kutta stage storage.
  !>
  !> \param[in,out] this Persistent IMEX tableau and Runge-Kutta stage workspace.
  !> \note Validates N_RK before allocation, including callers bypassing input parsing.

  SUBROUTINE initialize_time_integration( this )

    CLASS(time_integration_workspace_type), INTENT(INOUT) :: this
    REAL(wp) :: gamma, delta

    CALL validate_runge_kutta_stages(n_RK)

    ALLOCATE( this%a_tilde_ij(n_RK,n_RK) )
    ALLOCATE( this%a_dirk_ij(n_RK,n_RK) )
    ALLOCATE( this%omega_tilde(n_RK) )
    ALLOCATE( this%omega(n_RK) )

    ! Initialize the coefficients for the IMEX Runge-Kutta scheme
    ! Please note that with respect to the schemes described in Pareschi & Russo
    ! (2000) we do not have the coefficient vectors c_tilde and c, because the
    ! explicit and implicit terms do not depend explicitly on time.

    this%a_tilde_ij = 0.0_wp
    this%omega_tilde = 0.0_wp
    this%a_dirk_ij = 0.0_wp
    this%omega = 0.0_wp

    gamma = 1.0_wp - 1.0_wp / SQRT(2.0_wp)
    delta = 1.0_wp - 1.0_wp / ( 2.0_wp * gamma )

    IF ( n_RK .EQ. 2 ) THEN

       this%a_tilde_ij(2,1) = 1.0_wp
       this%omega_tilde(1) = 1.0_wp
       this%omega_tilde(2) = 0.0_wp
       this%a_dirk_ij(2,2) = 1.0_wp
       this%omega(1) = 0.0_wp
       this%omega(2) = 1.0_wp

    ELSEIF ( n_RK .EQ. 3 ) THEN

       ! Tableau for the IMEX-SSP(3,3,2) Stiffly Accurate Scheme
       ! from Pareschi & Russo (2005), Table IV
       this%a_tilde_ij(2,1) = 0.5_wp
       this%a_tilde_ij(3,1) = 0.5_wp
       this%a_tilde_ij(3,2) = 0.5_wp
       this%omega_tilde(1) = 1.0_wp / 3.0_wp
       this%omega_tilde(2) = 1.0_wp / 3.0_wp
       this%omega_tilde(3) = 1.0_wp / 3.0_wp
       this%a_dirk_ij(1,1) = 0.25_wp
       this%a_dirk_ij(2,2) = 0.25_wp
       this%a_dirk_ij(3,1) = 1.0_wp / 3.0_wp
       this%a_dirk_ij(3,2) = 1.0_wp / 3.0_wp
       this%a_dirk_ij(3,3) = 1.0_wp / 3.0_wp
       this%omega(1) = 1.0_wp / 3.0_wp
       this%omega(2) = 1.0_wp / 3.0_wp
       this%omega(3) = 1.0_wp / 3.0_wp


    ELSEIF ( n_RK .EQ. 4 ) THEN

       ! IMEX-SSP3(4,3,3), Pareschi & Russo (2005)
       ! Third-order, four-stage IMEX scheme.
       ! Explicit SSP(3,3) and L-stable implicit part.

       gamma = 0.24169426078821_wp
       delta = 0.06042356519705_wp
       ! eta is stored locally through the expression below.

       ! Explicit tableau
       this%a_tilde_ij(3,2) = 1.0_wp
       this%a_tilde_ij(4,2) = 0.25_wp
       this%a_tilde_ij(4,3) = 0.25_wp

       this%omega_tilde(1) = 0.0_wp
       this%omega_tilde(2) = 1.0_wp / 6.0_wp
       this%omega_tilde(3) = 1.0_wp / 6.0_wp
       this%omega_tilde(4) = 2.0_wp / 3.0_wp

       ! Implicit tableau
       this%a_dirk_ij(1,1) = gamma

       this%a_dirk_ij(2,1) = -gamma
       this%a_dirk_ij(2,2) = gamma

       this%a_dirk_ij(3,2) = 1.0_wp - gamma
       this%a_dirk_ij(3,3) = gamma

       this%a_dirk_ij(4,1) = delta
       this%a_dirk_ij(4,2) = 0.12915286960590_wp
       this%a_dirk_ij(4,3) = 0.5_wp - delta  &
            - 0.12915286960590_wp - gamma
       this%a_dirk_ij(4,4) = gamma

       this%omega(1) = 0.0_wp
       this%omega(2) = 1.0_wp / 6.0_wp
       this%omega(3) = 1.0_wp / 6.0_wp
       this%omega(4) = 2.0_wp / 3.0_wp


    END IF

    ALLOCATE( this%a_tilde(n_RK) )
    ALLOCATE( this%a_dirk(n_RK) )
    ALLOCATE( this%q_rk(n_vars,comp_cells_x,comp_cells_y) )
    ALLOCATE( this%qp_rk(n_vars+2,comp_cells_x,comp_cells_y) )
    ALLOCATE( this%spatial_terms(n_eqns,comp_cells_x,comp_cells_y,n_RK) )
    ALLOCATE( this%IM_terms(n_eqns,comp_cells_x,comp_cells_y,n_RK) )
    ALLOCATE( this%EX_terms(n_eqns,comp_cells_x,comp_cells_y,n_RK) )

  END SUBROUTINE initialize_time_integration

  !> \brief Release the IMEX tableau and stage work arrays.
  !>
  !> \param[in,out] this Persistent IMEX tableau and Runge-Kutta stage workspace.

  SUBROUTINE finalize_time_integration( this )

    CLASS(time_integration_workspace_type), INTENT(INOUT) :: this

    DEALLOCATE( this%a_tilde_ij )
    DEALLOCATE( this%a_dirk_ij )
    DEALLOCATE( this%omega_tilde )
    DEALLOCATE( this%omega )
    DEALLOCATE( this%a_tilde )
    DEALLOCATE( this%a_dirk )
    DEALLOCATE( this%q_rk )
    DEALLOCATE( this%qp_rk )
    DEALLOCATE( this%spatial_terms )
    DEALLOCATE( this%IM_terms )
    DEALLOCATE( this%EX_terms )

  END SUBROUTINE finalize_time_integration

  !******************************************************************************
  !> \brief Advance the conservative state through all explicit/implicit Runge-Kutta stages.
  !
  !> This subroutine integrates the explicit spatial and local source terms
  !> together with the implicit terms using an IMEX Runge-Kutta scheme.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !>
  !> \param[in,out] this Persistent IMEX tableau and Runge-Kutta stage workspace.
  !> \param[in,out] q Cell-centered conservative states, indexed as (variable,x-cell,y-cell).
  !> \param[in,out] qp Canonical physical states: h, hu, hv, T, component mass fractions, then
  !>                   appended u and v.
  !> \param[in] t Current simulation or stage time [s].
  !> \param[in] dt Time increment [s].
  !> \param[in] Z Coordinates of the three-point limiter stencil.
  !> \param[in] equation_partition Explicit/implicit equation mask and compact index maps.
  !> \param[in] source_registry Validated source metadata used to dispatch stage solves.
  !> \param[in] domain Active-cell/face lists and reconstruction halo for this simulation.
  !> \param[in,out] spatial_operator Spatial operator and its persistent reconstruction/KT
  !>                                 workspaces.
  !> \param[in] observer Optional thread-safe read-only stage audit callback. When absent,
  !>                    no diagnostic storage is allocated and no callbacks are made.
  !******************************************************************************

  SUBROUTINE imex_RK_solver(this, q, qp, t, dt, Z, equation_partition,      &
       source_registry, domain, spatial_operator, observer)

    USE state_conversion_2d, ONLY : qc_to_qp, project_dry_cell_momenta,      &
         clip_roundoff_negative_component_masses
    USE parameters_2d, ONLY : dry_thickness_tolerance, idx_u, idx_v

    USE equation_terms_2d, ONLY : eval_expl_terms

    USE constitutive_parameters_2d, ONLY : T_ambient

    USE geometry_2d, ONLY : B_nodata


    IMPLICIT NONE

    CLASS(time_integration_workspace_type), INTENT(INOUT) :: this
    REAL(wp), INTENT(INOUT) :: q(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(INOUT) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: t, dt
    REAL(wp), INTENT(IN) :: Z(comp_cells_x,comp_cells_y)
    TYPE(equation_partition_type), INTENT(IN) :: equation_partition
    TYPE(source_registry_type), INTENT(IN) :: source_registry
    CLASS(domain_type), INTENT(IN) :: domain
    CLASS(spatial_operator_type), INTENT(INOUT) :: spatial_operator
    PROCEDURE(imex_state_observer), OPTIONAL :: observer

    TYPE(source_group_type), ALLOCATABLE :: stage_groups(:)
    TYPE(source_solver_type), ALLOCATABLE :: source_solvers(:)
    REAL(wp) :: q_guess(n_vars) !< initial guess for the solution of the RK step
    REAL(wp) :: q_fv_cell(n_vars) !< finite-volume state for the current cell
    REAL(wp) :: residual_cell(n_vars) !< final RK residual for the current cell
    REAL(wp) :: q_old_cell(n_vars) !< state before the final RK assembly
    INTEGER :: i_RK
    INTEGER :: j,k,l            !< loop counter over the grid volumes
    REAL(wp) :: known_stage_residual(n_eqns)

    REAL(wp) :: a_diag
    REAL(wp) :: p_dyn
    CHARACTER(LEN=128) :: closure_context

    LOGICAL :: solid_excess_roundoff
    LOGICAL :: need_explicit_stage

    INTEGER :: newton_iterations
    INTEGER :: newton_linear_info
    LOGICAL :: newton_converged
    LOGICAL :: newton_line_search_failed
    INTEGER :: newton_calls_step
    INTEGER :: implicit_status
    INTEGER :: newton_iterations_step
    INTEGER :: newton_iterations_max_step
    INTEGER :: newton_failures_step
    INTEGER :: newton_linear_failures_step
    INTEGER :: newton_line_search_failures_step

    IF (.NOT.source_registry%is_initialized()) THEN
       CALL fatal_error('D17: IMEX advance requires an initialized source registry')
    END IF
    stage_groups=source_registry%stage_groups()
    source_solvers=source_registry%solvers()

    IF ( .NOT. equation_partition%is_initialized() ) THEN
       ERROR STOP 'imex_RK_solver: equation partition is not initialized'
    END IF

    IF ( equation_partition%n_equations .NE. n_eqns ) THEN
       ERROR STOP 'imex_RK_solver: equation partition size mismatch'
    END IF

    newton_calls_step = 0
    newton_iterations_step = 0
    newton_iterations_max_step = 0
    newton_failures_step = 0
    newton_linear_failures_step = 0
    newton_line_search_failures_step = 0
    
    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'solver, imex_RK_solver: beginning'

    !$OMP PARALLEL PRIVATE(j,k)

    ! The HP stencil reads two cells beyond solve_cells. Initialize its
    ! read-only halo to the physical dry state before active stages overwrite
    ! their own primitive variables.
    !$OMP DO
    DO l = 1,domain%reconstruction_cells
       j = domain%j_reconstruction(l)
       k = domain%k_reconstruction(l)
       this%qp_rk(:,j,k) = 0.0_wp
       this%qp_rk(4,j,k) = T_ambient
    END DO
    !$OMP END DO

    !$OMP DO
    DO l = 1,domain%solve_cells

       j = domain%j_cent(l)
       k = domain%k_cent(l)

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*) 'solver, imex_RK_solver: j,k',j,k

       END IF

       ! Initialization of the variables for the Runge-Kutta scheme
       this%q_rk( 1:n_vars , j , k ) = 0.0_wp
       this%spatial_terms(1:n_eqns,j,k,1:n_RK) = 0.0_wp
       this%IM_terms(1:n_eqns,j,k,1:n_RK) = 0.0_wp
       this%EX_terms(1:n_eqns,j,k,1:n_RK) = 0.0_wp
       
    END DO
    !$OMP END DO

    !$OMP END PARALLEL

    runge_kutta:DO i_RK = 1,n_RK

       IF ( verbose_level .GE. 1 ) WRITE(*,*) 'solver, imex_RK_solver: i_RK',i_RK

       ! An explicit stage is required not only when it contributes to the
       ! final RK assembly, but also when a later stage depends on it.
       need_explicit_stage = ( this%omega_tilde(i_RK) .NE. 0.0_wp )

       IF ( i_RK .LT. n_RK ) THEN
          need_explicit_stage = need_explicit_stage .OR.                       &
               ANY( this%a_tilde_ij(i_RK+1:n_RK,i_RK) .NE. 0.0_wp )
       END IF

       ! define the explicits coefficients for the i-th step of the Runge-Kutta
       this%a_tilde = 0.0_wp
       this%a_dirk = 0.0_wp

       ! in the first step of the RK scheme all the coefficients remain to 0
       this%a_tilde(1:i_RK-1) = this%a_tilde_ij(i_RK,1:i_RK-1)
       this%a_dirk(1:i_RK-1) = this%a_dirk_ij(i_RK,1:i_RK-1)

       ! define the implicit coefficient for the i-th step of the Runge-Kutta
       a_diag = this%a_dirk_ij(i_RK,i_RK)

       !$OMP PARALLEL 
       !$OMP DO schedule(guided)                                                &
       !$OMP & private(j,k,q_guess,q_fv_cell,known_stage_residual,p_dyn,       &
       !$OMP & closure_context,                                                &
       !$OMP & newton_iterations,newton_linear_info,newton_converged,          &
       !$OMP & newton_line_search_failed,implicit_status)

       solve_cells_loop:DO l = 1,domain%solve_cells

          j = domain%j_cent(l)
          k = domain%k_cent(l)
          implicit_status = 0

          IF ( verbose_level .GE. 2 ) THEN

             WRITE(*,*) 'solver, imex_RK_solver: j,k',j,k
          END IF

          ! initialize the RK step
          IF ( i_RK .EQ. 1 ) THEN

             ! solution from the previous time step
             q_guess(1:n_vars) = q( 1:n_vars , j , k)

          ELSE

             ! For stages after the first, q_guess is assembled below from
             ! q_fv_cell and the current implicit contribution.

          END IF

          ! Known part of the stage: assemble earlier spatial/explicit and
          ! implicit/semi-implicit contributions with their respective tableaux.
          ! The current diagonal implicit contribution is solved below.
          q_fv_cell(1:n_vars) = q( 1:n_vars , j , k )                            &
               - dt * (MATMUL( this%spatial_terms(1:n_eqns,j,k,1:i_RK)       &
               - this%EX_terms(1:n_eqns,j,k,1:i_RK) , this%a_tilde(1:i_RK) )  &
               - MATMUL( this%IM_terms(1:n_eqns,j,k,1:i_RK),                   &
               this%a_dirk(1:i_RK) ) )

          ! Capture the algebraic state before qc_to_qp can desingularize or
          ! map an inadmissible conservative state to a dry primitive state.
          IF ( PRESENT(observer) ) CALL observer(1, i_RK, j, k, q_fv_cell, 0)

          WRITE(closure_context,'(A,I0,A,I0,A,I0)')                           &
               'IMEX finite-volume stage=',i_RK,' cell=',j,',',k
          CALL clip_roundoff_negative_component_masses(q_fv_cell)
          CALL qc_to_qp(q_fv_cell , qp(1:n_vars+2,j,k) , p_dyn,               &
               TRIM(closure_context) )

          IF ( verbose_level .GE. 2 ) THEN

             WRITE(*,*) 'q_guess',q_guess
             IF ( q_guess(1) .GT. 0.0_wp  ) THEN 

                CALL qc_to_qp( q_guess , qp(1:n_vars+2,j,k) , p_dyn,         &
                     TRIM(closure_context)//' q_guess' )
                WRITE(*,*) 'q_guess: qp',qp(1:n_vars+2,j,k)

             END IF

          END IF

          adiag_pos:IF ( a_diag .NE. 0.0_wp ) THEN

             pos_thick:IF ( q_fv_cell(1) .GT.  0.0_wp )  THEN
                known_stage_residual = ( MATMUL(                             &
                     this%spatial_terms(1:n_eqns,j,k,1:i_RK-1) -            &
                     this%EX_terms(1:n_eqns,j,k,1:i_RK-1),                    &
                     this%a_tilde(1:i_RK-1) )                                 &
                     - MATMUL(this%IM_terms(1:n_eqns,j,k,1:i_RK-1),           &
                     this%a_dirk(1:i_RK-1) ) )

                ! The local stage coordinator applies SI and the ordered
                ! implicit-group sweep, returning one effective IM term.
                CALL solve_imex_stage(equation_partition, stage_groups,       &
                     source_solvers, q_fv_cell,                              &
                     q(1:n_vars,j,k), dt, a_diag,                             &
                     known_stage_residual, B_prime_x_geom(j,k),               &
                     B_prime_y_geom(j,k), B_second_xx_geom(j,k),              &
                     B_second_xy_geom(j,k), B_second_yy_geom(j,k),            &
                     grav_coeff(j,k), Z(j,k), comp_cells_x .EQ. 1,            &
                     comp_cells_y .EQ. 1, q_guess,                            &
                     this%IM_terms(1:n_eqns,j,k,i_RK),                        &
                     newton_iterations, newton_converged, newton_linear_info,   &
                     newton_line_search_failed )

                implicit_status = MERGE(1, -1, newton_converged)

                IF ( ( verbose_level .GE. 1 ) .OR.                             &
                     ( .NOT. newton_converged ) ) THEN

                   !$OMP CRITICAL(newton_diagnostics)

                   IF ( verbose_level .GE. 1 ) THEN
                      newton_calls_step = newton_calls_step + 1
                      newton_iterations_step = newton_iterations_step          &
                           + newton_iterations
                      newton_iterations_max_step = MAX(                        &
                           newton_iterations_max_step, newton_iterations )
                   END IF

                   IF ( .NOT. newton_converged ) THEN
                      newton_failures_step = newton_failures_step + 1
                      IF ( newton_linear_info .NE. 0 )                         &
                           newton_linear_failures_step =                       &
                           newton_linear_failures_step + 1
                      IF ( newton_line_search_failed )                         &
                           newton_line_search_failures_step =                  &
                           newton_line_search_failures_step + 1

                      IF ( verbose_level .GE. 1 ) THEN
                         WRITE(*,*) 'WARNING: Newton solve did not converge'
                         WRITE(*,*)                                           &
                              'cell, RK stage, iterations, linear info:',      &
                              j, k, i_RK, newton_iterations, newton_linear_info
                         WRITE(*,*) 'line search failed:',                     &
                              newton_line_search_failed
                      END IF
                   END IF

                   !$OMP END CRITICAL(newton_diagnostics)

                END IF
                
             ELSE

                ! If h=0 nothing has to be changed 
                q_guess(1:n_vars) = q_fv_cell
                this%IM_terms(1:n_eqns,j,k,i_RK) = 0.0_wp

             END IF pos_thick

          END IF adiag_pos

          ! Store the current stage. Previous stage states are no longer
          ! needed here: their evaluated terms are retained in spatial_terms,
          ! EX_terms and IM_terms.
          this%q_rk( 1:n_vars , j , k ) = q_guess
          IF ( PRESENT(observer) ) CALL observer(2, i_RK, j, k, q_guess, implicit_status)

          IF ( verbose_level .GE. 2 ) THEN

             WRITE(*,*) 'imex_RK_solver: qc',q_guess

             IF ( q_guess(1) .GT. 0.0_wp ) THEN

                CALL qc_to_qp( q_guess , qp(1:n_vars+2,j,k) , p_dyn,         &
                     TRIM(closure_context)//' stage output' )
                WRITE(*,*) 'imex_RK_solver: qp',qp(1:n_vars+2,j,k)

             END IF
             
             CALL debug_pause('IMEX Runge-Kutta stage state')

          END IF


          IF ( need_explicit_stage ) THEN
          
             IF ( this%q_rk(1,j,k) .GT. 0.0_wp ) THEN

                CALL qc_to_qp( this%q_rk(1:n_vars,j,k) ,                             &
                     this%qp_rk(1:n_vars+2,j,k) , p_dyn,                     &
                     TRIM(closure_context)//' explicit stage state' )

             ELSE

                this%qp_rk(1:n_vars+2,j,k) = 0.0_wp
                this%qp_rk(4,j,k) = T_ambient

             END IF

             ! Evaluate local explicit sources separately from the spatial
             ! hydrostatic path contribution handled by spatial_operator.
             CALL eval_expl_terms( B_prime_x_geom(j,k) , B_prime_y_geom(j,k) ,  &
                  B_second_xx_geom(j,k) , B_second_xy_geom(j,k) ,               &
                  B_second_yy_geom(j,k) , grav_coeff(j,k),                      &
                  this%qp_rk(1:n_vars+2,j,k), this%EX_terms(1:n_eqns,j,k,i_RK), t, &
                  cell_source_fractions(j,k), cell_fissure_fractions(j,k,:),    &
                  cell_arc_perim(j,k), cell_arc_n_x(j,k), cell_arc_n_y(j,k),    &
                  dx * dy )
  
          END IF

       END DO solve_cells_loop

       !$OMP END DO
       !$OMP END PARALLEL 

       IF ( need_explicit_stage ) THEN

          ! Evaluate the complete explicit spatial term. The time integrator
          ! is independent of the KT flux workspace.
          CALL spatial_operator%evaluate( this%qp_rk,                       &
               this%spatial_terms(:,:,:,i_RK), t, domain )

       END IF

    END DO runge_kutta

    !$OMP PARALLEL DO private(j,k,p_dyn,solid_excess_roundoff,                 &
    !$OMP & residual_cell,q_old_cell)

    assemble_sol:DO l = 1,domain%solve_cells

       j = domain%j_cent(l)
       k = domain%k_cent(l)

       ! q remains equal to Q^n throughout all RK stages. Preserve the old
       ! state locally before overwriting this cell during final assembly.
       q_old_cell = q(1:n_vars,j,k)

       residual_cell = MATMUL( this%spatial_terms(1:n_eqns,j,k,1:n_RK)       &
            - this%EX_terms(1:n_eqns,j,k,1:n_RK) , this%omega_tilde ) -      &
            MATMUL( this%IM_terms(1:n_eqns,j,k,1:n_RK), this%omega )


       IF ( ( SUM(ABS( this%omega_tilde(:)-this%a_tilde_ij(n_RK,:))) .EQ. 0.0_wp  )       &
            .AND. ( SUM(ABS(this%omega(:)-this%a_dirk_ij(n_RK,:))) .EQ. 0.0_wp ) ) THEN

          ! The assembling coeffs are equal to the last step of the RK scheme
          q(1:n_vars,j,k) = this%q_rk(1:n_vars,j,k)

       ELSE

          ! The assembling coeffs are different
          q(1:n_vars,j,k) = q_old_cell - dt*residual_cell

       END IF

       ! Record the final conservative assembly before any roundoff clipping,
       ! component repair or primitive-state conversion below.
       IF ( PRESENT(observer) ) CALL observer(3, 0, j, k, q(:,j,k), 0)

       IF ( ANY(ISNAN(q(:,j,k))) ) THEN
          CALL fatal_error('NaN detected after IMEX Runge-Kutta update')
          
       END IF

       CALL clip_roundoff_negative_component_masses(q(:,j,k))
       
       
       negative_thickness_check:IF ( q(1,j,k) .LT. EPSILON(1.0_wp) ) THEN

          IF ( q(1,j,k) .GT. -1.0E-7_wp ) THEN

             q(1,j,k) = 0.0_wp
             q(2:n_vars,j,k) = 0.0_wp

          ELSE
             CALL fatal_error('negative thickness after IMEX Runge-Kutta update')

          END IF

       END IF negative_thickness_check

       negative_alpha_check:IF ( ANY(q(5:4+n_solid,j,k) .LT. 0.0_wp ) ) THEN
          CALL fatal_error('negative solid mass after IMEX Runge-Kutta update')
       END IF negative_alpha_check

       CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars+2,j,k) , p_dyn )

       ! Do not carry unresolved momentum into a later rewetting event.
       ! The raw assembly has already been observed and checked above. Unlike
       ! the existing machine-epsilon cleanup, this projection retains mass,
       ! thermal energy and every transported component. Keep qp consistent
       ! without repeating the thermodynamic conversion (depth is unchanged).
       CALL project_dry_cell_momenta(q(:,j,k), qp(1,j,k))
       IF (qp(1,j,k) <= dry_thickness_tolerance) THEN
          qp(2:3,j,k) = 0.0_wp
          qp(idx_u:idx_v,j,k) = 0.0_wp
       END IF

       IF ( qp(4,j,k) .LT. 273.0_wp ) THEN
          CALL fatal_error('temperature below the admissible threshold')

       END IF
          
       IF ( SUM(q(5:4+n_solid,j,k)) .GT. q(1,j,k) ) THEN

          ! Fortran does not guarantee short-circuit evaluation of .OR.;
          ! evaluate the relative excess only when the denominator is safe.
          IF ( q(1,j,k) .LT. EPSILON(1.0_wp) ) THEN
             solid_excess_roundoff = .TRUE.
          ELSE
             solid_excess_roundoff = ( ( SUM(q(5:4+n_solid,j,k))             &
                  - q(1,j,k) ) / q(1,j,k) .LT. 1.0E-10_wp )
          END IF

          IF ( solid_excess_roundoff ) THEN

             CALL qc_to_qp(q_old_cell , qp(1:n_vars+2,j,k) , p_dyn )

             q(5:4+n_solid,j,k) = q(5:4+n_solid,j,k)                            &
                  / SUM(q(5:4+n_solid,j,k)) * q(1,j,k)

          ELSE
             CALL fatal_error('solid mass exceeds total mass after IMEX update')

          END IF

       END IF

       IF ( B_nodata(j,k) ) q(:,j,k) = 0.0_wp

    END DO assemble_sol

    !$OMP END PARALLEL DO

    IF ( ( verbose_level .GE. 1 ) .AND. ( newton_calls_step .GT. 0 ) ) THEN
       WRITE(*,*) 'Newton solves:',newton_calls_step
       WRITE(*,*) 'Newton iterations average/max:',                            &
            REAL(newton_iterations_step,wp) / REAL(newton_calls_step,wp),       &
            newton_iterations_max_step
       WRITE(*,*) 'Newton failures / linear failures:',                        &
            newton_failures_step,newton_linear_failures_step
       WRITE(*,*) 'Newton line-search failures:',                              &
            newton_line_search_failures_step
    ELSEIF ( newton_failures_step .GT. 0 ) THEN
       WRITE(*,*) 'WARNING: Newton failures / linear / line search:',           &
            newton_failures_step,newton_linear_failures_step,                  &
            newton_line_search_failures_step
    END IF
     
    RETURN

  END SUBROUTINE imex_RK_solver

END MODULE time_integration_2d
