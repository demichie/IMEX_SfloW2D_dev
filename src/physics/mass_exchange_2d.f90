!********************************************************************************
!> \brief Cell-local erosion, deposition, entrainment and evolving-bed updates.
!> \details Evaluates cell-local mass-exchange sources through
!>          equation_terms_2d and delegates conservative-state, inventory and
!>          admissibility updates to time_integration_2d. Limited bed-rate
!>          proposals are handed to geometry_2d for vertex projection and
!>          derived-geometry refresh.
!
!********************************************************************************

MODULE mass_exchange_2d

  USE constitutive_parameters_2d, ONLY : T_ambient

  USE geometry_2d, ONLY : B_prime_x_geom, B_prime_y_geom
  USE geometry_2d, ONLY : cell_source_fractions
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE geometry_2d, ONLY : begin_topography_rate_update,                      &
       record_topography_rate_cell, apply_topography_rate_update

  USE parameters_2d, ONLY : wp
  USE parameters_2d, ONLY : n_eqns, n_vars, n_solid
  USE parameters_2d, ONLY : bottom_radial_source_flag

  USE domain_2d, ONLY : domain_type

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: update_erosion_deposition_cell


CONTAINS

  !******************************************************************************
  !> \brief Advance cell-local mass exchange and optionally update the shared nodal bed.
  !> \details Runs as a post-IMEX fractional step. For each active solve cell,
  !>          reconstructs the physical state and evaluates mass-exchange
  !>          sources; time_integration_2d applies conservative and inventory
  !>          increments. If topography changes, geometry_2d projects the
  !>          limited cell bed rates after the loop, before the caller advances
  !>          the simulation clock.
  !> \param[in,out] q Cell-centered conservative state [variable,x-cell,y-cell]; updated by
  !>                  the mass-exchange equation source.
  !> \param[in,out] qp Cell-centered physical state with h, temperature, component fractions,
  !>                   and appended velocity components; refreshed for solved cells.
  !> \param[in] dt Physical time increment [s].
  !> \param[in] domain Domain descriptor containing the active solve-cell count and the
  !>                  corresponding j_cent/k_cent cell indices.
  !******************************************************************************

  SUBROUTINE update_erosion_deposition_cell(q, qp, dt, domain)

    USE constitutive_parameters_2d, ONLY : erosion_coeff, settling_flag,    &
         entrainment_flag, loss_rate
    
    USE geometry_2d, ONLY : deposit, erosion, erodible, B_zone

    USE equation_terms_2d, ONLY : eval_mass_exchange_terms
    USE time_integration_2d, ONLY : apply_mass_exchange_increment

    USE state_conversion_2d, ONLY : qc_to_qp
    USE parameters_2d, ONLY : topo_change_flag, bottom_radial_source_flag, &
         bottom_fissural_source_flag
    USE parameters_2d, ONLY : erodible_deposit_flag
    USE parameters_2d, ONLY : pore_pressure_flag
    USE parameters_2d, ONLY : liquid_flag

    IMPLICIT NONE

    REAL(wp), INTENT(INOUT) :: q(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(INOUT) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: dt
    CLASS(domain_type), INTENT(IN) :: domain

    REAL(wp) :: erosion_term(n_solid)
    REAL(wp) :: deposition_term(n_solid)
    REAL(wp) :: continuous_phase_erosion_term
    REAL(wp) :: continuous_phase_loss_term
    REAL(wp) :: eqns_term(n_eqns)
    REAL(wp) :: topo_term

    INTEGER :: j,k,l

    REAL(wp) :: out_of_source_fraction

    REAL(wp) :: p_dyn
    CHARACTER(LEN=96) :: closure_context

    LOGICAL :: liquid_loss_active



    ! Carrier loss is a mass exchange even without erosion or settling. Use
    ! nested guards: Fortran does not require short-circuit evaluation, and
    ! gas-only configurations need not allocate the liquid loss parameter.
    liquid_loss_active = .FALSE.
    IF ( liquid_flag ) THEN
       IF ( ALLOCATED(loss_rate) ) liquid_loss_active = loss_rate .GT. 0.0_wp
    END IF
    IF ( ( erosion_coeff .EQ. 0.0_wp ) .AND. ( .NOT.settling_flag ) &
         .AND. ( .NOT.pore_pressure_flag ) .AND. ( .NOT.entrainment_flag) &
         .AND. ( .NOT.liquid_loss_active ) ) RETURN

    IF ( topo_change_flag ) CALL begin_topography_rate_update

    !$OMP PARALLEL DO private(j,k,closure_context,erosion_term,deposition_term, &
    !$OMP & eqns_term,                                                          &
    !$OMP & topo_term,                                                          &
    !$OMP & continuous_phase_erosion_term,continuous_phase_loss_term,           &
    !$OMP & out_of_source_fraction,p_dyn)

    DO l = 1,domain%solve_cells

       j = domain%j_cent(l)
       k = domain%k_cent(l)
       WRITE(closure_context,'(A,I0,A,I0)') 'mass exchange cell=',j,',',k

       IF ( q(1,j,k) .GT. 0.0_wp ) THEN

          CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars+2,j,k) , p_dyn,        &
               TRIM(closure_context) )

       ELSE

          qp(1:n_vars+2,j,k) = 0.0_wp
          qp(4,j,k) = T_ambient

       END IF

       CALL eval_mass_exchange_terms( qp(1:n_vars+2,j,k) , B_zone(j,k) ,           &
            B_prime_x_geom(j,k) , B_prime_y_geom(j,k) , erodible(1:n_solid,j,k) ,  &
            dt , erosion_term , deposition_term , continuous_phase_erosion_term ,  &
            continuous_phase_loss_term , eqns_term , topo_term  )
          
       IF ( bottom_radial_source_flag .OR. bottom_fissural_source_flag ) THEN

          ! Both bottom-source geometries contribute to this shared coverage
          ! fraction. Apply the same mask to all increments before recording
          ! inventories or projecting the limited proposal onto the nodal bed.
          out_of_source_fraction = 1.0_wp - cell_source_fractions(j,k)
          deposition_term = deposition_term * out_of_source_fraction
          erosion_term = erosion_term * out_of_source_fraction
          eqns_term = eqns_term * out_of_source_fraction
          topo_term = topo_term * out_of_source_fraction

       END IF
       
       CALL apply_mass_exchange_increment(q(1:n_vars,j,k),                    &
            qp(1:n_vars+2,j,k), dt, eqns_term, deposition_term, erosion_term, &
            deposit(j,k,1:n_solid), erosion(j,k,1:n_solid),                   &
            erodible(1:n_solid,j,k), erodible_deposit_flag, closure_context)

       ! Store the final, already limited cell proposal. The nodal bed is
       ! assembled only after every cell-local flow/inventory update is done.
       IF ( topo_change_flag) CALL record_topography_rate_cell(j, k, topo_term)

    END DO

    !$OMP END PARALLEL DO

    IF ( topo_change_flag ) CALL apply_topography_rate_update(dt)

    RETURN

  END SUBROUTINE update_erosion_deposition_cell

END MODULE mass_exchange_2d
