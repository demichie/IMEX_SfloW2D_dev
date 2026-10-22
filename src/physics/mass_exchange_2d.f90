!********************************************************************************
!> \brief Cell-local erosion, deposition, entrainment and evolving-bed updates.
!> \details Evaluates cell-local mass-exchange sources through
!>          equation_terms_2d and delegates conservative-state, inventory and
!>          admissibility updates to time_integration_2d. When bed evolution is
!>          enabled, limited cell elevation-rate proposals are accumulated,
!>          projected to the nodal bed once per update, and followed by a
!>          refresh of derived geometry.
!>          The persistent work arrays belong to this module and can be
!>          released explicitly when a simulation is finalized.
!
!********************************************************************************

MODULE mass_exchange_2d

  USE constitutive_parameters_2d, ONLY : T_ambient

  USE geometry_2d, ONLY : B_cent, B_vertex
  USE geometry_2d, ONLY : B_prime_x_geom, B_prime_y_geom
  USE geometry_2d, ONLY : cell_source_fractions
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE geometry_2d, ONLY : comp_interfaces_x, comp_interfaces_y
  USE geometry_2d, ONLY : project_cell_field_to_vertices
  USE geometry_2d, ONLY : refresh_topography_geometry

  USE parameters_2d, ONLY : wp
  USE parameters_2d, ONLY : n_eqns, n_vars, n_solid
  USE parameters_2d, ONLY : bottom_radial_source_flag

  USE domain_2d, ONLY : domain_type

  USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: update_erosion_deposition_cell, release_topography_workspace

  ! Persistent work arrays for the vertex-first evolving-bed update. Inactive
  ! and dry cells retain a zero proposal but remain in the geometric stencil.
  REAL(wp), ALLOCATABLE :: topography_rate_cell(:,:)
  REAL(wp), ALLOCATABLE :: topography_rate_vertex(:,:)


CONTAINS

  !******************************************************************************
  !> \brief Advance cell-local mass exchange and optionally update the shared nodal bed.
  !> \details Runs as a post-IMEX fractional step. For each active solve cell,
  !>          reconstructs the physical state and evaluates mass-exchange
  !>          sources; time_integration_2d applies the conservative and
  !>          inventory increments and their admissibility checks. If
  !>          topography changes, each limited cell bed-rate proposal is
  !>          projected to vertices and integrated once after the cell loop,
  !>          before the caller advances the simulation clock.
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
    
    USE geometry_2d, ONLY : deposit , erosion , erodible
    USE geometry_2d, ONLY : B_zone

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

    IF ( topo_change_flag ) THEN
       CALL ensure_topography_workspace
       !$OMP PARALLEL DO COLLAPSE(2)
       DO k = 1, comp_cells_y
          DO j = 1, comp_cells_x
             topography_rate_cell(j,k) = 0.0_wp
          END DO
       END DO
       !$OMP END PARALLEL DO
    END IF

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
       IF ( topo_change_flag ) topography_rate_cell(j,k) = topo_term

    END DO

    !$OMP END PARALLEL DO

    IF ( topo_change_flag ) CALL apply_vertex_first_topography_update(dt)

    RETURN

  END SUBROUTINE update_erosion_deposition_cell

  !******************************************************************************
  !> \brief Ensure that persistent evolving-bed proposal arrays match the current grid.
  !> \details Allocates the cell-rate and vertex-rate arrays on first use. If
  !>          the computational grid dimensions changed, releases the old
  !>          arrays before allocating storage with the current dimensions.
  !>          Existing proposals are not preserved across a resize.
  !>
  !> \note Uses comp_cells_x/comp_cells_y and comp_interfaces_x/comp_interfaces_y;
  !>       this routine has no arguments or return value.
  !******************************************************************************

  SUBROUTINE ensure_topography_workspace

    IMPLICIT NONE

    LOGICAL :: grid_size_changed

    grid_size_changed = .FALSE.
    IF ( ALLOCATED(topography_rate_cell) ) THEN
       grid_size_changed = ( SIZE(topography_rate_cell,1) .NE. comp_cells_x ) &
            .OR. ( SIZE(topography_rate_cell,2) .NE. comp_cells_y )
    END IF

    IF ( grid_size_changed ) THEN
       CALL release_topography_workspace
    END IF

    IF ( .NOT.ALLOCATED(topography_rate_cell) ) THEN
       ALLOCATE(topography_rate_cell(comp_cells_x,comp_cells_y))
       ALLOCATE(topography_rate_vertex(comp_interfaces_x,comp_interfaces_y))
    END IF

  END SUBROUTINE ensure_topography_workspace

  !******************************************************************************
  !> \brief Project cell bed rates to vertices and refresh the derived geometry.
  !> \details Projects the already limited topography_rate_cell field onto
  !>          vertices, advances the authoritative nodal bed by
  !>          dt*topography_rate_vertex, and regenerates all dependent
  !>          center, face, slope and curvature fields.
  !> \param[in] dt Time increment [s].
  !> \note Requires the evolving-bed workspace to have been allocated and populated by
  !>       update_erosion_deposition_cell.
  !******************************************************************************

  SUBROUTINE apply_vertex_first_topography_update(dt)

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: dt
    INTEGER :: j, k

    CALL project_cell_field_to_vertices(topography_rate_cell,                 &
         topography_rate_vertex)

    !$OMP PARALLEL DO COLLAPSE(2)
    DO k = 1, comp_interfaces_y
       DO j = 1, comp_interfaces_x
          B_vertex(j,k) = B_vertex(j,k) + dt * topography_rate_vertex(j,k)
       END DO
    END DO
    !$OMP END PARALLEL DO

    CALL refresh_topography_geometry

  END SUBROUTINE apply_vertex_first_topography_update

  !******************************************************************************
  !> \brief Release the persistent evolving-bed proposal arrays.
  !> \details Deallocates the private cell- and vertex-rate work arrays when
  !>          allocated. It does not modify the authoritative bed or derived
  !>          geometry, and repeated calls are safe.
  !> \note This routine has no arguments or return value.
  !******************************************************************************

  SUBROUTINE release_topography_workspace
    IF (ALLOCATED(topography_rate_cell)) DEALLOCATE(topography_rate_cell)
    IF (ALLOCATED(topography_rate_vertex)) DEALLOCATE(topography_rate_vertex)
  END SUBROUTINE release_topography_workspace

END MODULE mass_exchange_2d
