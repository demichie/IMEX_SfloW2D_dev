!********************************************************************************
!> \brief Solver state ownership
!>
!> This module owns the conservative and physical state arrays together with
!> diagnostics that evolve during a simulation.
!********************************************************************************
MODULE state_2d

  USE constitutive_2d, ONLY : T_ambient

  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y

  USE parameters_2d, ONLY : wp
  USE parameters_2d, ONLY : n_vars
  USE parameters_2d, ONLY : n_thickness_levels, n_dyn_pres_levels

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: initialize_state
  PUBLIC :: finalize_state

  PUBLIC :: q, qp
  PUBLIC :: hpos, hpos_old
  PUBLIC :: hmax, pdynmax, mod_vel_max
  PUBLIC :: vuln_table, thck_table, pdyn_table

  !> Conservative variables
  REAL(wp), ALLOCATABLE :: q(:,:,:)

  !> Physical variables
  REAL(wp), ALLOCATABLE :: qp(:,:,:)

  !> Map of positive thickness
  LOGICAL, ALLOCATABLE :: hpos(:,:)

  !> Map of positive thickness at previous output step
  LOGICAL, ALLOCATABLE :: hpos_old(:,:)

  !> Maximum over time of thickness
  REAL(wp), ALLOCATABLE :: hmax(:,:)

  !> Maximum over time of dynamic pressure
  REAL(wp), ALLOCATABLE :: pdynmax(:,:)

  !> Maximum over time of velocity magnitude
  REAL(wp), ALLOCATABLE :: mod_vel_max(:,:)

  LOGICAL, ALLOCATABLE :: vuln_table(:,:,:)
  LOGICAL, ALLOCATABLE :: thck_table(:,:)
  LOGICAL, ALLOCATABLE :: pdyn_table(:,:)

CONTAINS

  !******************************************************************************
  !> \brief Allocate the prognostic state and diagnostic arrays.
  !******************************************************************************

  SUBROUTINE initialize_state

    ALLOCATE( q(n_vars,comp_cells_x,comp_cells_y) )

    ALLOCATE( hpos(comp_cells_x,comp_cells_y),                              &
         hpos_old(comp_cells_x,comp_cells_y) )

    ALLOCATE( qp(n_vars+2,comp_cells_x,comp_cells_y) )

    q(1:n_vars,1:comp_cells_x,1:comp_cells_y) = 0.0_wp
    qp(1:n_vars+2,1:comp_cells_x,1:comp_cells_y) = 0.0_wp
    qp(4,1:comp_cells_x,1:comp_cells_y) = T_ambient

    ALLOCATE( hmax(comp_cells_x,comp_cells_y) )
    ALLOCATE( pdynmax(comp_cells_x,comp_cells_y) )
    ALLOCATE( mod_vel_max(comp_cells_x,comp_cells_y) )

    ALLOCATE( vuln_table(n_thickness_levels*n_dyn_pres_levels,              &
         comp_cells_x,comp_cells_y) )

    ALLOCATE( thck_table(comp_cells_x,comp_cells_y) )
    ALLOCATE( pdyn_table(comp_cells_x,comp_cells_y) )

  END SUBROUTINE initialize_state

  !******************************************************************************
  !> \brief Deallocate the prognostic state and diagnostic arrays.
  !******************************************************************************

  SUBROUTINE finalize_state

    DEALLOCATE( q, hpos, hpos_old )

    DEALLOCATE( hmax, pdynmax, mod_vel_max )

    DEALLOCATE( vuln_table )

    DEALLOCATE( thck_table, pdyn_table )

    DEALLOCATE( qp )

  END SUBROUTINE finalize_state

END MODULE state_2d
