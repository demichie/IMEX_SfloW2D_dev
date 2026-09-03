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

  TYPE, PUBLIC :: state_type

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

     PROCEDURE :: initialize => initialize_state
     PROCEDURE :: finalize => finalize_state

  END TYPE state_type

  TYPE(state_type), PUBLIC, TARGET :: state

CONTAINS

  !******************************************************************************
  !> \brief Allocate the prognostic state and diagnostic arrays.
  !******************************************************************************

  SUBROUTINE initialize_state(this)

    CLASS(state_type), INTENT(INOUT) :: this

    ALLOCATE( this%q(n_vars,comp_cells_x,comp_cells_y) )

    ALLOCATE( this%hpos(comp_cells_x,comp_cells_y),                         &
         this%hpos_old(comp_cells_x,comp_cells_y) )

    ALLOCATE( this%qp(n_vars+2,comp_cells_x,comp_cells_y) )

    this%q(1:n_vars,1:comp_cells_x,1:comp_cells_y) = 0.0_wp
    this%qp(1:n_vars+2,1:comp_cells_x,1:comp_cells_y) = 0.0_wp
    this%qp(4,1:comp_cells_x,1:comp_cells_y) = T_ambient

    ALLOCATE( this%hmax(comp_cells_x,comp_cells_y) )
    ALLOCATE( this%pdynmax(comp_cells_x,comp_cells_y) )
    ALLOCATE( this%mod_vel_max(comp_cells_x,comp_cells_y) )

    ALLOCATE( this%vuln_table(n_thickness_levels*n_dyn_pres_levels,         &
         comp_cells_x,comp_cells_y) )

    ALLOCATE( this%thck_table(comp_cells_x,comp_cells_y) )
    ALLOCATE( this%pdyn_table(comp_cells_x,comp_cells_y) )

  END SUBROUTINE initialize_state

  !******************************************************************************
  !> \brief Deallocate the prognostic state and diagnostic arrays.
  !******************************************************************************

  SUBROUTINE finalize_state(this)

    CLASS(state_type), INTENT(INOUT) :: this

    DEALLOCATE( this%q, this%hpos, this%hpos_old )

    DEALLOCATE( this%hmax, this%pdynmax, this%mod_vel_max )

    DEALLOCATE( this%vuln_table )

    DEALLOCATE( this%thck_table, this%pdyn_table )

    DEALLOCATE( this%qp )

  END SUBROUTINE finalize_state

END MODULE state_2d
