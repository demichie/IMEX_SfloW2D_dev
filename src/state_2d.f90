!********************************************************************************
!> \brief Solver state ownership
!>
!> This module owns the conservative and physical state arrays together with
!> diagnostics that evolve during a simulation.
!>
!> Stores conservative/physical cell states and accumulated hazard/runout diagnostics. Configuration
!> and numerical scratch storage belong to other modules.
!********************************************************************************

MODULE state_2d

  USE constitutive_parameters_2d, ONLY : T_ambient

  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE model_layout_2d, ONLY : model_layout_type

  USE parameters_2d, ONLY : wp
  USE parameters_2d, ONLY : n_thickness_levels, n_dyn_pres_levels

  IMPLICIT NONE

  PRIVATE

  !> \brief Owned cell states and diagnostics for one simulation.
  !> \details State arrays use (component,x-cell,y-cell); diagnostic maps use
  !>          (x-cell,y-cell). The model layout fixes component counts and indices.
  TYPE, PUBLIC :: state_type

     !> Conservative variables: M=rho_m*h, M*u, M*v, M*cp*T, then transported masses.
     REAL(wp), ALLOCATABLE :: q(:,:,:)

     !> Physical variables: h, h*u, h*v, T, mass fractions, and appended u and v.
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

   CONTAINS

     PROCEDURE :: initialize => initialize_state
     PROCEDURE :: finalize => finalize_state

  END TYPE state_type

CONTAINS

  !******************************************************************************
  !> \brief Allocate and reset cell states and accumulated diagnostics for one model layout.
  !>
  !> \param[in,out] this Allocated conservative/physical cell state and accumulated diagnostics.
  !> \param[in] model_layout Validated variable counts and canonical component indices.
  !******************************************************************************

  SUBROUTINE initialize_state(this, model_layout)

    CLASS(state_type), INTENT(INOUT) :: this
    TYPE(model_layout_type), INTENT(IN) :: model_layout

    IF (.NOT. model_layout%is_initialized()) THEN
       ERROR STOP 'State allocation requires an initialized model layout'
    END IF

    ALLOCATE( this%q(model_layout%n_variables,comp_cells_x,comp_cells_y) )

    ALLOCATE( this%hpos(comp_cells_x,comp_cells_y),                         &
         this%hpos_old(comp_cells_x,comp_cells_y) )

    ALLOCATE( this%qp(model_layout%n_physical_variables,comp_cells_x,        &
         comp_cells_y) )

    this%q = 0.0_wp
    this%qp = 0.0_wp
    this%qp(4,1:comp_cells_x,1:comp_cells_y) = T_ambient

    ALLOCATE( this%hmax(comp_cells_x,comp_cells_y) )
    ALLOCATE( this%pdynmax(comp_cells_x,comp_cells_y) )
    ALLOCATE( this%mod_vel_max(comp_cells_x,comp_cells_y) )

    this%hpos = .FALSE.
    this%hpos_old = .FALSE.
    this%hmax = 0.0_wp
    this%pdynmax = 0.0_wp
    this%mod_vel_max = 0.0_wp

    ALLOCATE( this%vuln_table(n_thickness_levels*n_dyn_pres_levels,         &
         comp_cells_x,comp_cells_y) )


    this%vuln_table = .FALSE.

  END SUBROUTINE initialize_state

  !******************************************************************************
  !> \brief Release the physical/conservative state and diagnostic fields.
  !>
  !> \param[in,out] this Allocated conservative/physical cell state and accumulated diagnostics.
  !******************************************************************************

  SUBROUTINE finalize_state(this)

    CLASS(state_type), INTENT(INOUT) :: this

    DEALLOCATE( this%q, this%hpos, this%hpos_old )

    DEALLOCATE( this%hmax, this%pdynmax, this%mod_vel_max )

    DEALLOCATE( this%vuln_table )


    DEALLOCATE( this%qp )

  END SUBROUTINE finalize_state

END MODULE state_2d
