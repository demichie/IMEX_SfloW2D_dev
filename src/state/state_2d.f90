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
    USE parameters_2d, ONLY : idx_h, idx_hu, idx_hv, idx_T, &
       idx_totMassEqn, idx_solidEqn_first, idx_solidEqn_last

  USE constitutive_parameters_2d, ONLY : T_ambient

  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE model_layout_2d, ONLY : model_layout_type

  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
  USE parameters_2d, ONLY : n_eqns,idx_add_gas_first,idx_add_gas_last
  USE parameters_2d, ONLY : wp
  USE parameters_2d, ONLY : n_thickness_levels, n_dyn_pres_levels

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: initialize_dry_physical_cell
  PUBLIC :: conservative_mass_admissible
  PUBLIC :: physical_cell_state_finite

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

  !> Validate a primitive cell state without changing it.
  PURE LOGICAL FUNCTION physical_cell_state_finite(qp_cell) RESULT(valid)
    REAL(wp), INTENT(IN) :: qp_cell(:)
    valid=ALL(ieee_is_finite(qp_cell))
  END FUNCTION physical_cell_state_finite

  !> Check finite conservative mass budgets without modifying the state.
  PURE LOGICAL FUNCTION conservative_mass_admissible(q_cell) RESULT(valid)
    REAL(wp), INTENT(IN) :: q_cell(n_eqns)
    valid=.FALSE.
    IF (.NOT.ALL(ieee_is_finite(q_cell))) RETURN
    IF (q_cell(idx_totMassEqn)<0.0_wp) RETURN
    IF (ANY(q_cell(idx_solidEqn_first:idx_solidEqn_last)<0.0_wp)) RETURN
    IF (ANY(q_cell(idx_add_gas_first:idx_add_gas_last)<0.0_wp)) RETURN
    IF (SUM(q_cell(idx_solidEqn_first:idx_solidEqn_last))>q_cell(idx_totMassEqn)) RETURN
    valid=.TRUE.
  END FUNCTION conservative_mass_admissible



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
    this%qp(model_layout%idx_energy,1:comp_cells_x,1:comp_cells_y) = T_ambient

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

  !> \brief Set one primitive state (including reconstruction halo) to dry.
  !> \details Preserve the established zero state and ambient temperature.
  PURE SUBROUTINE initialize_dry_physical_cell(qp_cell)
    REAL(wp), INTENT(OUT) :: qp_cell(:)
    qp_cell = 0.0_wp
    qp_cell(idx_T) = T_ambient
  END SUBROUTINE initialize_dry_physical_cell

  !******************************************************************************
END MODULE state_2d
