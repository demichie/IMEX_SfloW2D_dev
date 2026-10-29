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

  USE parameters_2d, ONLY : wp
  USE parameters_2d, ONLY : n_thickness_levels, n_dyn_pres_levels

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: correct_friction_momentum, finalize_imex_cell_state
  PUBLIC :: initialize_dry_physical_cell

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

  !> \brief Constrain a momentum update from an isolated friction substep.
  !> \details Projects the candidate post-friction momentum onto the ray
  !> defined by the pre-friction momentum and clamps its scale factor to
  !> \f$[0,1]\f$. This preserves direction, prevents reversal, and ensures
  !> friction cannot increase momentum magnitude. If the pre-friction
  !> momentum is zero, the post-friction momentum is set to zero.
  !> \param[in] momentum_before_friction Momentum immediately before friction.
  !> \param[in,out] momentum_after_friction Candidate momentum immediately after friction.
  !******************************************************************************
  PURE SUBROUTINE correct_friction_momentum(momentum_before_friction,         &
       momentum_after_friction)

    USE parameters_2d, ONLY : wp
    REAL(wp), INTENT(IN) :: momentum_before_friction(2)
    REAL(wp), INTENT(INOUT) :: momentum_after_friction(2)

    REAL(wp) :: momentum_before_squared
    REAL(wp) :: friction_scale

    momentum_before_squared = DOT_PRODUCT(momentum_before_friction,           &
         momentum_before_friction)
    IF (momentum_before_squared .EQ. 0.0_wp) THEN
       momentum_after_friction = 0.0_wp
       RETURN
    END IF

    friction_scale = DOT_PRODUCT(momentum_before_friction,                   &
         momentum_after_friction) / momentum_before_squared
    friction_scale = MAX(0.0_wp, MIN(1.0_wp, friction_scale))
    momentum_after_friction = friction_scale * momentum_before_friction

  END SUBROUTINE correct_friction_momentum

  !> \brief Finalize a cell after the last IMEX Runge-Kutta assembly.
  !> \details Called after the observer sees the unmodified assembled state.
  !> Preserves legacy error thresholds, clipping, dry momentum projection,
  !> solid mass roundoff handling, and no-data masking in their original order.
  !> The previous conservative state is needed by the historical solid-mass
  !> roundoff branch, which refreshes qp before rescaling component masses.
  SUBROUTINE finalize_imex_cell_state(q_cell, qp_cell, q_previous, nodata)
    USE diagnostics_2d, ONLY : fatal_error
    USE state_conversion_2d, ONLY : qc_to_qp, project_dry_cell_momenta, &
         clip_roundoff_negative_component_masses
    USE parameters_2d, ONLY : n_vars, n_solid, dry_thickness_tolerance, idx_u, idx_v
    REAL(wp), INTENT(INOUT) :: q_cell(n_vars), qp_cell(n_vars+2)
    REAL(wp), INTENT(IN) :: q_previous(n_vars)
    LOGICAL, INTENT(IN) :: nodata
    REAL(wp) :: p_dyn
    LOGICAL :: solid_excess_roundoff

       IF ( ANY(ISNAN(q_cell)) ) THEN
          CALL fatal_error('NaN detected after IMEX Runge-Kutta update')
          
       END IF

       CALL clip_roundoff_negative_component_masses(q_cell)
       
       
       negative_thickness_check:IF ( q_cell(idx_totMassEqn) .LT. EPSILON(1.0_wp) ) THEN

          IF ( q_cell(idx_totMassEqn) .GT. -1.0E-7_wp ) THEN

             q_cell(idx_totMassEqn) = 0.0_wp
             q_cell(2:n_vars) = 0.0_wp

          ELSE
             CALL fatal_error('negative thickness after IMEX Runge-Kutta update')

          END IF

       END IF negative_thickness_check

       negative_alpha_check:IF ( ANY(q_cell(idx_solidEqn_first:idx_solidEqn_last) .LT. 0.0_wp ) ) THEN
          CALL fatal_error('negative solid mass after IMEX Runge-Kutta update')
       END IF negative_alpha_check

       CALL qc_to_qp(q_cell , qp_cell , p_dyn )

       ! Do not carry unresolved momentum into a later rewetting event.
       ! The raw assembly has already been observed and checked above. Unlike
       ! the existing machine-epsilon cleanup, this projection retains mass,
       ! thermal energy and every transported component. Keep qp consistent
       ! without repeating the thermodynamic conversion (depth is unchanged).
       CALL project_dry_cell_momenta(q_cell, qp_cell(idx_h))
       IF (qp_cell(idx_h) <= dry_thickness_tolerance) THEN
          qp_cell(idx_hu:idx_hv) = 0.0_wp
          qp_cell(idx_u:idx_v) = 0.0_wp
       END IF

       IF ( qp_cell(idx_T) .LT. 273.0_wp ) THEN
          CALL fatal_error('temperature below the admissible threshold')

       END IF
          
       IF ( SUM(q_cell(idx_solidEqn_first:idx_solidEqn_last)) .GT. q_cell(idx_totMassEqn) ) THEN

          ! Fortran does not guarantee short-circuit evaluation of .OR.;
          ! evaluate the relative excess only when the denominator is safe.
          IF ( q_cell(idx_totMassEqn) .LT. EPSILON(1.0_wp) ) THEN
             solid_excess_roundoff = .TRUE.
          ELSE
             solid_excess_roundoff = ( ( SUM(q_cell(idx_solidEqn_first:idx_solidEqn_last))             &
                  - q_cell(idx_totMassEqn) ) / q_cell(idx_totMassEqn) .LT. 1.0E-10_wp )
          END IF

          IF ( solid_excess_roundoff ) THEN

             CALL qc_to_qp(q_previous , qp_cell , p_dyn )

             q_cell(idx_solidEqn_first:idx_solidEqn_last) = q_cell(idx_solidEqn_first:idx_solidEqn_last) &
                  / SUM(q_cell(idx_solidEqn_first:idx_solidEqn_last)) * q_cell(idx_totMassEqn)

          ELSE
             CALL fatal_error('solid mass exceeds total mass after IMEX update')

          END IF

       END IF

       IF ( nodata ) q_cell = 0.0_wp
  END SUBROUTINE finalize_imex_cell_state

  !> \brief Set one primitive state (including reconstruction halo) to dry.
  !> \details Preserve the established zero state and ambient temperature.
  PURE SUBROUTINE initialize_dry_physical_cell(qp_cell)
    REAL(wp), INTENT(OUT) :: qp_cell(:)
    qp_cell = 0.0_wp
    qp_cell(idx_T) = T_ambient
  END SUBROUTINE initialize_dry_physical_cell

  !******************************************************************************
END MODULE state_2d
