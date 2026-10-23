!********************************************************************************
!> \brief Model-facing explicit spatial operator
!>
!> This module is the boundary between time integration and the numerical
!> discretization of spatial terms. It owns the KT backend and
!> exposes only complete spatial-term and CFL evaluations.  Time integration
!> therefore does not depend on flux, path-contribution, or wave-speed storage.
!********************************************************************************

MODULE spatial_operator_2d

  USE parameters_2d, ONLY : wp, n_eqns, n_vars, max_dt, cfl, cfl_policy
  USE diagnostics_2d, ONLY : fatal_error
  USE constitutive_parameters_2d, ONLY : T_ambient
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y, dx, dy
  USE state_conversion_2d, ONLY : qc_to_qp
  USE domain_2d, ONLY : domain_type
  USE reconstruction_2d, ONLY : reconstruction_workspace_type
  USE hyperbolic_2d, ONLY : hyperbolic_workspace_type

  IMPLICIT NONE

  PRIVATE

  !> \brief Numerical backend owned by one model-facing spatial operator.
  !> \details Its private reconstruction/flux workspaces are reused by both
  !>          spatial-term evaluation and the CFL calculation.
  TYPE, PUBLIC :: spatial_operator_type
     PRIVATE
     TYPE(reconstruction_workspace_type) :: reconstruction
     TYPE(hyperbolic_workspace_type) :: hp_pccu
   CONTAINS
     PROCEDURE, PUBLIC :: initialize => initialize_spatial_operator
     PROCEDURE, PUBLIC :: finalize => finalize_spatial_operator
     PROCEDURE, PUBLIC :: evaluate => evaluate_spatial_operator
     PROCEDURE, PUBLIC :: compute_timestep => compute_spatial_timestep
  END TYPE spatial_operator_type

CONTAINS

  !> \brief Initialize the reconstruction and KT flux workspace.
  !>
  !> \param[in,out] this Spatial operator owning reconstruction and HP-PCCU backend storage.

  SUBROUTINE initialize_spatial_operator( this )

    CLASS(spatial_operator_type), INTENT(INOUT) :: this

    CALL this%reconstruction%initialize
    CALL this%hp_pccu%initialize

  END SUBROUTINE initialize_spatial_operator

  !> \brief Release the reconstruction and KT flux workspace.
  !>
  !> \param[in,out] this Spatial operator owning reconstruction and HP-PCCU backend storage.

  SUBROUTINE finalize_spatial_operator( this )

    CLASS(spatial_operator_type), INTENT(INOUT) :: this

    CALL this%hp_pccu%finalize
    CALL this%reconstruction%finalize

  END SUBROUTINE finalize_spatial_operator

  !******************************************************************************
  !> \brief Evaluate the complete explicit spatial term on the active domain.
  !>
  !> The returned term includes conservative interface transport and the
  !> cell/interface path contributions assembled by HP-PCCU.  Its sign follows
  !> q_t + spatial_term = local_sources.
  !>
  !> \param[in,out] this Spatial operator owning reconstruction and HP-PCCU backend storage.
  !> \param[in] qp Canonical physical states: h, hu, hv, T, component mass fractions, then appended
  !>               u and v.
  !> \param[out] spatial_term Complete explicit spatial term with sign
  !>                          q_t+spatial_term=local_sources.
  !> \param[in] time Current simulation time [s].
  !> \param[in] domain Active-cell/face lists and reconstruction halo for this simulation.
  !>
  !> \note The sign convention is q_t+spatial_term=local_sources; hydrostatic pressure/topography
  !>       are already included through path terms.
  !******************************************************************************

  SUBROUTINE evaluate_spatial_operator( this, qp, spatial_term, time, domain, &
       boundary_source_enabled )

    CLASS(spatial_operator_type), INTENT(INOUT) :: this
    REAL(wp), INTENT(IN) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(OUT) :: spatial_term(n_eqns,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: time
    CLASS(domain_type), INTENT(IN) :: domain
    LOGICAL, INTENT(IN), OPTIONAL :: boundary_source_enabled

    CALL this%hp_pccu%evaluate_terms( this%reconstruction, qp, spatial_term, &
         time, domain%solve_cells, domain%j_cent, domain%k_cent,              &
         domain%solve_interfaces_x, domain%j_stag_x, domain%k_stag_x,        &
         domain%solve_interfaces_y, domain%j_stag_y, domain%k_stag_y,        &
         source_boundary_enabled=boundary_source_enabled )

  END SUBROUTINE evaluate_spatial_operator

  !******************************************************************************
  !> \brief Refresh the reconstruction states and compute the face-speed CFL timestep.
  !>
  !> \param[in,out] this Spatial operator owning reconstruction and HP-PCCU backend storage.
  !> \param[in] q Cell-centered conservative states, indexed as (variable,x-cell,y-cell).
  !> \param[in,out] qp Physical state cache refreshed from q on the reconstruction workset.
  !> \param[in] time Current simulation time [s].
  !> \param[out] dt Returned CFL-limited timestep, bounded by max_dt [s].
  !> \param[in] domain Active-cell/face lists and reconstruction halo for this simulation.
  !******************************************************************************

  SUBROUTINE compute_spatial_timestep( this, q, qp, time, dt, domain,        &
       boundary_source_enabled )

    CLASS(spatial_operator_type), INTENT(INOUT) :: this
    REAL(wp), INTENT(IN) :: q(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(INOUT) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: time
    REAL(wp), INTENT(OUT) :: dt
    CLASS(domain_type), INTENT(IN) :: domain
    LOGICAL, INTENT(IN), OPTIONAL :: boundary_source_enabled

    INTEGER :: j, k, l
    REAL(wp) :: max_a_x, max_a_y, cfl_spatial
    REAL(wp) :: dynamic_pressure, propagation_rate

    dt = max_dt

    ! The sentinel selects the prescribed max_dt without evaluating face speeds.
    IF ( cfl .EQ. -1.0_wp ) RETURN

    ! Refresh wet cells and the dry halo: the HP stencil reads beyond solve_cells.
    ! The end-of-loop barrier makes every cache entry available to reconstruction.
    !$OMP PARALLEL DO PRIVATE(j,k,dynamic_pressure)
    DO l = 1, domain%reconstruction_cells
       j = domain%j_reconstruction(l)
       k = domain%k_reconstruction(l)

       IF ( q(1,j,k) .GT. 0.0_wp ) THEN
          CALL qc_to_qp( q(:,j,k), qp(:,j,k), dynamic_pressure )
       ELSE
          qp(:,j,k) = 0.0_wp
          qp(4,j,k) = T_ambient
       END IF
    END DO
    !$OMP END PARALLEL DO

    ! CFL speeds must use the same selected advective traces as the spatial operator.
    CALL this%reconstruction%reconstruct( qp, time, domain%solve_cells,      &
         domain%j_cent, domain%k_cent,                                       &
         .TRUE., source_boundary_enabled=boundary_source_enabled )

    ! Match the KT face depths and desingularized velocities used by the
    ! spatial operator before evaluating the CFL wave speeds.
    CALL this%reconstruction%reconstruct_kt_eta( qp )

    CALL this%hp_pccu%evaluate_speeds( this%reconstruction,                  &
         domain%solve_interfaces_x, domain%j_stag_x, domain%k_stag_x,        &
         domain%solve_interfaces_y, domain%j_stag_y, domain%k_stag_y )

    max_a_x = 0.0_wp
    max_a_y = 0.0_wp

    ! Each active cell contributes all of its faces; a global maximum gives
    ! one timestep valid for every participating OpenMP thread.
    !$OMP PARALLEL DO PRIVATE(j,k) REDUCTION(MAX:max_a_x,max_a_y)
    DO l = 1, domain%solve_cells
       j = domain%j_cent(l)
       k = domain%k_cent(l)

       CALL accumulate_directional_face_speeds(this%hp_pccu,j,k,             &
            max_a_x,max_a_y)
    END DO
    !$OMP END PARALLEL DO

    ! Combine the SAME reconstructed face speeds, without another backend
    ! pass. Keep the legacy operation order exactly in the default branch.
    ! This is a spatial suggestion, not a general IMEX stability guarantee;
    ! the main loop subsequently caps accepted steps at events and growth.
    cfl_spatial = cfl
    cfl_spatial = MIN(cfl_spatial,0.25_wp)
    SELECT CASE (TRIM(cfl_policy))
    CASE ('directional_min')
       IF ( max_a_x .GT. 0.0_wp ) dt = MIN(dt,cfl_spatial*dx/max_a_x)
       IF ( max_a_y .GT. 0.0_wp ) dt = MIN(dt,cfl_spatial*dy/max_a_y)
    CASE ('combined')
       ! Use the legacy division order in 1D: mathematical equivalence alone
       ! would allow a one-ulp difference from cfl/(lambda/dx).
       IF ( comp_cells_y .EQ. 1 ) THEN
          IF ( comp_cells_x .GT. 1 .AND. max_a_x .GT. 0.0_wp ) dt = MIN(dt,cfl_spatial*dx/max_a_x)
       ELSEIF ( comp_cells_x .EQ. 1 ) THEN
          IF ( max_a_y .GT. 0.0_wp ) dt = MIN(dt,cfl_spatial*dy/max_a_y)
       ELSE
          propagation_rate = max_a_x/dx + max_a_y/dy
          IF ( propagation_rate .GT. 0.0_wp ) dt = MIN(dt,cfl_spatial/propagation_rate)
       END IF
    CASE DEFAULT
       ! Also protect callers that configure the module without input I/O.
       CALL fatal_error('CFL_POLICY must be directional_min or combined')
    END SELECT

  END SUBROUTINE compute_spatial_timestep

  !> \brief Include all four reconstructed face speeds for one active cell.
  !> \details The MAX operand order is kept identical to the original CFL loop.
  !> \param[in] speeds Reconstructed face-speed arrays.
  !> \param[in] j X index of the active cell.
  !> \param[in] k Y index of the active cell.
  !> \param[in,out] max_a_x Maximum outward wave speed in x.
  !> \param[in,out] max_a_y Maximum outward wave speed in y.
  SUBROUTINE accumulate_directional_face_speeds(speeds,j,k,max_a_x,max_a_y)
    CLASS(hyperbolic_workspace_type), INTENT(IN) :: speeds
    INTEGER, INTENT(IN) :: j,k
    REAL(wp), INTENT(INOUT) :: max_a_x,max_a_y

    max_a_x = MAX(max_a_x,                                                 &
         speeds%a_interface_xPos(j,k),                                     &
         -speeds%a_interface_xNeg(j,k),                                    &
         speeds%a_interface_xPos(j+1,k),                                   &
         -speeds%a_interface_xNeg(j+1,k))

    max_a_y = MAX(max_a_y,                                                 &
         speeds%a_interface_yPos(j,k),                                     &
         -speeds%a_interface_yNeg(j,k),                                    &
         speeds%a_interface_yPos(j,k+1),                                   &
         -speeds%a_interface_yNeg(j,k+1))
  END SUBROUTINE accumulate_directional_face_speeds

END MODULE spatial_operator_2d
