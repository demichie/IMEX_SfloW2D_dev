!********************************************************************************
!> \brief Model-facing explicit spatial operator
!>
!> This module is the boundary between time integration and the numerical
!> discretization of spatial terms.  It owns the active HP-PCCU backend and
!> exposes only complete spatial-term and CFL evaluations.  Time integration
!> therefore does not depend on flux, path-contribution, or wave-speed storage.
!********************************************************************************
MODULE spatial_operator_2d

  USE parameters_2d, ONLY : wp, n_eqns, n_vars, max_dt, cfl
  USE constitutive_parameters_2d, ONLY : T_ambient
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y, dx, dy
  USE state_conversion_2d, ONLY : qc_to_qp
  USE domain_2d, ONLY : domain_type
  USE reconstruction_2d, ONLY : reconstruction_workspace_type
  USE hyperbolic_2d, ONLY : hyperbolic_workspace_type

  IMPLICIT NONE

  PRIVATE

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

  SUBROUTINE initialize_spatial_operator( this )

    CLASS(spatial_operator_type), INTENT(INOUT) :: this

    CALL this%reconstruction%initialize
    CALL this%hp_pccu%initialize

  END SUBROUTINE initialize_spatial_operator

  SUBROUTINE finalize_spatial_operator( this )

    CLASS(spatial_operator_type), INTENT(INOUT) :: this

    CALL this%hp_pccu%finalize
    CALL this%reconstruction%finalize

  END SUBROUTINE finalize_spatial_operator

  !******************************************************************************
  !> \brief Evaluate the complete explicit spatial term on active cells
  !>
  !> The returned term includes conservative interface transport and the
  !> cell/interface path contributions assembled by HP-PCCU.  Its sign follows
  !> q_t + spatial_term = local_sources.
  !******************************************************************************
  SUBROUTINE evaluate_spatial_operator( this, qp, spatial_term, time, domain )

    CLASS(spatial_operator_type), INTENT(INOUT) :: this
    REAL(wp), INTENT(IN) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(OUT) :: spatial_term(n_eqns,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: time
    CLASS(domain_type), INTENT(IN) :: domain

    CALL this%hp_pccu%evaluate_terms( this%reconstruction, qp, spatial_term, &
         time, domain%solve_cells, domain%j_cent, domain%k_cent,              &
         domain%solve_interfaces_x, domain%j_stag_x, domain%k_stag_x,        &
         domain%solve_interfaces_y, domain%j_stag_y, domain%k_stag_y )

  END SUBROUTINE evaluate_spatial_operator

  !******************************************************************************
  !> \brief Compute the CFL timestep required by the active spatial operator
  !******************************************************************************
  SUBROUTINE compute_spatial_timestep( this, q, qp, time, dt, domain )

    CLASS(spatial_operator_type), INTENT(INOUT) :: this
    REAL(wp), INTENT(IN) :: q(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(INOUT) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: time
    REAL(wp), INTENT(OUT) :: dt
    CLASS(domain_type), INTENT(IN) :: domain

    INTEGER :: j, k, l
    REAL(wp) :: max_a_x, max_a_y
    REAL(wp) :: dynamic_pressure

    dt = max_dt

    IF ( cfl .EQ. -1.0_wp ) RETURN

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

    CALL this%reconstruction%reconstruct( qp, time, domain%solve_cells,      &
         domain%j_cent, domain%k_cent )

    CALL this%hp_pccu%evaluate_speeds( this%reconstruction,                  &
         domain%solve_interfaces_x, domain%j_stag_x, domain%k_stag_x,        &
         domain%solve_interfaces_y, domain%j_stag_y, domain%k_stag_y )

    max_a_x = 0.0_wp
    max_a_y = 0.0_wp

    !$OMP PARALLEL DO PRIVATE(j,k) REDUCTION(MAX:max_a_x,max_a_y)
    DO l = 1, domain%solve_cells
       j = domain%j_cent(l)
       k = domain%k_cent(l)

       max_a_x = MAX( max_a_x,                                               &
            MAXVAL(this%hp_pccu%a_interface_xPos(1:n_vars,j,k)),             &
            MAXVAL(-this%hp_pccu%a_interface_xNeg(1:n_vars,j,k)),            &
            MAXVAL(this%hp_pccu%a_interface_xPos(1:n_vars,j+1,k)),           &
            MAXVAL(-this%hp_pccu%a_interface_xNeg(1:n_vars,j+1,k)) )

       max_a_y = MAX( max_a_y,                                               &
            MAXVAL(this%hp_pccu%a_interface_yPos(1:n_vars,j,k)),             &
            MAXVAL(-this%hp_pccu%a_interface_yNeg(1:n_vars,j,k)),            &
            MAXVAL(this%hp_pccu%a_interface_yPos(1:n_vars,j,k+1)),           &
            MAXVAL(-this%hp_pccu%a_interface_yNeg(1:n_vars,j,k+1)) )
    END DO
    !$OMP END PARALLEL DO

    IF ( max_a_x .GT. 0.0_wp ) dt = MIN(dt,cfl*dx/max_a_x)
    IF ( max_a_y .GT. 0.0_wp ) dt = MIN(dt,cfl*dy/max_a_y)

  END SUBROUTINE compute_spatial_timestep

END MODULE spatial_operator_2d
