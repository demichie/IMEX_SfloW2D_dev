!********************************************************************************
!> \brief Reconstruction of cell states at computational interfaces
!>
!> This module owns the reconstructed conservative/physical interface states
!> together with free-surface traces and reusable HP scratch storage.
!>
!> West/east/south/north arrays are cell-owned traces; left/right/bottom/top arrays are
!> face-oriented traces. HP evaluation retains separate candidate and commit phases to avoid reading
!> partially updated neighbours.
!********************************************************************************

MODULE reconstruction_2d

  USE parameters_2d, ONLY : wp, n_vars, kt_trace_reconstruction
  USE parameters_2d, ONLY : lateral_source_flag, radial_source_flag
  USE parameters_2d, ONLY : bcW, bcE, bcS, bcN
  USE parameters_2d, ONLY : idx_u, idx_v

  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE geometry_2d, ONLY : comp_interfaces_x, comp_interfaces_y
  USE geometry_2d, ONLY : B_cent, B_face_x, B_face_y
  USE geometry_2d, ONLY : source_cell
  USE geometry_2d, ONLY : one_by_dx, one_by_dy, dx, dy
  USE geometry_2d, ONLY : limit

   USE hp_reconstruction_2d, ONLY : hp_dry_tolerance
   USE hp_reconstruction_2d, ONLY : hp_dynamic_residual_threshold

  IMPLICIT NONE

  PRIVATE

  ! Smooth transition from eta_x to h_x+B_x over this directional Froude band.
  REAL(wp), PARAMETER :: fr_blend_lower = 0.8_wp
  REAL(wp), PARAMETER :: fr_blend_upper = 1.2_wp

  !> \brief Persistent cell/face traces and thread-local storage for HP reconstruction.
  !> \details Arrays are indexed by component, x and y; hp_scratch uses its last
  !>          index to select the OpenMP thread. Allocation occurs at initialization,
  !>          not on each reconstruction call.
  TYPE, PUBLIC :: reconstruction_workspace_type
     REAL(wp), ALLOCATABLE :: q_interfaceL(:,:,:)
     REAL(wp), ALLOCATABLE :: q_interfaceR(:,:,:)
     REAL(wp), ALLOCATABLE :: q_interfaceB(:,:,:)
     REAL(wp), ALLOCATABLE :: q_interfaceT(:,:,:)

     REAL(wp), ALLOCATABLE :: qp_interfaceL(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_interfaceR(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_interfaceB(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_interfaceT(:,:,:)

     ! Direct h_x-based states for the advective part of the PCCU flux.
     ! The existing interface states retain the eta-based reconstruction and
     ! are used for central-upwind dissipation and hydrostatic path terms.
     REAL(wp), ALLOCATABLE :: q_adv_interfaceL(:,:,:)
     REAL(wp), ALLOCATABLE :: q_adv_interfaceR(:,:,:)
     REAL(wp), ALLOCATABLE :: q_adv_interfaceB(:,:,:)
     REAL(wp), ALLOCATABLE :: q_adv_interfaceT(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_adv_interfaceL(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_adv_interfaceR(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_adv_interfaceB(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_adv_interfaceT(:,:,:)

   ! Direct h_x-based cell-side candidates retained for the advective flux
   ! while the eta-based HP states are used by dissipation and path terms.
     REAL(wp), ALLOCATABLE :: qp_cellW(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_cellE(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_cellS(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_cellN(:,:,:)

     ! Final free-surface traces.  Both cell-side and face-oriented forms are
     ! stored because the following PCCU path discretization needs both views.
     REAL(wp), ALLOCATABLE :: eta_cellW(:,:), eta_cellE(:,:)
     REAL(wp), ALLOCATABLE :: eta_cellS(:,:), eta_cellN(:,:)
     REAL(wp), ALLOCATABLE :: eta_interfaceL(:,:), eta_interfaceR(:,:)
     REAL(wp), ALLOCATABLE :: eta_interfaceB(:,:), eta_interfaceT(:,:)
     REAL(wp), ALLOCATABLE :: h_minmod_interfaceL(:,:), h_minmod_interfaceR(:,:)
     REAL(wp), ALLOCATABLE :: h_minmod_interfaceB(:,:), h_minmod_interfaceT(:,:)
     REAL(wp), ALLOCATABLE :: h_froude_interfaceL(:,:), h_froude_interfaceR(:,:)
     REAL(wp), ALLOCATABLE :: h_froude_interfaceB(:,:), h_froude_interfaceT(:,:)

       LOGICAL, ALLOCATABLE :: hp_eta_mask(:,:)
       INTEGER, ALLOCATABLE :: hp_eta_j(:), hp_eta_k(:)
       INTEGER :: hp_eta_cells
       ! Per direction: h-/h+, normal momentum-/+, then tangential momentum-/+.
        
   CONTAINS
     PROCEDURE :: initialize => initialize_reconstruction
     PROCEDURE :: finalize => finalize_reconstruction
     PROCEDURE :: reconstruct => reconstruction
     PROCEDURE :: reconstruct_kt_eta => reconstruct_kt_eta_traces
  END TYPE reconstruction_workspace_type

CONTAINS

  !> \brief Allocate interface traces, HP candidates and per-thread scratch storage.
  !>
  !> \param[in,out] this Persistent reconstructed traces, HP candidates and per-thread scratch.

  SUBROUTINE initialize_reconstruction( this )

    CLASS(reconstruction_workspace_type), INTENT(INOUT) :: this

    ALLOCATE( this%q_interfaceL( n_vars, comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( this%q_interfaceR( n_vars, comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( this%q_interfaceB( n_vars, comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( this%q_interfaceT( n_vars, comp_cells_x, comp_interfaces_y ) )

    ALLOCATE( this%qp_interfaceL( n_vars+2, comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( this%qp_interfaceR( n_vars+2, comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( this%qp_interfaceB( n_vars+2, comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( this%qp_interfaceT( n_vars+2, comp_cells_x, comp_interfaces_y ) )

    ALLOCATE( this%q_adv_interfaceL( n_vars, comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( this%q_adv_interfaceR( n_vars, comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( this%q_adv_interfaceB( n_vars, comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( this%q_adv_interfaceT( n_vars, comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( this%qp_adv_interfaceL( n_vars+2, comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( this%qp_adv_interfaceR( n_vars+2, comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( this%qp_adv_interfaceB( n_vars+2, comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( this%qp_adv_interfaceT( n_vars+2, comp_cells_x, comp_interfaces_y ) )

    ALLOCATE( this%qp_cellW( n_vars+2, comp_cells_x, comp_cells_y ) )
    ALLOCATE( this%qp_cellE( n_vars+2, comp_cells_x, comp_cells_y ) )
    ALLOCATE( this%qp_cellS( n_vars+2, comp_cells_x, comp_cells_y ) )
    ALLOCATE( this%qp_cellN( n_vars+2, comp_cells_x, comp_cells_y ) )

    ALLOCATE( this%eta_cellW( comp_cells_x, comp_cells_y ) )
    ALLOCATE( this%eta_cellE( comp_cells_x, comp_cells_y ) )
    ALLOCATE( this%eta_cellS( comp_cells_x, comp_cells_y ) )
    ALLOCATE( this%eta_cellN( comp_cells_x, comp_cells_y ) )
    ALLOCATE( this%eta_interfaceL( comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( this%eta_interfaceR( comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( this%eta_interfaceB( comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( this%eta_interfaceT( comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( this%h_minmod_interfaceL( comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( this%h_minmod_interfaceR( comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( this%h_minmod_interfaceB( comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( this%h_minmod_interfaceT( comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( this%h_froude_interfaceL( comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( this%h_froude_interfaceR( comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( this%h_froude_interfaceB( comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( this%h_froude_interfaceT( comp_cells_x, comp_interfaces_y ) )

      ALLOCATE( this%hp_eta_mask(comp_cells_x,comp_cells_y) )
      ALLOCATE( this%hp_eta_j(comp_cells_x*comp_cells_y) )
      ALLOCATE( this%hp_eta_k(comp_cells_x*comp_cells_y) )
      this%hp_eta_mask = .FALSE.
      this%hp_eta_cells = 0


  END SUBROUTINE initialize_reconstruction

  !> \brief Release all reconstructed-state and HP work arrays.
  !>
  !> \param[in,out] this Persistent reconstructed traces, HP candidates and per-thread scratch.

  SUBROUTINE finalize_reconstruction( this )

    CLASS(reconstruction_workspace_type), INTENT(INOUT) :: this

    DEALLOCATE( this%q_interfaceL )
    DEALLOCATE( this%q_interfaceR )
    DEALLOCATE( this%q_interfaceB )
    DEALLOCATE( this%q_interfaceT )

    DEALLOCATE( this%qp_interfaceL )
    DEALLOCATE( this%qp_interfaceR )
    DEALLOCATE( this%qp_interfaceB )
    DEALLOCATE( this%qp_interfaceT )

    DEALLOCATE( this%q_adv_interfaceL )
    DEALLOCATE( this%q_adv_interfaceR )
    DEALLOCATE( this%q_adv_interfaceB )
    DEALLOCATE( this%q_adv_interfaceT )
    DEALLOCATE( this%qp_adv_interfaceL )
    DEALLOCATE( this%qp_adv_interfaceR )
    DEALLOCATE( this%qp_adv_interfaceB )
    DEALLOCATE( this%qp_adv_interfaceT )

    DEALLOCATE( this%qp_cellW )
    DEALLOCATE( this%qp_cellE )
    DEALLOCATE( this%qp_cellS )
    DEALLOCATE( this%qp_cellN )

    DEALLOCATE( this%eta_cellW )
    DEALLOCATE( this%eta_cellE )
    DEALLOCATE( this%eta_cellS )
    DEALLOCATE( this%eta_cellN )
    DEALLOCATE( this%eta_interfaceL )
    DEALLOCATE( this%eta_interfaceR )
    DEALLOCATE( this%eta_interfaceB )
    DEALLOCATE( this%eta_interfaceT )
    DEALLOCATE( this%h_minmod_interfaceL )
    DEALLOCATE( this%h_minmod_interfaceR )
    DEALLOCATE( this%h_minmod_interfaceB )
    DEALLOCATE( this%h_minmod_interfaceT )
    DEALLOCATE( this%h_froude_interfaceL )
    DEALLOCATE( this%h_froude_interfaceR )
    DEALLOCATE( this%h_froude_interfaceB )
    DEALLOCATE( this%h_froude_interfaceT )

      DEALLOCATE( this%hp_eta_mask )
      DEALLOCATE( this%hp_eta_j )
      DEALLOCATE( this%hp_eta_k )
    

  END SUBROUTINE finalize_reconstruction

   !> \brief Build a deduplicated eta workset containing active cells and immediate neighbours.
   !>
   !> \param[in,out] this Persistent reconstructed traces, HP candidates and per-thread scratch.
   !> \param[in] solve_cells Number of entries in the active-cell index lists.
   !> \param[in] j_cent X indices of active cells; only the first solve_cells entries are used.
   !> \param[in] k_cent Y indices of active cells; only the first solve_cells entries are used.
   !>
   !> \note Writes only the parent workspace eta mask/list; immediate neighbours provide the
   !>       additional continuity stencil needed by active cells.

   SUBROUTINE build_hp_workset(this, solve_cells, j_cent, k_cent)

      CLASS(reconstruction_workspace_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN) :: solve_cells
      INTEGER, INTENT(IN) :: j_cent(:), k_cent(:)
      INTEGER :: l, j, k

      DO l = 1, this%hp_eta_cells
          this%hp_eta_mask(this%hp_eta_j(l),this%hp_eta_k(l)) = .FALSE.
      END DO
      this%hp_eta_cells = 0

      DO l = 1, solve_cells
          j = j_cent(l)
          k = k_cent(l)
          CALL append_eta_cell(j,k)
          IF ( j .GT. 1 ) CALL append_eta_cell(j-1,k)
          IF ( j .LT. comp_cells_x ) CALL append_eta_cell(j+1,k)
          IF ( k .GT. 1 ) CALL append_eta_cell(j,k-1)
          IF ( k .LT. comp_cells_y ) CALL append_eta_cell(j,k+1)
      END DO

   CONTAINS

      !> \brief Append one cell to the parent eta workset only if it is not already present.
      !>
      !> \param[in] j_cell X index of the cell appended to the parent eta workset.
      !> \param[in] k_cell Y index of the cell appended to the parent eta workset.
      !>
      !> \note Updates the host-associated this workspace; the active-cell workset is not modified.

      SUBROUTINE append_eta_cell(j_cell,k_cell)

         INTEGER, INTENT(IN) :: j_cell, k_cell

         IF ( this%hp_eta_mask(j_cell,k_cell) ) RETURN

         this%hp_eta_mask(j_cell,k_cell) = .TRUE.
         this%hp_eta_cells = this%hp_eta_cells + 1
         this%hp_eta_j(this%hp_eta_cells) = j_cell
         this%hp_eta_k(this%hp_eta_cells) = k_cell

      END SUBROUTINE append_eta_cell

   END SUBROUTINE build_hp_workset

   !> \brief Reconstruct physical traces, apply boundary conditions and finalize HP states.
   !>
   !> \param[in,out] this Persistent reconstructed traces, HP candidates and per-thread scratch.
   !> \param[in] qp_expl Cell-centered physical states at the explicit stage, including the required
   !>                    stencil halo.
   !> \param[in] t Current simulation or stage time [s].
   !> \param[in] solve_cells Number of entries in the active-cell index lists.
   !> \param[in] j_cent X indices of active cells; only the first solve_cells entries are used.
   !> \param[in] k_cent Y indices of active cells; only the first solve_cells entries are used.
   !> \param[in] use_kt_direct Use direct KT traces instead of the HP dry-state projection.
   !> \param[in] source_boundary_enabled Optional override for boundary-source traces.


   SUBROUTINE reconstruction( this, qp_expl, t, solve_cells, j_cent, k_cent, &
        use_kt_direct, source_boundary_enabled )

    ! External procedures
      USE state_conversion_2d, ONLY : qp_to_qp2
    USE source_terms_2d, ONLY : eval_boundary_injection_explicit
    USE parameters_2d, ONLY : limiter

    USE geometry_2d, ONLY : x_comp , x_stag , y_comp , y_stag , dx2 , dy2

    USE geometry_2d, ONLY : sourceW , sourceE , sourceN , sourceS
    USE geometry_2d, ONLY : sourceW_vect_x , sourceW_vect_y
    USE geometry_2d, ONLY : sourceE_vect_x , sourceE_vect_y
    USE geometry_2d, ONLY : sourceN_vect_x , sourceN_vect_y
    USE geometry_2d, ONLY : sourceS_vect_x , sourceS_vect_y

    USE parameters_2d, ONLY : reconstr_coeff

    IMPLICIT NONE

    CLASS(reconstruction_workspace_type), INTENT(INOUT) :: this
    REAL(wp), INTENT(IN) :: qp_expl(:,:,:)
    REAL(wp), INTENT(IN) :: t
    INTEGER, INTENT(IN) :: solve_cells
    INTEGER, INTENT(IN) :: j_cent(:)
    INTEGER, INTENT(IN) :: k_cent(:)
    LOGICAL, INTENT(IN) :: use_kt_direct
    LOGICAL, INTENT(IN), OPTIONAL :: source_boundary_enabled

    REAL(wp) :: qrecW(n_vars+2) !< recons var at the west edge of the cells
    REAL(wp) :: qrecE(n_vars+2) !< recons var at the east edge of the cells
    REAL(wp) :: qrecS(n_vars+2) !< recons var at the south edge of the cells
    REAL(wp) :: qrecN(n_vars+2) !< recons var at the north edge of the cells

    REAL(wp) :: source_bdry(n_vars+2)
    REAL(wp) :: qrec_prime_x(n_vars+2)      !< recons variables slope
    REAL(wp) :: qrec_prime_y(n_vars+2)      !< recons variables slope

    REAL(wp) :: qp2recW(3) , qp2recE(3)
    REAL(wp) :: qp2recS(3) , qp2recN(3)

    REAL(wp) :: qrec_stencil(3) !< recons variables stencil for the limiter
    REAL(wp) :: x_stencil(3)    !< grid stencil for the limiter
    REAL(wp) :: y_stencil(3)    !< grid stencil for the limiter

    INTEGER :: l,j,k            !< loop counters (cells)
    INTEGER :: i                !< loop counter (variables)

    REAL(wp) :: dq

    LOGICAL :: regular_interior
    LOGICAL :: use_boundary_source

    ! Keep the same constant candidates for inactive neighbours, but initialize
    ! only cells that can contribute to an active HP trace.
    use_boundary_source=lateral_source_flag
    IF (PRESENT(source_boundary_enabled)) &
         use_boundary_source=source_boundary_enabled
    CALL build_hp_workset(this, solve_cells, j_cent, k_cent)

    !$OMP PARALLEL DO PRIVATE(l,j,k)
    DO l = 1, this%hp_eta_cells
       j = this%hp_eta_j(l)
       k = this%hp_eta_k(l)
       this%qp_cellW(:,j,k) = qp_expl(:,j,k)
       IF (qp_expl(1,j,k) .LE. hp_dry_tolerance) THEN
          this%qp_cellW(idx_u:idx_v,j,k) = 0.0_wp
       END IF
       this%qp_cellE(:,j,k) = this%qp_cellW(:,j,k)
       this%qp_cellS(:,j,k) = this%qp_cellW(:,j,k)
       this%qp_cellN(:,j,k) = this%qp_cellW(:,j,k)
    END DO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO private(j,k,i,qrecW,qrecE,qrecS,qrecN,x_stencil,y_stencil,&
    !$OMP & qrec_stencil,qrec_prime_x,qrec_prime_y,qp2recW,qp2recE,qp2recS,     &
    !$OMP & qp2recN,source_bdry,dq,regular_interior)

    DO l = 1,solve_cells

       j = j_cent(l)
       k = k_cent(l)

       qrecW(1:n_vars+2) = this%qp_cellW(1:n_vars+2,j,k)
       qrecE(1:n_vars+2) = this%qp_cellE(1:n_vars+2,j,k)
       qrecS(1:n_vars+2) = this%qp_cellS(1:n_vars+2,j,k)
       qrecN(1:n_vars+2) = this%qp_cellN(1:n_vars+2,j,k)

       x_stencil(2) = x_comp(j)
       y_stencil(2) = y_comp(k)

       ! Default source-side ghost state. For the conservative radial source,
       ! the ring is a wall and mass/momentum are injected through eval_expl_terms.
       ! The existing lateral boundary source keeps its Dirichlet treatment.
       source_bdry(1:n_vars+2) = qp_expl(1:n_vars+2,j,k)

       IF ( use_boundary_source .AND. lateral_source_flag .AND.             &
            ( source_cell(j,k) .EQ. 2 ) ) THEN

          IF ( sourceE(j,k) ) THEN

             CALL eval_boundary_injection_explicit( t, sourceE_vect_x(j,k), sourceE_vect_y(j,k), &
                  source_bdry )

          ELSEIF ( sourceW(j,k) ) THEN

             CALL eval_boundary_injection_explicit( t, sourceW_vect_x(j,k), sourceW_vect_y(j,k), &
                  source_bdry )

          ELSEIF ( sourceS(j,k) ) THEN

             CALL eval_boundary_injection_explicit( t, sourceS_vect_x(j,k), sourceS_vect_y(j,k), &
                  source_bdry )

          ELSEIF ( sourceN(j,k) ) THEN

             CALL eval_boundary_injection_explicit( t, sourceN_vect_x(j,k), sourceN_vect_y(j,k), &
                  source_bdry )

          END IF

       END IF

       regular_interior = ( j .GT. 1 ) .AND. ( j .LT. comp_cells_x ) .AND.     &
            ( k .GT. 1 ) .AND. ( k .LT. comp_cells_y ) .AND.                   &
            ( source_cell(j,k) .NE. 2 )

       IF ( regular_interior ) THEN

          x_stencil(1) = x_comp(j-1)
          x_stencil(3) = x_comp(j+1)
          y_stencil(1) = y_comp(k-1)
          y_stencil(3) = y_comp(k+1)

          fast_vars_loop:DO i=1,n_vars+2

             qrec_stencil(2) = reconstruction_cell_value(qp_expl,i,j,k)

             qrec_stencil(1) = reconstruction_cell_value(qp_expl,i,j-1,k)
             qrec_stencil(3) = reconstruction_cell_value(qp_expl,i,j+1,k)
             CALL limit( qrec_stencil , x_stencil , limiter(i) ,               &
                  qrec_prime_x(i) )

             dq = reconstr_coeff * dx2 * qrec_prime_x(i)
             qrecW(i) = qrec_stencil(2) - dq
             qrecE(i) = qrec_stencil(2) + dq

             qrec_stencil(1) = reconstruction_cell_value(qp_expl,i,j,k-1)
             qrec_stencil(3) = reconstruction_cell_value(qp_expl,i,j,k+1)
             CALL limit( qrec_stencil , y_stencil , limiter(i) ,               &
                  qrec_prime_y(i) )

             dq = reconstr_coeff * dy2 * qrec_prime_y(i)
             qrecS(i) = qrec_stencil(2) - dq
             qrecN(i) = qrec_stencil(2) + dq

          END DO fast_vars_loop

       ELSE

       vars_loop:DO i=1,n_vars

          qrec_stencil(2) = reconstruction_cell_value(qp_expl,i,j,k)

          ! x direction
          check_comp_cells_x:IF ( comp_cells_x .GT. 1 ) THEN

             ! west boundary
             check_x_boundary:IF ( j .EQ. 1 ) THEN

                x_stencil(1) = x_stag(1)
                x_stencil(3) = x_comp(j+1)

                IF ( source_cell(j,k).EQ.2 ) THEN

                   CALL boundary_reconstruction_slope(                         &
                        qrec_stencil(2), reconstruction_cell_value(qp_expl,i,j+1,k), &
                        x_stag(1), x_comp(j), x_comp(j+1), 0, source_bdry(i),   &
                        one_by_dx,                                              &
                        limiter(i), .TRUE., qrec_prime_x(i))

                ELSE

                   CALL boundary_reconstruction_slope(                         &
                        qrec_stencil(2), reconstruction_cell_value(qp_expl,i,j+1,k), &
                        x_stag(1), x_comp(j), x_comp(j+1), bcW(i)%flag,         &
                        bcW(i)%value, one_by_dx, limiter(i), .TRUE.,            &
                        qrec_prime_x(i))

                END IF

                !east boundary
             ELSEIF ( j .EQ. comp_cells_x ) THEN

                x_stencil(3) = x_stag(comp_interfaces_x)
                x_stencil(1) = x_comp(j-1)

                IF ( source_cell(j,k).EQ.2 ) THEN

                   CALL boundary_reconstruction_slope(                         &
                        qrec_stencil(2), reconstruction_cell_value(qp_expl,i,j-1,k), &
                        x_stag(comp_interfaces_x), x_comp(j), x_comp(j-1), 0,  &
                        source_bdry(i), one_by_dx, limiter(i), .FALSE.,         &
                        qrec_prime_x(i))

                ELSE

                   CALL boundary_reconstruction_slope(                         &
                        qrec_stencil(2), reconstruction_cell_value(qp_expl,i,j-1,k), &
                        x_stag(comp_interfaces_x), x_comp(j), x_comp(j-1),     &
                        bcE(i)%flag, bcE(i)%value, one_by_dx, limiter(i), .FALSE., &
                        qrec_prime_x(i))

                END IF

             ELSE

                ! internal x cells

                x_stencil(1) = x_comp(j-1)
                x_stencil(3) = x_comp(j+1)

                qrec_stencil(1) = reconstruction_cell_value(qp_expl,i,j-1,k)
                qrec_stencil(3) = reconstruction_cell_value(qp_expl,i,j+1,k)

                ! correction for radial source inlet x-interfaces values
                ! used for the linear reconstruction
                IF ( radial_source_flag .AND. ( source_cell(j,k).EQ.2 ) ) THEN

                   IF ( sourceE(j,k) ) THEN

                      x_stencil(3) = x_stag(j+1)
                      qrec_stencil(3) = source_bdry(i)

                   ELSEIF ( sourceW(j,k) ) THEN

                      x_stencil(1) = x_stag(j)
                      qrec_stencil(1) = source_bdry(i)

                   END IF

                END IF

                CALL limit( qrec_stencil , x_stencil , limiter(i) ,             &
                     qrec_prime_x(i) )

             ENDIF check_x_boundary

             dq = reconstr_coeff* dx2 * qrec_prime_x(i)

             qrecW(i) = qrec_stencil(2) - dq
             qrecE(i) = qrec_stencil(2) + dq


             IF ( j .EQ. 1 ) THEN

                IF ( source_cell(j,k).EQ.2 ) THEN

                   qrecW(i) = source_bdry(i)

                ELSE

                   ! Dirichelet boundary condition at the west of the domain
                   IF ( bcW(i)%flag .EQ. 0 ) THEN

                      qrecW(i) = bcW(i)%value

                   ELSE

                      IF ( i .EQ. 2 ) qrecW(i) = MIN( qrecW(i) , 0.0_wp )

                   END IF

                END IF

             END IF

             IF ( j .EQ. comp_cells_x ) THEN

                IF ( source_cell(j,k).EQ.2 ) THEN

                   qrecE(i) = source_bdry(i)

                ELSE

                   ! Dirichelet boundary condition at the east of the domain
                   IF ( bcE(i)%flag .EQ. 0 ) THEN

                      qrecE(i) = bcE(i)%value

                   ELSE

                      IF ( i .EQ. 2 ) qrecE(i) = MAX( qrecE(i) , 0.0_wp )

                   END IF

                END IF

             END IF

          END IF check_comp_cells_x

          ! y-direction
          check_comp_cells_y:IF ( comp_cells_y .GT. 1 ) THEN

             ! South boundary
             check_y_boundary:IF ( k .EQ. 1 ) THEN

                y_stencil(1) = y_stag(1)
                y_stencil(3) = y_comp(k+1)

                CALL boundary_reconstruction_slope(                            &
                     qrec_stencil(2), reconstruction_cell_value(qp_expl,i,j,k+1), &
                     y_stag(1), y_comp(k), y_comp(k+1), bcS(i)%flag,            &
                     bcS(i)%value,                                              &
                     one_by_dy, limiter(i), .TRUE., qrec_prime_y(i))

                ! North boundary
             ELSEIF ( k .EQ. comp_cells_y ) THEN

                   y_stencil(1) = y_comp(k-1)
                   y_stencil(3) = y_stag(comp_interfaces_y)

                   CALL boundary_reconstruction_slope(                            &
                     qrec_stencil(2), reconstruction_cell_value(qp_expl,i,j,k-1), &
                     y_stag(comp_interfaces_y), y_comp(k), y_comp(k-1),        &
                     bcN(i)%flag, bcN(i)%value, one_by_dy, limiter(i), .FALSE., &
                     qrec_prime_y(i))

             ELSE

                ! Internal y cells

                y_stencil(1) = y_comp(k-1)
                y_stencil(3) = y_comp(k+1)

                qrec_stencil(1) = reconstruction_cell_value(qp_expl,i,j,k-1)
                qrec_stencil(3) = reconstruction_cell_value(qp_expl,i,j,k+1)

                ! correction for radial source inlet y-interfaces
                ! used for the linear reconstruction
                IF ( radial_source_flag .AND. ( source_cell(j,k).EQ.2 ) ) THEN

                   IF ( sourceS(j,k) ) THEN

                      y_stencil(1) = y_stag(k)
                      qrec_stencil(1) = source_bdry(i)

                   ELSEIF ( sourceN(j,k) ) THEN

                      y_stencil(3) = y_stag(k+1)
                      qrec_stencil(3) = source_bdry(i)

                   END IF

                END IF

                CALL limit( qrec_stencil , y_stencil , limiter(i) ,             &
                     qrec_prime_y(i) )

             ENDIF check_y_boundary

             dq = reconstr_coeff * dy2 * qrec_prime_y(i)

             qrecS(i) = qrec_stencil(2) - dq
             qrecN(i) = qrec_stencil(2) + dq

             IF ( k .EQ. 1 ) THEN

                IF ( source_cell(j,k) .EQ. 2 ) THEN

                   qrecS(i) = source_bdry(i)

                ELSE

                   ! Dirichlet boundary condition at the south of the domain
                   IF ( bcS(i)%flag .EQ. 0 ) THEN
                      qrecS(i) = bcS(i)%value
                   ELSE
                      IF ( i .EQ. 3 ) qrecS(i) = MIN( qrecS(i) , 0.0_wp )
                   END IF

                END IF

             END IF

             IF ( k .EQ. comp_cells_y ) THEN

                IF ( source_cell(j,k) .EQ. 2 ) THEN

                   qrecN(i) = source_bdry(i)

                ELSE

                   ! Dirichlet boundary condition at the north of the domain
                   IF ( bcN(i)%flag .EQ. 0 ) THEN
                      qrecN(i) = bcN(i)%value
                   ELSE
                      IF ( i .EQ. 3 ) qrecN(i) = MAX( qrecN(i) , 0.0_wp )
                   END IF

                END IF

             END IF

          ENDIF check_comp_cells_y

       ENDDO vars_loop

       add_vars_loop:DO i=n_vars+1,n_vars+2
          ! reconstruction on u and v with same limiters of hu,hv

          ! x direction
          check_comp_cells_x2:IF ( comp_cells_x .GT. 1 ) THEN

             qrec_stencil(2) = reconstruction_cell_value(qp_expl,i,j,k)

             IF ( j .EQ. 1 ) THEN

                CALL qp_to_qp2( qrecW(1:n_vars+2) , B_cent(j,k) , qp2recW )
                qrec_stencil(1) = qp2recW(i-n_vars+1)
                qrec_stencil(3) = reconstruction_cell_value(qp_expl,i,j+1,k)

             ELSEIF ( j .EQ. comp_cells_x ) THEN

                CALL qp_to_qp2( qrecE(1:n_vars+2) , B_cent(j,k) , qp2recE )
                qrec_stencil(1) = reconstruction_cell_value(qp_expl,i,j-1,k)
                qrec_stencil(3) = qp2recE(i-n_vars+1)

             ELSE

                qrec_stencil(1) = reconstruction_cell_value(qp_expl,i,j-1,k)
                qrec_stencil(3) = reconstruction_cell_value(qp_expl,i,j+1,k)

                ! correction for radial source inlet x-interfaces values
                ! used for the linear reconstruction
                IF ( radial_source_flag .AND. ( source_cell(j,k).EQ.2 ) ) THEN

                   IF ( sourceE(j,k) ) THEN

                      x_stencil(3) = x_stag(j+1)
                      qrec_stencil(3) = source_bdry(i)

                   ELSEIF ( sourceW(j,k) ) THEN

                      x_stencil(1) = x_stag(j)
                      qrec_stencil(1) = source_bdry(i)

                   END IF

                END IF

             END IF

             CALL limit( qrec_stencil , x_stencil , limiter(i) ,                &
                  qrec_prime_x(i) )

             dq = reconstr_coeff*dx2*qrec_prime_x(i)

             qrecW(i) = qrec_stencil(2) - dq
             qrecE(i) = qrec_stencil(2) + dq

             IF ( j .EQ. 1 ) THEN

                CALL qp_to_qp2( qrecW(1:n_vars+2) , B_cent(j,k) , qp2recW )
                qrecW(i) = qp2recW(i-n_vars+1)

             ELSEIF ( j .EQ. comp_cells_x ) THEN

                CALL qp_to_qp2( qrecE(1:n_vars+2) , B_cent(j,k) , qp2recE )
                qrecE(i) = qp2recE(i-n_vars+1)

             ELSE

                ! correction for radial source inlet x-interfaces:
                ! the physical variables at the x-interfaces qrecW or
                ! qrecE are computed from the radial inlet values
                IF ( radial_source_flag .AND. ( source_cell(j,k).EQ.2 ) ) THEN

                   IF ( sourceE(j,k) ) THEN

                      qrecE(1:n_vars+2) = source_bdry(1:n_vars+2)

                   ELSEIF ( sourceW(j,k) ) THEN

                      qrecW(1:n_vars+2) = source_bdry(1:n_vars+2)

                   END IF

                END IF

             END IF

          END IF check_comp_cells_x2

          ! y-direction
          check_comp_cells_y2:IF ( comp_cells_y .GT. 1 ) THEN

             qrec_stencil(2) = reconstruction_cell_value(qp_expl,i,j,k)

             IF ( k .EQ. 1 ) THEN

                CALL qp_to_qp2( qrecS(1:n_vars+2) , B_cent(j,k) , qp2recS )
                qrec_stencil(1) = qp2recS(i-n_vars+1)
                qrec_stencil(3) = reconstruction_cell_value(qp_expl,i,j,k+1)

             ELSEIF ( k .EQ. comp_cells_y ) THEN

                CALL qp_to_qp2( qrecN(1:n_vars+2) , B_cent(j,k) , qp2recN )
                qrec_stencil(1) = reconstruction_cell_value(qp_expl,i,j,k-1)
                qrec_stencil(3) = qp2recN(i-n_vars+1)

             ELSE

                qrec_stencil(1) = reconstruction_cell_value(qp_expl,i,j,k-1)
                qrec_stencil(3) = reconstruction_cell_value(qp_expl,i,j,k+1)

                ! correction for radial source inlet y-interfaces
                ! used for the linear reconstruction
                IF ( radial_source_flag .AND. ( source_cell(j,k).EQ.2 ) ) THEN

                   IF ( sourceS(j,k) ) THEN

                      y_stencil(1) = y_stag(k)
                      qrec_stencil(1) = source_bdry(i)

                   ELSEIF ( sourceN(j,k) ) THEN

                      y_stencil(3) = y_stag(k+1)
                      qrec_stencil(3) = source_bdry(i)

                   END IF

                END IF

             ENDIF

             CALL limit( qrec_stencil , y_stencil , limiter(i) ,                &
                  qrec_prime_y(i) )

             dq = reconstr_coeff*dy2*qrec_prime_y(i)

             qrecS(i) = qrec_stencil(2) - dq
             qrecN(i) = qrec_stencil(2) + dq


             IF ( k .EQ. 1 ) THEN

                CALL qp_to_qp2( qrecS(1:n_vars+2) , B_cent(j,k) , qp2recS )
                qrecS(i) = qp2recS(i-n_vars+1)

             ELSEIF ( k .EQ. comp_cells_y ) THEN

                CALL qp_to_qp2( qrecN(1:n_vars+2) , B_cent(j,k) , qp2recN )
                qrecN(i) = qp2recN(i-n_vars+1)

             ELSE

                ! correction for radial source inlet y-interfaces:
                ! the physical variables at the y-interfaces qrecS or
                ! qrecN are computed from the radial inlet values
                IF ( radial_source_flag .AND. ( source_cell(j,k) .EQ. 2 ) ) THEN

                   IF ( sourceS(j,k) ) THEN

                      qrecS(1:n_vars+2) = source_bdry(1:n_vars+2)

                   ELSEIF ( sourceN(j,k) ) THEN

                      qrecN(1:n_vars+2) = source_bdry(1:n_vars+2)

                   END IF

                END IF

             END IF

          ENDIF check_comp_cells_y2

       ENDDO add_vars_loop

       END IF

       this%qp_cellW(:,j,k) = qrecW
       this%qp_cellE(:,j,k) = qrecE
       this%qp_cellS(:,j,k) = qrecS
       this%qp_cellN(:,j,k) = qrecN

    END DO

    !$OMP END PARALLEL DO

    CALL apply_hp_reconstruction( this, qp_expl, solve_cells, j_cent, k_cent, &
         use_kt_direct )

    RETURN

  END SUBROUTINE reconstruction

  !******************************************************************************
  !> \brief Reconstruct the KP free surface for KT dissipation and path terms.
  !> \param[in,out] this Workspace containing the active HP-eta cell list and traces.
  !> \param[in] qp_expl Cell-centered physical states at the explicit stage.
  !******************************************************************************
  SUBROUTINE reconstruct_kt_eta_traces(this,qp_expl)
    USE parameters_2d, ONLY : limiter,reconstr_coeff
    USE geometry_2d, ONLY : x_comp,y_comp,dx2,dy2
    USE geometry_2d, ONLY : sourceW,sourceE,sourceS,sourceN
    USE constitutive_parameters_2d, ONLY : grav

    CLASS(reconstruction_workspace_type), INTENT(INOUT) :: this
    REAL(wp), INTENT(IN) :: qp_expl(:,:,:)
    REAL(wp) :: sx(3),sy(3),hx(3),hy(3),xx(3),yy(3)
    REAL(wp) :: slope_x,slope_y,slope_hx,slope_hy,eta_center
    REAL(wp) :: eta_w,eta_e,eta_s,eta_n,slope_eta2,slope_minmod
    REAL(wp) :: eta_w_mm,eta_e_mm,eta_s_mm,eta_n_mm
    REAL(wp) :: slope_eta2_x,slope_eta2_y
    REAL(wp) :: slope_blend_x,slope_blend_y,eta_blend_w,eta_blend_e
    REAL(wp) :: eta_blend_s,eta_blend_n,fr_weight_x,fr_weight_y
    INTEGER :: l,j,k

    !$OMP PARALLEL DO PRIVATE(l,j,k,sx,sy,hx,hy,xx,yy,slope_x,slope_y,      &
    !$OMP& slope_hx,slope_hy,eta_center,eta_w,eta_e,eta_s,eta_n,            &
    !$OMP& slope_eta2,slope_minmod,eta_w_mm,eta_e_mm,eta_s_mm,eta_n_mm)
    DO l=1,this%hp_eta_cells
       j=this%hp_eta_j(l)
       k=this%hp_eta_k(l)
       eta_center=qp_expl(1,j,k)+B_cent(j,k)
       slope_x=0.0_wp
       slope_y=0.0_wp
       slope_hx=0.0_wp
       slope_hy=0.0_wp
       IF ((j.GT.1).AND.(j.LT.comp_cells_x)) THEN
          sx=[qp_expl(1,j-1,k)+B_cent(j-1,k),eta_center, &
               qp_expl(1,j+1,k)+B_cent(j+1,k)]
          hx=[qp_expl(1,j-1,k),qp_expl(1,j,k),qp_expl(1,j+1,k)]
          xx=[x_comp(j-1),x_comp(j),x_comp(j+1)]
          CALL limit(sx,xx,limiter(1),slope_x)
          CALL limit(hx,xx,limiter(1),slope_hx)
       END IF
       IF ((k.GT.1).AND.(k.LT.comp_cells_y)) THEN
          sy=[qp_expl(1,j,k-1)+B_cent(j,k-1),eta_center, &
               qp_expl(1,j,k+1)+B_cent(j,k+1)]
          hy=[qp_expl(1,j,k-1),qp_expl(1,j,k),qp_expl(1,j,k+1)]
          yy=[y_comp(k-1),y_comp(k),y_comp(k+1)]
          CALL limit(sy,yy,limiter(1),slope_y)
          CALL limit(hy,yy,limiter(1),slope_hy)
       END IF
       eta_w=eta_center-reconstr_coeff*dx2*slope_x
       eta_e=eta_center+reconstr_coeff*dx2*slope_x
       eta_s=eta_center-reconstr_coeff*dy2*slope_y
       eta_n=eta_center+reconstr_coeff*dy2*slope_y

       ! K&P positivity correction for the eta_x traces used by the
       ! hydrostatic/path discretization.
       IF ( eta_e .LT. B_face_x(j+1,k) ) THEN
          eta_e=B_face_x(j+1,k)
          eta_w=2.0_wp*eta_center-eta_e
       ELSEIF ( eta_w .LT. B_face_x(j,k) ) THEN
          eta_w=B_face_x(j,k)
          eta_e=2.0_wp*eta_center-eta_w
       END IF
       IF ( eta_n .LT. B_face_y(j,k+1) ) THEN
          eta_n=B_face_y(j,k+1)
          eta_s=2.0_wp*eta_center-eta_n
       ELSEIF ( eta_s .LT. B_face_y(j,k) ) THEN
          eta_s=B_face_y(j,k)
          eta_n=2.0_wp*eta_center-eta_s
       END IF
       this%eta_cellW(j,k)=eta_w
       this%eta_cellE(j,k)=eta_e
       this%eta_cellS(j,k)=eta_s
       this%eta_cellN(j,k)=eta_n

       IF (kt_trace_reconstruction .EQ. 'minmod') THEN
          ! Alternative transport traces: minmod(eta_x, h_x+B_x), followed
          ! by the face-bed positivity correction. Hydrostatic eta stays above.
          slope_eta2=reconstr_coeff*slope_hx + &
               (B_face_x(j+1,k)-B_face_x(j,k))*one_by_dx
          slope_minmod=minmod_slope_pair(reconstr_coeff*slope_x,slope_eta2)
          eta_w_mm=eta_center-dx2*slope_minmod
          eta_e_mm=eta_center+dx2*slope_minmod
          IF ( eta_e_mm .LT. B_face_x(j+1,k) ) THEN
             eta_e_mm=B_face_x(j+1,k)
             eta_w_mm=2.0_wp*eta_center-eta_e_mm
          ELSEIF ( eta_w_mm .LT. B_face_x(j,k) ) THEN
             eta_w_mm=B_face_x(j,k)
             eta_e_mm=2.0_wp*eta_center-eta_w_mm
          END IF
          this%h_minmod_interfaceR(j,k)= &
               MAX(eta_w_mm-B_face_x(j,k),0.0_wp)
          this%h_minmod_interfaceL(j+1,k)= &
               MAX(eta_e_mm-B_face_x(j+1,k),0.0_wp)

          slope_eta2=reconstr_coeff*slope_hy + &
               (B_face_y(j,k+1)-B_face_y(j,k))*one_by_dy
          slope_minmod=minmod_slope_pair(reconstr_coeff*slope_y,slope_eta2)
          eta_s_mm=eta_center-dy2*slope_minmod
          eta_n_mm=eta_center+dy2*slope_minmod
          IF ( eta_n_mm .LT. B_face_y(j,k+1) ) THEN
             eta_n_mm=B_face_y(j,k+1)
             eta_s_mm=2.0_wp*eta_center-eta_n_mm
          ELSEIF ( eta_s_mm .LT. B_face_y(j,k) ) THEN
             eta_s_mm=B_face_y(j,k)
             eta_n_mm=2.0_wp*eta_center-eta_s_mm
          END IF
          this%h_minmod_interfaceT(j,k)= &
               MAX(eta_s_mm-B_face_y(j,k),0.0_wp)
          this%h_minmod_interfaceB(j,k+1)= &
               MAX(eta_n_mm-B_face_y(j,k+1),0.0_wp)
       END IF
    END DO
    !$OMP END PARALLEL DO

    IF (kt_trace_reconstruction .EQ. 'froude_blend') THEN
       ! Cell-local directional blend: eta_x below Fr=0.8, h_x+B_x above
       ! Fr=1.2, with a smoothstep transition in between. The resulting
       ! h-face traces are shared by the advective and dissipative fluxes;
       ! eta_cell* above remains the independent hydrostatic/path state.
       !$OMP PARALLEL DO PRIVATE(l,j,k,sx,sy,hx,hy,xx,yy,slope_x,slope_y,  &
       !$OMP& slope_hx,slope_hy,eta_center,slope_eta2_x,slope_eta2_y,     &
       !$OMP& slope_blend_x,slope_blend_y,eta_blend_w,eta_blend_e,        &
       !$OMP& eta_blend_s,eta_blend_n,fr_weight_x,fr_weight_y)
       DO l=1,this%hp_eta_cells
          j=this%hp_eta_j(l)
          k=this%hp_eta_k(l)
          eta_center=qp_expl(1,j,k)+B_cent(j,k)
          slope_x=0.0_wp
          slope_y=0.0_wp
          slope_hx=0.0_wp
          slope_hy=0.0_wp
          IF ((j.GT.1).AND.(j.LT.comp_cells_x)) THEN
             sx=[qp_expl(1,j-1,k)+B_cent(j-1,k),eta_center, &
                  qp_expl(1,j+1,k)+B_cent(j+1,k)]
             hx=[qp_expl(1,j-1,k),qp_expl(1,j,k),qp_expl(1,j+1,k)]
             xx=[x_comp(j-1),x_comp(j),x_comp(j+1)]
             CALL limit(sx,xx,limiter(1),slope_x)
             CALL limit(hx,xx,limiter(1),slope_hx)
          END IF
          IF ((k.GT.1).AND.(k.LT.comp_cells_y)) THEN
             sy=[qp_expl(1,j,k-1)+B_cent(j,k-1),eta_center, &
                  qp_expl(1,j,k+1)+B_cent(j,k+1)]
             hy=[qp_expl(1,j,k-1),qp_expl(1,j,k),qp_expl(1,j,k+1)]
             yy=[y_comp(k-1),y_comp(k),y_comp(k+1)]
             CALL limit(sy,yy,limiter(1),slope_y)
             CALL limit(hy,yy,limiter(1),slope_hy)
          END IF

          fr_weight_x=froude_h_blend_weight(qp_expl(1,j,k), &
               qp_expl(idx_u,j,k),grav)
          fr_weight_y=froude_h_blend_weight(qp_expl(1,j,k), &
               qp_expl(idx_v,j,k),grav)
          slope_eta2_x=reconstr_coeff*slope_hx + &
               (B_face_x(j+1,k)-B_face_x(j,k))*one_by_dx
          slope_eta2_y=reconstr_coeff*slope_hy + &
               (B_face_y(j,k+1)-B_face_y(j,k))*one_by_dy
          slope_blend_x=(1.0_wp-fr_weight_x)*reconstr_coeff*slope_x + &
               fr_weight_x*slope_eta2_x
          slope_blend_y=(1.0_wp-fr_weight_y)*reconstr_coeff*slope_y + &
               fr_weight_y*slope_eta2_y

          eta_blend_w=eta_center-dx2*slope_blend_x
          eta_blend_e=eta_center+dx2*slope_blend_x
          IF ( eta_blend_e .LT. B_face_x(j+1,k) ) THEN
             eta_blend_e=B_face_x(j+1,k)
             eta_blend_w=2.0_wp*eta_center-eta_blend_e
          ELSEIF ( eta_blend_w .LT. B_face_x(j,k) ) THEN
             eta_blend_w=B_face_x(j,k)
             eta_blend_e=2.0_wp*eta_center-eta_blend_w
          END IF
          this%h_froude_interfaceR(j,k)= &
               MAX(eta_blend_w-B_face_x(j,k),0.0_wp)
          this%h_froude_interfaceL(j+1,k)= &
               MAX(eta_blend_e-B_face_x(j+1,k),0.0_wp)

          eta_blend_s=eta_center-dy2*slope_blend_y
          eta_blend_n=eta_center+dy2*slope_blend_y
          IF ( eta_blend_n .LT. B_face_y(j,k+1) ) THEN
             eta_blend_n=B_face_y(j,k+1)
             eta_blend_s=2.0_wp*eta_center-eta_blend_n
          ELSEIF ( eta_blend_s .LT. B_face_y(j,k) ) THEN
             eta_blend_s=B_face_y(j,k)
             eta_blend_n=2.0_wp*eta_center-eta_blend_s
          END IF
          this%h_froude_interfaceT(j,k)= &
               MAX(eta_blend_s-B_face_y(j,k),0.0_wp)
          this%h_froude_interfaceB(j,k+1)= &
               MAX(eta_blend_n-B_face_y(j,k+1),0.0_wp)
       END DO
       !$OMP END PARALLEL DO

       DO k=1,comp_cells_y
          this%h_froude_interfaceL(1,k)=this%h_froude_interfaceR(1,k)
          this%h_froude_interfaceR(comp_interfaces_x,k)= &
               this%h_froude_interfaceL(comp_interfaces_x,k)
       END DO
       DO j=1,comp_cells_x
          this%h_froude_interfaceB(j,1)=this%h_froude_interfaceT(j,1)
          this%h_froude_interfaceT(j,comp_interfaces_y)= &
               this%h_froude_interfaceB(j,comp_interfaces_y)
       END DO
    END IF

    !$OMP PARALLEL DO PRIVATE(l,j,k)
    DO l=1,this%hp_eta_cells
       j=this%hp_eta_j(l)
       k=this%hp_eta_k(l)
       this%eta_interfaceR(j,k)=this%eta_cellW(j,k)
       this%eta_interfaceL(j+1,k)=this%eta_cellE(j,k)
       this%eta_interfaceT(j,k)=this%eta_cellS(j,k)
       this%eta_interfaceB(j,k+1)=this%eta_cellN(j,k)
       IF (j.EQ.1) THEN
          this%eta_interfaceL(1,k)=this%eta_interfaceR(1,k)
          IF ( kt_trace_reconstruction .EQ. 'minmod' ) &
               this%h_minmod_interfaceL(1,k)=this%h_minmod_interfaceR(1,k)
       END IF
       IF (j.EQ.comp_cells_x) THEN
          this%eta_interfaceR(comp_interfaces_x,k)= &
               this%eta_interfaceL(comp_interfaces_x,k)
          IF ( kt_trace_reconstruction .EQ. 'minmod' ) &
               this%h_minmod_interfaceR(comp_interfaces_x,k)= &
               this%h_minmod_interfaceL(comp_interfaces_x,k)
       END IF
       IF (k.EQ.1) THEN
          this%eta_interfaceB(j,1)=this%eta_interfaceT(j,1)
          IF ( kt_trace_reconstruction .EQ. 'minmod' ) &
               this%h_minmod_interfaceB(j,1)=this%h_minmod_interfaceT(j,1)
       END IF
       IF (k.EQ.comp_cells_y) THEN
          this%eta_interfaceT(j,comp_interfaces_y)= &
               this%eta_interfaceB(j,comp_interfaces_y)
          IF ( kt_trace_reconstruction .EQ. 'minmod' ) &
               this%h_minmod_interfaceT(j,comp_interfaces_y)= &
               this%h_minmod_interfaceB(j,comp_interfaces_y)
       END IF
    END DO
    !$OMP END PARALLEL DO

    DO k=1,comp_cells_y
       this%q_adv_interfaceL(:,1,k)=this%q_adv_interfaceR(:,1,k)
       this%qp_adv_interfaceL(:,1,k)=this%qp_adv_interfaceR(:,1,k)
       this%q_adv_interfaceR(:,comp_interfaces_x,k)= &
            this%q_adv_interfaceL(:,comp_interfaces_x,k)
       this%qp_adv_interfaceR(:,comp_interfaces_x,k)= &
            this%qp_adv_interfaceL(:,comp_interfaces_x,k)
    END DO
    DO j=1,comp_cells_x
       this%q_adv_interfaceB(:,j,1)=this%q_adv_interfaceT(:,j,1)
       this%qp_adv_interfaceB(:,j,1)=this%qp_adv_interfaceT(:,j,1)
       this%q_adv_interfaceT(:,j,comp_interfaces_y)= &
            this%q_adv_interfaceB(:,j,comp_interfaces_y)
       this%qp_adv_interfaceT(:,j,comp_interfaces_y)= &
            this%qp_adv_interfaceB(:,j,comp_interfaces_y)
    END DO
    ! reconstruct() has already reflected the primitive/conservative states
    ! at the internal radial-source wall. The eta pass above overwrites its
    ! scalar ghost traces: restore equal depths on both sides, including the
    ! selected transport reconstruction, before evaluating speeds and fluxes.
    ! Injection is volumetric; no mass may cross this excluded inner boundary.
    IF ( radial_source_flag ) THEN
       DO l=1,this%hp_eta_cells
          j=this%hp_eta_j(l)
          k=this%hp_eta_k(l)
          IF ( source_cell(j,k) .NE. 2 ) CYCLE
          IF ( sourceE(j,k) ) THEN
             this%eta_interfaceR(j+1,k)=this%eta_interfaceL(j+1,k)
             IF ( kt_trace_reconstruction .EQ. 'minmod' ) &
                  this%h_minmod_interfaceR(j+1,k)=this%h_minmod_interfaceL(j+1,k)
             IF ( kt_trace_reconstruction .EQ. 'froude_blend' ) &
                  this%h_froude_interfaceR(j+1,k)=this%h_froude_interfaceL(j+1,k)
          ELSEIF ( sourceW(j,k) ) THEN
             this%eta_interfaceL(j,k)=this%eta_interfaceR(j,k)
             IF ( kt_trace_reconstruction .EQ. 'minmod' ) &
                  this%h_minmod_interfaceL(j,k)=this%h_minmod_interfaceR(j,k)
             IF ( kt_trace_reconstruction .EQ. 'froude_blend' ) &
                  this%h_froude_interfaceL(j,k)=this%h_froude_interfaceR(j,k)
          END IF
          IF ( sourceN(j,k) ) THEN
             this%eta_interfaceT(j,k+1)=this%eta_interfaceB(j,k+1)
             IF ( kt_trace_reconstruction .EQ. 'minmod' ) &
                  this%h_minmod_interfaceT(j,k+1)=this%h_minmod_interfaceB(j,k+1)
             IF ( kt_trace_reconstruction .EQ. 'froude_blend' ) &
                  this%h_froude_interfaceT(j,k+1)=this%h_froude_interfaceB(j,k+1)
          ELSEIF ( sourceS(j,k) ) THEN
             this%eta_interfaceB(j,k)=this%eta_interfaceT(j,k)
             IF ( kt_trace_reconstruction .EQ. 'minmod' ) &
                  this%h_minmod_interfaceB(j,k)=this%h_minmod_interfaceT(j,k)
             IF ( kt_trace_reconstruction .EQ. 'froude_blend' ) &
                  this%h_froude_interfaceB(j,k)=this%h_froude_interfaceT(j,k)
          END IF
       END DO
    END IF
  END SUBROUTINE reconstruct_kt_eta_traces

  !******************************************************************************
  !> \brief Compute a smooth blend weight from the local Froude number.
  !> \param[in] h Flow depth.
  !> \param[in] u_normal Normal velocity.
  !> \param[in] gravity Effective gravitational acceleration.
  !> \return Blend weight between zero and one.
  !******************************************************************************
  PURE FUNCTION froude_h_blend_weight(h,u_normal,gravity) RESULT(weight)
    REAL(wp), INTENT(IN) :: h,u_normal,gravity
    REAL(wp) :: weight,fr,scaled

    weight=0.0_wp
    IF ( (h .LE. hp_dry_tolerance) .OR. (gravity .LE. 0.0_wp) ) RETURN

    fr=ABS(u_normal)/SQRT(gravity*h)
    scaled=MIN(1.0_wp,MAX(0.0_wp, &
         (fr-fr_blend_lower)/(fr_blend_upper-fr_blend_lower)))
    weight=scaled*scaled*(3.0_wp-2.0_wp*scaled)
  END FUNCTION froude_h_blend_weight

  !******************************************************************************
  !> \brief Apply the two-argument minmod limiter to a pair of slopes.
  !> \param[in] slope_a First slope candidate.
  !> \param[in] slope_b Second slope candidate.
  !> \return Same-sign slope with the smaller magnitude, or zero.
  !******************************************************************************
  PURE FUNCTION minmod_slope_pair(slope_a,slope_b) RESULT(slope)
    REAL(wp), INTENT(IN) :: slope_a,slope_b
    REAL(wp) :: slope
    slope=0.0_wp
    IF ( slope_a*slope_b .GT. 0.0_wp ) &
         slope=SIGN(MIN(ABS(slope_a),ABS(slope_b)),slope_a)
  END FUNCTION minmod_slope_pair

  !> \brief Apply one side's boundary condition to the linear reconstruction slope.
  !> \details Flag 0 limits a Dirichlet boundary stencil, flag 1 prescribes the
  !>          Neumann slope, and flag 2 uses the adjacent-cell secant.
  !> \param[in] center_value Value in the boundary-adjacent cell.
  !> \param[in] adjacent_value Value in the next cell inward from the boundary.
  !> \param[in] boundary_coordinate Coordinate of the physical boundary.
  !> \param[in] center_coordinate Coordinate of the boundary-adjacent cell center.
  !> \param[in] adjacent_coordinate Coordinate of the next cell inward.
  !> \param[in] boundary_flag Boundary mode (Dirichlet, Neumann or secant).
  !> \param[in] boundary_value Dirichlet value or Neumann derivative.
  !> \param[in] inverse_spacing Inverse uniform grid spacing for secant mode.
  !> \param[in] limiter_id Limiter selector.
  !> \param[in] lower_boundary True for west/south, false for east/north.
  !> \param[out] slope Limited or prescribed derivative.
  SUBROUTINE boundary_reconstruction_slope(center_value,adjacent_value,       &
       boundary_coordinate,center_coordinate,adjacent_coordinate,             &
       boundary_flag,boundary_value,inverse_spacing,limiter_id,               &
       lower_boundary,slope)
    REAL(wp), INTENT(IN) :: center_value,adjacent_value
    REAL(wp), INTENT(IN) :: boundary_coordinate,center_coordinate
    REAL(wp), INTENT(IN) :: adjacent_coordinate
    REAL(wp), INTENT(IN) :: boundary_value,inverse_spacing
    INTEGER, INTENT(IN) :: boundary_flag,limiter_id
    LOGICAL, INTENT(IN) :: lower_boundary
    REAL(wp), INTENT(OUT) :: slope
    REAL(wp) :: stencil_values(3),stencil_coordinates(3)

    IF (boundary_flag .EQ. 0) THEN
       stencil_values(2)=center_value
       stencil_coordinates(2)=center_coordinate
       IF (lower_boundary) THEN
          stencil_values(1)=boundary_value
          stencil_values(3)=adjacent_value
          stencil_coordinates(1)=boundary_coordinate
          stencil_coordinates(3)=adjacent_coordinate
       ELSE
          stencil_values(1)=adjacent_value
          stencil_values(3)=boundary_value
          stencil_coordinates(1)=adjacent_coordinate
          stencil_coordinates(3)=boundary_coordinate
       END IF
       CALL limit(stencil_values,stencil_coordinates,limiter_id,slope)
    ELSEIF (boundary_flag .EQ. 1) THEN
       slope=boundary_value
    ELSEIF (boundary_flag .EQ. 2) THEN
       IF (lower_boundary) THEN
          slope=(adjacent_value-center_value)*inverse_spacing
       ELSE
          slope=(center_value-adjacent_value)*inverse_spacing
       END IF
    END IF
  END SUBROUTINE boundary_reconstruction_slope

  !> \brief Read one centre value for a direct reconstruction stencil.
  !> \details KT uses the cell-centred velocities already obtained by
  !>          qc_to_qp, which applies the v2 mass-based desingularization.
  !>          Auxiliary u/v are set to zero below the HP dry-depth threshold
  !>          so a dry cell cannot supply a velocity to a neighbouring limiter.
  !>          Depth, volumetric momenta, temperature and composition remain
  !>          unchanged; neither the input physical nor conservative state is
  !>          modified. The function is allocation-free and thread-safe.
  !> \param[in] qp_center Read-only physical centre fields, including halo.
  !> \param[in] component Physical variable index; idx_u/idx_v are auxiliary velocities.
  !> \param[in] j X index of the centre stencil entry.
  !> \param[in] k Y index of the centre stencil entry.
  !> \return Original centre value, or zero for a dry auxiliary velocity.
  PURE FUNCTION reconstruction_cell_value(qp_center,component,j,k) RESULT(value)
    REAL(wp), INTENT(IN) :: qp_center(:,:,:)
    INTEGER, INTENT(IN) :: component,j,k
    REAL(wp) :: value

    value = qp_center(component,j,k)
    IF ((component .EQ. idx_u) .OR. (component .EQ. idx_v)) THEN
       IF (qp_center(1,j,k) .LE. hp_dry_tolerance) value = 0.0_wp
    END IF
  END FUNCTION reconstruction_cell_value

  !******************************************************************************
  !> \brief Evaluate hydrostatic candidates, blend them and commit final cell/face traces.
  !>
   !> The first reconstruction pass supplies direct h_x-based candidates for
   !> advective fluxes. The HP phases build eta-based traces for dissipative
   !> jumps and hydrostatic paths, then apply the existing bounded momentum
   !> correction on the active target cells.
  !>
  !> \param[in,out] this Persistent reconstructed traces, HP candidates and per-thread scratch.
  !> \param[in] qp_center Unmodified cell-centered physical states used by all HP candidate phases.
  !> \param[in] solve_cells Number of entries in the active-cell index lists.
  !> \param[in] j_cent X indices of active cells; only the first solve_cells entries are used.
  !> \param[in] k_cent Y indices of active cells; only the first solve_cells entries are used.
  !> \param[in] use_kt_direct Use direct KT traces instead of the HP dry-state projection.
  !******************************************************************************

  SUBROUTINE apply_hp_reconstruction( this, qp_center, solve_cells, j_cent,    &
       k_cent, use_kt_direct )

    USE state_conversion_2d, ONLY : qp_to_qc
    USE state_conversion_2d, ONLY : enforce_primitive_mass_fraction_closure
    USE constitutive_parameters_2d, ONLY : T_ambient
    USE parameters_2d, ONLY : limiter, reconstr_coeff
    USE geometry_2d, ONLY : sourceW, sourceE, sourceS, sourceN

    IMPLICIT NONE

    CLASS(reconstruction_workspace_type), INTENT(INOUT) :: this
    REAL(wp), INTENT(IN) :: qp_center(:,:,:)
    INTEGER, INTENT(IN) :: solve_cells
    INTEGER, INTENT(IN) :: j_cent(:), k_cent(:)
    LOGICAL, INTENT(IN) :: use_kt_direct

    REAL(wp) :: q_final(n_vars), qp_final(n_vars+2)
    INTEGER :: j, k, l

    ! Reuse one thread team across dependent reconstruction phases. Each END DO keeps
    ! its implicit barrier because the following phase consumes its results.
    !$OMP PARALLEL PRIVATE(l,j,k,q_final,qp_final)

    !$OMP DO
    DO l = 1, this%hp_eta_cells
       j = this%hp_eta_j(l)
       k = this%hp_eta_k(l)
       CALL compute_eta_traces(j,k)
    END DO
    !$OMP END DO

    ! Reconstruct both face momenta from the directly limited velocity traces.
    ! The PCCU interface evaluator selects the transport depth and recomputes
    ! momentum using h_x, minmod, Froude blend or eta_x as configured.
    !$OMP DO
    DO l = 1, this%hp_eta_cells
       j = this%hp_eta_j(l)
       k = this%hp_eta_k(l)
       CALL form_direct_face_momentum(this%qp_cellW(:,j,k))
       CALL form_direct_face_momentum(this%qp_cellE(:,j,k))
       CALL form_direct_face_momentum(this%qp_cellS(:,j,k))
       CALL form_direct_face_momentum(this%qp_cellN(:,j,k))
    END DO
    !$OMP END DO

    ! Preserve the direct face traces for the KT advection path. The selected
    ! depth will be applied consistently to both momenta at each PCCU face.
    IF ( use_kt_direct ) THEN
       !$OMP DO
       DO l = 1, this%hp_eta_cells
          j = this%hp_eta_j(l)
          k = this%hp_eta_k(l)
          CALL final_kt_state( this%qp_cellW(:,j,k),                         &
               this%q_adv_interfaceR(:,j,k), this%qp_adv_interfaceR(:,j,k),  &
               'KT reconstruction west', &
               'KT reconstruction west unit-depth state' )
          CALL final_kt_state( this%qp_cellE(:,j,k),                         &
               this%q_adv_interfaceL(:,j+1,k), this%qp_adv_interfaceL(:,j+1,k), &
               'KT reconstruction east', &
               'KT reconstruction east unit-depth state' )
          CALL final_kt_state( this%qp_cellS(:,j,k),                         &
               this%q_adv_interfaceT(:,j,k), this%qp_adv_interfaceT(:,j,k),  &
               'KT reconstruction south', &
               'KT reconstruction south unit-depth state' )
          CALL final_kt_state( this%qp_cellN(:,j,k),                         &
               this%q_adv_interfaceB(:,j,k+1), this%qp_adv_interfaceB(:,j,k+1), &
               'KT reconstruction north', &
               'KT reconstruction north unit-depth state' )
       END DO
       !$OMP END DO
    END IF

    ! Map all eta cell traces to their oriented face storage.  Conservative
    ! only workset cells can contribute to an active interface.
    !$OMP DO
    DO l = 1, this%hp_eta_cells
       j = this%hp_eta_j(l)
       k = this%hp_eta_k(l)
       this%eta_interfaceR(j,k) = this%eta_cellW(j,k)
       this%eta_interfaceL(j+1,k) = this%eta_cellE(j,k)
       this%eta_interfaceT(j,k) = this%eta_cellS(j,k)
       this%eta_interfaceB(j,k+1) = this%eta_cellN(j,k)
       IF ( j .EQ. 1 ) this%eta_interfaceL(1,k) = this%eta_interfaceR(1,k)
       IF ( j .EQ. comp_cells_x ) this%eta_interfaceR(comp_interfaces_x,k) = &
            this%eta_interfaceL(comp_interfaces_x,k)
       IF ( k .EQ. 1 ) this%eta_interfaceB(j,1) = this%eta_interfaceT(j,1)
       IF ( k .EQ. comp_cells_y ) this%eta_interfaceT(j,comp_interfaces_y) = &
            this%eta_interfaceB(j,comp_interfaces_y)
    END DO
    !$OMP END DO

    !$OMP DO
    DO l = 1, this%hp_eta_cells

   j = this%hp_eta_j(l)
   k = this%hp_eta_k(l)

       IF ( .NOT. use_kt_direct ) THEN
          CALL final_hp_state( this%qp_cellW(:,j,k),                         &
               this%q_adv_interfaceR(:,j,k), this%qp_adv_interfaceR(:,j,k) )
          CALL final_hp_state( this%qp_cellE(:,j,k),                         &
               this%q_adv_interfaceL(:,j+1,k), this%qp_adv_interfaceL(:,j+1,k) )
          CALL final_hp_state( this%qp_cellS(:,j,k),                         &
               this%q_adv_interfaceT(:,j,k), this%qp_adv_interfaceT(:,j,k) )
          CALL final_hp_state( this%qp_cellN(:,j,k),                         &
               this%q_adv_interfaceB(:,j,k+1), this%qp_adv_interfaceB(:,j,k+1) )
       END IF

       CALL final_hp_state( this%qp_cellW(:,j,k), q_final, qp_final )
       this%q_interfaceR(:,j,k) = q_final
       this%qp_interfaceR(:,j,k) = qp_final

       CALL final_hp_state( this%qp_cellE(:,j,k), q_final, qp_final )
       this%q_interfaceL(:,j+1,k) = q_final
       this%qp_interfaceL(:,j+1,k) = qp_final

       CALL final_hp_state( this%qp_cellS(:,j,k), q_final, qp_final )
       this%q_interfaceT(:,j,k) = q_final
       this%qp_interfaceT(:,j,k) = qp_final

       CALL final_hp_state( this%qp_cellN(:,j,k), q_final, qp_final )
       this%q_interfaceB(:,j,k+1) = q_final
       this%qp_interfaceB(:,j,k+1) = qp_final

    END DO
    !$OMP END DO

    ! Preserve the solver's existing external ghost convention after replacing
    ! the cell-owned traces.
    !$OMP DO
    DO k = 1, comp_cells_y
       this%q_interfaceL(:,1,k) = this%q_interfaceR(:,1,k)
       this%qp_interfaceL(:,1,k) = this%qp_interfaceR(:,1,k)
       this%q_adv_interfaceL(:,1,k) = this%q_adv_interfaceR(:,1,k)
       this%qp_adv_interfaceL(:,1,k) = this%qp_adv_interfaceR(:,1,k)
       this%q_interfaceR(:,comp_interfaces_x,k) =                             &
            this%q_interfaceL(:,comp_interfaces_x,k)
       this%qp_interfaceR(:,comp_interfaces_x,k) =                            &
            this%qp_interfaceL(:,comp_interfaces_x,k)
       this%q_adv_interfaceR(:,comp_interfaces_x,k) =                         &
            this%q_adv_interfaceL(:,comp_interfaces_x,k)
       this%qp_adv_interfaceR(:,comp_interfaces_x,k) =                        &
            this%qp_adv_interfaceL(:,comp_interfaces_x,k)
    END DO
    !$OMP END DO

    !$OMP DO
    DO j = 1, comp_cells_x
       this%q_interfaceB(:,j,1) = this%q_interfaceT(:,j,1)
       this%qp_interfaceB(:,j,1) = this%qp_interfaceT(:,j,1)
       this%q_adv_interfaceB(:,j,1) = this%q_adv_interfaceT(:,j,1)
       this%qp_adv_interfaceB(:,j,1) = this%qp_adv_interfaceT(:,j,1)
       this%q_interfaceT(:,j,comp_interfaces_y) =                             &
            this%q_interfaceB(:,j,comp_interfaces_y)
       this%qp_interfaceT(:,j,comp_interfaces_y) =                            &
            this%qp_interfaceB(:,j,comp_interfaces_y)
       this%q_adv_interfaceT(:,j,comp_interfaces_y) =                        &
            this%q_adv_interfaceB(:,j,comp_interfaces_y)
       this%qp_adv_interfaceT(:,j,comp_interfaces_y) =                       &
            this%qp_adv_interfaceB(:,j,comp_interfaces_y)
    END DO
    !$OMP END DO

    ! Rebuild the internal reflecting side of radial-source cells from the new
    ! HP state.  The volumetric momentum and auxiliary normal velocity are both
    ! reflected so qp and q remain consistent.
    IF ( radial_source_flag ) THEN
       !$OMP DO
       DO l = 1, solve_cells
          j = j_cent(l)
          k = k_cent(l)
          IF ( source_cell(j,k) .NE. 2 ) CYCLE

          IF ( sourceE(j,k) ) THEN
             this%q_interfaceR(:,j+1,k) = this%q_interfaceL(:,j+1,k)
             this%q_interfaceR(2,j+1,k) = -this%q_interfaceR(2,j+1,k)
             this%qp_interfaceR(:,j+1,k) = this%qp_interfaceL(:,j+1,k)
             this%qp_interfaceR(2,j+1,k) = -this%qp_interfaceR(2,j+1,k)
             this%qp_interfaceR(idx_u,j+1,k) =                               &
                  -this%qp_interfaceR(idx_u,j+1,k)
             this%q_adv_interfaceR(:,j+1,k) = this%q_adv_interfaceL(:,j+1,k)
             this%q_adv_interfaceR(2,j+1,k) =                                &
                  -this%q_adv_interfaceR(2,j+1,k)
             this%qp_adv_interfaceR(:,j+1,k) =                                &
                  this%qp_adv_interfaceL(:,j+1,k)
             this%qp_adv_interfaceR(2,j+1,k) =                                &
                  -this%qp_adv_interfaceR(2,j+1,k)
             this%qp_adv_interfaceR(idx_u,j+1,k) =                            &
                  -this%qp_adv_interfaceR(idx_u,j+1,k)
             this%eta_interfaceR(j+1,k) = this%eta_interfaceL(j+1,k)
          ELSEIF ( sourceW(j,k) ) THEN
             this%q_interfaceL(:,j,k) = this%q_interfaceR(:,j,k)
             this%q_interfaceL(2,j,k) = -this%q_interfaceL(2,j,k)
             this%qp_interfaceL(:,j,k) = this%qp_interfaceR(:,j,k)
             this%qp_interfaceL(2,j,k) = -this%qp_interfaceL(2,j,k)
             this%qp_interfaceL(idx_u,j,k) = -this%qp_interfaceL(idx_u,j,k)
             this%q_adv_interfaceL(:,j,k) = this%q_adv_interfaceR(:,j,k)
             this%q_adv_interfaceL(2,j,k) = -this%q_adv_interfaceL(2,j,k)
             this%qp_adv_interfaceL(:,j,k) = this%qp_adv_interfaceR(:,j,k)
             this%qp_adv_interfaceL(2,j,k) = -this%qp_adv_interfaceL(2,j,k)
             this%qp_adv_interfaceL(idx_u,j,k) =                             &
                  -this%qp_adv_interfaceL(idx_u,j,k)
             this%eta_interfaceL(j,k) = this%eta_interfaceR(j,k)
          END IF

          IF ( sourceN(j,k) ) THEN
             this%q_interfaceT(:,j,k+1) = this%q_interfaceB(:,j,k+1)
             this%q_interfaceT(3,j,k+1) = -this%q_interfaceT(3,j,k+1)
             this%qp_interfaceT(:,j,k+1) = this%qp_interfaceB(:,j,k+1)
             this%qp_interfaceT(3,j,k+1) = -this%qp_interfaceT(3,j,k+1)
             this%qp_interfaceT(idx_v,j,k+1) =                               &
                  -this%qp_interfaceT(idx_v,j,k+1)
             this%q_adv_interfaceT(:,j,k+1) =                                 &
                  this%q_adv_interfaceB(:,j,k+1)
             this%q_adv_interfaceT(3,j,k+1) =                                 &
                  -this%q_adv_interfaceT(3,j,k+1)
             this%qp_adv_interfaceT(:,j,k+1) =                                &
                  this%qp_adv_interfaceB(:,j,k+1)
             this%qp_adv_interfaceT(3,j,k+1) =                                &
                  -this%qp_adv_interfaceT(3,j,k+1)
             this%qp_adv_interfaceT(idx_v,j,k+1) =                            &
                  -this%qp_adv_interfaceT(idx_v,j,k+1)
             this%eta_interfaceT(j,k+1) = this%eta_interfaceB(j,k+1)
          ELSEIF ( sourceS(j,k) ) THEN
             this%q_interfaceB(:,j,k) = this%q_interfaceT(:,j,k)
             this%q_interfaceB(3,j,k) = -this%q_interfaceB(3,j,k)
             this%qp_interfaceB(:,j,k) = this%qp_interfaceT(:,j,k)
             this%qp_interfaceB(3,j,k) = -this%qp_interfaceB(3,j,k)
             this%qp_interfaceB(idx_v,j,k) = -this%qp_interfaceB(idx_v,j,k)
             this%q_adv_interfaceB(:,j,k) = this%q_adv_interfaceT(:,j,k)
             this%q_adv_interfaceB(3,j,k) = -this%q_adv_interfaceB(3,j,k)
             this%qp_adv_interfaceB(:,j,k) = this%qp_adv_interfaceT(:,j,k)
             this%qp_adv_interfaceB(3,j,k) = -this%qp_adv_interfaceB(3,j,k)
             this%qp_adv_interfaceB(idx_v,j,k) =                             &
                  -this%qp_adv_interfaceB(idx_v,j,k)
             this%eta_interfaceB(j,k) = this%eta_interfaceT(j,k)
          END IF
       END DO
       !$OMP END DO
    END IF

    !$OMP END PARALLEL

  CONTAINS

    !******************************************************************************
    !> \brief Convert direct KT traces without the HP dry-state projection.
    !> Python clamps reconstructed h to zero but keeps independently
    !> reconstructed hu/hv for the inertial momentum flux. The unit-depth
    !> conversion below retains that conservative momentum even when h is
    !> clipped to zero, while mass-like quantities still scale with h.
    !> \param[in] qp_candidate Reconstructed physical state at the face.
    !> \param[out] q_conservative Conservative state converted from the trace.
    !> \param[out] qp_reconstructed Physical state after depth and composition closure.
    !> \param[in] j_cell X index of the adjacent cell.
    !> \param[in] k_cell Y index of the adjacent cell.
    !> \param[in] side Label identifying the face side for diagnostics.
    !******************************************************************************
    SUBROUTINE final_kt_state(qp_candidate,q_conservative,qp_reconstructed,  &
         context,unit_context)

      REAL(wp), INTENT(IN) :: qp_candidate(n_vars+2)
      REAL(wp), INTENT(OUT) :: q_conservative(n_vars)
      REAL(wp), INTENT(OUT) :: qp_reconstructed(n_vars+2)
      CHARACTER(LEN=*), INTENT(IN) :: context,unit_context

      REAL(wp) :: qp_unit(n_vars+2),q_unit(n_vars),h_trace

      ! The caller passes constant, side-specific diagnostic strings.
      ! No formatting, trimming or concatenation is needed per KT face.

      qp_reconstructed = qp_candidate
      h_trace = MAX(qp_candidate(1),0.0_wp)
      qp_reconstructed(1) = h_trace

      IF ( h_trace .GT. EPSILON(1.0_wp) ) THEN
         CALL enforce_primitive_mass_fraction_closure(qp_reconstructed,      &
              context)
         CALL qp_to_qc(qp_reconstructed,q_conservative,context)
      ELSE
         qp_unit = qp_reconstructed
         qp_unit(1) = 1.0_wp
         CALL enforce_primitive_mass_fraction_closure(qp_unit,               &
              unit_context)
         CALL qp_to_qc(qp_unit,q_unit,unit_context)
         q_conservative = h_trace*q_unit
         q_conservative(2) = q_unit(2)
         q_conservative(3) = q_unit(3)
      END IF

    END SUBROUTINE final_kt_state

    !> \brief Compute positivity-constrained free-surface traces for one parent-workset cell.
    !>
    !> \param[in] jc X index of the current cell in the parent reconstruction.
    !> \param[in] kc Y index of the current cell in the parent reconstruction.
    !>
    !> \note Reads the host qp_center and shared face bed; writes only this cell eta candidates in
    !>       the parent workspace.

    SUBROUTINE compute_eta_traces(jc,kc)

    INTEGER, INTENT(IN) :: jc, kc
    REAL(wp) :: eta_left, eta_center, eta_right

    eta_center = qp_center(1,jc,kc) +                                     &
       0.5_wp * ( B_face_x(jc,kc) + B_face_x(jc+1,kc) )
    eta_left = eta_center
    eta_right = eta_center
    IF ( ( jc .GT. 1 ) .AND. ( jc .LT. comp_cells_x ) ) THEN
       eta_left = qp_center(1,jc-1,kc) +                                  &
          0.5_wp * ( B_face_x(jc-1,kc) + B_face_x(jc,kc) )
       eta_right = qp_center(1,jc+1,kc) +                                 &
          0.5_wp * ( B_face_x(jc+1,kc) + B_face_x(jc+2,kc) )
    END IF
    CALL hp_eta_face_pair( eta_left, eta_center, eta_right,                &
       B_face_x(jc,kc), B_face_x(jc+1,kc),                              &
       ( jc .EQ. 1 ) .OR. ( jc .EQ. comp_cells_x ),                     &
       this%eta_cellW(jc,kc), this%eta_cellE(jc,kc) )

    eta_center = qp_center(1,jc,kc) +                                     &
       0.5_wp * ( B_face_y(jc,kc) + B_face_y(jc,kc+1) )
    eta_left = eta_center
    eta_right = eta_center
    IF ( ( kc .GT. 1 ) .AND. ( kc .LT. comp_cells_y ) ) THEN
       eta_left = qp_center(1,jc,kc-1) +                                  &
          0.5_wp * ( B_face_y(jc,kc-1) + B_face_y(jc,kc) )
       eta_right = qp_center(1,jc,kc+1) +                                 &
          0.5_wp * ( B_face_y(jc,kc+1) + B_face_y(jc,kc+2) )
    END IF
    CALL hp_eta_face_pair( eta_left, eta_center, eta_right,                &
       B_face_y(jc,kc), B_face_y(jc,kc+1),                              &
       ( kc .EQ. 1 ) .OR. ( kc .EQ. comp_cells_y ),                     &
       this%eta_cellS(jc,kc), this%eta_cellN(jc,kc) )

    END SUBROUTINE compute_eta_traces

    !> \brief Limit a three-cell free-surface slope against both endpoint bed elevations.
    !>
    !> \param[in] eta_left Free-surface elevation at the negative/left endpoint or stencil neighbour
    !>                     [m].
    !> \param[in] eta_center Free-surface elevation at the central stencil cell [m].
    !> \param[in] eta_right Free-surface elevation at the positive/right endpoint or stencil
    !>                      neighbour [m].
    !> \param[in] B_minus Bed elevations at the negative-side cell faces [m].
    !> \param[in] B_plus Bed elevations at the positive-side cell faces [m].
    !> \param[in] at_boundary True to use the zero-gradient boundary eta-slope convention.
    !> \param[out] eta_minus Negative-side hydrostatic free-surface candidate [m].
    !> \param[out] eta_plus Positive-side hydrostatic free-surface candidate [m].

    SUBROUTINE hp_eta_face_pair( eta_left, eta_center, eta_right,            &
       B_minus, B_plus, at_boundary, eta_minus, eta_plus )

    REAL(wp), INTENT(IN) :: eta_left, eta_center, eta_right
    REAL(wp), INTENT(IN) :: B_minus, B_plus
    LOGICAL, INTENT(IN) :: at_boundary
    REAL(wp), INTENT(OUT) :: eta_minus, eta_plus
    REAL(wp) :: eta_stencil(3), coordinate_stencil(3)
    REAL(wp) :: eta_slope, slope_min, slope_max

    eta_slope = 0.0_wp
    IF ( .NOT. at_boundary ) THEN
       eta_stencil = [ eta_left, eta_center, eta_right ]
       coordinate_stencil = [ -1.0_wp, 0.0_wp, 1.0_wp ]
       CALL limit( eta_stencil, coordinate_stencil, limiter(1), eta_slope )
       eta_slope = reconstr_coeff * eta_slope
    END IF

    slope_min = 2.0_wp * ( B_plus - eta_center )
    slope_max = 2.0_wp * ( eta_center - B_minus )
    eta_slope = MIN( MAX( eta_slope, slope_min ), slope_max )

    eta_minus = eta_center - 0.5_wp * eta_slope
    eta_plus = eta_center + 0.5_wp * eta_slope

    END SUBROUTINE hp_eta_face_pair

    !> Combine directly limited velocity traces with a nonnegative depth.
    !> The selected transport depth is applied consistently to both moments
    !> by the downstream PCCU face evaluator.
    SUBROUTINE form_direct_face_momentum(trace)
      REAL(wp), INTENT(INOUT) :: trace(n_vars+2)
      trace(1) = MAX(trace(1),0.0_wp)
      IF (trace(1) .LE. hp_dry_tolerance) THEN
         trace(1:3) = 0.0_wp
         trace(idx_u:idx_v) = 0.0_wp
      ELSE
         trace(2) = trace(1)*trace(idx_u)
         trace(3) = trace(1)*trace(idx_v)
      END IF
    END SUBROUTINE form_direct_face_momentum

    !> \brief Close component fractions and convert a final HP candidate into conservative form.
    !>
    !> \param[in] qp_candidate Final HP physical candidate before composition closure.
    !> \param[out] q_conservative Conservative state obtained from the closed physical candidate.
    !> \param[out] qp_reconstructed Closed physical state retaining independently reconstructed
    !>                              auxiliary velocity traces.

    SUBROUTINE final_hp_state(qp_candidate,q_conservative,qp_reconstructed)

      REAL(wp), INTENT(IN) :: qp_candidate(n_vars+2)
      REAL(wp), INTENT(OUT) :: q_conservative(n_vars)
      REAL(wp), INTENT(OUT) :: qp_reconstructed(n_vars+2)

      qp_reconstructed = qp_candidate

      IF ( qp_reconstructed(1) .LE. hp_dry_tolerance ) THEN
         qp_reconstructed = 0.0_wp
         qp_reconstructed(4) = T_ambient
         q_conservative = 0.0_wp
         RETURN
      END IF

      CALL enforce_primitive_mass_fraction_closure(qp_reconstructed)
      CALL qp_to_qc(qp_reconstructed,q_conservative)

    END SUBROUTINE final_hp_state

  END SUBROUTINE apply_hp_reconstruction

END MODULE reconstruction_2d
