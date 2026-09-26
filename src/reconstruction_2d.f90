!********************************************************************************
!> \brief Reconstruction of cell states at computational interfaces
!>
!> This module owns the reconstructed conservative/physical interface states
!> and the corresponding divergence flags.
!********************************************************************************
MODULE reconstruction_2d

  USE parameters_2d, ONLY : wp, n_vars
  USE parameters_2d, ONLY : lateral_source_flag, radial_source_flag
  USE parameters_2d, ONLY : bcW, bcE, bcS, bcN
  USE parameters_2d, ONLY : idx_u, idx_v

  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE geometry_2d, ONLY : comp_interfaces_x, comp_interfaces_y
  USE geometry_2d, ONLY : B_cent, B_face_x, B_face_y
  USE geometry_2d, ONLY : source_cell
  USE geometry_2d, ONLY : one_by_dx, one_by_dy
  USE geometry_2d, ONLY : limit

   USE hp_reconstruction_2d, ONLY : hp_dry_tolerance
   USE hp_reconstruction_2d, ONLY : reconstruct_hp_line
   USE hp_reconstruction_2d, ONLY : hp_dynamic_residual_threshold
   USE omp_lib, ONLY : omp_get_max_threads, omp_get_thread_num

  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: reconstruction_workspace_type
     REAL(wp), ALLOCATABLE :: q_interfaceL(:,:,:)
     REAL(wp), ALLOCATABLE :: q_interfaceR(:,:,:)
     REAL(wp), ALLOCATABLE :: q_interfaceB(:,:,:)
     REAL(wp), ALLOCATABLE :: q_interfaceT(:,:,:)

     REAL(wp), ALLOCATABLE :: qp_interfaceL(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_interfaceR(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_interfaceB(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_interfaceT(:,:,:)

   ! Direct cell-side candidates retained until the local HP blend can
   ! compare continuity errors across adjacent interfaces.
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

     REAL(wp), ALLOCATABLE :: w_eta_x(:,:), w_eta_y(:,:)

       REAL(wp), ALLOCATABLE :: hydrostatic_residual_2d(:,:)
       REAL(wp), ALLOCATABLE :: topographic_relief_ratio_2d(:,:)
       LOGICAL, ALLOCATABLE :: hp_eta_mask(:,:)
       INTEGER, ALLOCATABLE :: hp_eta_j(:), hp_eta_k(:)
       INTEGER :: hp_eta_cells
       REAL(wp), ALLOCATABLE :: hp_blended(:,:,:,:)
       REAL(wp), ALLOCATABLE :: hp_scratch(:,:,:)

     LOGICAL, ALLOCATABLE :: diverg_interfaceL(:,:)
     LOGICAL, ALLOCATABLE :: diverg_interfaceR(:,:)
     LOGICAL, ALLOCATABLE :: diverg_interfaceB(:,:)
     LOGICAL, ALLOCATABLE :: diverg_interfaceT(:,:)
   CONTAINS
     PROCEDURE :: initialize => initialize_reconstruction
     PROCEDURE :: finalize => finalize_reconstruction
     PROCEDURE :: reconstruct => reconstruction
  END TYPE reconstruction_workspace_type

CONTAINS

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
    ALLOCATE( this%w_eta_x( comp_cells_x, comp_cells_y ) )
    ALLOCATE( this%w_eta_y( comp_cells_x, comp_cells_y ) )

    ALLOCATE( this%hydrostatic_residual_2d(comp_cells_x,comp_cells_y) )
    ALLOCATE( this%topographic_relief_ratio_2d(comp_cells_x,comp_cells_y) )
      ALLOCATE( this%hp_eta_mask(comp_cells_x,comp_cells_y) )
      ALLOCATE( this%hp_eta_j(comp_cells_x*comp_cells_y) )
      ALLOCATE( this%hp_eta_k(comp_cells_x*comp_cells_y) )
      ALLOCATE( this%hp_blended(4,comp_cells_x,comp_cells_y,2) )
      ALLOCATE( this%hp_scratch(5,8,MAX(1,omp_get_max_threads())) )
      this%hp_eta_mask = .FALSE.
      this%hp_eta_cells = 0

    ALLOCATE( this%diverg_interfaceL( comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( this%diverg_interfaceR( comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( this%diverg_interfaceB( comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( this%diverg_interfaceT( comp_cells_x, comp_interfaces_y ) )

  END SUBROUTINE initialize_reconstruction

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
    DEALLOCATE( this%w_eta_x )
    DEALLOCATE( this%w_eta_y )

      DEALLOCATE( this%hydrostatic_residual_2d )
      DEALLOCATE( this%topographic_relief_ratio_2d )
      DEALLOCATE( this%hp_eta_mask )
      DEALLOCATE( this%hp_eta_j )
      DEALLOCATE( this%hp_eta_k )
      DEALLOCATE( this%hp_blended )
      DEALLOCATE( this%hp_scratch )

    DEALLOCATE( this%diverg_interfaceL )
    DEALLOCATE( this%diverg_interfaceR )
    DEALLOCATE( this%diverg_interfaceB )
    DEALLOCATE( this%diverg_interfaceT )

  END SUBROUTINE finalize_reconstruction

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

      SUBROUTINE append_eta_cell(j_cell,k_cell)

         INTEGER, INTENT(IN) :: j_cell, k_cell

         IF ( this%hp_eta_mask(j_cell,k_cell) ) RETURN

         this%hp_eta_mask(j_cell,k_cell) = .TRUE.
         this%hp_eta_cells = this%hp_eta_cells + 1
         this%hp_eta_j(this%hp_eta_cells) = j_cell
         this%hp_eta_k(this%hp_eta_cells) = k_cell

      END SUBROUTINE append_eta_cell

   END SUBROUTINE build_hp_workset

   SUBROUTINE reconstruction( this, qp_expl, t, solve_cells, j_cent, k_cent )

    ! External procedures
      USE state_conversion_2d, ONLY : qp_to_qp2
    USE equation_terms_2d, ONLY : eval_source_bdry
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

    LOGICAL :: diverging_flag
    LOGICAL :: regular_interior

    ! Keep the same constant candidates for inactive neighbours, but initialize
    ! only cells that can contribute to an active HP trace.
    CALL build_hp_workset(this, solve_cells, j_cent, k_cent)

    !$OMP PARALLEL DO PRIVATE(l,j,k)
    DO l = 1, this%hp_eta_cells
       j = this%hp_eta_j(l)
       k = this%hp_eta_k(l)
       this%qp_cellW(:,j,k) = qp_expl(:,j,k)
       this%qp_cellE(:,j,k) = qp_expl(:,j,k)
       this%qp_cellS(:,j,k) = qp_expl(:,j,k)
       this%qp_cellN(:,j,k) = qp_expl(:,j,k)
    END DO
    !$OMP END PARALLEL DO

    !WRITE(*,*) 'recontruction 0'
    !WRITE(*,*) 'nvars',n_vars
    !WRITE(*,*) 'qp_expl(:,1,1)',qp_expl(:,1,1)

    !$OMP PARALLEL DO private(j,k,i,qrecW,qrecE,qrecS,qrecN,x_stencil,y_stencil,&
    !$OMP & qrec_stencil,qrec_prime_x,qrec_prime_y,qp2recW,qp2recE,qp2recS,     &
    !$OMP & qp2recN,source_bdry,dq,diverging_flag,regular_interior)

    DO l = 1,solve_cells

       j = j_cent(l)
       k = k_cent(l)

       qrecW(1:n_vars+2) = qp_expl(1:n_vars+2,j,k)
       qrecE(1:n_vars+2) = qp_expl(1:n_vars+2,j,k)
       qrecS(1:n_vars+2) = qp_expl(1:n_vars+2,j,k)
       qrecN(1:n_vars+2) = qp_expl(1:n_vars+2,j,k)

       x_stencil(2) = x_comp(j)
       y_stencil(2) = y_comp(k)

       ! Default source-side ghost state. For the conservative radial source,
       ! the ring is a wall and mass/momentum are injected through eval_expl_terms.
       ! The existing lateral boundary source keeps its Dirichlet treatment.
       source_bdry(1:n_vars+2) = qp_expl(1:n_vars+2,j,k)

       IF ( lateral_source_flag .AND. ( source_cell(j,k) .EQ. 2 ) ) THEN

          IF ( sourceE(j,k) ) THEN

             CALL eval_source_bdry( t, sourceE_vect_x(j,k), sourceE_vect_y(j,k), &
                  source_bdry )

          ELSEIF ( sourceW(j,k) ) THEN

             CALL eval_source_bdry( t, sourceW_vect_x(j,k), sourceW_vect_y(j,k), &
                  source_bdry )

          ELSEIF ( sourceS(j,k) ) THEN

             CALL eval_source_bdry( t, sourceS_vect_x(j,k), sourceS_vect_y(j,k), &
                  source_bdry )

          ELSEIF ( sourceN(j,k) ) THEN

             CALL eval_source_bdry( t, sourceN_vect_x(j,k), sourceN_vect_y(j,k), &
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

             qrec_stencil(2) = qp_expl(i,j,k)

             qrec_stencil(1) = qp_expl(i,j-1,k)
             qrec_stencil(3) = qp_expl(i,j+1,k)
             CALL limit( qrec_stencil , x_stencil , limiter(i) ,               &
                  qrec_prime_x(i) )

             dq = reconstr_coeff * dx2 * qrec_prime_x(i)
             qrecW(i) = qrec_stencil(2) - dq
             qrecE(i) = qrec_stencil(2) + dq

             qrec_stencil(1) = qp_expl(i,j,k-1)
             qrec_stencil(3) = qp_expl(i,j,k+1)
             CALL limit( qrec_stencil , y_stencil , limiter(i) ,               &
                  qrec_prime_y(i) )

             dq = reconstr_coeff * dy2 * qrec_prime_y(i)
             qrecS(i) = qrec_stencil(2) - dq
             qrecN(i) = qrec_stencil(2) + dq

          END DO fast_vars_loop

       ELSE

       vars_loop:DO i=1,n_vars

          qrec_stencil(2) = qp_expl(i,j,k)

          ! x direction
          check_comp_cells_x:IF ( comp_cells_x .GT. 1 ) THEN

             ! west boundary
             check_x_boundary:IF ( j .EQ. 1 ) THEN

                x_stencil(1) = x_stag(1)
                x_stencil(3) = x_comp(j+1)

                IF ( source_cell(j,k).EQ.2 ) THEN

                   ! Dirichlet boundary condition
                   qrec_stencil(1) = source_bdry(i)
                   qrec_stencil(3) = qp_expl(i,j+1,k)

                   CALL limit( qrec_stencil , x_stencil , limiter(i) ,          &
                        qrec_prime_x(i) )

                ELSE

                   IF ( bcW(i)%flag .EQ. 0 ) THEN

                      ! Dirichlet boundary condition
                      qrec_stencil(1) = bcW(i)%value
                      qrec_stencil(3) = qp_expl(i,j+1,k)

                      CALL limit( qrec_stencil , x_stencil , limiter(i) ,          &
                           qrec_prime_x(i) )

                   ELSEIF ( bcW(i)%flag .EQ. 1 ) THEN

                      ! Neumann boundary condition
                      qrec_prime_x(i) = bcW(i)%value

                   ELSEIF ( bcW(i)%flag .EQ. 2 ) THEN

                      qrec_prime_x(i) = ( qp_expl(i,2,k) - qp_expl(i,1,k) )        &
                           * one_by_dx

                   END IF

                END IF

                !east boundary
             ELSEIF ( j .EQ. comp_cells_x ) THEN

                x_stencil(3) = x_stag(comp_interfaces_x)
                x_stencil(1) = x_comp(j-1)

                IF ( source_cell(j,k).EQ.2 ) THEN

                   ! Dirichlet boundary condition
                   qrec_stencil(3) = source_bdry(i)
                   qrec_stencil(1)= qp_expl(i,j-1,k)

                   CALL limit( qrec_stencil , x_stencil , limiter(i) ,          &
                        qrec_prime_x(i) )

                ELSE

                   IF ( bcE(i)%flag .EQ. 0 ) THEN

                      ! Dirichlet boundary condition
                      qrec_stencil(3) = bcE(i)%value
                      qrec_stencil(1)= qp_expl(i,j-1,k)

                      CALL limit( qrec_stencil , x_stencil , limiter(i) ,          &
                           qrec_prime_x(i) )

                   ELSEIF ( bcE(i)%flag .EQ. 1 ) THEN

                      ! Neumann boundary condition
                      qrec_prime_x(i) = bcE(i)%value

                   ELSEIF ( bcE(i)%flag .EQ. 2 ) THEN

                      qrec_prime_x(i) = ( qp_expl(i,comp_cells_x,k) -              &
                           qp_expl(i,comp_cells_x-1,k) ) * one_by_dx

                   END IF

                END IF

             ELSE

                ! internal x cells

                x_stencil(1) = x_comp(j-1)
                x_stencil(3) = x_comp(j+1)

                qrec_stencil(1) = qp_expl(i,j-1,k)
                qrec_stencil(3) = qp_expl(i,j+1,k)

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

                IF ( bcS(i)%flag .EQ. 0 ) THEN

                   ! Dirichlet boundary condition
                   qrec_stencil(1) = bcS(i)%value
                   qrec_stencil(3) = qp_expl(i,j,k+1)

                   CALL limit( qrec_stencil , y_stencil , limiter(i) ,          &
                        qrec_prime_y(i) )

                ELSEIF ( bcS(i)%flag .EQ. 1 ) THEN

                   ! Neumann boundary condition
                   qrec_prime_y(i) = bcS(i)%value

                ELSEIF ( bcS(i)%flag .EQ. 2 ) THEN

                   qrec_prime_y(i) = ( qp_expl(i,j,2) - qp_expl(i,j,1) )        &
                        * one_by_dy

                END IF

                ! North boundary
             ELSEIF ( k .EQ. comp_cells_y ) THEN

                y_stencil(1) = y_comp(k-1)
                y_stencil(3) = y_stag(comp_interfaces_y)

                IF ( bcN(i)%flag .EQ. 0 ) THEN

                   ! Dirichlet boundary condition
                   qrec_stencil(1)= qp_expl(i,j,k-1)
                   qrec_stencil(3) = bcN(i)%value

                   CALL limit( qrec_stencil , y_stencil , limiter(i) ,          &
                        qrec_prime_y(i) )

                ELSEIF ( bcN(i)%flag .EQ. 1 ) THEN

                   ! Neumann boundary condition
                   qrec_prime_y(i) = bcN(i)%value

                ELSEIF ( bcN(i)%flag .EQ. 2 ) THEN

                   qrec_prime_y(i) = ( qp_expl(i,j,comp_cells_y) -              &
                        qp_expl(i,j,comp_cells_y-1) ) * one_by_dy

                END IF

             ELSE

                ! Internal y cells

                y_stencil(1) = y_comp(k-1)
                y_stencil(3) = y_comp(k+1)

                qrec_stencil(1) = qp_expl(i,j,k-1)
                qrec_stencil(3) = qp_expl(i,j,k+1)

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

                ! Dirichelet boundary condition at the south of the domain
                IF ( bcS(i)%flag .EQ. 0 ) THEN

                   qrecS(i) = bcS(i)%value

                ELSE

                   IF ( i .EQ. 3 ) qrecS(i) = MIN( qrecS(i) , 0.0_wp )

                END IF

             END IF

             IF ( k .EQ. comp_cells_y ) THEN

                ! Dirichelet boundary condition at the north of the domain
                IF ( bcN(i)%flag .EQ. 0 ) THEN

                   qrecN(i) = bcN(i)%value

                ELSE

                   IF ( i .EQ. 3 ) qrecN(i) = MAX( qrecN(i) , 0.0_wp )

                END IF

             END IF

          ENDIF check_comp_cells_y

       ENDDO vars_loop

       add_vars_loop:DO i=n_vars+1,n_vars+2
          ! reconstruction on u and v with same limiters of hu,hv

          ! x direction
          check_comp_cells_x2:IF ( comp_cells_x .GT. 1 ) THEN

             qrec_stencil(2) = qp_expl(i,j,k)

             IF ( j .EQ. 1 ) THEN

                CALL qp_to_qp2( qrecW(1:n_vars+2) , B_cent(j,k) , qp2recW )
                qrec_stencil(1) = qp2recW(i-n_vars+1)
                qrec_stencil(3) = qp_expl(i,j+1,k)

             ELSEIF ( j .EQ. comp_cells_x ) THEN

                CALL qp_to_qp2( qrecE(1:n_vars+2) , B_cent(j,k) , qp2recE )
                qrec_stencil(1) = qp_expl(i,j-1,k)
                qrec_stencil(3) = qp2recE(i-n_vars+1)

             ELSE

                qrec_stencil(1) = qp_expl(i,j-1,k)
                qrec_stencil(3) = qp_expl(i,j+1,k)

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

             qrec_stencil(2) = qp_expl(i,j,k)

             IF ( k .EQ. 1 ) THEN

                CALL qp_to_qp2( qrecS(1:n_vars+2) , B_cent(j,k) , qp2recS )
                qrec_stencil(1) = qp2recS(i-n_vars+1)
                qrec_stencil(3) = qp_expl(i,j,k+1)

             ELSEIF ( k .EQ. comp_cells_y ) THEN

                CALL qp_to_qp2( qrecN(1:n_vars+2) , B_cent(j,k) , qp2recN )
                qrec_stencil(1) = qp_expl(i,j,k-1)
                qrec_stencil(3) = qp2recN(i-n_vars+1)

             ELSE

                qrec_stencil(1) = qp_expl(i,j,k-1)
                qrec_stencil(3) = qp_expl(i,j,k+1)

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

       ! check if du/dx + dv/dy > 0 (flow locally diverges)
       diverging_flag = ( ( qrec_prime_x(n_vars+1) + qrec_prime_y(n_vars+2) )   &
            .GT. 0.0_wp )

       this%qp_cellW(:,j,k) = qrecW
       this%qp_cellE(:,j,k) = qrecE
       this%qp_cellS(:,j,k) = qrecS
       this%qp_cellN(:,j,k) = qrecN

       IF ( comp_cells_x .GT. 1 ) THEN

          this%diverg_interfaceR(j,k) = diverging_flag
          this%diverg_interfaceL(j+1,k) = diverging_flag

          IF ( j.EQ.1 ) THEN

             ! Interface value at the left of first x-interface (external)
             this%diverg_interfaceR(j,k) = this%diverg_interfaceL(j,k)

          ELSEIF ( j.EQ.comp_cells_x ) THEN

             ! Interface value at the right of last x-interface (external)
             this%diverg_interfaceR(j+1,k) = this%diverg_interfaceL(j+1,k)

          ELSE

          END IF

       ELSE

          this%diverg_interfaceR(j,k) = diverging_flag
          this%diverg_interfaceL(j+1,k) = diverging_flag

       END IF

       IF ( comp_cells_y .GT. 1 ) THEN

          this%diverg_interfaceT(j,k) = diverging_flag
          this%diverg_interfaceB(j,k+1) = diverging_flag

          IF ( k .EQ. 1 ) THEN

             ! Interface value at the bottom of first y-interface (external)
             this%diverg_interfaceB(j,k) = this%diverg_interfaceT(j,k)

          ELSEIF ( k .EQ. comp_cells_y ) THEN

             ! Interface value at the top of last y-interface (external)
             this%diverg_interfaceT(j,k+1) = this%diverg_interfaceB(j,k+1)

          ELSE

          END IF

       ELSE

       END IF

    END DO

    !$OMP END PARALLEL DO

    CALL apply_hp_reconstruction( this, qp_expl, solve_cells, j_cent, k_cent )

    RETURN

  END SUBROUTINE reconstruction

  !******************************************************************************
  !> \brief Replace direct cell traces with the final HP face states
  !>
   !> The first reconstruction pass supplies direct limited candidates. The HP
   !> phases build local eta traces, evaluate continuity weights, then apply
   !> blended thickness and momentum on the active target cells.
  !******************************************************************************
  SUBROUTINE apply_hp_reconstruction( this, qp_center, solve_cells, j_cent,    &
       k_cent )

    USE state_conversion_2d, ONLY : qp_to_qc
    USE state_conversion_2d, ONLY : enforce_primitive_mass_fraction_closure
    USE state_conversion_2d, ONLY : velocity_from_conservative
    USE constitutive_parameters_2d, ONLY : T_ambient
    USE parameters_2d, ONLY : limiter, reconstr_coeff
    USE geometry_2d, ONLY : sourceW, sourceE, sourceS, sourceN

    IMPLICIT NONE

    CLASS(reconstruction_workspace_type), INTENT(INOUT) :: this
    REAL(wp), INTENT(IN) :: qp_center(:,:,:)
    INTEGER, INTENT(IN) :: solve_cells
    INTEGER, INTENT(IN) :: j_cent(:), k_cent(:)

    REAL(wp) :: q_final(n_vars), qp_final(n_vars+2)
    REAL(wp) :: relief2d
    INTEGER :: j, k, l

    ! Phase 1: compute the cell-local indicators and eta traces for cells
    ! adjacent to the solve mask.
    !$OMP PARALLEL DO PRIVATE(l,j,k,relief2d)
    DO l = 1, solve_cells
       j = j_cent(l)
       k = k_cent(l)
       this%hydrostatic_residual_2d(j,k) = local_hydrostatic_residual(j,k)
       relief2d = MAX( B_face_x(j,k), B_face_x(j+1,k),                      &
            B_face_y(j,k), B_face_y(j,k+1) ) -                              &
            MIN( B_face_x(j,k), B_face_x(j+1,k),                            &
            B_face_y(j,k), B_face_y(j,k+1) )
       this%topographic_relief_ratio_2d(j,k) = relief2d /                    &
            MAX(qp_center(1,j,k),hp_dry_tolerance)
    END DO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(l,j,k)
    DO l = 1, this%hp_eta_cells
       j = this%hp_eta_j(l)
       k = this%hp_eta_k(l)
       CALL compute_eta_traces(j,k)
    END DO
    !$OMP END PARALLEL DO

    ! Phase 2: compare the direct and eta thickness jumps at each target cell.
    !$OMP PARALLEL DO PRIVATE(l,j,k)
    DO l = 1, solve_cells
       j = j_cent(l)
       k = k_cent(l)
       CALL compute_cell_weights(j,k)
    END DO
    !$OMP END PARALLEL DO

    ! Phase 3: blend thickness and normal momentum only on solve_cells.
    !$OMP PARALLEL DO PRIVATE(l,j,k)
    DO l = 1, solve_cells
       j = j_cent(l)
       k = k_cent(l)
       CALL evaluate_hp_cell_x(j,k)
       CALL evaluate_hp_cell_y(j,k)
    END DO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(l,j,k)
    DO l = 1, solve_cells
       j = j_cent(l)
       k = k_cent(l)
       CALL commit_hp_cell(j,k)
    END DO
    !$OMP END PARALLEL DO

    ! Map all eta cell traces to their oriented face storage.  Conservative
    ! only workset cells can contribute to an active interface.
    !$OMP PARALLEL DO PRIVATE(l,j,k)
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
    !$OMP END PARALLEL DO

    DO l = 1, solve_cells

       j = j_cent(l)
       k = k_cent(l)

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

    ! Preserve the solver's existing external ghost convention after replacing
    ! the cell-owned traces.
    DO k = 1, comp_cells_y
       this%q_interfaceL(:,1,k) = this%q_interfaceR(:,1,k)
       this%qp_interfaceL(:,1,k) = this%qp_interfaceR(:,1,k)
       this%q_interfaceR(:,comp_interfaces_x,k) =                             &
            this%q_interfaceL(:,comp_interfaces_x,k)
       this%qp_interfaceR(:,comp_interfaces_x,k) =                            &
            this%qp_interfaceL(:,comp_interfaces_x,k)
    END DO

    DO j = 1, comp_cells_x
       this%q_interfaceB(:,j,1) = this%q_interfaceT(:,j,1)
       this%qp_interfaceB(:,j,1) = this%qp_interfaceT(:,j,1)
       this%q_interfaceT(:,j,comp_interfaces_y) =                             &
            this%q_interfaceB(:,j,comp_interfaces_y)
       this%qp_interfaceT(:,j,comp_interfaces_y) =                            &
            this%qp_interfaceB(:,j,comp_interfaces_y)
    END DO

    ! Rebuild the internal reflecting side of radial-source cells from the new
    ! HP state.  The volumetric momentum and auxiliary normal velocity are both
    ! reflected so qp and q remain consistent.
    IF ( radial_source_flag ) THEN
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
             this%eta_interfaceR(j+1,k) = this%eta_interfaceL(j+1,k)
          ELSEIF ( sourceW(j,k) ) THEN
             this%q_interfaceL(:,j,k) = this%q_interfaceR(:,j,k)
             this%q_interfaceL(2,j,k) = -this%q_interfaceL(2,j,k)
             this%qp_interfaceL(:,j,k) = this%qp_interfaceR(:,j,k)
             this%qp_interfaceL(2,j,k) = -this%qp_interfaceL(2,j,k)
             this%qp_interfaceL(idx_u,j,k) = -this%qp_interfaceL(idx_u,j,k)
             this%eta_interfaceL(j,k) = this%eta_interfaceR(j,k)
          END IF

          IF ( sourceN(j,k) ) THEN
             this%q_interfaceT(:,j,k+1) = this%q_interfaceB(:,j,k+1)
             this%q_interfaceT(3,j,k+1) = -this%q_interfaceT(3,j,k+1)
             this%qp_interfaceT(:,j,k+1) = this%qp_interfaceB(:,j,k+1)
             this%qp_interfaceT(3,j,k+1) = -this%qp_interfaceT(3,j,k+1)
             this%qp_interfaceT(idx_v,j,k+1) =                               &
                  -this%qp_interfaceT(idx_v,j,k+1)
             this%eta_interfaceT(j,k+1) = this%eta_interfaceB(j,k+1)
          ELSEIF ( sourceS(j,k) ) THEN
             this%q_interfaceB(:,j,k) = this%q_interfaceT(:,j,k)
             this%q_interfaceB(3,j,k) = -this%q_interfaceB(3,j,k)
             this%qp_interfaceB(:,j,k) = this%qp_interfaceT(:,j,k)
             this%qp_interfaceB(3,j,k) = -this%qp_interfaceB(3,j,k)
             this%qp_interfaceB(idx_v,j,k) = -this%qp_interfaceB(idx_v,j,k)
             this%eta_interfaceB(j,k) = this%eta_interfaceT(j,k)
          END IF
       END DO
    END IF

  CONTAINS

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

    SUBROUTINE compute_cell_weights(jc,kc)

    INTEGER, INTENT(IN) :: jc, kc
    REAL(wp) :: Eh, Eeta, h_minus_eta, h_plus_eta

    Eh = 0.0_wp
    Eeta = 0.0_wp
    h_minus_eta = MAX( this%eta_cellW(jc,kc) - B_face_x(jc,kc), 0.0_wp )
    h_plus_eta = MAX( this%eta_cellE(jc,kc) - B_face_x(jc+1,kc), 0.0_wp )

    IF ( jc .GT. 1 ) THEN
       Eh = Eh + ABS( MAX(this%qp_cellE(1,jc-1,kc),0.0_wp) -              &
          MAX(this%qp_cellW(1,jc,kc),0.0_wp) )
       Eeta = Eeta + ABS(                                               &
          MAX(this%eta_cellE(jc-1,kc)-B_face_x(jc,kc),0.0_wp) -       &
          MAX(this%eta_cellW(jc,kc)-B_face_x(jc,kc),0.0_wp) )
    END IF
    IF ( jc .LT. comp_cells_x ) THEN
       Eh = Eh + ABS( MAX(this%qp_cellE(1,jc,kc),0.0_wp) -              &
          MAX(this%qp_cellW(1,jc+1,kc),0.0_wp) )
       Eeta = Eeta + ABS(                                               &
          MAX(this%eta_cellE(jc,kc)-B_face_x(jc+1,kc),0.0_wp) -       &
          MAX(this%eta_cellW(jc+1,kc)-B_face_x(jc+1,kc),0.0_wp) )
    END IF
    this%w_eta_x(jc,kc) = continuity_weight(                             &
       qp_center(1,jc,kc), Eh, Eeta, h_minus_eta, h_plus_eta,          &
       this%hydrostatic_residual_2d(jc,kc),                            &
       this%topographic_relief_ratio_2d(jc,kc) )

    Eh = 0.0_wp
    Eeta = 0.0_wp
    h_minus_eta = MAX( this%eta_cellS(jc,kc) - B_face_y(jc,kc), 0.0_wp )
    h_plus_eta = MAX( this%eta_cellN(jc,kc) - B_face_y(jc,kc+1), 0.0_wp )

    IF ( kc .GT. 1 ) THEN
       Eh = Eh + ABS( MAX(this%qp_cellN(1,jc,kc-1),0.0_wp) -            &
          MAX(this%qp_cellS(1,jc,kc),0.0_wp) )
       Eeta = Eeta + ABS(                                               &
          MAX(this%eta_cellN(jc,kc-1)-B_face_y(jc,kc),0.0_wp) -       &
          MAX(this%eta_cellS(jc,kc)-B_face_y(jc,kc),0.0_wp) )
    END IF
    IF ( kc .LT. comp_cells_y ) THEN
       Eh = Eh + ABS( MAX(this%qp_cellN(1,jc,kc),0.0_wp) -              &
          MAX(this%qp_cellS(1,jc,kc+1),0.0_wp) )
       Eeta = Eeta + ABS(                                               &
          MAX(this%eta_cellN(jc,kc)-B_face_y(jc,kc+1),0.0_wp) -       &
          MAX(this%eta_cellS(jc,kc+1)-B_face_y(jc,kc+1),0.0_wp) )
    END IF
    this%w_eta_y(jc,kc) = continuity_weight(                             &
       qp_center(1,jc,kc), Eh, Eeta, h_minus_eta, h_plus_eta,          &
       this%hydrostatic_residual_2d(jc,kc),                            &
       this%topographic_relief_ratio_2d(jc,kc) )

    END SUBROUTINE compute_cell_weights

    FUNCTION continuity_weight(h_center,Eh,Eeta,h_minus_eta,h_plus_eta,    &
       hydrostatic_residual,topographic_relief_ratio) RESULT(weight)

    REAL(wp), INTENT(IN) :: h_center, Eh, Eeta
    REAL(wp), INTENT(IN) :: h_minus_eta, h_plus_eta
    REAL(wp), INTENT(IN) :: hydrostatic_residual, topographic_relief_ratio
    REAL(wp) :: weight, denominator, distribution_tolerance

    denominator = Eh + Eeta
    distribution_tolerance = 1.0E-14_wp * MAX(1.0_wp,h_center)
    IF ( denominator .GT. distribution_tolerance ) THEN
       weight = Eh / denominator
    ELSE
       weight = 0.5_wp
    END IF

    IF ( h_center .LE. hp_dry_tolerance ) THEN
       weight = 1.0_wp
    ELSEIF ( ( h_minus_eta .LE. hp_dry_tolerance ) .OR.                   &
       ( h_plus_eta .LE. hp_dry_tolerance ) ) THEN
       IF ( hydrostatic_residual .GT. hp_dynamic_residual_threshold ) THEN
        weight = 0.0_wp
       ELSE
        weight = 1.0_wp
       END IF
    ELSEIF ( ( hydrostatic_residual .GT. hp_dynamic_residual_threshold ) .AND. &
       ( topographic_relief_ratio .GT. 1.0_wp ) ) THEN
       weight = 0.0_wp
    END IF

    END FUNCTION continuity_weight

       SUBROUTINE evaluate_hp_cell_x(jc,kc)

       INTEGER, INTENT(IN) :: jc, kc
       REAL(wp) :: h_line(5), u_line(5), Bm_line(5), Bp_line(5)
       REAL(wp) :: hm_direct(5), hp_direct(5), hum_direct(5), hup_direct(5)
       REAL(wp) :: um_candidate(5), up_candidate(5)
       REAL(wp) :: residual_line(5), relief_line(5)
       REAL(wp) :: hm_line(5), hp_line(5), hum_line(5), hup_line(5)
       REAL(wp) :: eta_m_line(5), eta_p_line(5), weight_line(5)
       INTEGER :: j_start, j_end, line_size, center, offset, jj, thread_id

       j_start = MAX(1,jc-2)
       j_end = MIN(comp_cells_x,jc+2)
       line_size = j_end-j_start+1
       center = jc-j_start+1
       residual_line = 0.0_wp
       relief_line = 0.0_wp

       DO offset = 1, line_size
          jj = j_start+offset-1
          h_line(offset) = qp_center(1,jj,kc)
          u_line(offset) = qp_center(idx_u,jj,kc)
          Bm_line(offset) = B_face_x(jj,kc)
          Bp_line(offset) = B_face_x(jj+1,kc)
          IF ( this%hp_eta_mask(jj,kc) ) THEN
           hm_direct(offset) = this%qp_cellW(1,jj,kc)
           hp_direct(offset) = this%qp_cellE(1,jj,kc)
           hum_direct(offset) = this%qp_cellW(2,jj,kc)
           hup_direct(offset) = this%qp_cellE(2,jj,kc)
           um_candidate(offset) = this%qp_cellW(idx_u,jj,kc)
           up_candidate(offset) = this%qp_cellE(idx_u,jj,kc)
          ELSE
           hm_direct(offset) = qp_center(1,jj,kc)
           hp_direct(offset) = qp_center(1,jj,kc)
           hum_direct(offset) = qp_center(2,jj,kc)
           hup_direct(offset) = qp_center(2,jj,kc)
           um_candidate(offset) = qp_center(idx_u,jj,kc)
           up_candidate(offset) = qp_center(idx_u,jj,kc)
          END IF
       END DO
       residual_line(center) = this%hydrostatic_residual_2d(jc,kc)
       relief_line(center) = this%topographic_relief_ratio_2d(jc,kc)
       thread_id = omp_get_thread_num()+1

       CALL reconstruct_hp_line( h_line(1:line_size),u_line(1:line_size),    &
          Bm_line(1:line_size),Bp_line(1:line_size),                       &
          hm_direct(1:line_size),hp_direct(1:line_size),                   &
          hum_direct(1:line_size),hup_direct(1:line_size),                 &
          um_candidate(1:line_size),up_candidate(1:line_size),             &
          residual_line(1:line_size),relief_line(1:line_size),             &
          limiter(1),reconstr_coeff,hm_line(1:line_size),                 &
          hp_line(1:line_size),hum_line(1:line_size),hup_line(1:line_size),&
          eta_m_line(1:line_size),eta_p_line(1:line_size),                 &
          weight_line(1:line_size),this%hp_scratch(1:line_size,:,thread_id) )

       this%hp_blended(:,jc,kc,1) = [ hm_line(center),hp_line(center),      &
          hum_line(center),hup_line(center) ]

       END SUBROUTINE evaluate_hp_cell_x

       SUBROUTINE evaluate_hp_cell_y(jc,kc)

       INTEGER, INTENT(IN) :: jc, kc
       REAL(wp) :: h_line(5), u_line(5), Bm_line(5), Bp_line(5)
       REAL(wp) :: hm_direct(5), hp_direct(5), hum_direct(5), hup_direct(5)
       REAL(wp) :: um_candidate(5), up_candidate(5)
       REAL(wp) :: residual_line(5), relief_line(5)
       REAL(wp) :: hm_line(5), hp_line(5), hum_line(5), hup_line(5)
       REAL(wp) :: eta_m_line(5), eta_p_line(5), weight_line(5)
       INTEGER :: k_start, k_end, line_size, center, offset, kk, thread_id

       k_start = MAX(1,kc-2)
       k_end = MIN(comp_cells_y,kc+2)
       line_size = k_end-k_start+1
       center = kc-k_start+1
       residual_line = 0.0_wp
       relief_line = 0.0_wp

       DO offset = 1, line_size
          kk = k_start+offset-1
          h_line(offset) = qp_center(1,jc,kk)
          u_line(offset) = qp_center(idx_v,jc,kk)
          Bm_line(offset) = B_face_y(jc,kk)
          Bp_line(offset) = B_face_y(jc,kk+1)
          IF ( this%hp_eta_mask(jc,kk) ) THEN
           hm_direct(offset) = this%qp_cellS(1,jc,kk)
           hp_direct(offset) = this%qp_cellN(1,jc,kk)
           hum_direct(offset) = this%qp_cellS(3,jc,kk)
           hup_direct(offset) = this%qp_cellN(3,jc,kk)
           um_candidate(offset) = this%qp_cellS(idx_v,jc,kk)
           up_candidate(offset) = this%qp_cellN(idx_v,jc,kk)
          ELSE
           hm_direct(offset) = qp_center(1,jc,kk)
           hp_direct(offset) = qp_center(1,jc,kk)
           hum_direct(offset) = qp_center(3,jc,kk)
           hup_direct(offset) = qp_center(3,jc,kk)
           um_candidate(offset) = qp_center(idx_v,jc,kk)
           up_candidate(offset) = qp_center(idx_v,jc,kk)
          END IF
       END DO
       residual_line(center) = this%hydrostatic_residual_2d(jc,kc)
       relief_line(center) = this%topographic_relief_ratio_2d(jc,kc)
       thread_id = omp_get_thread_num()+1

       CALL reconstruct_hp_line( h_line(1:line_size),u_line(1:line_size),    &
          Bm_line(1:line_size),Bp_line(1:line_size),                       &
          hm_direct(1:line_size),hp_direct(1:line_size),                   &
          hum_direct(1:line_size),hup_direct(1:line_size),                 &
          um_candidate(1:line_size),up_candidate(1:line_size),             &
          residual_line(1:line_size),relief_line(1:line_size),             &
          limiter(1),reconstr_coeff,hm_line(1:line_size),                 &
          hp_line(1:line_size),hum_line(1:line_size),hup_line(1:line_size),&
          eta_m_line(1:line_size),eta_p_line(1:line_size),                 &
          weight_line(1:line_size),this%hp_scratch(1:line_size,:,thread_id) )

       this%hp_blended(:,jc,kc,2) = [ hm_line(center),hp_line(center),      &
          hum_line(center),hup_line(center) ]

       END SUBROUTINE evaluate_hp_cell_y

       SUBROUTINE commit_hp_cell(jc,kc)

       INTEGER, INTENT(IN) :: jc, kc

       this%qp_cellW(1,jc,kc) = this%hp_blended(1,jc,kc,1)
       this%qp_cellE(1,jc,kc) = this%hp_blended(2,jc,kc,1)
       this%qp_cellW(2,jc,kc) = this%hp_blended(3,jc,kc,1)
       this%qp_cellE(2,jc,kc) = this%hp_blended(4,jc,kc,1)
       this%qp_cellW(3,jc,kc) = this%qp_cellW(1,jc,kc)*this%qp_cellW(idx_v,jc,kc)
       this%qp_cellE(3,jc,kc) = this%qp_cellE(1,jc,kc)*this%qp_cellE(idx_v,jc,kc)
       this%qp_cellW(idx_u,jc,kc) = safe_velocity(this%qp_cellW(2,jc,kc),   &
          this%qp_cellW(1,jc,kc))
       this%qp_cellE(idx_u,jc,kc) = safe_velocity(this%qp_cellE(2,jc,kc),   &
          this%qp_cellE(1,jc,kc))

       this%qp_cellS(1,jc,kc) = this%hp_blended(1,jc,kc,2)
       this%qp_cellN(1,jc,kc) = this%hp_blended(2,jc,kc,2)
       this%qp_cellS(3,jc,kc) = this%hp_blended(3,jc,kc,2)
       this%qp_cellN(3,jc,kc) = this%hp_blended(4,jc,kc,2)
       this%qp_cellS(2,jc,kc) = this%qp_cellS(1,jc,kc)*this%qp_cellS(idx_u,jc,kc)
       this%qp_cellN(2,jc,kc) = this%qp_cellN(1,jc,kc)*this%qp_cellN(idx_u,jc,kc)
       this%qp_cellS(idx_v,jc,kc) = safe_velocity(this%qp_cellS(3,jc,kc),   &
          this%qp_cellS(1,jc,kc))
       this%qp_cellN(idx_v,jc,kc) = safe_velocity(this%qp_cellN(3,jc,kc),   &
          this%qp_cellN(1,jc,kc))

       END SUBROUTINE commit_hp_cell

    FUNCTION local_hydrostatic_residual(jc,kc) RESULT(residual)

      INTEGER, INTENT(IN) :: jc, kc
      REAL(wp) :: residual
      REAL(wp) :: h0, eta0, numerator, denominator, scale
      REAL(wp) :: hn, etan
      INTEGER :: jj, kk

      h0 = qp_center(1,jc,kc)
      eta0 = h0 + B_cent(jc,kc)
      numerator = 0.0_wp
      denominator = 0.0_wp

      IF ( jc .GT. 1 ) THEN
         jj = jc-1; kk = kc
         hn = qp_center(1,jj,kk); etan = hn + B_cent(jj,kk)
         numerator = numerator + ABS(etan-eta0)
         denominator = denominator + ABS(hn-h0) + ABS(B_cent(jj,kk)-B_cent(jc,kc))
      END IF
      IF ( jc .LT. comp_cells_x ) THEN
         jj = jc+1; kk = kc
         hn = qp_center(1,jj,kk); etan = hn + B_cent(jj,kk)
         numerator = numerator + ABS(etan-eta0)
         denominator = denominator + ABS(hn-h0) + ABS(B_cent(jj,kk)-B_cent(jc,kc))
      END IF
      IF ( kc .GT. 1 ) THEN
         jj = jc; kk = kc-1
         hn = qp_center(1,jj,kk); etan = hn + B_cent(jj,kk)
         numerator = numerator + ABS(etan-eta0)
         denominator = denominator + ABS(hn-h0) + ABS(B_cent(jj,kk)-B_cent(jc,kc))
      END IF
      IF ( kc .LT. comp_cells_y ) THEN
         jj = jc; kk = kc+1
         hn = qp_center(1,jj,kk); etan = hn + B_cent(jj,kk)
         numerator = numerator + ABS(etan-eta0)
         denominator = denominator + ABS(hn-h0) + ABS(B_cent(jj,kk)-B_cent(jc,kc))
      END IF

      scale = 128.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,ABS(eta0),ABS(h0),ABS(B_cent(jc,kc)))
      IF ( denominator .GT. scale ) THEN
         residual = MIN(1.0_wp,MAX(0.0_wp,numerator/denominator))
      ELSE
         residual = 0.0_wp
      END IF

    END FUNCTION local_hydrostatic_residual

   PURE ELEMENTAL FUNCTION safe_velocity(momentum,thickness) RESULT(velocity)

      REAL(wp), INTENT(IN) :: momentum, thickness
      REAL(wp) :: velocity

      IF ( thickness .GT. hp_dry_tolerance ) THEN
         velocity = momentum / thickness
      ELSE
         velocity = 0.0_wp
      END IF

    END FUNCTION safe_velocity

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
      CALL velocity_from_conservative(q_conservative,                        &
           qp_reconstructed(idx_u),qp_reconstructed(idx_v))

    END SUBROUTINE final_hp_state

  END SUBROUTINE apply_hp_reconstruction

END MODULE reconstruction_2d
