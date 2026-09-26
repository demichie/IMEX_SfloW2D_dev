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

  USE hp_reconstruction_2d, ONLY : reconstruct_hp_line
  USE hp_reconstruction_2d, ONLY : hp_dry_tolerance
   USE omp_lib, ONLY : omp_get_max_threads, omp_get_thread_num

  IMPLICIT NONE

  PRIVATE

   INTEGER, PARAMETER :: hp_h_center = 1, hp_u_center = 2
   INTEGER, PARAMETER :: hp_B_minus = 3, hp_B_plus = 4
   INTEGER, PARAMETER :: hp_h_minus_direct = 5, hp_h_plus_direct = 6
   INTEGER, PARAMETER :: hp_hu_minus_direct = 7, hp_hu_plus_direct = 8
   INTEGER, PARAMETER :: hp_u_minus_candidate = 9, hp_u_plus_candidate = 10
   INTEGER, PARAMETER :: hp_h_minus = 11, hp_h_plus = 12
   INTEGER, PARAMETER :: hp_hu_minus = 13, hp_hu_plus = 14
   INTEGER, PARAMETER :: hp_eta_minus = 15, hp_eta_plus = 16, hp_weight = 17
   INTEGER, PARAMETER :: hp_scratch_first = 18, hp_scratch_last = 25

  TYPE, PUBLIC :: reconstruction_workspace_type
     REAL(wp), ALLOCATABLE :: q_interfaceL(:,:,:)
     REAL(wp), ALLOCATABLE :: q_interfaceR(:,:,:)
     REAL(wp), ALLOCATABLE :: q_interfaceB(:,:,:)
     REAL(wp), ALLOCATABLE :: q_interfaceT(:,:,:)

     REAL(wp), ALLOCATABLE :: qp_interfaceL(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_interfaceR(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_interfaceB(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_interfaceT(:,:,:)

     ! Direct cell-side candidates retained until the line-wise HP blend can
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
    ALLOCATE( this%hp_scratch(MAX(comp_cells_x,comp_cells_y),                &
         hp_scratch_last,MAX(1,omp_get_max_threads())) )

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
      DEALLOCATE( this%hp_scratch )

    DEALLOCATE( this%diverg_interfaceL )
    DEALLOCATE( this%diverg_interfaceR )
    DEALLOCATE( this%diverg_interfaceB )
    DEALLOCATE( this%diverg_interfaceT )

  END SUBROUTINE finalize_reconstruction

  SUBROUTINE reconstruction( this, q_expl, qp_expl, t, solve_cells, j_cent, k_cent )

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
    REAL(wp), INTENT(IN) :: q_expl(:,:,:)
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

    ! Every cell has a defined constant candidate.  Active cells overwrite it
    ! below with the ordinary limited reconstruction.  This also gives the HP
    ! continuity indicator a well-defined neighbour at solve-mask boundaries.
    this%qp_cellW = qp_expl
    this%qp_cellE = qp_expl
    this%qp_cellS = qp_expl
    this%qp_cellN = qp_expl

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
  !> The first reconstruction pass supplies direct limited candidates.  This
  !> second, line-wise pass can evaluate the neighbouring face mismatches used
  !> by the parameter-free h/eta blend.  It then rebuilds thermodynamics and
  !> conservative states from the final thickness and momenta.
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
   INTEGER :: j, k, l, line_size, thread_id

   !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(j,k,relief2d)
    DO k = 1, comp_cells_y
       DO j = 1, comp_cells_x
          this%hydrostatic_residual_2d(j,k) = local_hydrostatic_residual(j,k)
          relief2d = MAX( B_face_x(j,k), B_face_x(j+1,k),                    &
               B_face_y(j,k), B_face_y(j,k+1) ) -                           &
               MIN( B_face_x(j,k), B_face_x(j+1,k),                         &
               B_face_y(j,k), B_face_y(j,k+1) )
          this%topographic_relief_ratio_2d(j,k) = relief2d /                 &
               MAX(qp_center(1,j,k),hp_dry_tolerance)
       END DO
    END DO
    !$OMP END PARALLEL DO

    line_size = comp_cells_x
    !$OMP PARALLEL DO PRIVATE(thread_id)
    DO k = 1, comp_cells_y

      thread_id = omp_get_thread_num() + 1
      this%hp_scratch(1:line_size,hp_h_center,thread_id) = qp_center(1,:,k)
      this%hp_scratch(1:line_size,hp_u_center,thread_id) = qp_center(idx_u,:,k)
      this%hp_scratch(1:line_size,hp_B_minus,thread_id) = B_face_x(1:comp_cells_x,k)
      this%hp_scratch(1:line_size,hp_B_plus,thread_id) = B_face_x(2:comp_interfaces_x,k)
      this%hp_scratch(1:line_size,hp_h_minus_direct,thread_id) = this%qp_cellW(1,:,k)
      this%hp_scratch(1:line_size,hp_h_plus_direct,thread_id) = this%qp_cellE(1,:,k)
      this%hp_scratch(1:line_size,hp_hu_minus_direct,thread_id) = this%qp_cellW(2,:,k)
      this%hp_scratch(1:line_size,hp_hu_plus_direct,thread_id) = this%qp_cellE(2,:,k)
      this%hp_scratch(1:line_size,hp_u_minus_candidate,thread_id) = this%qp_cellW(idx_u,:,k)
      this%hp_scratch(1:line_size,hp_u_plus_candidate,thread_id) = this%qp_cellE(idx_u,:,k)

       CALL reconstruct_hp_line(                                                &
           this%hp_scratch(1:line_size,hp_h_center,thread_id),              &
           this%hp_scratch(1:line_size,hp_u_center,thread_id),               &
           this%hp_scratch(1:line_size,hp_B_minus,thread_id),                &
           this%hp_scratch(1:line_size,hp_B_plus,thread_id),                 &
           this%hp_scratch(1:line_size,hp_h_minus_direct,thread_id),         &
           this%hp_scratch(1:line_size,hp_h_plus_direct,thread_id),          &
           this%hp_scratch(1:line_size,hp_hu_minus_direct,thread_id),        &
           this%hp_scratch(1:line_size,hp_hu_plus_direct,thread_id),         &
           this%hp_scratch(1:line_size,hp_u_minus_candidate,thread_id),      &
           this%hp_scratch(1:line_size,hp_u_plus_candidate,thread_id),       &
           this%hydrostatic_residual_2d(:,k),                               &
           this%topographic_relief_ratio_2d(:,k), limiter(1),                &
           reconstr_coeff, this%hp_scratch(1:line_size,hp_h_minus,thread_id),&
           this%hp_scratch(1:line_size,hp_h_plus,thread_id),                 &
           this%hp_scratch(1:line_size,hp_hu_minus,thread_id),               &
           this%hp_scratch(1:line_size,hp_hu_plus,thread_id),                &
           this%hp_scratch(1:line_size,hp_eta_minus,thread_id),              &
           this%hp_scratch(1:line_size,hp_eta_plus,thread_id),               &
           this%hp_scratch(1:line_size,hp_weight,thread_id),                 &
           this%hp_scratch(1:line_size,hp_scratch_first:hp_scratch_last,     &
           thread_id) )

      this%eta_cellW(:,k) = this%hp_scratch(1:line_size,hp_eta_minus,thread_id)
      this%eta_cellE(:,k) = this%hp_scratch(1:line_size,hp_eta_plus,thread_id)
      this%w_eta_x(:,k) = this%hp_scratch(1:line_size,hp_weight,thread_id)

        this%qp_cellW(1,:,k) = this%hp_scratch(1:line_size,hp_h_minus,thread_id)
        this%qp_cellE(1,:,k) = this%hp_scratch(1:line_size,hp_h_plus,thread_id)
        this%qp_cellW(2,:,k) = this%hp_scratch(1:line_size,hp_hu_minus,thread_id)
        this%qp_cellE(2,:,k) = this%hp_scratch(1:line_size,hp_hu_plus,thread_id)
        this%qp_cellW(3,:,k) = this%hp_scratch(1:line_size,hp_h_minus,thread_id) &
           * this%qp_cellW(idx_v,:,k)
        this%qp_cellE(3,:,k) = this%hp_scratch(1:line_size,hp_h_plus,thread_id)  &
           * this%qp_cellE(idx_v,:,k)
        this%qp_cellW(idx_u,:,k) = safe_velocity(                            &
           this%hp_scratch(1:line_size,hp_hu_minus,thread_id),              &
           this%hp_scratch(1:line_size,hp_h_minus,thread_id))
        this%qp_cellE(idx_u,:,k) = safe_velocity(                            &
           this%hp_scratch(1:line_size,hp_hu_plus,thread_id),               &
           this%hp_scratch(1:line_size,hp_h_plus,thread_id))

    END DO
    !$OMP END PARALLEL DO

    line_size = comp_cells_y
    !$OMP PARALLEL DO PRIVATE(thread_id)
    DO j = 1, comp_cells_x

      thread_id = omp_get_thread_num() + 1
      this%hp_scratch(1:line_size,hp_h_center,thread_id) = qp_center(1,j,:)
      this%hp_scratch(1:line_size,hp_u_center,thread_id) = qp_center(idx_v,j,:)
      this%hp_scratch(1:line_size,hp_B_minus,thread_id) = B_face_y(j,1:comp_cells_y)
      this%hp_scratch(1:line_size,hp_B_plus,thread_id) = B_face_y(j,2:comp_interfaces_y)
      this%hp_scratch(1:line_size,hp_h_minus_direct,thread_id) = this%qp_cellS(1,j,:)
      this%hp_scratch(1:line_size,hp_h_plus_direct,thread_id) = this%qp_cellN(1,j,:)
      this%hp_scratch(1:line_size,hp_hu_minus_direct,thread_id) = this%qp_cellS(3,j,:)
      this%hp_scratch(1:line_size,hp_hu_plus_direct,thread_id) = this%qp_cellN(3,j,:)
      this%hp_scratch(1:line_size,hp_u_minus_candidate,thread_id) = this%qp_cellS(idx_v,j,:)
      this%hp_scratch(1:line_size,hp_u_plus_candidate,thread_id) = this%qp_cellN(idx_v,j,:)

       CALL reconstruct_hp_line(                                                &
           this%hp_scratch(1:line_size,hp_h_center,thread_id),              &
           this%hp_scratch(1:line_size,hp_u_center,thread_id),               &
           this%hp_scratch(1:line_size,hp_B_minus,thread_id),                &
           this%hp_scratch(1:line_size,hp_B_plus,thread_id),                 &
           this%hp_scratch(1:line_size,hp_h_minus_direct,thread_id),         &
           this%hp_scratch(1:line_size,hp_h_plus_direct,thread_id),          &
           this%hp_scratch(1:line_size,hp_hu_minus_direct,thread_id),        &
           this%hp_scratch(1:line_size,hp_hu_plus_direct,thread_id),         &
           this%hp_scratch(1:line_size,hp_u_minus_candidate,thread_id),      &
           this%hp_scratch(1:line_size,hp_u_plus_candidate,thread_id),       &
           this%hydrostatic_residual_2d(j,:),                               &
           this%topographic_relief_ratio_2d(j,:), limiter(1),                &
           reconstr_coeff, this%hp_scratch(1:line_size,hp_h_minus,thread_id),&
           this%hp_scratch(1:line_size,hp_h_plus,thread_id),                 &
           this%hp_scratch(1:line_size,hp_hu_minus,thread_id),               &
           this%hp_scratch(1:line_size,hp_hu_plus,thread_id),                &
           this%hp_scratch(1:line_size,hp_eta_minus,thread_id),              &
           this%hp_scratch(1:line_size,hp_eta_plus,thread_id),               &
           this%hp_scratch(1:line_size,hp_weight,thread_id),                 &
           this%hp_scratch(1:line_size,hp_scratch_first:hp_scratch_last,     &
           thread_id) )

      this%eta_cellS(j,:) = this%hp_scratch(1:line_size,hp_eta_minus,thread_id)
      this%eta_cellN(j,:) = this%hp_scratch(1:line_size,hp_eta_plus,thread_id)
      this%w_eta_y(j,:) = this%hp_scratch(1:line_size,hp_weight,thread_id)

        this%qp_cellS(1,j,:) = this%hp_scratch(1:line_size,hp_h_minus,thread_id)
        this%qp_cellN(1,j,:) = this%hp_scratch(1:line_size,hp_h_plus,thread_id)
        this%qp_cellS(2,j,:) = this%hp_scratch(1:line_size,hp_h_minus,thread_id) &
           * this%qp_cellS(idx_u,j,:)
        this%qp_cellN(2,j,:) = this%hp_scratch(1:line_size,hp_h_plus,thread_id)  &
           * this%qp_cellN(idx_u,j,:)
        this%qp_cellS(3,j,:) = this%hp_scratch(1:line_size,hp_hu_minus,thread_id)
        this%qp_cellN(3,j,:) = this%hp_scratch(1:line_size,hp_hu_plus,thread_id)
        this%qp_cellS(idx_v,j,:) = safe_velocity(                            &
           this%hp_scratch(1:line_size,hp_hu_minus,thread_id),              &
           this%hp_scratch(1:line_size,hp_h_minus,thread_id))
        this%qp_cellN(idx_v,j,:) = safe_velocity(                            &
           this%hp_scratch(1:line_size,hp_hu_plus,thread_id),               &
           this%hp_scratch(1:line_size,hp_h_plus,thread_id))

   END DO
   !$OMP END PARALLEL DO

    ! Map all eta cell traces to their oriented face storage.  Conservative
    ! and primitive states below are restricted to the active solve mask.
    !$OMP PARALLEL DO
    DO k = 1, comp_cells_y
       this%eta_interfaceR(1:comp_cells_x,k) = this%eta_cellW(:,k)
       this%eta_interfaceL(2:comp_interfaces_x,k) = this%eta_cellE(:,k)
       this%eta_interfaceL(1,k) = this%eta_interfaceR(1,k)
       this%eta_interfaceR(comp_interfaces_x,k) =                             &
            this%eta_interfaceL(comp_interfaces_x,k)
    END DO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO
    DO j = 1, comp_cells_x
       this%eta_interfaceT(j,1:comp_cells_y) = this%eta_cellS(j,:)
       this%eta_interfaceB(j,2:comp_interfaces_y) = this%eta_cellN(j,:)
       this%eta_interfaceB(j,1) = this%eta_interfaceT(j,1)
       this%eta_interfaceT(j,comp_interfaces_y) =                             &
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
