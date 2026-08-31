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
  USE geometry_2d, ONLY : B_cent
  USE geometry_2d, ONLY : source_cell
  USE geometry_2d, ONLY : one_by_dx, one_by_dy
  USE geometry_2d, ONLY : limit

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: initialize_reconstruction
  PUBLIC :: finalize_reconstruction
  PUBLIC :: reconstruction

  PUBLIC :: q_interfaceL, q_interfaceR, q_interfaceB, q_interfaceT
  PUBLIC :: qp_interfaceL, qp_interfaceR, qp_interfaceB, qp_interfaceT
  PUBLIC :: diverg_interfaceL, diverg_interfaceR
  PUBLIC :: diverg_interfaceB, diverg_interfaceT

  REAL(wp), ALLOCATABLE :: q_interfaceL(:,:,:)
  REAL(wp), ALLOCATABLE :: q_interfaceR(:,:,:)
  REAL(wp), ALLOCATABLE :: q_interfaceB(:,:,:)
  REAL(wp), ALLOCATABLE :: q_interfaceT(:,:,:)

  REAL(wp), ALLOCATABLE :: qp_interfaceL(:,:,:)
  REAL(wp), ALLOCATABLE :: qp_interfaceR(:,:,:)
  REAL(wp), ALLOCATABLE :: qp_interfaceB(:,:,:)
  REAL(wp), ALLOCATABLE :: qp_interfaceT(:,:,:)

  LOGICAL, ALLOCATABLE :: diverg_interfaceL(:,:)
  LOGICAL, ALLOCATABLE :: diverg_interfaceR(:,:)
  LOGICAL, ALLOCATABLE :: diverg_interfaceB(:,:)
  LOGICAL, ALLOCATABLE :: diverg_interfaceT(:,:)

CONTAINS

  SUBROUTINE initialize_reconstruction

    ALLOCATE( q_interfaceL( n_vars, comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( q_interfaceR( n_vars, comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( q_interfaceB( n_vars, comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( q_interfaceT( n_vars, comp_cells_x, comp_interfaces_y ) )

    ALLOCATE( qp_interfaceL( n_vars+2, comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( qp_interfaceR( n_vars+2, comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( qp_interfaceB( n_vars+2, comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( qp_interfaceT( n_vars+2, comp_cells_x, comp_interfaces_y ) )

    ALLOCATE( diverg_interfaceL( comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( diverg_interfaceR( comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( diverg_interfaceB( comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( diverg_interfaceT( comp_cells_x, comp_interfaces_y ) )

  END SUBROUTINE initialize_reconstruction

  SUBROUTINE finalize_reconstruction

    DEALLOCATE( q_interfaceL )
    DEALLOCATE( q_interfaceR )
    DEALLOCATE( q_interfaceB )
    DEALLOCATE( q_interfaceT )

    DEALLOCATE( qp_interfaceL )
    DEALLOCATE( qp_interfaceR )
    DEALLOCATE( qp_interfaceB )
    DEALLOCATE( qp_interfaceT )

    DEALLOCATE( diverg_interfaceL )
    DEALLOCATE( diverg_interfaceR )
    DEALLOCATE( diverg_interfaceB )
    DEALLOCATE( diverg_interfaceT )

  END SUBROUTINE finalize_reconstruction

  SUBROUTINE reconstruction( q_expl, qp_expl, t, solve_cells, j_cent, k_cent )

    ! External procedures
    USE constitutive_2d, ONLY : qp_to_qc, qp_to_qp2
    USE constitutive_2d, ONLY : eval_source_bdry
    USE parameters_2d, ONLY : limiter

    USE geometry_2d, ONLY : x_comp , x_stag , y_comp , y_stag , dx2 , dy2

    USE geometry_2d, ONLY : sourceW , sourceE , sourceN , sourceS
    USE geometry_2d, ONLY : sourceW_vect_x , sourceW_vect_y
    USE geometry_2d, ONLY : sourceE_vect_x , sourceE_vect_y
    USE geometry_2d, ONLY : sourceN_vect_x , sourceN_vect_y
    USE geometry_2d, ONLY : sourceS_vect_x , sourceS_vect_y

    USE parameters_2d, ONLY : alpha_flag
    USE parameters_2d, ONLY : reconstr_coeff

    IMPLICIT NONE

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

                CALL qp_to_qp2( qrecW(1:n_vars) , B_cent(j,k) , qp2recW )
                qrec_stencil(1) = qp2recW(i-n_vars+1)
                qrec_stencil(3) = qp_expl(i,j+1,k)

             ELSEIF ( j .EQ. comp_cells_x ) THEN

                CALL qp_to_qp2( qrecE(1:n_vars) , B_cent(j,k) , qp2recE )
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

                CALL qp_to_qp2( qrecW(1:n_vars) , B_cent(j,k) , qp2recW )
                qrecW(i) = qp2recW(i-n_vars+1)

             ELSEIF ( j .EQ. comp_cells_x ) THEN

                CALL qp_to_qp2( qrecE(1:n_vars) , B_cent(j,k) , qp2recE )
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

                CALL qp_to_qp2( qrecS(1:n_vars) , B_cent(j,k) , qp2recS )
                qrec_stencil(1) = qp2recS(i-n_vars+1)
                qrec_stencil(3) = qp_expl(i,j,k+1)

             ELSEIF ( k .EQ. comp_cells_y ) THEN

                CALL qp_to_qp2( qrecN(1:n_vars) , B_cent(j,k) , qp2recN )
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

                CALL qp_to_qp2( qrecS(1:n_vars) , B_cent(j,k) , qp2recS )
                qrecS(i) = qp2recS(i-n_vars+1)

             ELSEIF ( k .EQ. comp_cells_y ) THEN

                CALL qp_to_qp2( qrecN(1:n_vars) , B_cent(j,k) , qp2recN )
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

       IF ( comp_cells_x .GT. 1 ) THEN

          IF ( ( j .GT. 1 ) .AND. ( j .LT. comp_cells_x ) ) THEN

             IF ( q_expl(1,j,k) .EQ. 0.0_wp ) THEN

                IF ( ( .NOT. radial_source_flag ) .OR.                          &
                     ( ( radial_source_flag ) .AND.                             &
                     ( source_cell(j,k) .EQ. 0 ) ) ) THEN

                   ! In the internal cell, if thickness h is 0 at the center
                   ! of the cell, then all the variables are 0 at the center
                   ! and at the interfaces (no conversion back is needed from
                   ! reconstructed to conservative)
                   q_interfaceR(:,j,k) = 0.0_wp
                   q_interfaceL(:,j+1,k) = 0.0_wp

                   qp_interfaceR(1:3,j,k) = 0.0_wp
                   qp_interfaceR(4:n_vars,j,k) = qrecW(4:n_vars)
                   qp_interfaceR(n_vars+1:n_vars+2,j,k) = 0.0_wp

                   qp_interfaceL(1:3,j+1,k) = 0.0_wp
                   qp_interfaceL(4:n_vars,j+1,k) = qrecE(4:n_vars)
                   qp_interfaceL(n_vars+1:n_vars+2,j+1,k) = 0.0_wp

                   diverg_interfaceR(j,k) = .FALSE.
                   diverg_interfaceL(j+1,k) = .FALSE.

                END IF

             END IF

          END IF

          ! Correction for residual volume fraction of continuous phase
          IF ( alpha_flag ) THEN

             !qrecW(5:4+n_solid) = qrecW(5:4+n_solid) *                          &
             !     MIN( 1.0_wp , maximum_solid_packing /                         &
             !     SUM( qrecW(5:4+n_solid) ) )

             !qrecE(5:4+n_solid) = qrecE(5:4+n_solid) *                          &
             !     MIN( 1.0_wp , maximum_solid_packing /                         &
             !     SUM( qrecE(5:4+n_solid) ) )

          ELSE

             !qrecW(5:4+n_solid) = qrecW(5:4+n_solid) *                          &
             !     MIN( 1.0_wp , maximum_solid_packing * qrecW(1) /              &
             !     SUM( qrecW(5:4+n_solid) ) )

             !qrecE(5:4+n_solid) = qrecE(5:4+n_solid) *                          &
             !     MIN( 1.0_wp , maximum_solid_packing * qrecE(1) /              &
             !     SUM( qrecE(5:4+n_solid) ) )

          END IF

          CALL qp_to_qc( qrecW,q_interfaceR(:,j,k) )
          CALL qp_to_qc( qrecE,q_interfaceL(:,j+1,k) )

          qp_interfaceR(1:n_vars+2,j,k) = qrecW(1:n_vars+2)
          qp_interfaceL(1:n_vars+2,j+1,k) = qrecE(1:n_vars+2)

          diverg_interfaceR(j,k) = diverging_flag
          diverg_interfaceL(j+1,k) = diverging_flag

          IF ( j.EQ.1 ) THEN

             ! Interface value at the left of first x-interface (external)
             q_interfaceL(:,j,k) = q_interfaceR(:,j,k)
             qp_interfaceL(:,j,k) = qp_interfaceR(:,j,k)

             !WRITE(*,*) 'j,k',j,k
             !WRITE(*,*) 'qp_interfaceL(:,j,k)',qp_interfaceL(:,j,k)
             !READ(*,*)

             diverg_interfaceR(j,k) = diverg_interfaceL(j,k)

          ELSEIF ( j.EQ.comp_cells_x ) THEN

             ! Interface value at the right of last x-interface (external)
             q_interfaceR(:,j+1,k) = q_interfaceL(:,j+1,k)
             qp_interfaceR(:,j+1,k) = qp_interfaceL(:,j+1,k)

             diverg_interfaceR(j+1,k) = diverg_interfaceL(j+1,k)

          ELSE

             IF ( radial_source_flag .AND. ( source_cell(j,k) .EQ. 2 ) ) THEN

                IF ( sourceE(j,k) ) THEN

                   q_interfaceR(:,j+1,k) = q_interfaceL(:,j+1,k)
                   q_interfaceR(2,j+1,k) = -q_interfaceL(2,j+1,k)
                   qp_interfaceR(:,j+1,k) = qp_interfaceL(:,j+1,k)
                   qp_interfaceR(idx_u,j+1,k) = -qp_interfaceL(idx_u,j+1,k)

                ELSEIF ( sourceW(j,k) ) THEN

                   q_interfaceL(:,j,k) = q_interfaceR(:,j,k)
                   q_interfaceL(2,j,k) = -q_interfaceR(2,j,k)
                   qp_interfaceL(:,j,k) = qp_interfaceR(:,j,k)
                   qp_interfaceL(idx_u,j,k) = -qp_interfaceR(idx_u,j,k)

                END IF

             END IF

          END IF

       ELSE

          ! for case comp_cells_x = 1
          q_interfaceR(1:n_vars,j,k) = q_expl(1:n_vars,j,k)
          q_interfaceL(1:n_vars,j+1,k) = q_expl(1:n_vars,j,k)

          qp_interfaceR(1:n_vars+2,j,k) = qp_expl(1:n_vars+2,j,k)
          qp_interfaceL(1:n_vars+2,j+1,k) = qp_expl(1:n_vars+2,j,k)

          diverg_interfaceR(j,k) = diverging_flag
          diverg_interfaceL(j+1,k) = diverging_flag

       END IF

       IF ( comp_cells_y .GT. 1 ) THEN

          IF ( ( k .GT. 1 ) .AND. ( k .LT. comp_cells_y ) ) THEN

             IF ( q_expl(1,j,k) .EQ. 0.0_wp ) THEN

                IF ( ( .NOT. radial_source_flag ) .OR.                          &
                     ( ( radial_source_flag ) .AND.                             &
                     ( source_cell(j,k) .EQ. 0 ) ) ) THEN

                   ! In the internal cell, if thickness h is 0 at the center
                   ! of the cell, then all the variables are 0 at the center
                   ! and at the interfaces (no conversion back is needed from
                   ! reconstructed to conservative)

                   q_interfaceT(:,j,k) = 0.0_wp
                   q_interfaceB(:,j,k+1) = 0.0_wp

                   qp_interfaceT(1:3,j,k) = 0.0_wp
                   qp_interfaceT(4:n_vars,j,k) = qrecS(4:n_vars)
                   qp_interfaceT(n_vars+1:n_vars+2,j,k) = 0.0_wp

                   qp_interfaceB(1:3,j,k+1) = 0.0_wp
                   qp_interfaceB(4:n_vars,j,k+1) = qrecN(4:n_vars)
                   qp_interfaceB(n_vars+1:n_vars+2,j,k+1) = 0.0_wp

                   diverg_interfaceT(j,k) = .FALSE.
                   diverg_interfaceB(j,k+1) = .FALSE.

                END IF

             END IF

          END IF

          ! Correction for maximum solid packing
          IF ( alpha_flag ) THEN

             !qrecS(5:4+n_solid) = qrecS(5:4+n_solid) *                          &
             !     MIN( 1.0_wp , maximum_solid_packing /                         &
             !     SUM( qrecS(5:4+n_solid) ) )

             !qrecN(5:4+n_solid) = qrecN(5:4+n_solid) *                          &
             !     MIN( 1.0_wp , maximum_solid_packing /                         &
             !     SUM( qrecN(5:4+n_solid) ) )

          ELSE

             !qrecS(5:4+n_solid) = qrecS(5:4+n_solid) *                          &
             !     MIN( 1.0_wp , maximum_solid_packing * qrecS(1) /              &
             !     SUM( qrecS(5:4+n_solid) ) )

             !qrecN(5:4+n_solid) = qrecN(5:4+n_solid) *                          &
             !     MIN( 1.0_wp , maximum_solid_packing * qrecN(1) /              &
             !     SUM( qrecN(5:4+n_solid) ) )

          END IF

          CALL qp_to_qc( qrecS, q_interfaceT(:,j,k) )
          CALL qp_to_qc( qrecN, q_interfaceB(:,j,k+1) )

          qp_interfaceT(1:n_vars+2,j,k) = qrecS(1:n_vars+2)
          qp_interfaceB(1:n_vars+2,j,k+1) = qrecN(1:n_vars+2)

          diverg_interfaceT(j,k) = diverging_flag
          diverg_interfaceB(j,k+1) = diverging_flag

          IF ( k .EQ. 1 ) THEN

             ! Interface value at the bottom of first y-interface (external)
             q_interfaceB(:,j,k) = q_interfaceT(:,j,k)
             qp_interfaceB(:,j,k) = qp_interfaceT(:,j,k)

             diverg_interfaceB(j,k) = diverg_interfaceT(j,k)

          ELSEIF ( k .EQ. comp_cells_y ) THEN

             ! Interface value at the top of last y-interface (external)
             q_interfaceT(:,j,k+1) = q_interfaceB(:,j,k+1)
             qp_interfaceT(:,j,k+1) = qp_interfaceB(:,j,k+1)

             diverg_interfaceT(j,k+1) = diverg_interfaceB(j,k+1)

          ELSE

             IF ( radial_source_flag .AND. ( source_cell(j,k) .EQ. 2 ) ) THEN

                IF ( sourceS(j,k) ) THEN

                   q_interfaceB(:,j,k) = q_interfaceT(:,j,k)
                   q_interfaceB(3,j,k) = -q_interfaceT(3,j,k)
                   qp_interfaceB(:,j,k) = qp_interfaceT(:,j,k)
                   qp_interfaceB(idx_v,j,k) = -qp_interfaceT(idx_v,j,k)

                ELSEIF ( sourceN(j,k) ) THEN

                   q_interfaceT(:,j,k+1) = q_interfaceB(:,j,k+1)
                   q_interfaceT(3,j,k+1) = -q_interfaceB(3,j,k+1)
                   qp_interfaceT(:,j,k+1) = qp_interfaceB(:,j,k+1)
                   qp_interfaceT(idx_v,j,k+1) = -qp_interfaceB(idx_v,j,k+1)

                END IF

             END IF

          END IF

       ELSE

          ! case comp_cells_y = 1

          q_interfaceB(:,j,k) = q_expl(:,j,k)
          q_interfaceT(:,j,k) = q_expl(:,j,k)
          q_interfaceB(:,j,k+1) = q_expl(:,j,k)
          q_interfaceT(:,j,k+1) = q_expl(:,j,k)

          qp_interfaceB(:,j,k) = qp_expl(:,j,k)
          qp_interfaceT(:,j,k) = qp_expl(:,j,k)
          qp_interfaceB(:,j,k+1) = qp_expl(:,j,k)
          qp_interfaceT(:,j,k+1) = qp_expl(:,j,k)

       END IF

    END DO

    !$OMP END PARALLEL DO

    RETURN

  END SUBROUTINE reconstruction

END MODULE reconstruction_2d
