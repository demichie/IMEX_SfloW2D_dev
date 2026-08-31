!********************************************************************************
!> \brief Hyperbolic flux evaluation
!>
!> This module owns characteristic speeds and semidiscrete numerical interface
!> fluxes. Active-cell and active-interface lists are supplied explicitly by
!> the solver driver.
!********************************************************************************
MODULE hyperbolic_2d

  USE parameters_2d, ONLY : wp, n_eqns, n_vars
  USE parameters_2d, ONLY : idx_solidEqn_first, idx_solidEqn_last
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE geometry_2d, ONLY : comp_interfaces_x, comp_interfaces_y
  USE geometry_2d, ONLY : B_prime_x_geom, B_prime_y_geom
  USE geometry_2d, ONLY : grav_coeff_stag_x, grav_coeff_stag_y
  USE geometry_2d, ONLY : one_by_dx, one_by_dy

  USE reconstruction_2d, ONLY : reconstruction
  USE reconstruction_2d, ONLY : q_interfaceL, q_interfaceR
  USE reconstruction_2d, ONLY : q_interfaceB, q_interfaceT
  USE reconstruction_2d, ONLY : qp_interfaceL, qp_interfaceR
  USE reconstruction_2d, ONLY : qp_interfaceB, qp_interfaceT

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: initialize_hyperbolic
  PUBLIC :: finalize_hyperbolic
  PUBLIC :: eval_hyperbolic_terms
  PUBLIC :: eval_speeds

  ! Transitional public workspace: timestep CFL evaluation and existing solver
  ! diagnostics still consume these arrays directly.
  PUBLIC :: a_interface_xNeg, a_interface_xPos
  PUBLIC :: a_interface_yNeg, a_interface_yPos
  PUBLIC :: H_interface_x, H_interface_y

  REAL(wp), ALLOCATABLE :: a_interface_xNeg(:,:,:)
  REAL(wp), ALLOCATABLE :: a_interface_xPos(:,:,:)
  REAL(wp), ALLOCATABLE :: a_interface_yNeg(:,:,:)
  REAL(wp), ALLOCATABLE :: a_interface_yPos(:,:,:)
  REAL(wp), ALLOCATABLE :: H_interface_x(:,:,:)
  REAL(wp), ALLOCATABLE :: H_interface_y(:,:,:)

CONTAINS

  SUBROUTINE initialize_hyperbolic

    ALLOCATE( a_interface_xNeg(n_eqns,comp_interfaces_x,comp_cells_y) )
    ALLOCATE( a_interface_xPos(n_eqns,comp_interfaces_x,comp_cells_y) )
    ALLOCATE( a_interface_yNeg(n_eqns,comp_cells_x,comp_interfaces_y) )
    ALLOCATE( a_interface_yPos(n_eqns,comp_cells_x,comp_interfaces_y) )

    ALLOCATE( H_interface_x(n_eqns,comp_interfaces_x,comp_cells_y) )
    ALLOCATE( H_interface_y(n_eqns,comp_cells_x,comp_interfaces_y) )

    a_interface_xNeg = 0.0_wp
    a_interface_xPos = 0.0_wp
    a_interface_yNeg = 0.0_wp
    a_interface_yPos = 0.0_wp

  END SUBROUTINE initialize_hyperbolic

  SUBROUTINE finalize_hyperbolic

    DEALLOCATE( a_interface_xNeg )
    DEALLOCATE( a_interface_xPos )
    DEALLOCATE( a_interface_yNeg )
    DEALLOCATE( a_interface_yPos )
    DEALLOCATE( H_interface_x )
    DEALLOCATE( H_interface_y )

  END SUBROUTINE finalize_hyperbolic

  SUBROUTINE eval_hyperbolic_terms( q_expl, qp_expl, divFlux_iRK, t,         &
       solve_cells, j_cent, k_cent, solve_interfaces_x, j_stag_x, k_stag_x,   &
       solve_interfaces_y, j_stag_y, k_stag_y )

    ! External variables
    USE parameters_2d, ONLY : solver_scheme

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: q_expl(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: qp_expl(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(OUT) :: divFlux_iRK(n_eqns,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: t
    INTEGER, INTENT(IN) :: solve_cells
    INTEGER, INTENT(IN) :: j_cent(:), k_cent(:)
    INTEGER, INTENT(IN) :: solve_interfaces_x, solve_interfaces_y
    INTEGER, INTENT(IN) :: j_stag_x(:), k_stag_x(:)
    INTEGER, INTENT(IN) :: j_stag_y(:), k_stag_y(:)

    INTEGER :: l , i, j, k      !< loop counters

    !WRITE(*,*) 'SUBROUTINE eval_hyperbolic_terms'
    !WRITE(*,*) 'qp_expl(4,1,1)',qp_expl(4,1,1)
    !WRITE(*,*)
    
    ! Linear reconstruction of the physical variables at the interfaces
    CALL reconstruction( q_expl, qp_expl, t, solve_cells, j_cent, k_cent )

    ! Evaluation of the maximum local speeds at the interfaces
    CALL eval_speeds( solve_interfaces_x, j_stag_x, k_stag_x,               &
         solve_interfaces_y, j_stag_y, k_stag_y )

    ! Evaluation of the numerical fluxes
    SELECT CASE ( solver_scheme )

    CASE ("LxF")

       CALL eval_flux_LxF

    CASE ("GFORCE")

       CALL eval_flux_GFORCE

    CASE ("KT")

       CALL eval_flux_KT( solve_interfaces_x, j_stag_x, k_stag_x,           &
            solve_interfaces_y, j_stag_y, k_stag_y )

    CASE ("UP")

       CALL eval_flux_UP( solve_interfaces_x, j_stag_x, k_stag_x,           &
            solve_interfaces_y, j_stag_y, k_stag_y )

    END SELECT

    !$OMP PARALLEL DO private(l,j,k,i)

    cells_loop:DO l = 1,solve_cells

       j = j_cent(l)
       k = k_cent(l)

       DO i=1,n_eqns

          divFlux_iRK(i,j,k) = 0.0_wp

          IF ( comp_cells_x .GT. 1 ) THEN

             divFlux_iRK(i,j,k) = divFlux_iRK(i,j,k) +                          &
                  ( H_interface_x(i,j+1,k) - H_interface_x(i,j,k) ) * one_by_dx

          END IF

          IF ( comp_cells_y .GT. 1 ) THEN

             divFlux_iRK(i,j,k) = divFlux_iRK(i,j,k) +                          &
                  ( H_interface_y(i,j,k+1) - H_interface_y(i,j,k) ) * one_by_dy

          END IF

       END DO

    END DO cells_loop

    !$OMP END PARALLEL DO

    RETURN

  END SUBROUTINE eval_hyperbolic_terms

  !******************************************************************************
  !> \brief Upwind numerical fluxes
  !
  !> This subroutine evaluates the numerical fluxes H at the 
  !> cells interfaces with an upwind discretization.
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 2019/11/16
  !******************************************************************************
  
  SUBROUTINE eval_flux_UP( solve_interfaces_x, j_stag_x, k_stag_x,          &
       solve_interfaces_y, j_stag_y, k_stag_y )

    ! External procedures
    USE constitutive_2d, ONLY : eval_fluxes
    USE geometry_2d, ONLY : grav_coeff_stag_x , grav_coeff_stag_y

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: solve_interfaces_x, solve_interfaces_y
    INTEGER, INTENT(IN) :: j_stag_x(:), k_stag_x(:)
    INTEGER, INTENT(IN) :: j_stag_y(:), k_stag_y(:)

    REAL(wp) :: fluxL(n_eqns)           !< Numerical fluxes from the eqns 
    REAL(wp) :: fluxR(n_eqns)           !< Numerical fluxes from the eqns
    REAL(wp) :: fluxB(n_eqns)           !< Numerical fluxes from the eqns 
    REAL(wp) :: fluxT(n_eqns)           !< Numerical fluxes from the eqns

    INTEGER :: j,k,l                  !< Loop counters

    H_interface_x = 0.0_wp
    H_interface_y = 0.0_wp

    IF ( comp_cells_x .GT. 1 ) THEN

       !$OMP PARALLEL DO private(l,j,k,fluxL,fluxR)

       DO l = 1,solve_interfaces_x

          j = j_stag_x(l)
          k = k_stag_x(l)

          CALL eval_fluxes( q_interfaceL(1:n_vars,j,k) ,                        &
               qp_interfaceL(1:n_vars+2,j,k) ,                                  &
               B_prime_x_geom(MAX(1,j-1),MIN(k,comp_cells_y)) ,                 &
               B_prime_y_geom(MAX(1,j-1),MIN(k,comp_cells_y)) ,                 &
               grav_coeff_stag_x(j,k) , 1 , fluxL )

          CALL eval_fluxes( q_interfaceR(1:n_vars,j,k) ,                        &
               qp_interfaceR(1:n_vars+2,j,k) ,                                  &
               B_prime_x_geom(MIN(j,comp_cells_x),MIN(k,comp_cells_y)) ,        &
               B_prime_y_geom(MIN(j,comp_cells_x),MIN(k,comp_cells_y)) ,        &
               grav_coeff_stag_x(j,k) , 1 , fluxR )

          IF ( ( qp_interfaceL(n_vars+1,j,k) .GT. 0.0_wp ) .AND.                &
               ( qp_interfaceR(n_vars+1,j,k) .GE. 0.0_wp ) ) THEN

             H_interface_x(:,j,k) = fluxL

          ELSEIF ( ( qp_interfaceL(n_vars+1,j,k) .LE. 0.0_wp ) .AND.            &
               ( qp_interfaceR(n_vars+1,j,k) .LT. 0.0_wp ) ) THEN

             H_interface_x(:,j,k) = fluxR

          ELSE

             H_interface_x(:,j,k) = 0.5_wp * ( fluxL + fluxR )

          END IF

          IF ( (  qp_interfaceL(n_vars+1,j,k) .EQ. 0.0_wp ) .AND.               &
               (  qp_interfaceR(n_vars+1,j,k) .EQ. 0.0_wp ) ) THEN

             H_interface_x(1,j,k) = 0.0_wp
             H_interface_x(4:n_vars,j,k) = 0.0_wp

          END IF
               
       END DO

       !$OMP END PARALLEL DO

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN

       !$OMP PARALLEL DO private(l,j,k,fluxB,fluxT)
       
       DO l = 1,solve_interfaces_y

          j = j_stag_y(l)
          k = k_stag_y(l)

          CALL eval_fluxes( q_interfaceB(1:n_vars,j,k) ,                        &
               qp_interfaceB(1:n_vars+2,j,k) ,                                  &
               B_prime_x_geom(MIN(j,comp_cells_x),MAX(1,k-1)) ,                 &
               B_prime_y_geom(MIN(j,comp_cells_x),MAX(1,k-1)) ,                 &
               grav_coeff_stag_y(j,k) , 2 , fluxB )

          CALL eval_fluxes( q_interfaceT(1:n_vars,j,k) ,                        &
               qp_interfaceT(1:n_vars+2,j,k) ,                                  &
               B_prime_x_geom(MIN(j,comp_cells_x),MIN(k,comp_cells_y)) ,        &
               B_prime_y_geom(MIN(j,comp_cells_x),MIN(k,comp_cells_y)) ,        &
               grav_coeff_stag_y(j,k) , 2 , fluxT )

          IF ( ( q_interfaceB(3,j,k) .GT. 0.0_wp ) .AND.                        &
               ( q_interfaceT(3,j,k) .GE. 0.0_wp ) ) THEN

             H_interface_y(:,j,k) = fluxB

          ELSEIF ( ( q_interfaceB(3,j,k) .LE. 0.0_wp ) .AND.                    &
               ( q_interfaceT(3,j,k) .LT. 0.0_wp ) ) THEN

             H_interface_y(:,j,k) = fluxT

          ELSE

             H_interface_y(:,j,k) = 0.5_wp * ( fluxB + fluxT )

          END IF

          ! In the equation for mass and for trasnport (T,alphas) if the 
          ! velocities at the interfaces are null, then the flux is null
          IF ( (  qp_interfaceB(n_vars+2,j,k) .EQ. 0.0_wp ) .AND.               &
               (  qp_interfaceT(n_vars+2,j,k) .EQ. 0.0_wp ) ) THEN

             H_interface_y(1,j,k) = 0.0_wp
             H_interface_y(4:n_vars,j,k) = 0.0_wp

          END IF
          
       END DO

       !$OMP END PARALLEL DO
       
    END IF

    RETURN

  END SUBROUTINE eval_flux_UP


  !******************************************************************************
  !> \brief Semidiscrete numerical fluxes
  !
  !> This subroutine evaluates the numerical fluxes H at the 
  !> cells interfaces according to Kurganov et al. 2001. 
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 16/08/2011
  !******************************************************************************

  SUBROUTINE eval_flux_KT( solve_interfaces_x, j_stag_x, k_stag_x,          &
       solve_interfaces_y, j_stag_y, k_stag_y )

    ! External procedures
    USE constitutive_2d, ONLY : eval_fluxes
    USE geometry_2d, ONLY : grav_coeff_stag_x , grav_coeff_stag_y

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: solve_interfaces_x, solve_interfaces_y
    INTEGER, INTENT(IN) :: j_stag_x(:), k_stag_x(:)
    INTEGER, INTENT(IN) :: j_stag_y(:), k_stag_y(:)

    REAL(wp) :: fluxL(n_eqns)           !< Numerical fluxes from the eqns 
    REAL(wp) :: fluxR(n_eqns)           !< Numerical fluxes from the eqns
    REAL(wp) :: fluxB(n_eqns)           !< Numerical fluxes from the eqns 
    REAL(wp) :: fluxT(n_eqns)           !< Numerical fluxes from the eqns

    REAL(wp) :: flux_avg_x(n_eqns)   
    REAL(wp) :: flux_avg_y(n_eqns)   

    INTEGER :: i,j,k,l                  !< Loop counters

    ! WRITE(*,*) 'eval_flux_KT: qp_interfaceR(1,1,1)',qp_interfaceR(1,1,1)


    !H_interface_x = 0.0_wp
    !H_interface_y = 0.0_wp

    !$OMP PARALLEL

    IF ( comp_cells_x .GT. 1 ) THEN

       !$OMP DO private(j,k,i,fluxL,fluxR,flux_avg_x)

       x_interfaces_loop:DO l = 1,solve_interfaces_x

          j = j_stag_x(l)
          k = k_stag_x(l)

          CALL eval_fluxes( q_interfaceL(1:n_vars,j,k) ,                        &
               qp_interfaceL(1:n_vars+2,j,k) ,                                  &
               B_prime_x_geom(MAX(1,j-1),MIN(k,comp_cells_y)) ,                 &
               B_prime_y_geom(MAX(1,j-1),MIN(k,comp_cells_y)) ,                 &
               grav_coeff_stag_x(j,k) , 1 , fluxL )

          CALL eval_fluxes( q_interfaceR(1:n_vars,j,k) ,                        &
               qp_interfaceR(1:n_vars+2,j,k) ,                                  &
               B_prime_x_geom(MIN(j,comp_cells_x),MIN(k,comp_cells_y)) ,        &
               B_prime_y_geom(MIN(j,comp_cells_x),MIN(k,comp_cells_y)) ,        &
               grav_coeff_stag_x(j,k) , 1 , fluxR )

          ! First term in Eq. 25 GMD paper
          CALL average_KT( a_interface_xNeg(:,j,k), a_interface_xPos(:,j,k) ,   &
               fluxL , fluxR , flux_avg_x )

          eqns_loop:DO i=1,n_eqns

             IF ( a_interface_xNeg(i,j,k) .EQ. a_interface_xPos(i,j,k) ) THEN

                H_interface_x(i,j,k) = 0.0_wp

             ELSE

                ! Eq. 25 from GMD paper
                H_interface_x(i,j,k) = flux_avg_x(i)                            &
                     + ( a_interface_xPos(i,j,k) * a_interface_xNeg(i,j,k) )    &
                     / ( a_interface_xPos(i,j,k) - a_interface_xNeg(i,j,k) )    &
                     * ( q_interfaceR(i,j,k) - q_interfaceL(i,j,k) )             

             END IF

          ENDDO eqns_loop

          ! Fix to avoid sum of solid fluxes larger tham flux for mixture.
          ! Guarded: H_interface_x(1,j,k) is zero at dry or zero-flux
          ! interfaces and the test used to divide by it unconditionally.
          IF ( H_interface_x(1,j,k) .GT. 0.0_wp ) THEN

             IF ( SUM(H_interface_x(idx_solidEqn_first:idx_solidEqn_last,j,k))  &
                  .GE. H_interface_x(1,j,k) ) THEN

                H_interface_x(idx_solidEqn_first:idx_solidEqn_last,j,k) =       &
                     H_interface_x(idx_solidEqn_first:idx_solidEqn_last,j,k) /  &
                     ( SUM(H_interface_x(idx_solidEqn_first:idx_solidEqn_last,  &
                     j,k)) / H_interface_x(1,j,k) )

             END IF

          END IF
          
          ! In the equation for mass and for trasnport (T,alphas) if the 
          ! velocities at the interfaces are null, then the flux is null
          IF ( (  qp_interfaceL(2,j,k) .EQ. 0.0_wp ) .AND.                      &
               (  qp_interfaceR(2,j,k) .EQ. 0.0_wp ) ) THEN

             H_interface_x(1,j,k) = 0.0_wp
             H_interface_x(4:n_vars,j,k) = 0.0_wp

          END IF
          
       END DO x_interfaces_loop
       
       !$OMP END DO NOWAIT

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN

       !$OMP DO private(j,k,i,fluxB,fluxT,flux_avg_y)
       
       y_interfaces_loop:DO l = 1,solve_interfaces_y

          j = j_stag_y(l)
          k = k_stag_y(l)

          CALL eval_fluxes( q_interfaceB(1:n_vars,j,k) ,                        &
               qp_interfaceB(1:n_vars+2,j,k) ,                                  &
               B_prime_x_geom(MIN(j,comp_cells_x),MAX(1,k-1)) ,                 &
               B_prime_y_geom(MIN(j,comp_cells_x),MAX(1,k-1)) ,                 &
               grav_coeff_stag_y(j,k) , 2 , fluxB )

          CALL eval_fluxes( q_interfaceT(1:n_vars,j,k) ,                        &
               qp_interfaceT(1:n_vars+2,j,k) ,                                  &
               B_prime_x_geom(MIN(j,comp_cells_x),MIN(k,comp_cells_y)) ,        &
               B_prime_y_geom(MIN(j,comp_cells_x),MIN(k,comp_cells_y)) ,        &
               grav_coeff_stag_y(j,k) , 2 , fluxT )
          
          CALL average_KT( a_interface_yNeg(:,j,k) ,                            &
               a_interface_yPos(:,j,k) , fluxB , fluxT , flux_avg_y )

          DO i=1,n_eqns

             IF ( a_interface_yNeg(i,j,k) .EQ. a_interface_yPos(i,j,k) ) THEN

                H_interface_y(i,j,k) = 0.0_wp

             ELSE

                H_interface_y(i,j,k) = flux_avg_y(i)                            &
                     + ( a_interface_yPos(i,j,k) * a_interface_yNeg(i,j,k) )    &
                     / ( a_interface_yPos(i,j,k) - a_interface_yNeg(i,j,k) )    &
                     * ( q_interfaceT(i,j,k) - q_interfaceB(i,j,k) )             

             END IF

          END DO

          ! Fix to avoid sum of solid fluxes larger tham flux for mixture.
          ! Guarded: see the x-interface limiter above.
          IF ( H_interface_y(1,j,k) .GT. 0.0_wp ) THEN

             IF ( SUM(H_interface_y(idx_solidEqn_first:idx_solidEqn_last,j,k))  &
                  .GT. H_interface_y(1,j,k) ) THEN

                H_interface_y(idx_solidEqn_first:idx_solidEqn_last,j,k) =       &
                     H_interface_y(idx_solidEqn_first:idx_solidEqn_last,j,k) /  &
                     ( SUM(H_interface_y(idx_solidEqn_first:idx_solidEqn_last,  &
                     j,k)) / H_interface_y(1,j,k) )

             END IF

          END IF
          
          ! In the equation for mass and for trasnport (T,alphas) if the 
          ! velocities at the interfaces are null, then the flux is null
          IF ( (  q_interfaceB(3,j,k) .EQ. 0.0_wp ) .AND.                       &
               (  q_interfaceT(3,j,k) .EQ. 0.0_wp ) ) THEN

             H_interface_y(1,j,k) = 0.0_wp
             H_interface_y(4:n_vars,j,k) = 0.0_wp

          END IF

       END DO y_interfaces_loop
       
       !$OMP END DO

    END IF

    !$OMP END PARALLEL

    RETURN
    
  END SUBROUTINE eval_flux_KT

  !******************************************************************************
  !> \brief averaged KT flux
  !
  !> This subroutine compute n averaged flux from the fluxes at the two sides of
  !> a cell interface and the max an min speed at the two sides.
  !> \param[in]     a1            speed at one side of the interface
  !> \param[in]     a2            speed at the other side of the interface
  !> \param[in]     w1            fluxes at one side of the interface
  !> \param[in]     w2            fluxes at the other side of the interface
  !> \param[out]    w_avg         array of averaged fluxes
  !> \date 2019/12/13
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE average_KT( a1 , a2 , w1 , w2 , w_avg )

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: a1(:) , a2(:)
    REAL(wp), INTENT(IN) :: w1(:) , w2(:)
    REAL(wp), INTENT(OUT) :: w_avg(:)

    INTEGER :: n
    INTEGER :: i 

    n = SIZE( a1 )

    DO i=1,n

       IF ( a1(i) .EQ. a2(i) ) THEN

          w_avg(i) = 0.5_wp * ( w1(i) + w2(i) )
          w_avg(i) = 0.0_wp

       ELSE

          w_avg(i) = ( a2(i) * w1(i) - a1(i) * w2(i) ) / ( a2(i) - a1(i) )  

       END IF

    END DO

    RETURN
    
  END SUBROUTINE average_KT

  !******************************************************************************
  !> \brief Numerical fluxes GFORCE
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_flux_GFORCE

    ! to be implemented
    WRITE(*,*) 'method not yet implemented in 2-d case'

  END SUBROUTINE eval_flux_GFORCE

  !******************************************************************************
  !> \brief Numerical fluxes Lax-Friedrichs
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_flux_LxF

    ! to be implemented
    WRITE(*,*) 'method not yet implemented in 2-d case'

  END SUBROUTINE eval_flux_LxF


  !******************************************************************************
  !> \brief Characteristic speeds
  !
  !> This subroutine evaluates the largest characteristic speed at the
  !> cells interfaces from the reconstructed states.
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 2019/11/11
  !******************************************************************************

  SUBROUTINE eval_speeds( solve_interfaces_x, j_stag_x, k_stag_x,          &
       solve_interfaces_y, j_stag_y, k_stag_y )

    ! External procedures
    USE constitutive_2d, ONLY : eval_local_speeds_x, eval_local_speeds_y 

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: solve_interfaces_x, solve_interfaces_y
    INTEGER, INTENT(IN) :: j_stag_x(:), k_stag_x(:)
    INTEGER, INTENT(IN) :: j_stag_y(:), k_stag_y(:)

    REAL(wp) :: abslambdaL_min(n_vars) , abslambdaL_max(n_vars)
    REAL(wp) :: abslambdaR_min(n_vars) , abslambdaR_max(n_vars)
    REAL(wp) :: abslambdaB_min(n_vars) , abslambdaB_max(n_vars)
    REAL(wp) :: abslambdaT_min(n_vars) , abslambdaT_max(n_vars)
    REAL(wp) :: min_r(n_vars) , max_r(n_vars)

    INTEGER :: j,k,l

    !$OMP PARALLEL

    IF ( comp_cells_x .GT. 1 ) THEN

       !$OMP DO private(j , k , abslambdaL_min , abslambdaL_max ,               &
       !$OMP & abslambdaR_min , abslambdaR_max , min_r , max_r )

       x_interfaces_loop:DO l = 1,solve_interfaces_x

          j = j_stag_x(l)
          k = k_stag_x(l)

          CALL eval_local_speeds_x( qp_interfaceL(:,j,k) ,                      &
               grav_coeff_stag_x(j,k) , abslambdaL_min , abslambdaL_max )

          CALL eval_local_speeds_x( qp_interfaceR(:,j,k) ,                      &
               grav_coeff_stag_x(j,k) , abslambdaR_min , abslambdaR_max )

          min_r = MIN(abslambdaL_min , abslambdaR_min , 0.0_wp)
          max_r = MAX(abslambdaL_max , abslambdaR_max , 0.0_wp)

          a_interface_xNeg(:,j,k) = min_r
          a_interface_xPos(:,j,k) = max_r

       END DO x_interfaces_loop

       !$OMP END DO NOWAIT

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN

       !$OMP DO private(j , k , abslambdaB_min , abslambdaB_max ,               &
       !$OMP & abslambdaT_min , abslambdaT_max , min_r , max_r )

       y_interfaces_loop:DO l = 1,solve_interfaces_y

          j = j_stag_y(l)
          k = k_stag_y(l)

          CALL eval_local_speeds_y( qp_interfaceB(:,j,k) ,                      &
               grav_coeff_stag_y(j,k) , abslambdaB_min , abslambdaB_max )
          
          CALL eval_local_speeds_y( qp_interfaceT(:,j,k) ,                      &
               grav_coeff_stag_y(j,k) , abslambdaT_min , abslambdaT_max )

          min_r = MIN(abslambdaB_min , abslambdaT_min , 0.0_wp)
          max_r = MAX(abslambdaB_max , abslambdaT_max , 0.0_wp)

          a_interface_yNeg(:,j,k) = min_r
          a_interface_yPos(:,j,k) = max_r

       END DO y_interfaces_loop

       !$OMP END DO

    END IF

    !$OMP END PARALLEL

    RETURN
    
  END SUBROUTINE eval_speeds

END MODULE hyperbolic_2d

