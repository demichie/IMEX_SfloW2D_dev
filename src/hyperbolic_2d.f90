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

  USE reconstruction_2d, ONLY : recon => reconstruction_workspace

  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: hyperbolic_workspace_type
     REAL(wp), ALLOCATABLE :: a_interface_xNeg(:,:,:)
     REAL(wp), ALLOCATABLE :: a_interface_xPos(:,:,:)
     REAL(wp), ALLOCATABLE :: a_interface_yNeg(:,:,:)
     REAL(wp), ALLOCATABLE :: a_interface_yPos(:,:,:)
     REAL(wp), ALLOCATABLE :: H_interface_x(:,:,:)
     REAL(wp), ALLOCATABLE :: H_interface_y(:,:,:)
   CONTAINS
     PROCEDURE :: initialize => initialize_hyperbolic
     PROCEDURE :: finalize => finalize_hyperbolic
     PROCEDURE :: evaluate_terms => eval_hyperbolic_terms
     PROCEDURE :: evaluate_speeds => eval_speeds
  END TYPE hyperbolic_workspace_type

  TYPE(hyperbolic_workspace_type), PUBLIC, TARGET :: hyperbolic_workspace

CONTAINS

  SUBROUTINE initialize_hyperbolic( this )

    CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this

    ALLOCATE( this%a_interface_xNeg(n_eqns,comp_interfaces_x,comp_cells_y) )
    ALLOCATE( this%a_interface_xPos(n_eqns,comp_interfaces_x,comp_cells_y) )
    ALLOCATE( this%a_interface_yNeg(n_eqns,comp_cells_x,comp_interfaces_y) )
    ALLOCATE( this%a_interface_yPos(n_eqns,comp_cells_x,comp_interfaces_y) )

    ALLOCATE( this%H_interface_x(n_eqns,comp_interfaces_x,comp_cells_y) )
    ALLOCATE( this%H_interface_y(n_eqns,comp_cells_x,comp_interfaces_y) )

    this%a_interface_xNeg = 0.0_wp
    this%a_interface_xPos = 0.0_wp
    this%a_interface_yNeg = 0.0_wp
    this%a_interface_yPos = 0.0_wp

  END SUBROUTINE initialize_hyperbolic

  SUBROUTINE finalize_hyperbolic( this )

    CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this

    DEALLOCATE( this%a_interface_xNeg )
    DEALLOCATE( this%a_interface_xPos )
    DEALLOCATE( this%a_interface_yNeg )
    DEALLOCATE( this%a_interface_yPos )
    DEALLOCATE( this%H_interface_x )
    DEALLOCATE( this%H_interface_y )

  END SUBROUTINE finalize_hyperbolic

  SUBROUTINE eval_hyperbolic_terms( this, q_expl, qp_expl, divFlux_iRK, t,   &
       solve_cells, j_cent, k_cent, solve_interfaces_x, j_stag_x, k_stag_x,   &
       solve_interfaces_y, j_stag_y, k_stag_y )

    ! External variables
    USE parameters_2d, ONLY : solver_scheme

    IMPLICIT NONE

    CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this
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
    CALL recon%reconstruct( q_expl, qp_expl, t, solve_cells, &
         j_cent, k_cent )

    ! Evaluation of the maximum local speeds at the interfaces
    CALL this%evaluate_speeds( solve_interfaces_x, j_stag_x, k_stag_x,      &
         solve_interfaces_y, j_stag_y, k_stag_y )

    ! Evaluation of the numerical fluxes
    SELECT CASE ( solver_scheme )

    CASE ("LxF")

       CALL eval_flux_LxF

    CASE ("GFORCE")

       CALL eval_flux_GFORCE

    CASE ("KT")

       CALL eval_flux_KT( this, solve_interfaces_x, j_stag_x, k_stag_x,     &
            solve_interfaces_y, j_stag_y, k_stag_y )

    CASE ("UP")

       CALL eval_flux_UP( this, solve_interfaces_x, j_stag_x, k_stag_x,     &
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
                  ( this%H_interface_x(i,j+1,k) - this%H_interface_x(i,j,k) ) * one_by_dx

          END IF

          IF ( comp_cells_y .GT. 1 ) THEN

             divFlux_iRK(i,j,k) = divFlux_iRK(i,j,k) +                          &
                  ( this%H_interface_y(i,j,k+1) - this%H_interface_y(i,j,k) ) * one_by_dy

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
  
  SUBROUTINE eval_flux_UP( this, solve_interfaces_x, j_stag_x, k_stag_x,    &
       solve_interfaces_y, j_stag_y, k_stag_y )

    ! External procedures
    USE constitutive_2d, ONLY : eval_fluxes
    USE geometry_2d, ONLY : grav_coeff_stag_x , grav_coeff_stag_y

    IMPLICIT NONE

    CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: solve_interfaces_x, solve_interfaces_y
    INTEGER, INTENT(IN) :: j_stag_x(:), k_stag_x(:)
    INTEGER, INTENT(IN) :: j_stag_y(:), k_stag_y(:)

    REAL(wp) :: fluxL(n_eqns)           !< Numerical fluxes from the eqns 
    REAL(wp) :: fluxR(n_eqns)           !< Numerical fluxes from the eqns
    REAL(wp) :: fluxB(n_eqns)           !< Numerical fluxes from the eqns 
    REAL(wp) :: fluxT(n_eqns)           !< Numerical fluxes from the eqns

    INTEGER :: j,k,l                  !< Loop counters

    this%H_interface_x = 0.0_wp
    this%H_interface_y = 0.0_wp

    IF ( comp_cells_x .GT. 1 ) THEN

       !$OMP PARALLEL DO private(l,j,k,fluxL,fluxR)

       DO l = 1,solve_interfaces_x

          j = j_stag_x(l)
          k = k_stag_x(l)

          CALL eval_fluxes( recon%q_interfaceL(1:n_vars,j,k) ,                        &
               recon%qp_interfaceL(1:n_vars+2,j,k) ,                                  &
               B_prime_x_geom(MAX(1,j-1),MIN(k,comp_cells_y)) ,                 &
               B_prime_y_geom(MAX(1,j-1),MIN(k,comp_cells_y)) ,                 &
               grav_coeff_stag_x(j,k) , 1 , fluxL )

          CALL eval_fluxes( recon%q_interfaceR(1:n_vars,j,k) ,                        &
               recon%qp_interfaceR(1:n_vars+2,j,k) ,                                  &
               B_prime_x_geom(MIN(j,comp_cells_x),MIN(k,comp_cells_y)) ,        &
               B_prime_y_geom(MIN(j,comp_cells_x),MIN(k,comp_cells_y)) ,        &
               grav_coeff_stag_x(j,k) , 1 , fluxR )

          IF ( ( recon%qp_interfaceL(n_vars+1,j,k) .GT. 0.0_wp ) .AND.                &
               ( recon%qp_interfaceR(n_vars+1,j,k) .GE. 0.0_wp ) ) THEN

             this%H_interface_x(:,j,k) = fluxL

          ELSEIF ( ( recon%qp_interfaceL(n_vars+1,j,k) .LE. 0.0_wp ) .AND.            &
               ( recon%qp_interfaceR(n_vars+1,j,k) .LT. 0.0_wp ) ) THEN

             this%H_interface_x(:,j,k) = fluxR

          ELSE

             this%H_interface_x(:,j,k) = 0.5_wp * ( fluxL + fluxR )

          END IF

          IF ( (  recon%qp_interfaceL(n_vars+1,j,k) .EQ. 0.0_wp ) .AND.               &
               (  recon%qp_interfaceR(n_vars+1,j,k) .EQ. 0.0_wp ) ) THEN

             this%H_interface_x(1,j,k) = 0.0_wp
             this%H_interface_x(4:n_vars,j,k) = 0.0_wp

          END IF
               
       END DO

       !$OMP END PARALLEL DO

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN

       !$OMP PARALLEL DO private(l,j,k,fluxB,fluxT)
       
       DO l = 1,solve_interfaces_y

          j = j_stag_y(l)
          k = k_stag_y(l)

          CALL eval_fluxes( recon%q_interfaceB(1:n_vars,j,k) ,                        &
               recon%qp_interfaceB(1:n_vars+2,j,k) ,                                  &
               B_prime_x_geom(MIN(j,comp_cells_x),MAX(1,k-1)) ,                 &
               B_prime_y_geom(MIN(j,comp_cells_x),MAX(1,k-1)) ,                 &
               grav_coeff_stag_y(j,k) , 2 , fluxB )

          CALL eval_fluxes( recon%q_interfaceT(1:n_vars,j,k) ,                        &
               recon%qp_interfaceT(1:n_vars+2,j,k) ,                                  &
               B_prime_x_geom(MIN(j,comp_cells_x),MIN(k,comp_cells_y)) ,        &
               B_prime_y_geom(MIN(j,comp_cells_x),MIN(k,comp_cells_y)) ,        &
               grav_coeff_stag_y(j,k) , 2 , fluxT )

          IF ( ( recon%q_interfaceB(3,j,k) .GT. 0.0_wp ) .AND.                        &
               ( recon%q_interfaceT(3,j,k) .GE. 0.0_wp ) ) THEN

             this%H_interface_y(:,j,k) = fluxB

          ELSEIF ( ( recon%q_interfaceB(3,j,k) .LE. 0.0_wp ) .AND.                    &
               ( recon%q_interfaceT(3,j,k) .LT. 0.0_wp ) ) THEN

             this%H_interface_y(:,j,k) = fluxT

          ELSE

             this%H_interface_y(:,j,k) = 0.5_wp * ( fluxB + fluxT )

          END IF

          ! In the equation for mass and for trasnport (T,alphas) if the 
          ! velocities at the interfaces are null, then the flux is null
          IF ( (  recon%qp_interfaceB(n_vars+2,j,k) .EQ. 0.0_wp ) .AND.               &
               (  recon%qp_interfaceT(n_vars+2,j,k) .EQ. 0.0_wp ) ) THEN

             this%H_interface_y(1,j,k) = 0.0_wp
             this%H_interface_y(4:n_vars,j,k) = 0.0_wp

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

  SUBROUTINE eval_flux_KT( this, solve_interfaces_x, j_stag_x, k_stag_x,    &
       solve_interfaces_y, j_stag_y, k_stag_y )

    ! External procedures
    USE constitutive_2d, ONLY : eval_fluxes
    USE geometry_2d, ONLY : grav_coeff_stag_x , grav_coeff_stag_y

    IMPLICIT NONE

    CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this
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

          CALL eval_fluxes( recon%q_interfaceL(1:n_vars,j,k) ,                        &
               recon%qp_interfaceL(1:n_vars+2,j,k) ,                                  &
               B_prime_x_geom(MAX(1,j-1),MIN(k,comp_cells_y)) ,                 &
               B_prime_y_geom(MAX(1,j-1),MIN(k,comp_cells_y)) ,                 &
               grav_coeff_stag_x(j,k) , 1 , fluxL )

          CALL eval_fluxes( recon%q_interfaceR(1:n_vars,j,k) ,                        &
               recon%qp_interfaceR(1:n_vars+2,j,k) ,                                  &
               B_prime_x_geom(MIN(j,comp_cells_x),MIN(k,comp_cells_y)) ,        &
               B_prime_y_geom(MIN(j,comp_cells_x),MIN(k,comp_cells_y)) ,        &
               grav_coeff_stag_x(j,k) , 1 , fluxR )

          ! First term in Eq. 25 GMD paper
          CALL average_KT( this%a_interface_xNeg(:,j,k), this%a_interface_xPos(:,j,k) ,   &
               fluxL , fluxR , flux_avg_x )

          eqns_loop:DO i=1,n_eqns

             IF ( this%a_interface_xNeg(i,j,k) .EQ. this%a_interface_xPos(i,j,k) ) THEN

                this%H_interface_x(i,j,k) = 0.0_wp

             ELSE

                ! Eq. 25 from GMD paper
                this%H_interface_x(i,j,k) = flux_avg_x(i)                            &
                     + ( this%a_interface_xPos(i,j,k) * this%a_interface_xNeg(i,j,k) )    &
                     / ( this%a_interface_xPos(i,j,k) - this%a_interface_xNeg(i,j,k) )    &
                     * ( recon%q_interfaceR(i,j,k) - recon%q_interfaceL(i,j,k) )

             END IF

          ENDDO eqns_loop

          ! Fix to avoid sum of solid fluxes larger tham flux for mixture.
          ! Guarded: H_interface_x(1,j,k) is zero at dry or zero-flux
          ! interfaces and the test used to divide by it unconditionally.
          IF ( this%H_interface_x(1,j,k) .GT. 0.0_wp ) THEN

             IF ( SUM(this%H_interface_x(idx_solidEqn_first:idx_solidEqn_last,j,k))  &
                  .GE. this%H_interface_x(1,j,k) ) THEN

                this%H_interface_x(idx_solidEqn_first:idx_solidEqn_last,j,k) =       &
                     this%H_interface_x(idx_solidEqn_first:idx_solidEqn_last,j,k) /  &
                     ( SUM(this%H_interface_x(idx_solidEqn_first:idx_solidEqn_last,  &
                     j,k)) / this%H_interface_x(1,j,k) )

             END IF

          END IF
          
          ! In the equation for mass and for trasnport (T,alphas) if the 
          ! velocities at the interfaces are null, then the flux is null
          IF ( (  recon%qp_interfaceL(2,j,k) .EQ. 0.0_wp ) .AND.                      &
               (  recon%qp_interfaceR(2,j,k) .EQ. 0.0_wp ) ) THEN

             this%H_interface_x(1,j,k) = 0.0_wp
             this%H_interface_x(4:n_vars,j,k) = 0.0_wp

          END IF
          
       END DO x_interfaces_loop
       
       !$OMP END DO NOWAIT

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN

       !$OMP DO private(j,k,i,fluxB,fluxT,flux_avg_y)
       
       y_interfaces_loop:DO l = 1,solve_interfaces_y

          j = j_stag_y(l)
          k = k_stag_y(l)

          CALL eval_fluxes( recon%q_interfaceB(1:n_vars,j,k) ,                        &
               recon%qp_interfaceB(1:n_vars+2,j,k) ,                                  &
               B_prime_x_geom(MIN(j,comp_cells_x),MAX(1,k-1)) ,                 &
               B_prime_y_geom(MIN(j,comp_cells_x),MAX(1,k-1)) ,                 &
               grav_coeff_stag_y(j,k) , 2 , fluxB )

          CALL eval_fluxes( recon%q_interfaceT(1:n_vars,j,k) ,                        &
               recon%qp_interfaceT(1:n_vars+2,j,k) ,                                  &
               B_prime_x_geom(MIN(j,comp_cells_x),MIN(k,comp_cells_y)) ,        &
               B_prime_y_geom(MIN(j,comp_cells_x),MIN(k,comp_cells_y)) ,        &
               grav_coeff_stag_y(j,k) , 2 , fluxT )
          
          CALL average_KT( this%a_interface_yNeg(:,j,k) ,                            &
               this%a_interface_yPos(:,j,k) , fluxB , fluxT , flux_avg_y )

          DO i=1,n_eqns

             IF ( this%a_interface_yNeg(i,j,k) .EQ. this%a_interface_yPos(i,j,k) ) THEN

                this%H_interface_y(i,j,k) = 0.0_wp

             ELSE

                this%H_interface_y(i,j,k) = flux_avg_y(i)                            &
                     + ( this%a_interface_yPos(i,j,k) * this%a_interface_yNeg(i,j,k) )    &
                     / ( this%a_interface_yPos(i,j,k) - this%a_interface_yNeg(i,j,k) )    &
                     * ( recon%q_interfaceT(i,j,k) - recon%q_interfaceB(i,j,k) )

             END IF

          END DO

          ! Fix to avoid sum of solid fluxes larger tham flux for mixture.
          ! Guarded: see the x-interface limiter above.
          IF ( this%H_interface_y(1,j,k) .GT. 0.0_wp ) THEN

             IF ( SUM(this%H_interface_y(idx_solidEqn_first:idx_solidEqn_last,j,k))  &
                  .GT. this%H_interface_y(1,j,k) ) THEN

                this%H_interface_y(idx_solidEqn_first:idx_solidEqn_last,j,k) =       &
                     this%H_interface_y(idx_solidEqn_first:idx_solidEqn_last,j,k) /  &
                     ( SUM(this%H_interface_y(idx_solidEqn_first:idx_solidEqn_last,  &
                     j,k)) / this%H_interface_y(1,j,k) )

             END IF

          END IF
          
          ! In the equation for mass and for trasnport (T,alphas) if the 
          ! velocities at the interfaces are null, then the flux is null
          IF ( (  recon%q_interfaceB(3,j,k) .EQ. 0.0_wp ) .AND.                       &
               (  recon%q_interfaceT(3,j,k) .EQ. 0.0_wp ) ) THEN

             this%H_interface_y(1,j,k) = 0.0_wp
             this%H_interface_y(4:n_vars,j,k) = 0.0_wp

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

  SUBROUTINE eval_speeds( this, solve_interfaces_x, j_stag_x, k_stag_x,    &
       solve_interfaces_y, j_stag_y, k_stag_y )

    ! External procedures
    USE constitutive_2d, ONLY : eval_local_speeds_x, eval_local_speeds_y 

    IMPLICIT NONE

    CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this
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

          CALL eval_local_speeds_x( recon%qp_interfaceL(:,j,k) ,                      &
               grav_coeff_stag_x(j,k) , abslambdaL_min , abslambdaL_max )

          CALL eval_local_speeds_x( recon%qp_interfaceR(:,j,k) ,                      &
               grav_coeff_stag_x(j,k) , abslambdaR_min , abslambdaR_max )

          min_r = MIN(abslambdaL_min , abslambdaR_min , 0.0_wp)
          max_r = MAX(abslambdaL_max , abslambdaR_max , 0.0_wp)

          this%a_interface_xNeg(:,j,k) = min_r
          this%a_interface_xPos(:,j,k) = max_r

       END DO x_interfaces_loop

       !$OMP END DO NOWAIT

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN

       !$OMP DO private(j , k , abslambdaB_min , abslambdaB_max ,               &
       !$OMP & abslambdaT_min , abslambdaT_max , min_r , max_r )

       y_interfaces_loop:DO l = 1,solve_interfaces_y

          j = j_stag_y(l)
          k = k_stag_y(l)

          CALL eval_local_speeds_y( recon%qp_interfaceB(:,j,k) ,                      &
               grav_coeff_stag_y(j,k) , abslambdaB_min , abslambdaB_max )
          
          CALL eval_local_speeds_y( recon%qp_interfaceT(:,j,k) ,                      &
               grav_coeff_stag_y(j,k) , abslambdaT_min , abslambdaT_max )

          min_r = MIN(abslambdaB_min , abslambdaT_min , 0.0_wp)
          max_r = MAX(abslambdaB_max , abslambdaT_max , 0.0_wp)

          this%a_interface_yNeg(:,j,k) = min_r
          this%a_interface_yPos(:,j,k) = max_r

       END DO y_interfaces_loop

       !$OMP END DO

    END IF

    !$OMP END PARALLEL

    RETURN
    
  END SUBROUTINE eval_speeds

END MODULE hyperbolic_2d
