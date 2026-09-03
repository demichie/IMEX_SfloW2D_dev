!********************************************************************************
!> \brief IMEX Runge-Kutta time integration
!>
!> This module owns the IMEX tableau and Runge-Kutta stage workspace and
!> orchestrates timestep selection and explicit/implicit stage advancement.
!********************************************************************************
MODULE time_integration_2d

  USE parameters_2d, ONLY : wp
  USE parameters_2d, ONLY : n_eqns, n_vars, n_RK, n_solid
  USE parameters_2d, ONLY : verbose_level

  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE geometry_2d, ONLY : dx, dy
  USE geometry_2d, ONLY : B_cent
  USE geometry_2d, ONLY : B_prime_x_geom, B_prime_y_geom
  USE geometry_2d, ONLY : B_second_xx_geom, B_second_xy_geom, B_second_yy_geom
  USE geometry_2d, ONLY : grav_coeff, d_grav_coeff_dx, d_grav_coeff_dy
  USE geometry_2d, ONLY : cell_source_fractions
  USE geometry_2d, ONLY : cell_arc_perim, cell_arc_n_x, cell_arc_n_y

  USE constitutive_2d, ONLY : rheology_model, T_ambient

  USE nonlinear_solver_2d, ONLY : solve_rk_step

  USE reconstruction_2d, ONLY : recon => reconstruction_workspace

  USE hyperbolic_2d, ONLY : hyper => hyperbolic_workspace

  USE domain_2d, ONLY : domain

  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: time_integration_workspace_type
     PRIVATE

     REAL(wp), ALLOCATABLE :: a_tilde_ij(:,:)
     REAL(wp), ALLOCATABLE :: a_dirk_ij(:,:)
     REAL(wp), ALLOCATABLE :: omega_tilde(:)
     REAL(wp), ALLOCATABLE :: omega(:)
     REAL(wp), ALLOCATABLE :: a_tilde(:)
     REAL(wp), ALLOCATABLE :: a_dirk(:)
     REAL(wp), ALLOCATABLE :: q_rk(:,:,:)
     REAL(wp), ALLOCATABLE :: qp_rk(:,:,:)
     REAL(wp), ALLOCATABLE :: divFlux(:,:,:,:)
     REAL(wp), ALLOCATABLE :: NH(:,:,:,:)
     REAL(wp), ALLOCATABLE :: SI_NH(:,:,:,:)
     REAL(wp), ALLOCATABLE :: expl_terms(:,:,:,:)
   CONTAINS
     PROCEDURE, PUBLIC :: initialize => initialize_time_integration
     PROCEDURE, PUBLIC :: finalize => finalize_time_integration
     PROCEDURE, NOPASS, PUBLIC :: compute_timestep => timestep
     PROCEDURE, PUBLIC :: advance => imex_RK_solver
  END TYPE time_integration_workspace_type

  TYPE(time_integration_workspace_type), PUBLIC, TARGET :: time_integration_workspace

CONTAINS

  SUBROUTINE initialize_time_integration( this )

    CLASS(time_integration_workspace_type), INTENT(INOUT) :: this
    REAL(wp) :: gamma, delta

    ALLOCATE( this%a_tilde_ij(n_RK,n_RK) )
    ALLOCATE( this%a_dirk_ij(n_RK,n_RK) )
    ALLOCATE( this%omega_tilde(n_RK) )
    ALLOCATE( this%omega(n_RK) )

    ! Initialize the coefficients for the IMEX Runge-Kutta scheme
    ! Please note that with respect to the schemes described in Pareschi & Russo
    ! (2000) we do not have the coefficient vectors c_tilde and c, because the
    ! explicit and implicit terms do not depend explicitly on time.

    this%a_tilde_ij = 0.0_wp
    this%omega_tilde = 0.0_wp
    this%a_dirk_ij = 0.0_wp
    this%omega = 0.0_wp

    gamma = 1.0_wp - 1.0_wp / SQRT(2.0_wp)
    delta = 1.0_wp - 1.0_wp / ( 2.0_wp * gamma )

    IF ( n_RK .EQ. 1 ) THEN

       this%a_tilde_ij(1,1) = 1.0_wp
       this%omega_tilde(1) = 1.0_wp
       this%a_dirk_ij(1,1) = 0.0_wp
       this%omega(1) = 0.0_wp

    ELSEIF ( n_RK .EQ. 2 ) THEN

       this%a_tilde_ij(2,1) = 1.0_wp
       this%omega_tilde(1) = 1.0_wp
       this%omega_tilde(2) = 0.0_wp
       this%a_dirk_ij(2,2) = 1.0_wp
       this%omega(1) = 0.0_wp
       this%omega(2) = 1.0_wp

    ELSEIF ( n_RK .EQ. 3 ) THEN

       ! Tableau for the IMEX-SSP(3,3,2) Stiffly Accurate Scheme
       ! from Pareschi & Russo (2005), Table IV
       this%a_tilde_ij(2,1) = 0.5_wp
       this%a_tilde_ij(3,1) = 0.5_wp
       this%a_tilde_ij(3,2) = 0.5_wp
       this%omega_tilde(1) = 1.0_wp / 3.0_wp
       this%omega_tilde(2) = 1.0_wp / 3.0_wp
       this%omega_tilde(3) = 1.0_wp / 3.0_wp
       this%a_dirk_ij(1,1) = 0.25_wp
       this%a_dirk_ij(2,2) = 0.25_wp
       this%a_dirk_ij(3,1) = 1.0_wp / 3.0_wp
       this%a_dirk_ij(3,2) = 1.0_wp / 3.0_wp
       this%a_dirk_ij(3,3) = 1.0_wp / 3.0_wp
       this%omega(1) = 1.0_wp / 3.0_wp
       this%omega(2) = 1.0_wp / 3.0_wp
       this%omega(3) = 1.0_wp / 3.0_wp

    ELSEIF ( n_RK .EQ. 4 ) THEN

       ! LRR(3,2,2) from Table 3 in Pareschi & Russo (2000)
       this%a_tilde_ij(2,1) = 0.5_wp
       this%a_tilde_ij(3,1) = 1.0_wp / 3.0_wp
       this%a_tilde_ij(4,2) = 1.0_wp
       this%omega_tilde(1) = 0.0_wp
       this%omega_tilde(2) = 1.0_wp
       this%omega_tilde(3) = 0.0_wp
       this%omega_tilde(4) = 0.0_wp
       this%a_dirk_ij(2,2) = 0.5_wp
       this%a_dirk_ij(3,3) = 1.0_wp / 3.0_wp
       this%a_dirk_ij(4,3) = 0.75_wp
       this%a_dirk_ij(4,4) = 0.25_wp
       this%omega(1) = 0.0_wp
       this%omega(2) = 0.0_wp
       this%omega(3) = 0.75_wp
       this%omega(4) = 0.25_wp

    END IF

    ALLOCATE( this%a_tilde(n_RK) )
    ALLOCATE( this%a_dirk(n_RK) )
    ALLOCATE( this%q_rk(n_vars,comp_cells_x,comp_cells_y) )
    ALLOCATE( this%qp_rk(n_vars+2,comp_cells_x,comp_cells_y) )
    ALLOCATE( this%divFlux(n_eqns,comp_cells_x,comp_cells_y,n_RK) )
    ALLOCATE( this%NH(n_eqns,comp_cells_x,comp_cells_y,n_RK) )
    ALLOCATE( this%SI_NH(n_eqns,comp_cells_x,comp_cells_y,n_RK) )
    ALLOCATE( this%expl_terms(n_eqns,comp_cells_x,comp_cells_y,n_RK) )

  END SUBROUTINE initialize_time_integration

  SUBROUTINE finalize_time_integration( this )

    CLASS(time_integration_workspace_type), INTENT(INOUT) :: this

    DEALLOCATE( this%a_tilde_ij )
    DEALLOCATE( this%a_dirk_ij )
    DEALLOCATE( this%omega_tilde )
    DEALLOCATE( this%omega )
    DEALLOCATE( this%a_tilde )
    DEALLOCATE( this%a_dirk )
    DEALLOCATE( this%q_rk )
    DEALLOCATE( this%qp_rk )
    DEALLOCATE( this%divFlux )
    DEALLOCATE( this%NH )
    DEALLOCATE( this%SI_NH )
    DEALLOCATE( this%expl_terms )

  END SUBROUTINE finalize_time_integration

  SUBROUTINE timestep(q, qp, t, dt)

    ! External variables
    USE geometry_2d, ONLY : dx,dy
    USE parameters_2d, ONLY : max_dt , cfl

    USE constitutive_2d, ONLY : qc_to_qp

    IMPLICIT none

    REAL(wp), INTENT(IN) :: q(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(INOUT) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: t
    REAL(wp), INTENT(OUT) :: dt

    INTEGER :: j,k,l          !< loop counter

    REAL(wp) :: max_a_x
    REAL(wp) :: max_a_y
    REAL(wp) p_dyn

    dt = max_dt

    IF ( cfl .NE. -1.0_wp ) THEN

       !$OMP PARALLEL DO private(j,k,p_dyn)

       DO l = 1,domain%solve_cells

          j = domain%j_cent(l)
          k = domain%k_cent(l)

          IF ( q(1,j,k) .GT. 0.0_wp ) THEN

             CALL qc_to_qp( q(1:n_vars,j,k) , qp(1:n_vars+2,j,k) , p_dyn )

          ELSE

             qp(1:n_vars+2,j,k) = 0.0_wp
             qp(4,j,k) = T_ambient

          END IF

       END DO

       !$OMP END PARALLEL DO

       !WRITE(*,*) 'qp(1:n_vars+2,1,1)',qp(1:n_vars+2,1,1)
       !READ(*,*)

       ! Compute the physical and conservative variables at the interfaces
       CALL recon%reconstruct( q, qp, t, domain%solve_cells, &
            domain%j_cent, domain%k_cent )

       ! Compute the max/min eigenvalues at the interfaces
       CALL hyper%evaluate_speeds( domain%solve_interfaces_x, domain%j_stag_x, domain%k_stag_x,    &
            domain%solve_interfaces_y, domain%j_stag_y, domain%k_stag_y )

       max_a_x = 0.0_wp
       max_a_y = 0.0_wp

       ! The minimum CFL step over all active cells is determined by the
       ! maximum characteristic speed on their adjacent interfaces.  Compute
       ! those two maxima directly, avoiding full-domain scratch arrays and
       ! an atomic update of dt for every cell.
       !$OMP PARALLEL DO private(j,k) reduction(max:max_a_x,max_a_y)
       DO l = 1,domain%solve_cells

          j = domain%j_cent(l)
          k = domain%k_cent(l)

          max_a_x = MAX( max_a_x,                                               &
               MAXVAL(hyper%a_interface_xPos(1:n_vars,j,k)),                         &
               MAXVAL(-hyper%a_interface_xNeg(1:n_vars,j,k)),                        &
               MAXVAL(hyper%a_interface_xPos(1:n_vars,j+1,k)),                       &
               MAXVAL(-hyper%a_interface_xNeg(1:n_vars,j+1,k)) )

          max_a_y = MAX( max_a_y,                                               &
               MAXVAL(hyper%a_interface_yPos(1:n_vars,j,k)),                         &
               MAXVAL(-hyper%a_interface_yNeg(1:n_vars,j,k)),                        &
               MAXVAL(hyper%a_interface_yPos(1:n_vars,j,k+1)),                       &
               MAXVAL(-hyper%a_interface_yNeg(1:n_vars,j,k+1)) )

       END DO
       !$OMP END PARALLEL DO

       IF ( max_a_x .GT. 0.0_wp ) dt = MIN(dt,cfl*dx/max_a_x)
       IF ( max_a_y .GT. 0.0_wp ) dt = MIN(dt,cfl*dy/max_a_y)

    END IF

    RETURN

  END SUBROUTINE timestep

  !******************************************************************************
  !> \brief Runge-Kutta integration
  !
  !> This subroutine integrate the hyperbolic conservation law with
  !> non-hyperbolic terms using an implicit-explicit runge-kutta scheme.
  !> The fluxes are integrated explicitely while the non-hyperbolic terms
  !> are integrated implicitely.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE imex_RK_solver(this, q, qp, t, dt, Z)

    USE constitutive_2d, ONLY : maximum_solid_packing
    
    USE constitutive_2d, ONLY : eval_implicit_terms

    USE constitutive_2d, ONLY : eval_nh_semi_impl_terms

    USE constitutive_2d, ONLY : qc_to_qp

    USE constitutive_2d, ONLY : eval_expl_terms

    USE constitutive_2d, ONLY : T_ambient

    USE geometry_2d, ONLY : B_nodata

    USE parameters_2d, ONLY : alpha_flag
    
!!$    USE parameters_2d, ONLY : time_param , bottom_radial_source_flag
    
    IMPLICIT NONE

    CLASS(time_integration_workspace_type), INTENT(INOUT) :: this
    REAL(wp), INTENT(INOUT) :: q(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(INOUT) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: t, dt
    REAL(wp), INTENT(IN) :: Z(comp_cells_x,comp_cells_y)

    REAL(wp) :: q_si(n_vars) !< solution after the semi-implicit step
    REAL(wp) :: q_guess(n_vars) !< initial guess for the solution of the RK step
    REAL(wp) :: q_fv_cell(n_vars) !< finite-volume state for the current cell
    REAL(wp) :: residual_cell(n_vars) !< final RK residual for the current cell
    REAL(wp) :: q_old_cell(n_vars) !< state before the final RK assembly
    INTEGER :: i_RK
    INTEGER :: j,k,l            !< loop counter over the grid volumes
    REAL(wp) :: Rj_not_impl(n_eqns)

    REAL(wp) :: a_diag
    REAL(wp) :: p_dyn

    REAL(wp) :: alpha_s
    LOGICAL :: solid_excess_roundoff
    LOGICAL :: need_explicit_stage

    INTEGER :: newton_iterations
    INTEGER :: newton_linear_info
    LOGICAL :: newton_converged
    LOGICAL :: newton_line_search_failed
    INTEGER :: newton_calls_step
    INTEGER :: newton_iterations_step
    INTEGER :: newton_iterations_max_step
    INTEGER :: newton_failures_step
    INTEGER :: newton_linear_failures_step
    INTEGER :: newton_line_search_failures_step

    newton_calls_step = 0
    newton_iterations_step = 0
    newton_iterations_max_step = 0
    newton_failures_step = 0
    newton_linear_failures_step = 0
    newton_line_search_failures_step = 0
    
    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'solver, imex_RK_solver: beginning'

    !$OMP PARALLEL
 
    !$OMP DO private(j,k)
    DO l = 1,domain%solve_cells

       j = domain%j_cent(l)
       k = domain%k_cent(l)

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*) 'solver, imex_RK_solver: j,k',j,k
          !READ(*,*)
          
       END IF

       ! Initialization of the variables for the Runge-Kutta scheme
       this%q_rk( 1:n_vars , j , k ) = 0.0_wp
       this%qp_rk( 1:n_vars+2 , j , k ) = 0.0_wp
       this%qp_rk( 4 , j , k ) = T_ambient
       

       this%divFlux(1:n_eqns , j , k , 1:n_RK ) = 0.0_wp
       this%NH( 1:n_eqns, j , k , 1:n_RK ) = 0.0_wp
       this%SI_NH( 1:n_eqns , j , k , 1:n_RK ) = 0.0_wp
       this%expl_terms(1:n_eqns , j , k , 1:n_RK) = 0.0_wp
       
    END DO
    !$OMP END DO

    !$OMP END PARALLEL

    runge_kutta:DO i_RK = 1,n_RK

       IF ( verbose_level .GE. 1 ) WRITE(*,*) 'solver, imex_RK_solver: i_RK',i_RK

       ! An explicit stage is required not only when it contributes to the
       ! final RK assembly, but also when a later stage depends on it.
       need_explicit_stage = ( this%omega_tilde(i_RK) .NE. 0.0_wp )

       IF ( i_RK .LT. n_RK ) THEN
          need_explicit_stage = need_explicit_stage .OR.                       &
               ANY( this%a_tilde_ij(i_RK+1:n_RK,i_RK) .NE. 0.0_wp )
       END IF

       ! define the explicits coefficients for the i-th step of the Runge-Kutta
       this%a_tilde = 0.0_wp
       this%a_dirk = 0.0_wp

       ! in the first step of the RK scheme all the coefficients remain to 0
       this%a_tilde(1:i_RK-1) = this%a_tilde_ij(i_RK,1:i_RK-1)
       this%a_dirk(1:i_RK-1) = this%a_dirk_ij(i_RK,1:i_RK-1)

       ! define the implicit coefficient for the i-th step of the Runge-Kutta
       a_diag = this%a_dirk_ij(i_RK,i_RK)

       !$OMP PARALLEL 
       !$OMP DO schedule(guided)                                                &
       !$OMP & private(j,k,q_guess,q_si,q_fv_cell,Rj_not_impl,p_dyn,           &
       !$OMP & newton_iterations,newton_linear_info,newton_converged,          &
       !$OMP & newton_line_search_failed)

       solve_cells_loop:DO l = 1,domain%solve_cells

          j = domain%j_cent(l)
          k = domain%k_cent(l)

          IF ( verbose_level .GE. 2 ) THEN

             WRITE(*,*) 'solver, imex_RK_solver: j,k',j,k
             ! READ(*,*)

          END IF

          ! initialize the RK step
          IF ( i_RK .EQ. 1 ) THEN

             ! solution from the previous time step
             q_guess(1:n_vars) = q( 1:n_vars , j , k)

          ELSE

             ! For stages after the first, q_guess is assembled below from
             ! q_fv_cell and the current implicit contribution.

          END IF

          ! New solution at the i_RK step without the implicit  and
          ! semi-implicit term
          q_fv_cell(1:n_vars) = q( 1:n_vars , j , k )                            &
               - dt * (MATMUL( this%divFlux(1:n_eqns,j,k,1:i_RK)                     &
               - this%expl_terms(1:n_eqns,j,k,1:i_RK) , this%a_tilde(1:i_RK) )            &
               - MATMUL( this%NH(1:n_eqns,j,k,1:i_RK) + this%SI_NH(1:n_eqns,j,k,1:i_RK) , &
               this%a_dirk(1:i_RK) ) )

          CALL qc_to_qp(q_fv_cell , qp(1:n_vars+2,j,k) , p_dyn )

          IF ( verbose_level .GE. 2 ) THEN

             WRITE(*,*) 'q_guess',q_guess
             IF ( q_guess(1) .GT. 0.0_wp  ) THEN 

                CALL qc_to_qp( q_guess , qp(1:n_vars+2,j,k) , p_dyn )
                WRITE(*,*) 'q_guess: qp',qp(1:n_vars+2,j,k)

             END IF

          END IF

          adiag_pos:IF ( a_diag .NE. 0.0_wp ) THEN

             pos_thick:IF ( q_fv_cell(1) .GT.  0.0_wp )  THEN

                ! Eval the semi-implicit terms
                ! (terms which non depend on velocity magnitude)
                CALL eval_nh_semi_impl_terms( B_prime_x_geom(j,k) ,             &
                     B_prime_y_geom(j,k) , B_second_xx_geom(j,k) ,              &
                     B_second_xy_geom(j,k) , B_second_yy_geom(j,k) ,            &
                     grav_coeff(j,k) , q_fv_cell ,                              &
                     qp( 1:n_vars , j , k ) , this%SI_NH(1:n_eqns,j,k,i_RK) ,        &
                     Z(j,k) )

                ! Assemble the initial guess for the implicit solver
                q_si(1:n_vars) = q_fv_cell + dt * a_diag *                     &
                     this%SI_NH(1:n_eqns,j,k,i_RK)

                IF ( ( q_fv_cell(2)**2 + q_fv_cell(3)**2 ) .EQ. 0.0_wp ) THEN

                   !Case 1: if the velocity was null, then it must stay null
                   q_si(2:3) = 0.0_wp 

                ELSEIF ( ( q_si(2)*q_fv_cell(2) .LT. 0.0_wp ) .OR.              &
                     ( q_si(3)*q_fv_cell(3) .LT. 0.0_wp ) ) THEN

                   ! If the semi-impl. friction term changed the sign of the 
                   ! velocity then set it to zero
                   q_si(2:3) = 0.0_wp 

                ELSE

                   ! Align the velocity vector with previous one
                   q_si(2:3) = SQRT( q_si(2)**2 + q_si(3)**2 ) *                &
                        q_fv_cell(2:3) / SQRT( q_fv_cell(2)**2                  &
                        + q_fv_cell(3)**2 )

                END IF

                ! Update the semi-implicit term accordingly with the
                ! corrections above
                this%SI_NH(1:n_eqns,j,k,i_RK) = ( q_si(1:n_vars) -                   &
                     q_fv_cell ) / ( dt*a_diag )

                ! Initialize the guess for the NR solver
                q_guess(1:n_vars) = q_si(1:n_vars)


                Rj_not_impl =  ( MATMUL( this%divFlux(1:n_eqns,j,k,1:i_RK-1) -       &
                     this%expl_terms(1:n_eqns,j,k,1:i_RK-1), this%a_tilde(1:i_RK-1) )     &
                     - MATMUL( this%NH(1:n_eqns,j,k,1:i_RK-1)                        &
                     + this%SI_NH(1:n_eqns,j,k,1:i_RK-1) , this%a_dirk(1:i_RK-1) ) )      &
                     - a_diag * this%SI_NH(1:n_eqns,j,k,i_RK)

                ! Solve the implicit system to find the solution at the 
                ! i_RK step of the IMEX RK procedure
                CALL solve_rk_step( q_guess(1:n_vars) , q(1:n_vars,j,k ) ,      &
                     dt, a_diag , Rj_not_impl , B_prime_x_geom(j,k) ,            &
                     B_prime_y_geom(j,k), Z(j,k),                              &
                     newton_iterations, newton_converged, newton_linear_info,   &
                     newton_line_search_failed )

                IF ( ( verbose_level .GE. 1 ) .OR.                             &
                     ( .NOT. newton_converged ) ) THEN

                   !$OMP CRITICAL(newton_diagnostics)

                   IF ( verbose_level .GE. 1 ) THEN
                      newton_calls_step = newton_calls_step + 1
                      newton_iterations_step = newton_iterations_step          &
                           + newton_iterations
                      newton_iterations_max_step = MAX(                        &
                           newton_iterations_max_step, newton_iterations )
                   END IF

                   IF ( .NOT. newton_converged ) THEN
                      newton_failures_step = newton_failures_step + 1
                      IF ( newton_linear_info .NE. 0 )                         &
                           newton_linear_failures_step =                       &
                           newton_linear_failures_step + 1
                      IF ( newton_line_search_failed )                         &
                           newton_line_search_failures_step =                  &
                           newton_line_search_failures_step + 1

                      IF ( verbose_level .GE. 1 ) THEN
                         WRITE(*,*) 'WARNING: Newton solve did not converge'
                         WRITE(*,*)                                           &
                              'cell, RK stage, iterations, linear info:',      &
                              j, k, i_RK, newton_iterations, newton_linear_info
                         WRITE(*,*) 'line search failed:',                     &
                              newton_line_search_failed
                      END IF
                   END IF

                   !$OMP END CRITICAL(newton_diagnostics)

                END IF
                
                IF ( comp_cells_y .EQ. 1 ) THEN

                   q_guess(3) = 0.0_wp

                END IF

                IF ( comp_cells_x .EQ. 1 ) THEN

                   q_guess(2) = 0.0_wp

                END IF

                IF ( rheology_model .EQ. 8 ) THEN
                   
                   this%NH(1:n_eqns,j,k,i_RK) = ( q_guess(1:n_vars)                  &
                        - q_si(1:n_vars) ) / ( dt*a_diag )
                   
                ELSE
                   
                   ! Eval and store the implicit term at the i_RK step
                   CALL eval_implicit_terms( B_prime_x_geom(j,k) ,              &
                        B_prime_y_geom(j,k) , Z(j,k),                          &
                        r_qj = q_guess , r_nh_term_impl = this%NH(1:n_eqns,j,k,i_RK) )
                   
                   IF ( q_si(2)**2 + q_si(3)**2 .EQ. 0.0_wp ) THEN
                      
                      q_guess(2:3) = 0.0_wp 
                      
                   ELSEIF ( ( q_guess(2)*q_si(2) .LE. 0.0_wp ) .AND.            &
                        ( q_guess(3)*q_si(3) .LE. 0.0_wp ) ) THEN
                      
                   ! If the impl. friction term changed the sign of the 
                      ! velocity then set it to zero
                      q_guess(2:3) = 0.0_wp 
                      
                   ELSE
                      
                      ! Align the velocity vector with previous one
                      q_guess(2:3) = SQRT( q_guess(2)**2 + q_guess(3)**2 ) *      &
                           q_si(2:3) / SQRT( q_si(2)**2 + q_si(3)**2 ) 
                      
                   END IF
                   
                END IF
                
             ELSE

                ! If h=0 nothing has to be changed 
                q_guess(1:n_vars) = q_fv_cell
                q_si(1:n_vars) = q_fv_cell
                this%SI_NH(1:n_eqns,j,k,i_RK) = 0.0_wp
                this%NH(1:n_eqns,j,k,i_RK) = 0.0_wp

             END IF pos_thick

          END IF adiag_pos

          IF ( a_diag .NE. 0.0_wp ) THEN

             ! Update the implicit term with correction on the new velocity
             this%NH(1:n_vars,j,k,i_RK) = ( q_guess(1:n_vars) - q_si(1:n_vars))      &
                  / ( dt*a_diag ) 

          END IF

          ! Store the current stage. Previous stage states are no longer
          ! needed here: their evaluated terms are retained in divFlux, NH,
          ! SI_NH and expl_terms.
          this%q_rk( 1:n_vars , j , k ) = q_guess

          IF ( verbose_level .GE. 2 ) THEN

             WRITE(*,*) 'imex_RK_solver: qc',q_guess

             IF ( q_guess(1) .GT. 0.0_wp ) THEN

                CALL qc_to_qp( q_guess , qp(1:n_vars+2,j,k) , p_dyn )
                WRITE(*,*) 'imex_RK_solver: qp',qp(1:n_vars+2,j,k)

             END IF
             
             READ(*,*)

          END IF


          IF ( need_explicit_stage ) THEN
          
             IF ( this%q_rk(1,j,k) .GT. 0.0_wp ) THEN

                CALL qc_to_qp( this%q_rk(1:n_vars,j,k) ,                             &
                     this%qp_rk(1:n_vars+2,j,k) , p_dyn )

             ELSE

                this%qp_rk(1:n_vars+2,j,k) = 0.0_wp
                this%qp_rk(4,j,k) = T_ambient

             END IF

             ! Eval gravity term and radial bottom source terms
             CALL eval_expl_terms( B_prime_x_geom(j,k) , B_prime_y_geom(j,k) ,  &
                  B_second_xx_geom(j,k) , B_second_xy_geom(j,k) ,               &
                  B_second_yy_geom(j,k) , grav_coeff(j,k), d_grav_coeff_dx(j,k),&
                  d_grav_coeff_dy(j,k) ,                                       &
                  this%qp_rk(1:n_vars+2,j,k), this%expl_terms(1:n_eqns,j,k,i_RK), t,      &
                  cell_source_fractions(j,k),                                   &
                  cell_arc_perim(j,k), cell_arc_n_x(j,k), cell_arc_n_y(j,k),    &
                  dx * dy )
  
          END IF

       END DO solve_cells_loop

       !$OMP END DO
       !$OMP END PARALLEL 

       IF ( need_explicit_stage ) THEN

          ! Eval and store the explicit hyperbolic (fluxes) terms
          CALL hyper%evaluate_terms(                                            &
               this%q_rk , this%qp_rk ,                                                   &
               this%divFlux(1:n_eqns,1:comp_cells_x,1:comp_cells_y,i_RK), t,        &
               domain%solve_cells, domain%j_cent, domain%k_cent,                                    &
               domain%solve_interfaces_x, domain%j_stag_x, domain%k_stag_x,                        &
               domain%solve_interfaces_y, domain%j_stag_y, domain%k_stag_y )

       END IF

    END DO runge_kutta

    !$OMP PARALLEL DO private(j,k,p_dyn,alpha_s,solid_excess_roundoff,          &
    !$OMP & residual_cell,q_old_cell)

    assemble_sol:DO l = 1,domain%solve_cells

       j = domain%j_cent(l)
       k = domain%k_cent(l)

       ! q remains equal to Q^n throughout all RK stages. Preserve the old
       ! state locally before overwriting this cell during final assembly.
       q_old_cell = q(1:n_vars,j,k)

       residual_cell = MATMUL( this%divFlux(1:n_eqns,j,k,1:n_RK)                     &
            - this%expl_terms(1:n_eqns,j,k,1:n_RK) , this%omega_tilde ) -                 &
            MATMUL( this%NH(1:n_eqns,j,k,1:n_RK) + this%SI_NH(1:n_eqns,j,k,1:n_RK) ,      &
            this%omega )


       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'cell jk =',j,k
          WRITE(*,*) 'before imex_RK_solver: qc',q_old_cell

          IF ( q_old_cell(1) .GT. 0.0_wp ) THEN

             CALL qc_to_qp(q_old_cell , qp(1:n_vars+2,j,k) , p_dyn )
             WRITE(*,*) 'before imex_RK_solver: qp',qp(1:n_vars+2,j,k)
 
          END IF

       END IF

       IF ( ( SUM(ABS( this%omega_tilde(:)-this%a_tilde_ij(n_RK,:))) .EQ. 0.0_wp  )       &
            .AND. ( SUM(ABS(this%omega(:)-this%a_dirk_ij(n_RK,:))) .EQ. 0.0_wp ) ) THEN

          ! The assembling coeffs are equal to the last step of the RK scheme
          q(1:n_vars,j,k) = this%q_rk(1:n_vars,j,k)

       ELSE

          ! The assembling coeffs are different
          q(1:n_vars,j,k) = q_old_cell - dt*residual_cell

       END IF

       IF ( ANY(ISNAN(q(:,j,k))) ) THEN
          
          WRITE(*,*) 'j,k,n_RK',j,k,n_RK
          WRITE(*,*) 'dt',dt
          WRITE(*,*) 'before imex_RK_solver: qc',q_old_cell
          IF ( q_old_cell(1) .GT. 0.0_wp ) THEN

             CALL qc_to_qp(q_old_cell , qp(1:n_vars+2,j,k) , p_dyn )
             WRITE(*,*) 'before imex_RK_solver: qp',qp(1:n_vars+2,j,k)
             
          END IF
          WRITE(*,*) 'after imex_RK_solver: qc',q(1:n_vars,j,k)
          
          READ(*,*)
          
       END IF
       
       
       ! negative_thickness_check:IF ( q(1,j,k) .LT. 0.0_wp ) THEN
       negative_thickness_check:IF ( q(1,j,k) .LT. EPSILON(1.0_wp) ) THEN

          IF ( q(1,j,k) .GT. -1.0E-7_wp ) THEN

             q(1,j,k) = 0.0_wp
             q(2:n_vars,j,k) = 0.0_wp

          ELSE

             WRITE(*,*) 'j,k,n_RK',j,k,n_RK
             WRITE(*,*) 'dt',dt
             WRITE(*,*) 'before imex_RK_solver: qc',q_old_cell
             IF ( q_old_cell(1) .GT. 0.0_wp ) THEN

                CALL qc_to_qp(q_old_cell , qp(1:n_vars+2,j,k) , p_dyn )
                WRITE(*,*) 'before imex_RK_solver: qp',qp(1:n_vars+2,j,k)

             END IF
             WRITE(*,*) 'after imex_RK_solver: qc',q(1:n_vars,j,k)

             WRITE(*,*) 'divFlux(1,j,k,1:n_RK)',this%divFlux(1,j,k,1:n_RK)

             WRITE(*,*) hyper%H_interface_x(1,j+1,k), hyper%H_interface_x(1,j,k)
             WRITE(*,*) recon%qp_interfaceR(1:n_vars,j,k)
             WRITE(*,*) qp(1:n_vars,j,k)
             WRITE(*,*) recon%qp_interfaceL(1:n_vars,j+1,k)

             WRITE(*,*) 'expl_terms(1,j,k,1:n_RK)',this%expl_terms(1,j,k,1:n_RK)
             WRITE(*,*) 'NH(1,j,k,1:n_RK)',this%NH(1,j,k,1:n_RK)
             WRITE(*,*) 'SI_NH(1,j,k,1:n_RK)',this%SI_NH(1,j,k,1:n_RK)

             WRITE(*,*) 'B_cent(j,k)',B_cent(j,k)

             READ(*,*)

          END IF

       END IF negative_thickness_check

       negative_alpha_check:IF ( ANY(q(5:4+n_solid,j,k) .LT. 0.0_wp ) ) THEN

          IF ( ANY(q(5:4+n_solid,j,k) .LE. -1.0E-7_wp ) ) THEN
             
             WRITE(*,*) 'WARNINIG: negative solid mass'
             WRITE(*,*) 'j,k,n_RK',j,k,n_RK
             WRITE(*,*) 'dt',dt
             WRITE(*,*) 'before imex_RK_solver: qc',q_old_cell
             IF ( q_old_cell(1) .GT. 0.0_wp ) THEN

                CALL qc_to_qp(q_old_cell , qp(1:n_vars+2,j,k) , p_dyn )
                WRITE(*,*) 'before imex_RK_solver: qp',qp(1:n_vars+2,j,k)
                
             END IF
             WRITE(*,*) 'after imex_RK_solver: qc',q(1:n_vars,j,k)

             WRITE(*,*) 'H_interface(1)'
             WRITE(*,*) hyper%H_interface_x(1,j+1,k)/dx*dt, hyper%H_interface_x(1,j,k)/dx*dt
             WRITE(*,*) hyper%H_interface_y(1,j,k+1)/dy*dt, hyper%H_interface_y(1,j,k)/dy*dt
             
             WRITE(*,*) 'H_interface(5)'
             WRITE(*,*) hyper%H_interface_x(5,j+1,k)/dx*dt, hyper%H_interface_x(5,j,k)/dx*dt
             WRITE(*,*) hyper%H_interface_y(5,j,k+1)/dy*dt, hyper%H_interface_y(5,j,k)/dy*dt
             

             READ(*,*)

          ELSE

             WHERE ( ( q(5:4+n_solid,j,k) .LT. 0.0_wp ) .AND.                 &
                  ( q(5:4+n_solid,j,k) .GT. -1.0E-7_wp ) )                   &
                  q(5:4+n_solid,j,k) = 0.0_wp

          END IF
             
       END IF negative_alpha_check

       CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars+2,j,k) , p_dyn )

       IF ( qp(1,j,k) .GT. 1.e-10_wp ) THEN
       
          IF ( alpha_flag ) THEN
             
             alpha_s = SUM(qp(5:4+n_solid,j,k))
             
          ELSE
             
             alpha_s = SUM(qp(5:4+n_solid,j,k)) / qp(1,j,k)
             
          END IF

       ELSE

          alpha_s = 0.0_wp

       END IF
          
       IF ( alpha_s .GT. maximum_solid_packing ) THEN

          IF ( ( alpha_s - maximum_solid_packing ) .LT. 1.e-4_wp ) THEN 

             !q(5:4+n_solid,j,k) = q(5:4+n_solid,j,k) * maximum_solid_packing / alpha_s
                
          ELSE
          
             !WRITE(*,*) 'j,k',j,k
             !WRITE(*,*) 'alpha_s',alpha_s
             
             !WRITE(*,*) 'before imex_RK_solver: qc',q_old_cell
             !WRITE(*,*) 'after imex_RK_solver: qc',q(1:n_vars,j,k)
             
             !CALL qc_to_qp(q_old_cell , qp(1:n_vars+2,j,k) , p_dyn )
             !WRITE(*,*) 'before imex_RK_solver: qp',qp(1:n_vars+2,j,k)
             
             !CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars+2,j,k) , p_dyn )
             !WRITE(*,*) 'after imex_RK_solver: qp',qp(1:n_vars+2,j,k)
             
             !WRITE(*,*) 'H_interface(1)'
             !WRITE(*,*) H_interface_x(1,j+1,k)/dx*dt, H_interface_x(1,j,k)/dx*dt
             !WRITE(*,*) H_interface_y(1,j,k+1)/dy*dt, H_interface_y(1,j,k)/dy*dt
             
             !WRITE(*,*) 'H_interface(5)'
             !WRITE(*,*) H_interface_x(5,j+1,k)/dx*dt, H_interface_x(5,j,k)/dx*dt
             !WRITE(*,*) H_interface_y(5,j,k+1)/dy*dt, H_interface_y(5,j,k)/dy*dt
             
             !WRITE(*,*) 'divFlux(1)',divFlux(1,j,k,1:n_RK)
             !WRITE(*,*) 'expl_terms(1)', expl_terms(1,j,k,1:n_RK)
             !WRITE(*,*) 'NH(1)', NH(1,j,k,1:n_RK)
             !WRITE(*,*) 'SI(1)', SI_NH(1,j,k,1:n_RK) 
             
             !WRITE(*,*) 'divFlux(5)',divFlux(5,j,k,1:n_RK)
             !WRITE(*,*) 'expl_terms(5)', expl_terms(5,j,k,1:n_RK)
             !WRITE(*,*) 'NH(5)', NH(5,j,k,1:n_RK)
             !WRITE(*,*) 'SI(5)', SI_NH(5,j,k,1:n_RK)
             
             !READ(*,*)

          END IF
             
       END IF
          
       
       IF ( qp(4,j,k) .LT. 273.0_wp ) THEN
          
          WRITE(*,*) 'temperature check'
          WRITE(*,*) j,k
          WRITE(*,*) 'qp new',qp(1:n_vars+2,j,k)
          WRITE(*,*) 'qc new',q(1:n_vars,j,k)

          CALL qc_to_qp(q_old_cell , qp(1:n_vars+2,j,k) , p_dyn )
          WRITE(*,*) j,k
          WRITE(*,*) 'qp old',qp(1:n_vars+2,j,k)
          WRITE(*,*) 'qc old',q_old_cell

          WRITE(*,*) 'H_interface(4)'
          WRITE(*,*) hyper%H_interface_x(4,j+1,k)/dx*dt, hyper%H_interface_x(4,j,k)/dx*dt
          WRITE(*,*) hyper%H_interface_y(4,j,k+1)/dy*dt, hyper%H_interface_y(4,j,k)/dy*dt

          WRITE(*,*) hyper%H_interface_y(:,j,k)/dy*dt
          READ(*,*)

       END IF
          
       IF ( SUM(q(5:4+n_solid,j,k)) .GT. q(1,j,k) ) THEN

          ! Fortran does not guarantee short-circuit evaluation of .OR.;
          ! evaluate the relative excess only when the denominator is safe.
          IF ( q(1,j,k) .LT. EPSILON(1.0_wp) ) THEN
             solid_excess_roundoff = .TRUE.
          ELSE
             solid_excess_roundoff = ( ( SUM(q(5:4+n_solid,j,k))             &
                  - q(1,j,k) ) / q(1,j,k) .LT. 1.0E-10_wp )
          END IF

          IF ( solid_excess_roundoff ) THEN

             CALL qc_to_qp(q_old_cell , qp(1:n_vars+2,j,k) , p_dyn )

             q(5:4+n_solid,j,k) = q(5:4+n_solid,j,k)                            &
                  / SUM(q(5:4+n_solid,j,k)) * q(1,j,k)

          ELSE

             WRITE(*,*) 'WARNING:SUM(qsolid)>q1',SUM(q(5:4+n_solid,j,k))-q(1,j,k)
             
             WRITE(*,*) 'j,k,n_RK',j,k,n_RK
             WRITE(*,*) 'dt',dt
             WRITE(*,*) ' B_cent(j,k)', B_cent(j,k)
             WRITE(*,*) 'before imex_RK_solver: qc',q_old_cell
             IF ( q_old_cell(1) .GT. 0.0_wp ) THEN

                CALL qc_to_qp(q_old_cell , qp(1:n_vars+2,j,k) , p_dyn )
                WRITE(*,*) 'before imex_RK_solver: qp',qp(1:n_vars+2,j,k)

             END IF
             WRITE(*,*) 'after imex_RK_solver: qc',q(1:n_vars,j,k)
             
             IF ( q(1,j,k) .GT. 0.0_wp ) THEN

                CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars+2,j,k) , p_dyn )
                WRITE(*,*) 'after imex_RK_solver: qp',qp(1:n_vars+2,j,k)

             END IF

             WRITE(*,*) 'H_interface(1)'
             WRITE(*,*) hyper%H_interface_x(1,j+1,k)/dx*dt, hyper%H_interface_x(1,j,k)/dx*dt
             WRITE(*,*) hyper%H_interface_y(1,j,k+1)/dy*dt, hyper%H_interface_y(1,j,k)/dy*dt
             
             WRITE(*,*) 'H_interface(5)'
             WRITE(*,*) hyper%H_interface_x(5,j+1,k)/dx*dt, hyper%H_interface_x(5,j,k)/dx*dt
             WRITE(*,*) hyper%H_interface_y(5,j,k+1)/dy*dt, hyper%H_interface_y(5,j,k)/dy*dt

             WRITE(*,*) 'divFlux(1)',this%divFlux(1,j,k,1:n_RK)
             WRITE(*,*) 'expl_terms(1)', this%expl_terms(1,j,k,1:n_RK)
             WRITE(*,*) 'NH(1)', this%NH(1,j,k,1:n_RK)
             WRITE(*,*) 'SI(1)', this%SI_NH(1,j,k,1:n_RK)

             WRITE(*,*) 'divFlux(5)',this%divFlux(5,j,k,1:n_RK)
             WRITE(*,*) 'expl_terms(5)', this%expl_terms(5,j,k,1:n_RK)
             WRITE(*,*) 'NH(5)', this%NH(5,j,k,1:n_RK)
             WRITE(*,*) 'SI(5)', this%SI_NH(5,j,k,1:n_RK)
             

             READ(*,*)

          END IF

          IF ( verbose_level .GE. 1 ) THEN

             WRITE(*,*) 'h new',q(1,j,k) 
             READ(*,*)

          END IF

       END IF

       IF ( B_nodata(j,k) ) q(:,j,k) = 0.0_wp

    END DO assemble_sol

    !$OMP END PARALLEL DO

    IF ( ( verbose_level .GE. 1 ) .AND. ( newton_calls_step .GT. 0 ) ) THEN
       WRITE(*,*) 'Newton solves:',newton_calls_step
       WRITE(*,*) 'Newton iterations average/max:',                            &
            REAL(newton_iterations_step,wp) / REAL(newton_calls_step,wp),       &
            newton_iterations_max_step
       WRITE(*,*) 'Newton failures / linear failures:',                        &
            newton_failures_step,newton_linear_failures_step
       WRITE(*,*) 'Newton line-search failures:',                              &
            newton_line_search_failures_step
    ELSEIF ( newton_failures_step .GT. 0 ) THEN
       WRITE(*,*) 'WARNING: Newton failures / linear / line search:',           &
            newton_failures_step,newton_linear_failures_step,                  &
            newton_line_search_failures_step
    END IF
     
    RETURN

  END SUBROUTINE imex_RK_solver

END MODULE time_integration_2d
