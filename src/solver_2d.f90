!********************************************************************************
!> \brief Numerical solver
!
!> This module contains the variables and the subroutines for the 
!> numerical solution of the equations.  
!
!> \date 07/10/2016
!> @author 
!> Mattia de' Michieli Vitturi
!
!********************************************************************************
MODULE solver_2d

  ! external variables

  USE constitutive_2d, ONLY : implicit_flag, implicit_map, rheology_model
  USE constitutive_2d, ONLY : T_ambient
    
  USE geometry_2d, ONLY : comp_cells_x,comp_cells_y,comp_cells_xy
  USE geometry_2d, ONLY : comp_interfaces_x,comp_interfaces_y

  USE geometry_2d, ONLY : B_cent

  USE geometry_2d, ONLY : B_prime_x , B_prime_y  
  USE geometry_2d, ONLY : B_second_xx , B_second_xy , B_second_yy

  USE geometry_2d, ONLY : B_prime_x_geom , B_prime_y_geom  
  USE geometry_2d, ONLY : B_second_xx_geom , B_second_xy_geom , B_second_yy_geom

  USE geometry_2d, ONLY : grav_coeff
  USE geometry_2d, ONLY : grav_coeff_stag_x , grav_coeff_stag_y
  
  USE geometry_2d, ONLY : d_grav_coeff_dx , d_grav_coeff_dy
  USE geometry_2d, ONLY : source_cell
  USE geometry_2d, ONLY : cell_source_fractions
  USE geometry_2d, ONLY : cell_arc_perim , cell_arc_n_x , cell_arc_n_y

  USE parameters_2d, ONLY : wp , sp

  USE parameters_2d, ONLY : n_eqns , n_vars , n_nh , n_solid
  USE parameters_2d, ONLY : n_RK
  USE parameters_2d, ONLY : verbose_level
  USE parameters_2d, ONLY : radial_source_flag , bottom_radial_source_flag
  USE parameters_2d, ONLY : stochastic_flag

  USE parameters_2d, ONLY : idx_h, idx_hu, idx_hv, idx_T, idx_alfas_first,      &
       idx_alfas_last, idx_addGas_first, idx_addGas_last, idx_stoch, idx_pore

  USE parameters_2d, ONLY : idx_totMassEqn, idx_uEqn, idx_vEqn, idx_engyEqn,    &
       idx_solidEqn_first, idx_solidEqn_last, idx_addGasEqn_first,              &
       idx_addGasEqn_last, idx_stochEqn, idx_poreEqn
    
  ! external procedures
  USE geometry_2d, ONLY : dx , dy , one_by_dx , one_by_dy 

  USE OMP_LIB

  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver,                 &
       finalize_nonlinear_solver, solve_rk_step

  USE reconstruction_2d, ONLY : initialize_reconstruction,                    &
       finalize_reconstruction, reconstruction
  USE reconstruction_2d, ONLY : q_interfaceL, q_interfaceR,                   &
       q_interfaceB, q_interfaceT
  USE reconstruction_2d, ONLY : qp_interfaceL, qp_interfaceR,                 &
       qp_interfaceB, qp_interfaceT
  USE reconstruction_2d, ONLY : diverg_interfaceL, diverg_interfaceR,         &
       diverg_interfaceB, diverg_interfaceT

  IMPLICIT none

  !> time
  REAL(wp) :: t

  !> Conservative variables
  REAL(wp), ALLOCATABLE :: q(:,:,:)        
  !> Map of positive thickness 
  LOGICAL, ALLOCATABLE :: hpos(:,:)        
  !> Map of positive thickness at previous output step
  LOGICAL, ALLOCATABLE :: hpos_old(:,:)        


  !> Maximum over time of thickness
  REAL(wp), ALLOCATABLE :: hmax(:,:)

  !> Maximum over time of dynamic pressure
  REAL(wp), ALLOCATABLE :: pdynmax(:,:)

  !> Maximum over time of dynamic velocity
  REAL(wp), ALLOCATABLE :: mod_vel_max(:,:)

  !> Maximum over time of thickness
  LOGICAL, ALLOCATABLE :: vuln_table(:,:,:)

  LOGICAL, ALLOCATABLE :: thck_table(:,:)

  LOGICAL, ALLOCATABLE :: pdyn_table(:,:)

  !> Local speeds at the left of the x-interface
  REAL(wp), ALLOCATABLE :: a_interface_xNeg(:,:,:)
  !> Local speeds at the right of the x-interface
  REAL(wp), ALLOCATABLE :: a_interface_xPos(:,:,:)
  !> Local speeds at the bottom of the y-interface
  REAL(wp), ALLOCATABLE :: a_interface_yNeg(:,:,:)
  !> Local speeds at the top of the y-interface
  REAL(wp), ALLOCATABLE :: a_interface_yPos(:,:,:)
  !> Semidiscrete numerical interface fluxes 
  REAL(wp), ALLOCATABLE :: H_interface_x(:,:,:)
  !> Semidiscrete numerical interface fluxes 
  REAL(wp), ALLOCATABLE :: H_interface_y(:,:,:)
  !> Physical variables (\f$\alpha_1, p_1, p_2, \rho u, w, T\f$)
  REAL(wp), ALLOCATABLE :: qp(:,:,:)

  !> Array defining fraction of cells affected by source term
  REAL(wp), ALLOCATABLE :: source_xy(:,:)

  REAL(wp), ALLOCATABLE :: solve_mask_time(:,:)


  LOGICAL, ALLOCATABLE :: solve_mask(:,:)
  LOGICAL, ALLOCATABLE :: solve_mask_temp(:,:)
  LOGICAL, ALLOCATABLE :: solve_mask_x(:,:)
  LOGICAL, ALLOCATABLE :: solve_mask_y(:,:)

  INTEGER :: solve_cells
  INTEGER :: solve_interfaces_x
  INTEGER :: solve_interfaces_y

  !> Time step
  REAL(wp) :: dt

  INTEGER :: i_RK           !< loop counter for the RK iteration

  !> Butcher Tableau for the explicit part of the Runge-Kutta scheme
  REAL(wp), ALLOCATABLE :: a_tilde_ij(:,:)
  !> Butcher Tableau for the implicit part of the Runge-Kutta scheme
  REAL(wp), ALLOCATABLE :: a_dirk_ij(:,:)

  !> Coefficients for the explicit part of the Runge-Kutta scheme
  REAL(wp), ALLOCATABLE :: omega_tilde(:)

  !> Coefficients for the implicit part of the Runge-Kutta scheme
  REAL(wp), ALLOCATABLE :: omega(:)

  !> Explicit coeff. for the hyperbolic part for a single step of the R-K scheme
  REAL(wp), ALLOCATABLE :: a_tilde(:)

  !> Explicit coeff. for the non-hyp. part for a single step of the R-K scheme
  REAL(wp), ALLOCATABLE :: a_dirk(:)

  !> Implicit coeff. for the non-hyp. part for a single step of the R-K scheme
  REAL(wp) :: a_diag

  !> Conservative solution at the current Runge-Kutta stage
  REAL(wp), ALLOCATABLE :: q_rk(:,:,:)

  !> Physical solution at the current Runge-Kutta stage
  REAL(wp), ALLOCATABLE :: qp_rk(:,:,:)

  !> Intermediate hyperbolic terms of the Runge-Kutta scheme
  REAL(wp), ALLOCATABLE :: divFlux(:,:,:,:)

  !> Intermediate non-hyperbolic terms of the Runge-Kutta scheme
  REAL(wp), ALLOCATABLE :: NH(:,:,:,:)

  !> Intermediate semi-implicit non-hyperbolic terms of the Runge-Kutta scheme
  REAL(wp), ALLOCATABLE :: SI_NH(:,:,:,:)

  !> Intermediate explicit terms of the Runge-Kutta scheme
  REAL(wp), ALLOCATABLE :: expl_terms(:,:,:,:)

  INTEGER, ALLOCATABLE :: j_cent(:)
  INTEGER, ALLOCATABLE :: k_cent(:)

  INTEGER, ALLOCATABLE :: j_stag_x(:)
  INTEGER, ALLOCATABLE :: k_stag_x(:)

  INTEGER, ALLOCATABLE :: j_stag_y(:)
  INTEGER, ALLOCATABLE :: k_stag_y(:)

  !> Stochastic Noise
  REAL(wp), ALLOCATABLE :: Z(:,:)
  !> Array for kernel
  REAL(wp), ALLOCATABLE :: conv_kernel(:,:)
  !> Friction values at cells
  REAL(wp), ALLOCATABLE :: fric_array(:,:)
  
CONTAINS

  !******************************************************************************
  !> \brief Memory allocation
  !
  !> This subroutine allocate the memory for the variables of the 
  !> solver module.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE allocate_solver_variables

    USE parameters_2d, ONLY : n_thickness_levels , n_dyn_pres_levels
    
    IMPLICIT NONE

    REAL(wp) :: gamma , delta

    INTEGER :: i,j

    ALLOCATE( q( n_vars , comp_cells_x , comp_cells_y ) )

    ALLOCATE( hpos( comp_cells_x , comp_cells_y ) , hpos_old ( comp_cells_x ,   &
         comp_cells_y ) )

    ALLOCATE( qp( n_vars+2 , comp_cells_x , comp_cells_y ) )

    q(1:n_vars,1:comp_cells_x,1:comp_cells_y) = 0.0_wp
    qp(1:n_vars+2,1:comp_cells_x,1:comp_cells_y) = 0.0_wp
    qp(4,1:comp_cells_x,1:comp_cells_y) = T_ambient

    ALLOCATE( hmax( comp_cells_x , comp_cells_y ) )
    ALLOCATE( pdynmax( comp_cells_x , comp_cells_y ) )
    ALLOCATE( mod_vel_max( comp_cells_x , comp_cells_y ) )

    ALLOCATE( vuln_table( n_thickness_levels * n_dyn_pres_levels ,              &
         comp_cells_x , comp_cells_y ) )

    ALLOCATE( thck_table(comp_cells_x , comp_cells_y) )

    ALLOCATE( pdyn_table(comp_cells_x , comp_cells_y) )

    CALL initialize_reconstruction

    ALLOCATE( a_interface_xNeg( n_eqns , comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( a_interface_xPos( n_eqns , comp_interfaces_x, comp_cells_y ) )

    a_interface_xNeg = 0.0_wp
    a_interface_xPos = 0.0_wp
    
    ALLOCATE( H_interface_x( n_eqns , comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( H_interface_y( n_eqns , comp_cells_x, comp_interfaces_y ) )


    ALLOCATE( a_interface_yNeg( n_eqns , comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( a_interface_yPos( n_eqns , comp_cells_x, comp_interfaces_y ) )

    a_interface_yNeg = 0.0_wp
    a_interface_yPos = 0.0_wp
   
    ALLOCATE( solve_mask_time( comp_cells_x , comp_cells_y ) )
    ALLOCATE( solve_mask( comp_cells_x , comp_cells_y ) )
    ALLOCATE( solve_mask_temp( comp_cells_x , comp_cells_y ) )

    solve_mask_time(1:comp_cells_x,1:comp_cells_y) = 0.0_wp

    solve_mask(1,1:comp_cells_y) = .TRUE.
    solve_mask(comp_cells_x,1:comp_cells_y) = .TRUE.
    solve_mask(1:comp_cells_x,1) = .TRUE.
    solve_mask(1:comp_cells_x,comp_cells_y) = .TRUE.

    ALLOCATE( solve_mask_x( comp_interfaces_x , comp_cells_y ) )
    ALLOCATE( solve_mask_y( comp_cells_x , comp_interfaces_y ) )

    ALLOCATE( source_xy( comp_cells_x , comp_cells_y ) )


    ALLOCATE( a_tilde_ij(n_RK,n_RK) )
    ALLOCATE( a_dirk_ij(n_RK,n_RK) )
    ALLOCATE( omega_tilde(n_RK) )
    ALLOCATE( omega(n_RK) )


    CALL initialize_nonlinear_solver

    ! Initialize the coefficients for the IMEX Runge-Kutta scheme
    ! Please note that with respect to the schemes described in Pareschi & Russo 
    ! (2000) we do not have the coefficient vectors c_tilde and c, because the 
    ! explicit and implicit terms do not depend explicitly on time.

    ! Explicit part coefficients (a_tilde_ij=0 for j>=i)
    a_tilde_ij = 0.0_wp

    ! Weight coefficients of the explicit part in the final assemblage
    omega_tilde = 0.0_wp

    ! Implicit part coefficients (a_dirk_ij=0 for j>i)
    a_dirk_ij = 0.0_wp

    ! Weight coefficients of the explicit part in the final assemblage
    omega = 0.0_wp

    gamma = 1.0_wp - 1.0_wp / SQRT(2.0_wp)
    delta = 1.0_wp - 1.0_wp / ( 2.0_wp * gamma )

    IF ( n_RK .EQ. 1 ) THEN

       a_tilde_ij(1,1) = 1.0_wp

       omega_tilde(1) = 1.0_wp

       a_dirk_ij(1,1) = 0.0_wp

       omega(1) = 0.0_wp

    ELSEIF ( n_RK .EQ. 2 ) THEN

       a_tilde_ij(2,1) = 1.0_wp

       omega_tilde(1) = 1.0_wp
       omega_tilde(2) = 0.0_wp

       a_dirk_ij(2,2) = 1.0_wp

       omega(1) = 0.0_wp
       omega(2) = 1.0_wp

    ELSEIF ( n_RK .EQ. 3 ) THEN

       ! Tableau for the IMEX-SSP(3,3,2) Stiffly Accurate Scheme
       ! from Pareschi & Russo (2005), Table IV

       a_tilde_ij(2,1) = 0.5_wp
       a_tilde_ij(3,1) = 0.5_wp
       a_tilde_ij(3,2) = 0.5_wp

       omega_tilde(1) =  1.0_wp / 3.0_wp
       omega_tilde(2) =  1.0_wp / 3.0_wp
       omega_tilde(3) =  1.0_wp / 3.0_wp

       a_dirk_ij(1,1) = 0.25_wp
       a_dirk_ij(2,2) = 0.25_wp
       a_dirk_ij(3,1) = 1.0_wp / 3.0_wp
       a_dirk_ij(3,2) = 1.0_wp / 3.0_wp
       a_dirk_ij(3,3) = 1.0_wp / 3.0_wp

       omega(1) =  1.0_wp / 3.0_wp
       omega(2) =  1.0_wp / 3.0_wp
       omega(3) =  1.0_wp / 3.0_wp

    ELSEIF ( n_RK .EQ. 4 ) THEN

       ! LRR(3,2,2) from Table 3 in Pareschi & Russo (2000)

       a_tilde_ij(2,1) = 0.5_wp
       a_tilde_ij(3,1) = 1.0_wp / 3.0_wp
       a_tilde_ij(4,2) = 1.0_wp

       omega_tilde(1) = 0.0_wp
       omega_tilde(2) = 1.0_wp
       omega_tilde(3) = 0.0_wp
       omega_tilde(4) = 0.0_wp

       a_dirk_ij(2,2) = 0.5_wp
       a_dirk_ij(3,3) = 1.0_wp / 3.0_wp
       a_dirk_ij(4,3) = 0.75_wp
       a_dirk_ij(4,4) = 0.25_wp

       omega(1) = 0.0_wp
       omega(2) = 0.0_wp
       omega(3) = 0.75_wp
       omega(4) = 0.25_wp

    END IF

    ALLOCATE( a_tilde(n_RK) )
    ALLOCATE( a_dirk(n_RK) )

    ALLOCATE( q_rk( n_vars , comp_cells_x , comp_cells_y ) )
    ALLOCATE( qp_rk( n_vars+2 , comp_cells_x , comp_cells_y ) )
    ALLOCATE( divFlux( n_eqns , comp_cells_x , comp_cells_y , n_RK ) )
    ALLOCATE( NH( n_eqns , comp_cells_x , comp_cells_y , n_RK ) )
    ALLOCATE( SI_NH( n_eqns , comp_cells_x , comp_cells_y , n_RK ) )
    ALLOCATE( expl_terms( n_eqns , comp_cells_x , comp_cells_y , n_RK ) )

    comp_cells_xy = comp_cells_x * comp_cells_y

    ALLOCATE( j_cent( comp_cells_xy ) )
    ALLOCATE( k_cent( comp_cells_xy ) )

    ALLOCATE( j_stag_x( comp_interfaces_x * comp_cells_y ) )
    ALLOCATE( k_stag_x( comp_interfaces_x * comp_cells_y ) )

    ALLOCATE( j_stag_y( comp_cells_x * comp_interfaces_y ) )
    ALLOCATE( k_stag_y( comp_cells_x * comp_interfaces_y ) )

    ! Allocate array containing the stochastic noise.
    ! Allocated unconditionally: Z(j,k) is passed as a scalar actual
    ! argument to the semi-implicit and implicit routines whatever
    ! stochastic_flag is, so leaving it unallocated is invalid. With the
    ! flag off the values stay zero and have no effect.
    ALLOCATE ( Z(comp_cells_x , comp_cells_y) )
    Z(1:comp_cells_x,1:comp_cells_y) = 0.0_wp
    
    ! Allocate array containing the friction values.
    ! Allocated unconditionally for the same reason as Z above:
    ! fric_array(j,k) is passed as a scalar actual argument for every
    ! rheology model, not only those that populate it.
    ALLOCATE (fric_array(comp_cells_x , comp_cells_y))
    fric_array(1:comp_cells_x,1:comp_cells_y) = 0.0_wp
    
    WRITE(*,*) 'ALLOCATION OF ARRAYS COMPLETED'
    
    RETURN
    
  END SUBROUTINE allocate_solver_variables

  !******************************************************************************
  !> \brief Memory deallocation
  !
  !> This subroutine de-allocate the memory for the variables of the 
  !> solver module.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE deallocate_solver_variables

    DEALLOCATE( q , hpos , hpos_old )

    DEALLOCATE( hmax , pdynmax , mod_vel_max )

    DEALLOCATE( vuln_table )

    DEALLOCATE( thck_table ,  pdyn_table )

    CALL finalize_reconstruction

    DEALLOCATE( a_interface_xNeg )
    DEALLOCATE( a_interface_xPos )
    DEALLOCATE( a_interface_yNeg )
    DEALLOCATE( a_interface_yPos )

    DEALLOCATE( H_interface_x )
    DEALLOCATE( H_interface_y )

    DEALLOCATE( solve_mask_time )
    DEALLOCATE( solve_mask )
    DEALLOCATE( solve_mask_temp )
    DEALLOCATE( solve_mask_x )
    DEALLOCATE( solve_mask_y )

    DEALLOCATE( qp )

    DEALLOCATE( source_xy )

    DEALLOCATE( a_tilde_ij )
    DEALLOCATE( a_dirk_ij )
    DEALLOCATE( omega_tilde )
    DEALLOCATE( omega )

    DEALLOCATE( a_tilde )
    DEALLOCATE( a_dirk )

    DEALLOCATE( q_rk )
    DEALLOCATE( qp_rk )
    DEALLOCATE( divFlux )
    DEALLOCATE( NH )
    DEALLOCATE( SI_NH )
    DEALLOCATE( expl_terms )

    CALL finalize_nonlinear_solver

    DEALLOCATE( j_cent , k_cent )
    DEALLOCATE ( j_stag_x , k_stag_x )
    DEALLOCATE ( j_stag_y , k_stag_y )

    IF ( ALLOCATED(implicit_flag) ) DEALLOCATE(implicit_flag)
    IF ( ALLOCATED(implicit_map) ) DEALLOCATE(implicit_map)
    IF ( ALLOCATED(Z) ) DEALLOCATE(Z)
    IF ( ALLOCATED(conv_kernel) ) DEALLOCATE(conv_kernel)
    IF ( ALLOCATED(fric_array) ) DEALLOCATE(fric_array)

    
    RETURN
    
  END SUBROUTINE deallocate_solver_variables


  !******************************************************************************
  !> \brief Masking of cells to solve
  !
  !> This subroutine compute a 2D array of logicals defining the cells where the
  !> systems of equations have to be solved. It is defined according to the 
  !> positive thickness in the cell and in the neighbour cells
  !
  !> \date 20/04/2017
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE check_solve(solve_all)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: solve_all
    
    INTEGER :: i,j,k

    !$OMP PARALLEL
         
    IF ( solve_all ) THEN

       !$OMP WORKSHARE
       solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) = .TRUE. 
       !$OMP END WORKSHARE

    ELSE
       
       !$OMP WORKSHARE
       solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) = .FALSE.
       !$OMP END WORKSHARE

    END IF
    !$OMP BARRIER
    
    !$OMP WORKSHARE
    WHERE ( ( q(1,2:comp_cells_x-1,2:comp_cells_y-1) .GT. 0.0_wp ) .AND.        &
         ( solve_mask_time(2:comp_cells_x-1,2:comp_cells_y-1) .LE. t ) )       & 
         solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) = .TRUE.
    !$OMP END WORKSHARE
    
    !$OMP BARRIER

    IF ( bottom_radial_source_flag ) THEN

       !$OMP WORKSHARE
       WHERE ( cell_source_fractions .GT. 0.0_wp ) solve_mask = .TRUE.
       !$OMP END WORKSHARE
       
       !$OMP BARRIER

    END IF

    IF ( radial_source_flag ) THEN
             
       !$OMP DO private(j,k)
    
       DO k = 2,comp_cells_y-1
   
          DO j = 2,comp_cells_x-1

             IF ( source_cell(j,k) .EQ. 2 ) THEN
                
                solve_mask(j,k) = .TRUE.
                
             END IF

          END DO
          
       END DO

       !$OMP END DO

    END IF

    !$OMP BARRIER
    !$OMP MASTER

    DO i = 1,n_RK

       solve_mask_temp = solve_mask 
       
       ! solution domain is extended to neighbours of positive-mass cells
       solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) =                          &
            solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) .OR.                  &
            solve_mask_temp(1:comp_cells_x-2,2:comp_cells_y-1)
       
       solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) =                          &
            solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) .OR.                  &
            solve_mask_temp(3:comp_cells_x,2:comp_cells_y-1)
       
       solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) =                          &
            solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) .OR.                  &
            solve_mask_temp(2:comp_cells_x-1,1:comp_cells_y-2)
       
       solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) =                          &
            solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) .OR.                  &
            solve_mask_temp(2:comp_cells_x-1,3:comp_cells_y) 
       
    END DO

    !$OMP END MASTER
    !$OMP BARRIER

    !$OMP DO private(j,k)
    
    DO k = 1,comp_cells_y
       
       DO j = 1,comp_cells_x
          
          IF ( radial_source_flag ) THEN
             
             IF ( source_cell(j,k) .EQ. 1 ) solve_mask(j,k) = .FALSE.
             
          END IF
          
       END DO
       
    END DO
    
    !$OMP END DO
    !$OMP WORKSHARE

    solve_mask_x(1:comp_interfaces_x,1:comp_cells_y) = .FALSE.
    solve_mask_y(1:comp_cells_x,1:comp_interfaces_y) = .FALSE.

    !$OMP END WORKSHARE

    !$OMP END PARALLEL

    !----- check for cells where computation is needed
    i = 0
        
    DO k = 1,comp_cells_y

       DO j = 1,comp_cells_x

          IF ( solve_mask(j,k) ) THEN

             i = i+1
             j_cent(i) = j
             k_cent(i) = k

             solve_mask_x(j,k) = .TRUE.
             solve_mask_x(j+1,k) = .TRUE.
             solve_mask_y(j,k) = .TRUE.
             solve_mask_y(j,k+1) = .TRUE.

          END IF

       END DO

    END DO

    solve_cells = i
    
    !----- check for y-interfaces where computation is needed
    i = 0
        
    DO k = 1,comp_cells_y

       DO j = 1,comp_interfaces_x

          IF ( solve_mask_x(j,k) ) THEN

             i = i+1
             j_stag_x(i) = j
             k_stag_x(i) = k

          END IF

       END DO

    END DO

    solve_interfaces_x = i
   
    !----- check for y-interfaces where computation is needed

    i = 0
        
    DO k = 1,comp_interfaces_y

       DO j = 1,comp_cells_x

          IF ( solve_mask_y(j,k) ) THEN

             i = i+1
             j_stag_y(i) = j
             k_stag_y(i) = k

          END IF

       END DO

    END DO

    solve_interfaces_y = i

    RETURN

  END SUBROUTINE check_solve

  !*****************************************************************************
  !> \brief Time-step computation
  !
  !> This subroutine evaluate the maximum time step according to the CFL
  !> condition. The local speed are evaluated with the characteristic
  !> polynomial of the Jacobian of the fluxes.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !*****************************************************************************

  SUBROUTINE timestep

    ! External variables
    USE geometry_2d, ONLY : dx,dy
    USE parameters_2d, ONLY : max_dt , cfl

    USE constitutive_2d, ONLY : qc_to_qp

    IMPLICIT none

    INTEGER :: j,k,l          !< loop counter

    REAL(wp) :: max_a_x
    REAL(wp) :: max_a_y
    REAL(wp) p_dyn

    dt = max_dt

    IF ( cfl .NE. -1.0_wp ) THEN

       !$OMP PARALLEL DO private(j,k,p_dyn)

       DO l = 1,solve_cells

          j = j_cent(l)
          k = k_cent(l)

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
       CALL reconstruction( q, qp, t, solve_cells, j_cent, k_cent )

       ! Compute the max/min eigenvalues at the interfaces
       CALL eval_speeds

       max_a_x = 0.0_wp
       max_a_y = 0.0_wp

       ! The minimum CFL step over all active cells is determined by the
       ! maximum characteristic speed on their adjacent interfaces.  Compute
       ! those two maxima directly, avoiding full-domain scratch arrays and
       ! an atomic update of dt for every cell.
       !$OMP PARALLEL DO private(j,k) reduction(max:max_a_x,max_a_y)
       DO l = 1,solve_cells

          j = j_cent(l)
          k = k_cent(l)

          max_a_x = MAX( max_a_x,                                               &
               MAXVAL(a_interface_xPos(1:n_vars,j,k)),                         &
               MAXVAL(-a_interface_xNeg(1:n_vars,j,k)),                        &
               MAXVAL(a_interface_xPos(1:n_vars,j+1,k)),                       &
               MAXVAL(-a_interface_xNeg(1:n_vars,j+1,k)) )

          max_a_y = MAX( max_a_y,                                               &
               MAXVAL(a_interface_yPos(1:n_vars,j,k)),                         &
               MAXVAL(-a_interface_yNeg(1:n_vars,j,k)),                        &
               MAXVAL(a_interface_yPos(1:n_vars,j,k+1)),                       &
               MAXVAL(-a_interface_yNeg(1:n_vars,j,k+1)) )

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

  SUBROUTINE imex_RK_solver

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

    REAL(wp) :: q_si(n_vars) !< solution after the semi-implicit step
    REAL(wp) :: q_guess(n_vars) !< initial guess for the solution of the RK step
    REAL(wp) :: q_fv_cell(n_vars) !< finite-volume state for the current cell
    REAL(wp) :: residual_cell(n_vars) !< final RK residual for the current cell
    REAL(wp) :: q_old_cell(n_vars) !< state before the final RK assembly
    INTEGER :: j,k,l            !< loop counter over the grid volumes
    REAL(wp) :: Rj_not_impl(n_eqns)

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
    DO l = 1,solve_cells

       j = j_cent(l)
       k = k_cent(l)

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*) 'solver, imex_RK_solver: j,k',j,k
          !READ(*,*)
          
       END IF

       ! Initialization of the variables for the Runge-Kutta scheme
       q_rk( 1:n_vars , j , k ) = 0.0_wp
       qp_rk( 1:n_vars+2 , j , k ) = 0.0_wp
       qp_rk( 4 , j , k ) = T_ambient
       

       divFlux(1:n_eqns , j , k , 1:n_RK ) = 0.0_wp
       NH( 1:n_eqns, j , k , 1:n_RK ) = 0.0_wp
       SI_NH( 1:n_eqns , j , k , 1:n_RK ) = 0.0_wp
       expl_terms(1:n_eqns , j , k , 1:n_RK) = 0.0_wp
       
    END DO
    !$OMP END DO

    !$OMP END PARALLEL

    runge_kutta:DO i_RK = 1,n_RK

       IF ( verbose_level .GE. 1 ) WRITE(*,*) 'solver, imex_RK_solver: i_RK',i_RK

       ! An explicit stage is required not only when it contributes to the
       ! final RK assembly, but also when a later stage depends on it.
       need_explicit_stage = ( omega_tilde(i_RK) .NE. 0.0_wp )

       IF ( i_RK .LT. n_RK ) THEN
          need_explicit_stage = need_explicit_stage .OR.                       &
               ANY( a_tilde_ij(i_RK+1:n_RK,i_RK) .NE. 0.0_wp )
       END IF

       ! define the explicits coefficients for the i-th step of the Runge-Kutta
       a_tilde = 0.0_wp
       a_dirk = 0.0_wp

       ! in the first step of the RK scheme all the coefficients remain to 0
       a_tilde(1:i_RK-1) = a_tilde_ij(i_RK,1:i_RK-1)
       a_dirk(1:i_RK-1) = a_dirk_ij(i_RK,1:i_RK-1)

       ! define the implicit coefficient for the i-th step of the Runge-Kutta
       a_diag = a_dirk_ij(i_RK,i_RK)

       !$OMP PARALLEL 
       !$OMP DO schedule(guided)                                                &
       !$OMP & private(j,k,q_guess,q_si,q_fv_cell,Rj_not_impl,p_dyn,           &
       !$OMP & newton_iterations,newton_linear_info,newton_converged,          &
       !$OMP & newton_line_search_failed)

       solve_cells_loop:DO l = 1,solve_cells

          j = j_cent(l)
          k = k_cent(l)

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
               - dt * (MATMUL( divFlux(1:n_eqns,j,k,1:i_RK)                     &
               - expl_terms(1:n_eqns,j,k,1:i_RK) , a_tilde(1:i_RK) )            &
               - MATMUL( NH(1:n_eqns,j,k,1:i_RK) + SI_NH(1:n_eqns,j,k,1:i_RK) , &
               a_dirk(1:i_RK) ) )

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
                     qp( 1:n_vars , j , k ) , SI_NH(1:n_eqns,j,k,i_RK) ,        &
                     Z(j,k), fric_array(j,k) )

                ! Assemble the initial guess for the implicit solver
                q_si(1:n_vars) = q_fv_cell + dt * a_diag *                     &
                     SI_NH(1:n_eqns,j,k,i_RK)

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
                SI_NH(1:n_eqns,j,k,i_RK) = ( q_si(1:n_vars) -                   &
                     q_fv_cell ) / ( dt*a_diag )

                ! Initialize the guess for the NR solver
                q_guess(1:n_vars) = q_si(1:n_vars)


                Rj_not_impl =  ( MATMUL( divFlux(1:n_eqns,j,k,1:i_RK-1) -       &
                     expl_terms(1:n_eqns,j,k,1:i_RK-1), a_tilde(1:i_RK-1) )     &
                     - MATMUL( NH(1:n_eqns,j,k,1:i_RK-1)                        &
                     + SI_NH(1:n_eqns,j,k,1:i_RK-1) , a_dirk(1:i_RK-1) ) )      &
                     - a_diag * SI_NH(1:n_eqns,j,k,i_RK)

                ! Solve the implicit system to find the solution at the 
                ! i_RK step of the IMEX RK procedure
                CALL solve_rk_step( q_guess(1:n_vars) , q(1:n_vars,j,k ) ,      &
                     dt, a_diag , Rj_not_impl , B_prime_x_geom(j,k) ,            &
                     B_prime_y_geom(j,k), Z(j,k), fric_array(j,k),              &
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
                   
                   NH(1:n_eqns,j,k,i_RK) = ( q_guess(1:n_vars)                  &
                        - q_si(1:n_vars) ) / ( dt*a_diag )
                   
                ELSE
                   
                   ! Eval and store the implicit term at the i_RK step
                   CALL eval_implicit_terms( B_prime_x_geom(j,k) ,              &
                        B_prime_y_geom(j,k) , Z(j,k), fric_array(j,k),          &
                        r_qj = q_guess , r_nh_term_impl = NH(1:n_eqns,j,k,i_RK) )
                   
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
                SI_NH(1:n_eqns,j,k,i_RK) = 0.0_wp
                NH(1:n_eqns,j,k,i_RK) = 0.0_wp

             END IF pos_thick

          END IF adiag_pos

          IF ( a_diag .NE. 0.0_wp ) THEN

             ! Update the implicit term with correction on the new velocity
             NH(1:n_vars,j,k,i_RK) = ( q_guess(1:n_vars) - q_si(1:n_vars))      &
                  / ( dt*a_diag ) 

          END IF

          ! Store the current stage. Previous stage states are no longer
          ! needed here: their evaluated terms are retained in divFlux, NH,
          ! SI_NH and expl_terms.
          q_rk( 1:n_vars , j , k ) = q_guess

          IF ( verbose_level .GE. 2 ) THEN

             WRITE(*,*) 'imex_RK_solver: qc',q_guess

             IF ( q_guess(1) .GT. 0.0_wp ) THEN

                CALL qc_to_qp( q_guess , qp(1:n_vars+2,j,k) , p_dyn )
                WRITE(*,*) 'imex_RK_solver: qp',qp(1:n_vars+2,j,k)

             END IF
             
             READ(*,*)

          END IF


          IF ( need_explicit_stage ) THEN
          
             IF ( q_rk(1,j,k) .GT. 0.0_wp ) THEN

                CALL qc_to_qp( q_rk(1:n_vars,j,k) ,                             &
                     qp_rk(1:n_vars+2,j,k) , p_dyn )

             ELSE

                qp_rk(1:n_vars+2,j,k) = 0.0_wp
                qp_rk(4,j,k) = T_ambient

             END IF

             ! Eval gravity term and radial bottom source terms
             CALL eval_expl_terms( B_prime_x_geom(j,k) , B_prime_y_geom(j,k) ,  &
                  B_second_xx_geom(j,k) , B_second_xy_geom(j,k) ,               &
                  B_second_yy_geom(j,k) , grav_coeff(j,k), d_grav_coeff_dx(j,k),&
                  d_grav_coeff_dy(j,k) , source_xy(j,k),                        &
                  qp_rk(1:n_vars+2,j,k), expl_terms(1:n_eqns,j,k,i_RK), t,      &
                  cell_source_fractions(j,k),                                   &
                  cell_arc_perim(j,k), cell_arc_n_x(j,k), cell_arc_n_y(j,k),    &
                  dx * dy )
  
          END IF

       END DO solve_cells_loop

       !$OMP END DO
       !$OMP END PARALLEL 

       IF ( need_explicit_stage ) THEN

          ! Eval and store the explicit hyperbolic (fluxes) terms
          CALL eval_hyperbolic_terms(                                           &
               q_rk , qp_rk ,                                                   &
               divFlux(1:n_eqns,1:comp_cells_x,1:comp_cells_y,i_RK) )

       END IF

    END DO runge_kutta

    !$OMP PARALLEL DO private(j,k,p_dyn,alpha_s,solid_excess_roundoff,          &
    !$OMP & residual_cell,q_old_cell)

    assemble_sol:DO l = 1,solve_cells

       j = j_cent(l)
       k = k_cent(l)

       ! q remains equal to Q^n throughout all RK stages. Preserve the old
       ! state locally before overwriting this cell during final assembly.
       q_old_cell = q(1:n_vars,j,k)

       residual_cell = MATMUL( divFlux(1:n_eqns,j,k,1:n_RK)                     &
            - expl_terms(1:n_eqns,j,k,1:n_RK) , omega_tilde ) -                 &
            MATMUL( NH(1:n_eqns,j,k,1:n_RK) + SI_NH(1:n_eqns,j,k,1:n_RK) ,      &
            omega )


       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'cell jk =',j,k
          WRITE(*,*) 'before imex_RK_solver: qc',q_old_cell

          IF ( q_old_cell(1) .GT. 0.0_wp ) THEN

             CALL qc_to_qp(q_old_cell , qp(1:n_vars+2,j,k) , p_dyn )
             WRITE(*,*) 'before imex_RK_solver: qp',qp(1:n_vars+2,j,k)
 
          END IF

       END IF

       IF ( ( SUM(ABS( omega_tilde(:)-a_tilde_ij(n_RK,:))) .EQ. 0.0_wp  )       &
            .AND. ( SUM(ABS(omega(:)-a_dirk_ij(n_RK,:))) .EQ. 0.0_wp ) ) THEN

          ! The assembling coeffs are equal to the last step of the RK scheme
          q(1:n_vars,j,k) = q_rk(1:n_vars,j,k)

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

             WRITE(*,*) 'divFlux(1,j,k,1:n_RK)',divFlux(1,j,k,1:n_RK) 

             WRITE(*,*) H_interface_x(1,j+1,k), H_interface_x(1,j,k)
             WRITE(*,*) qp_interfaceR(1:n_vars,j,k)
             WRITE(*,*) qp(1:n_vars,j,k)
             WRITE(*,*) qp_interfaceL(1:n_vars,j+1,k)

             WRITE(*,*) 'expl_terms(1,j,k,1:n_RK)',expl_terms(1,j,k,1:n_RK) 
             WRITE(*,*) 'NH(1,j,k,1:n_RK)',NH(1,j,k,1:n_RK) 
             WRITE(*,*) 'SI_NH(1,j,k,1:n_RK)',SI_NH(1,j,k,1:n_RK) 

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
             WRITE(*,*) H_interface_x(1,j+1,k)/dx*dt, H_interface_x(1,j,k)/dx*dt
             WRITE(*,*) H_interface_y(1,j,k+1)/dy*dt, H_interface_y(1,j,k)/dy*dt
             
             WRITE(*,*) 'H_interface(5)'
             WRITE(*,*) H_interface_x(5,j+1,k)/dx*dt, H_interface_x(5,j,k)/dx*dt
             WRITE(*,*) H_interface_y(5,j,k+1)/dy*dt, H_interface_y(5,j,k)/dy*dt
             

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
          WRITE(*,*) H_interface_x(4,j+1,k)/dx*dt, H_interface_x(4,j,k)/dx*dt
          WRITE(*,*) H_interface_y(4,j,k+1)/dy*dt, H_interface_y(4,j,k)/dy*dt

          WRITE(*,*) H_interface_y(:,j,k)/dy*dt
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
             WRITE(*,*) H_interface_x(1,j+1,k)/dx*dt, H_interface_x(1,j,k)/dx*dt
             WRITE(*,*) H_interface_y(1,j,k+1)/dy*dt, H_interface_y(1,j,k)/dy*dt
             
             WRITE(*,*) 'H_interface(5)'
             WRITE(*,*) H_interface_x(5,j+1,k)/dx*dt, H_interface_x(5,j,k)/dx*dt
             WRITE(*,*) H_interface_y(5,j,k+1)/dy*dt, H_interface_y(5,j,k)/dy*dt

             WRITE(*,*) 'divFlux(1)',divFlux(1,j,k,1:n_RK)
             WRITE(*,*) 'expl_terms(1)', expl_terms(1,j,k,1:n_RK)
             WRITE(*,*) 'NH(1)', NH(1,j,k,1:n_RK)
             WRITE(*,*) 'SI(1)', SI_NH(1,j,k,1:n_RK) 

             WRITE(*,*) 'divFlux(5)',divFlux(5,j,k,1:n_RK)
             WRITE(*,*) 'expl_terms(5)', expl_terms(5,j,k,1:n_RK)
             WRITE(*,*) 'NH(5)', NH(5,j,k,1:n_RK)
             WRITE(*,*) 'SI(5)', SI_NH(5,j,k,1:n_RK) 
             

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

  !******************************************************************************
  !> \brief Runge-Kutta single step integration
  !
  !> This subroutine find the solution of the non-linear system 
  !> given the a step of the implicit-explicit Runge-Kutta scheme for a
  !> cell:\n
  !> \f$ Q^{(i)} = Q^n - dt \sum_{j=1}^{i-1}\tilde{a}_{j}\partial_x 
  !> F(Q^{(j)}) +  dt \sum_{j=1}^{i-1} a_j  NH(Q^{(j)}) 
  !> + dt a_{diag} NH(Q^{(i)}) \f$\n
  !
  !> \param[in,out] qj        conservative variables 
  !> \param[in]     qj_old    conservative variables at the old time step
  !> \param[in]     a_diag    implicit coefficient for the non-hyperbolic terms 
  !> \param[in]     Rj_not_impl
  !> \param[out]    line_search_failed  true when no acceptable step is found
  !
  !> \date 2019/12/16
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************
  !> \brief Evaluate the eroion/deposition terms
  !
  !> This subroutine update the solution and the topography computing the 
  !> erosion and deposition terms and the solution only because of entrainment.
  !
  !> \param[in]    dt      time step
  !
  !> \date 2019/11/08
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE update_erosion_deposition_cell(dt)

    USE constitutive_2d, ONLY : erosion_coeff , settling_flag
    USE constitutive_2d, ONLY : maximum_solid_packing
    
    USE geometry_2d, ONLY : deposit , erosion , erodible
    USE geometry_2d, ONLY : B_zone

    USE constitutive_2d, ONLY : eval_mass_exchange_terms

    USE constitutive_2d, ONLY : qc_to_qp , mixt_var
    USE constitutive_2d, ONLY : entrainment_flag
    USE parameters_2d, ONLY : topo_change_flag , bottom_radial_source_flag
    USE parameters_2d, ONLY : erodible_deposit_flag
    USE parameters_2d, ONLY : pore_pressure_flag

    IMPLICIT NONE
    
    REAL(wp), INTENT(IN) :: dt

    REAL(wp) :: erosion_term(n_solid)
    REAL(wp) :: deposition_term(n_solid)
    REAL(wp) :: continuous_phase_erosion_term
    REAL(wp) :: continuous_phase_loss_term
    REAL(wp) :: eqns_term(n_eqns)
    REAL(wp) :: topo_term

    REAL(wp) :: r_Ri , r_rho_m
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity

    INTEGER :: j,k,l

    REAL(wp) :: out_of_source_fraction

    REAL(wp) :: p_dyn

    LOGICAL :: sp_flag
    REAL(wp) :: r_sp_heat_c
    REAL(wp) :: r_sp_heat_mix

    sp_flag = .FALSE.


    IF ( ( erosion_coeff .EQ. 0.0_wp ) .AND. ( .NOT.settling_flag ) &
         .AND. ( .NOT.pore_pressure_flag ) .AND. ( .NOT.entrainment_flag) ) RETURN

    !$OMP PARALLEL DO private(j,k,erosion_term,deposition_term,eqns_term,       &
    !$OMP & topo_term,r_Ri,r_rho_m,r_rho_c,r_red_grav,                          &
    !$OMP & continuous_phase_erosion_term,continuous_phase_loss_term,           &
    !$OMP & out_of_source_fraction,p_dyn,r_sp_heat_c,r_sp_heat_mix)

    DO l = 1,solve_cells

       j = j_cent(l)
       k = k_cent(l)

       IF ( q(1,j,k) .GT. 0.0_wp ) THEN

          CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars+2,j,k) , p_dyn )

       ELSE

          qp(1:n_vars+2,j,k) = 0.0_wp
          qp(4,j,k) = T_ambient

       END IF

       CALL eval_mass_exchange_terms( qp(1:n_vars+2,j,k) , B_zone(j,k) ,           &
            B_prime_x_geom(j,k) , B_prime_y_geom(j,k) , erodible(1:n_solid,j,k) ,  &
            dt , erosion_term , deposition_term , continuous_phase_erosion_term ,  &
            continuous_phase_loss_term , eqns_term , topo_term  )
          
       IF ( bottom_radial_source_flag ) THEN

          ! entrainment, erosion and deposition occurs only outside source
          out_of_source_fraction = 1.0_wp - cell_source_fractions(j,k)
          deposition_term = deposition_term * out_of_source_fraction
          erosion_term = erosion_term * out_of_source_fraction
          eqns_term = eqns_term * out_of_source_fraction
          topo_term = topo_term * out_of_source_fraction

       END IF
       
       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*) 'before update erosion/deposition: j,k,q(:,j,k),B(j,k)',   &
               j,k,q(:,j,k),B_cent(j,k)

       END IF

       ! Update the solution with erosion/deposition terms
       q(1:n_eqns,j,k) = q(1:n_eqns,j,k) + dt * eqns_term(1:n_eqns)
       q(5:4+n_solid,j,k) = MAX( 0.0_wp , q(5:4+n_solid,j,k) )
       
       deposit(j,k,1:n_solid) = deposit(j,k,1:n_solid)                          &
            + dt * deposition_term(1:n_solid)

       erosion(j,k,1:n_solid) = erosion(j,k,1:n_solid)                          &
            + dt * erosion_term(1:n_solid)

       erodible(1:n_solid,j,k) = erodible(1:n_solid,j,k)                        &
            - dt * erosion_term(1:n_solid)

       IF ( erodible_deposit_flag ) THEN

          erodible(1:n_solid,j,k) = erodible(1:n_solid,j,k)                     &
               + dt * deposition_term(1:n_solid)

       END IF
       
       ! Update the topography with erosion/deposition terms
       IF ( topo_change_flag ) THEN

          B_cent(j,k) = B_cent(j,k) + dt * topo_term

       END IF

       negative_alpha_check:IF ( ANY(q(5:4+n_solid,j,k) .LT. 0.0_wp ) ) THEN

          WRITE(*,*) 'WARNINIG: negative solid mass'
          WRITE(*,*) 'j,k',j,k
          WRITE(*,*) 'dt',dt
          WRITE(*,*) 'before erosion: qc',q(1:n_vars,j,k) - dt * eqns_term(1:n_eqns)
          WRITE(*,*) 'deposition_term',deposition_term
          WRITE(*,*) 'erosion_term',erosion_term
          WRITE(*,*) 'after erosion: qc',q(1:n_vars,j,k)

          READ(*,*)
          
       END IF negative_alpha_check
       
       ! Check for negative thickness
       IF ( q(1,j,k) .LE. 0.0_wp ) THEN

          IF ( q(1,j,k) .GT. -1.0E-10_wp ) THEN

             q(1:n_vars,j,k) = 0.0_wp

          ELSE

             WRITE(*,*) 'j,k',j,k
             WRITE(*,*) 'dt',dt
             WRITE(*,*) 'before erosion'
             WRITE(*,*) 'qp',qp(1:n_eqns+2,j,k)
             WRITE(*,*) 'q',q(1:n_eqns,j,k) - dt * eqns_term(1:n_eqns)
             WRITE(*,*) 'deposition_term',deposition_term
             WRITE(*,*) 'erosion_term',erosion_term
             WRITE(*,*) 'continuous_phase_loss_term',continuous_phase_loss_term
             WRITE(*,*) 'eqns_term',eqns_term
             WRITE(*,*) 'after erosion'
             CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars+2,j,k) , p_dyn )
             WRITE(*,*) 'q',q(1:n_eqns,j,k)
             WRITE(*,*) 'qp',qp(1:n_eqns+2,j,k)
                
             READ(*,*)

          END IF

       END IF

       IF ( SUM(q(5:4+n_solid,j,k)) .GT. q(1,j,k) ) THEN

          IF ( q(1,j,k) .LT. 1.0e-10_wp ) THEN

             q(5:4+n_solid,j,k) = q(5:4+n_solid,j,k)                            &
                  / SUM(q(5:4+n_solid,j,k)) * q(1,j,k)
             
          ELSE

             WRITE(*,*) 'SUM SOLID > TOT'
             WRITE(*,*) 'j,k',j,k
             WRITE(*,*) 'dt',dt
             WRITE(*,*) 'before erosion'
             WRITE(*,*) 'qp',qp(1:n_eqns+2,j,k)
             WRITE(*,*) 'q',q(1:n_eqns,j,k) - dt * eqns_term(1:n_eqns)
             WRITE(*,*) 'deposition_term',deposition_term
             WRITE(*,*) 'erosion_term',erosion_term
             WRITE(*,*) 'continuous_phase_loss_term',continuous_phase_loss_term
             WRITE(*,*) 'after erosion'
             CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars+2,j,k) , p_dyn )
             WRITE(*,*) 'qp',qp(1:n_eqns+2,j,k)
             WRITE(*,*) 'q',q(1:n_eqns,j,k)          
             READ(*,*)
             
          END IF

       END IF


       IF ( q(1,j,k) .GT. 0.0_wp ) THEN

          CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars+2,j,k) , p_dyn )
          CALL mixt_var(qp(1:n_vars+2,j,k),r_Ri,r_rho_m,r_rho_c,r_red_grav,     &
               sp_flag,r_sp_heat_c,r_sp_heat_mix)

       ELSE

          qp(1:n_vars+2,j,k) = 0.0_wp
          qp(4,j,k) = T_ambient
          r_red_grav = 0.0_wp

       END IF

       IF ( r_red_grav .LE. 0.0_wp ) THEN

          q(1:n_vars,j,k) = 0.0_wp

       END IF

    END DO

    !$OMP END PARALLEL DO

    RETURN

  END SUBROUTINE update_erosion_deposition_cell

  !******************************************************************************
  !> \brief Semidiscrete finite volume central scheme
  !
  !> This subroutine compute the divergence part of the system of the eqns,
  !> with a modified version of the finite volume scheme from Kurganov et al.  
  !> 2001, where the reconstruction at the cells interfaces is applied to a
  !> set of physical variables derived from the conservative vriables.
  !
  !> \param[in]     q_expl         conservative variables
  !> \param[in]     qp_expl        conservative variables
  !> \param[out]    divFlux_iRK    divergence term
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_hyperbolic_terms( q_expl , qp_expl , divFlux_iRK )

    ! External variables
    USE parameters_2d, ONLY : solver_scheme

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: q_expl(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: qp_expl(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(OUT) :: divFlux_iRK(n_eqns,comp_cells_x,comp_cells_y)

    INTEGER :: l , i, j, k      !< loop counters

    !WRITE(*,*) 'SUBROUTINE eval_hyperbolic_terms'
    !WRITE(*,*) 'qp_expl(4,1,1)',qp_expl(4,1,1)
    !WRITE(*,*)
    
    ! Linear reconstruction of the physical variables at the interfaces
    CALL reconstruction( q_expl, qp_expl, t, solve_cells, j_cent, k_cent )

    ! Evaluation of the maximum local speeds at the interfaces
    CALL eval_speeds

    ! Evaluation of the numerical fluxes
    SELECT CASE ( solver_scheme )

    CASE ("LxF")

       CALL eval_flux_LxF

    CASE ("GFORCE")

       CALL eval_flux_GFORCE

    CASE ("KT")

       CALL eval_flux_KT

    CASE ("UP")

       CALL eval_flux_UP

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
  
  SUBROUTINE eval_flux_UP

    ! External procedures
    USE constitutive_2d, ONLY : eval_fluxes
    USE geometry_2d, ONLY : grav_coeff_stag_x , grav_coeff_stag_y

    IMPLICIT NONE

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

  SUBROUTINE eval_flux_KT

    ! External procedures
    USE constitutive_2d, ONLY : eval_fluxes
    USE geometry_2d, ONLY : grav_coeff_stag_x , grav_coeff_stag_y

    IMPLICIT NONE

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

  SUBROUTINE eval_speeds

    ! External procedures
    USE constitutive_2d, ONLY : eval_local_speeds_x, eval_local_speeds_y 

    IMPLICIT NONE

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

END MODULE solver_2d
