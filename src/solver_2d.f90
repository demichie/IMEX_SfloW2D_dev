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

  USE constitutive_2d, ONLY : implicit_flag, implicit_map
  USE constitutive_2d, ONLY : T_ambient
    
  USE geometry_2d, ONLY : comp_cells_x,comp_cells_y
  USE parameters_2d, ONLY : wp

  USE parameters_2d, ONLY : n_vars
  USE parameters_2d, ONLY : verbose_level

  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver,                 &
       finalize_nonlinear_solver

  USE reconstruction_2d, ONLY : initialize_reconstruction,                    &
       finalize_reconstruction

  USE hyperbolic_2d, ONLY : initialize_hyperbolic, finalize_hyperbolic

  USE domain_2d, ONLY : initialize_domain, finalize_domain
  USE time_integration_2d, ONLY : initialize_time_integration,                &
       finalize_time_integration

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

  !> Physical variables (\f$\alpha_1, p_1, p_2, \rho u, w, T\f$)
  REAL(wp), ALLOCATABLE :: qp(:,:,:)

  !> Array defining fraction of cells affected by source term
  REAL(wp), ALLOCATABLE :: source_xy(:,:)

  !> Time step
  REAL(wp) :: dt

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
    CALL initialize_hyperbolic
    CALL initialize_domain

    ALLOCATE( source_xy( comp_cells_x , comp_cells_y ) )


    CALL initialize_nonlinear_solver
    CALL initialize_time_integration

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

    CALL finalize_domain
    CALL finalize_time_integration
    CALL finalize_hyperbolic
    CALL finalize_reconstruction

    DEALLOCATE( qp )

    DEALLOCATE( source_xy )

    CALL finalize_nonlinear_solver

    IF ( ALLOCATED(implicit_flag) ) DEALLOCATE(implicit_flag)
    IF ( ALLOCATED(implicit_map) ) DEALLOCATE(implicit_map)
    IF ( ALLOCATED(Z) ) DEALLOCATE(Z)
    IF ( ALLOCATED(conv_kernel) ) DEALLOCATE(conv_kernel)
    IF ( ALLOCATED(fric_array) ) DEALLOCATE(fric_array)

    
    RETURN
    
  END SUBROUTINE deallocate_solver_variables


END MODULE solver_2d
