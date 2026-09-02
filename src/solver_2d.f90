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
  USE geometry_2d, ONLY : B_cent
  USE geometry_2d, ONLY : B_prime_x_geom , B_prime_y_geom  
  USE geometry_2d, ONLY : cell_source_fractions

  USE parameters_2d, ONLY : wp

  USE parameters_2d, ONLY : n_eqns , n_vars , n_solid
  USE parameters_2d, ONLY : verbose_level
  USE parameters_2d, ONLY : bottom_radial_source_flag

  USE OMP_LIB

  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver,                 &
       finalize_nonlinear_solver

  USE reconstruction_2d, ONLY : initialize_reconstruction,                    &
       finalize_reconstruction

  USE hyperbolic_2d, ONLY : initialize_hyperbolic, finalize_hyperbolic

  USE domain_2d, ONLY : initialize_domain, finalize_domain
  USE domain_2d, ONLY : solve_cells, j_cent, k_cent

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

END MODULE solver_2d
