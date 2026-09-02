!********************************************************************************
!> \brief Erosion, deposition and entrainment updates
!
!> This module applies mass-exchange terms to the conservative state and
!> updates the associated deposits, erodible material and topography.
!
!********************************************************************************
MODULE mass_exchange_2d

  USE constitutive_2d, ONLY : T_ambient

  USE geometry_2d, ONLY : B_cent
  USE geometry_2d, ONLY : B_prime_x_geom, B_prime_y_geom
  USE geometry_2d, ONLY : cell_source_fractions
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y

  USE parameters_2d, ONLY : wp
  USE parameters_2d, ONLY : n_eqns, n_vars, n_solid
  USE parameters_2d, ONLY : verbose_level
  USE parameters_2d, ONLY : bottom_radial_source_flag

  USE domain_2d, ONLY : solve_cells, j_cent, k_cent

  USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: update_erosion_deposition_cell

CONTAINS

  !******************************************************************************
  !> \brief Apply erosion, deposition and entrainment over the active cells.
  !>
  !> \param[in,out] q   conservative variables
  !> \param[in,out] qp  physical variables
  !> \param[in]     dt  time step
  !******************************************************************************

  SUBROUTINE update_erosion_deposition_cell(q, qp, dt)

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

    REAL(wp), INTENT(INOUT) :: q(n_vars,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(INOUT) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
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

END MODULE mass_exchange_2d
