!> \brief Physical treatment interfaces for bottom_circular_injection.
MODULE mass_bottom_circular_injection_2d
    USE parameters_2d, ONLY : idx_h, idx_hu, idx_hv, idx_T
  USE parameters_2d
  USE constitutive_parameters_2d
  USE state_conversion_2d, ONLY : mixt_var
  USE thermal_laws_2d, ONLY : eval_inlet_thermal_energy_source
  USE pore_pressure_laws_2d, ONLY : eval_pore_pressure_inlet_source
  USE geometry_2d, ONLY : pi_g, bottom_source_area_ratio
  USE mass_injection_coefficients_2d
  USE diagnostics_2d, ONLY : fatal_error
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_bottom_circular_injection_explicit
  PUBLIC :: eval_bottom_circular_injection_semi_implicit
  PUBLIC :: eval_bottom_circular_injection_implicit
  PUBLIC :: solve_bottom_circular_injection_analytic
CONTAINS
  SUBROUTINE eval_bottom_circular_injection_explicit(qpj, time, cell_fract_jk, q1, expl_term)
    USE parameters_2d, ONLY : vel_source, T_source, xs_source, xg_source, &
         time_param, bottom_radial_source_flag, pore_pressure_flag, &
         pore_pres_fract, n_intervals, t_intervals, vel_intervals, &
         radial_source_flag, h_source, bottom_fissural_source_flag, &
         n_fissures, linear_vel_fissures, T_fissures, time_param_fissures, &
         xs_fissures, xg_fissures
    USE geometry_2d, ONLY : pi_g, bottom_source_area_ratio
    REAL(wp), INTENT(IN) :: qpj(n_vars+2), time, cell_fract_jk, q1
    REAL(wp), INTENT(INOUT) :: expl_term(n_eqns)
    REAL(wp) :: t_rem, t_coeff, h_dot, vel_local
    REAL(wp) :: qp_source(n_vars+2)
    REAL(wp) :: r_Ri, r_rho_m, r_rho_c, r_red_grav
    REAL(wp) :: r_sp_heat_c, r_sp_heat_mix, exc_pore_pres
    INTEGER :: i_int, i_fissure

    ! ----------- ADDITIONAL EXPLICIT TERMS FOR BOTTOM RADIAL SOURCE ------------

    IF ( bottom_radial_source_flag .AND. ( cell_fract_jk .GT. 0.0_wp ) ) THEN

    IF ( n_intervals .GT. 0 ) THEN

       ! --- Variable velocity intervals mode ---
       ! Find which interval the current time falls in (reverse search)
       vel_local = 0.0_wp

       DO i_int = n_intervals, 1, -1

          IF ( time .GE. t_intervals(i_int) ) THEN

             vel_local = vel_intervals(i_int)
             EXIT

          END IF

       END DO

       IF ( bottom_vel_step .GE. 0.0_wp ) vel_local = bottom_vel_step

       ! vel_local = 0 means source is off
       IF ( vel_local .GT. 0.0_wp ) THEN

          t_coeff = 1.0_wp

       ELSE

          t_coeff = 0.0_wp

       END IF

       h_dot = cell_fract_jk * vel_local * bottom_source_area_ratio

    ELSE

       ! --- Original pulse mode (backward compatible) ---
       t_rem = MOD( time + time_param(4) , time_param(1) )

       IF ( time_param(3) .EQ. 0.0_wp ) THEN

          IF ( t_rem .LE. time_param(2) ) THEN

             t_coeff = 1.0_wp

          ELSE

             t_coeff = 0.0_wp

          END IF

       ELSE

          IF ( t_rem .LE. time_param(3) ) THEN

             t_coeff = ( t_rem / time_param(3) )

          ELSEIF ( t_rem .LE. time_param(2) - time_param(3) ) THEN

             t_coeff = 1.0_wp

          ELSEIF ( t_rem .LE. time_param(2) ) THEN

             t_coeff = 1.0_wp - ( t_rem - time_param(2) + time_param(3) ) /     &
                  time_param(3)

          ELSE

             t_coeff = 0.0_wp

          END IF

       END IF

       IF ( bottom_coeff_step .GE. 0.0_wp ) t_coeff = bottom_coeff_step

       h_dot = cell_fract_jk * vel_source * bottom_source_area_ratio

    END IF

    qp_source = 0.0_wp
    qp_source(idx_h) = 1.0_wp
    qp_source(idx_hu) = 0.0_wp
    qp_source(idx_hv) = 0.0_wp
    qp_source(idx_T) = t_source

    qp_source(idx_solid_first:idx_solid_last) = xs_source(1:n_solid)
    qp_source(idx_add_gas_first:idx_add_gas_last) = xg_source(1:n_add_gas)

    ! Source term transport stochastic equation
    IF ( stoch_transport_flag) qp_source(idx_stoch) = 0.0_wp

    IF ( pore_pressure_flag ) qp_source(idx_poreEqn) = 0.0_wp

    qp_source(idx_u) = 0.0_wp
    qp_source(idx_v) = 0.0_wp


    CALL mixt_var(qp_source,r_Ri,r_rho_m,r_rho_c,r_red_grav,                   &
         r_sp_heat_c,r_sp_heat_mix)

    expl_term(1) = expl_term(1) + t_coeff * h_dot * r_rho_m
    expl_term(2) = expl_term(2) + 0.0_wp
    expl_term(3) = expl_term(3) + 0.0_wp

    expl_term(4) = expl_term(4) + eval_inlet_thermal_energy_source(t_coeff,    &
         h_dot, r_rho_m, r_sp_heat_mix, t_source)

    ! source terms for the solid equations
    expl_term(idx_solidEqn_first:idx_solidEqn_last) =                          &
         expl_term(idx_solidEqn_first:idx_solidEqn_last) + t_coeff             &
         * h_dot * r_rho_m * xs_source(1:n_solid)

    ! source terms for the additional gas equations
    expl_term(idx_addGasEqn_first:idx_addGasEqn_last) =                         &
         expl_term(idx_addGasEqn_first:idx_addGasEqn_last) + t_coeff            &
         * h_dot * r_rho_m * xg_source(1:n_add_gas)


    IF ( pore_pressure_flag ) THEN

       exc_pore_pres = qpj(idx_pore)

       ! we multiply the pore pressure inlet rate by the rate for q1
       ! the units of this source term are: kg^2 m^-3 s^-3
       expl_term(idx_poreEqn) = expl_term(idx_poreEqn) +                       &
            eval_pore_pressure_inlet_source(t_coeff, q1, pore_pres_fract,      &
            h_dot, r_rho_m, r_red_grav, exc_pore_pres)

    END IF

    END IF   ! bottom_radial_source_flag
  END SUBROUTINE eval_bottom_circular_injection_explicit

  SUBROUTINE eval_bottom_circular_injection_semi_implicit(qpj, time, cell_fract_jk, q1, expl_term)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2), time, cell_fract_jk, q1
    REAL(wp), INTENT(INOUT) :: expl_term(n_eqns)
    CALL fatal_error('bottom_circular_injection: semi_implicit treatment is unsupported')
  END SUBROUTINE eval_bottom_circular_injection_semi_implicit

  SUBROUTINE eval_bottom_circular_injection_implicit(qpj, time, cell_fract_jk, q1, expl_term)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2), time, cell_fract_jk, q1
    REAL(wp), INTENT(INOUT) :: expl_term(n_eqns)
    CALL fatal_error('bottom_circular_injection: implicit treatment is unsupported')
  END SUBROUTINE eval_bottom_circular_injection_implicit

  SUBROUTINE solve_bottom_circular_injection_analytic(qpj, time, cell_fract_jk, q1, expl_term)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2), time, cell_fract_jk, q1
    REAL(wp), INTENT(INOUT) :: expl_term(n_eqns)
    CALL fatal_error('bottom_circular_injection: analytic treatment is unsupported')
  END SUBROUTINE solve_bottom_circular_injection_analytic
END MODULE mass_bottom_circular_injection_2d
