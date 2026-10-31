!> \brief Physical treatment interfaces for bottom_fissural_injection.
MODULE mass_bottom_fissural_injection_2d
    USE parameters_2d, ONLY : idx_h, idx_T
  USE parameters_2d
  USE constitutive_parameters_2d
  USE state_conversion_2d, ONLY : mixt_var, physical_state_type, complex_physical_state_type

  USE geometry_2d, ONLY : pi_g, bottom_source_area_ratio
  USE mass_injection_coefficients_2d
  USE diagnostics_2d, ONLY : fatal_error
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_bottom_fissural_injection_explicit
  PUBLIC :: eval_bottom_fissural_injection_semi_implicit
  PUBLIC :: eval_bottom_fissural_injection_implicit
  PUBLIC :: eval_bottom_fissural_injection_analytic
CONTAINS
  SUBROUTINE eval_bottom_fissural_injection_explicit(state, time, cell_fissure_fract, rate)
    USE parameters_2d, ONLY : vel_source, T_source, xs_source, xg_source, &
         time_param, bottom_radial_source_flag, pore_pressure_flag, &
         pore_pres_fract, n_intervals, t_intervals, vel_intervals, &
         radial_source_flag, h_source, bottom_fissural_source_flag, &
         n_fissures, linear_vel_fissures, T_fissures, time_param_fissures, &
         xs_fissures, xg_fissures
    USE geometry_2d, ONLY : pi_g, bottom_source_area_ratio
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: time, cell_fissure_fract(:)
    REAL(wp), INTENT(OUT) :: rate(n_eqns)
    REAL(wp) :: t_rem, t_coeff, h_dot
    REAL(wp) :: qp_source(n_vars+2)
    REAL(wp) :: r_Ri, r_rho_m, r_rho_c, r_red_grav
    REAL(wp) :: r_sp_heat_c, r_sp_heat_mix, exc_pore_pres
    INTEGER :: i_fissure

    rate=0.0_wp

    ! Each fissure uses its own clipped cell coverage, pulse and temperature,
    ! and the source composition converted at that temperature.
    IF ( bottom_fissural_source_flag ) THEN

       DO i_fissure = 1, n_fissures

          IF ( cell_fissure_fract(i_fissure) .LE. 0.0_wp ) CYCLE

          t_rem = MOD(time + time_param_fissures(4,i_fissure),                &
               time_param_fissures(1,i_fissure))
          IF (time_param_fissures(3,i_fissure) .EQ. 0.0_wp) THEN
             IF (t_rem .LE. time_param_fissures(2,i_fissure)) THEN
                t_coeff = 1.0_wp
             ELSE
                t_coeff = 0.0_wp
             END IF
          ELSEIF (t_rem .LT. time_param_fissures(3,i_fissure)) THEN
             t_coeff = t_rem / time_param_fissures(3,i_fissure)
          ELSEIF (t_rem .LE. time_param_fissures(2,i_fissure) -              &
               time_param_fissures(3,i_fissure)) THEN
             t_coeff = 1.0_wp
          ELSEIF (t_rem .LE. time_param_fissures(2,i_fissure)) THEN
             t_coeff = 1.0_wp - (t_rem - time_param_fissures(2,i_fissure) +  &
                  time_param_fissures(3,i_fissure)) /                        &
                  time_param_fissures(3,i_fissure)
          ELSE
             t_coeff = 0.0_wp
          END IF

          IF ( ALLOCATED( fissure_coeff_step ) )                              &
               t_coeff = fissure_coeff_step(i_fissure)

          IF (t_coeff .LE. 0.0_wp) CYCLE

          h_dot = cell_fissure_fract(i_fissure) * linear_vel_fissures(i_fissure)
          qp_source = 0.0_wp
          qp_source(idx_h) = 1.0_wp
          qp_source(idx_T) = T_fissures(i_fissure)
          qp_source(idx_solid_first:idx_solid_last) =                          &
               xs_fissures(1:n_solid,i_fissure)
          qp_source(idx_add_gas_first:idx_add_gas_last) =                      &
               xg_fissures(1:n_add_gas,i_fissure)
          IF (stoch_transport_flag) qp_source(idx_stoch) = 0.0_wp
          IF (pore_pressure_flag) qp_source(idx_poreEqn) = 0.0_wp
          qp_source(idx_u) = 0.0_wp
          qp_source(idx_v) = 0.0_wp

          CALL mixt_var(qp_source, r_Ri, r_rho_m, r_rho_c, r_red_grav,       &
               r_sp_heat_c, r_sp_heat_mix)

          rate(1) = rate(1) + t_coeff * h_dot * r_rho_m
          rate(4) = rate(4) + (t_coeff * &
               h_dot * r_rho_m * r_sp_heat_mix * T_fissures(i_fissure))
          rate(idx_solidEqn_first:idx_solidEqn_last) =                   &
               rate(idx_solidEqn_first:idx_solidEqn_last) +              &
               t_coeff * h_dot * r_rho_m * xs_fissures(1:n_solid,i_fissure)
          rate(idx_addGasEqn_first:idx_addGasEqn_last) =                  &
               rate(idx_addGasEqn_first:idx_addGasEqn_last) +             &
               t_coeff * h_dot * r_rho_m * xg_fissures(1:n_add_gas,i_fissure)


          IF (pore_pressure_flag) THEN
             exc_pore_pres = state%exc_pore_pres
             rate(idx_poreEqn) = rate(idx_poreEqn) +                &
                  t_coeff * ((state%h*state%rho_m) * pore_pres_fract * h_dot * &
                  r_rho_m * r_red_grav + exc_pore_pres * h_dot * r_rho_m)
          END IF

       END DO

    END IF
  END SUBROUTINE eval_bottom_fissural_injection_explicit

  SUBROUTINE eval_bottom_fissural_injection_semi_implicit(state, time, cell_fissure_fract, rate)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: time, cell_fissure_fract(:)
    REAL(wp), INTENT(OUT) :: rate(n_eqns)
    CALL fatal_error('bottom_fissural_injection: semi_implicit treatment is unsupported')
  END SUBROUTINE eval_bottom_fissural_injection_semi_implicit

  SUBROUTINE eval_bottom_fissural_injection_implicit(state, time, cell_fissure_fract, rate)
    TYPE(complex_physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: time, cell_fissure_fract(:)
    COMPLEX(wp), INTENT(OUT) :: rate(n_eqns)
    CALL fatal_error('bottom_fissural_injection: implicit treatment is unsupported')
  END SUBROUTINE eval_bottom_fissural_injection_implicit

  SUBROUTINE eval_bottom_fissural_injection_analytic(state, time, cell_fissure_fract, rate)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: time, cell_fissure_fract(:)
    REAL(wp), INTENT(OUT) :: rate(n_eqns)
    CALL fatal_error('bottom_fissural_injection: analytic treatment is unsupported')
  END SUBROUTINE eval_bottom_fissural_injection_analytic
END MODULE mass_bottom_fissural_injection_2d
