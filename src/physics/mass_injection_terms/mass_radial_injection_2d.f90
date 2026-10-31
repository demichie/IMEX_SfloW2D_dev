!> \brief Physical treatment interfaces for radial_injection.
MODULE mass_radial_injection_2d
    USE parameters_2d, ONLY : idx_h, idx_T
  USE parameters_2d
  USE constitutive_parameters_2d
  USE state_conversion_2d, ONLY : mixt_var, physical_state_type, complex_physical_state_type

  USE geometry_2d, ONLY : pi_g, bottom_source_area_ratio
  USE mass_injection_coefficients_2d
  USE diagnostics_2d, ONLY : fatal_error
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_radial_injection_explicit
  PUBLIC :: eval_radial_injection_semi_implicit
  PUBLIC :: eval_radial_injection_implicit
  PUBLIC :: eval_radial_injection_analytic
CONTAINS
  SUBROUTINE eval_radial_injection_explicit(state, time, lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, &
       cell_area_jk, rate)
    USE parameters_2d, ONLY : vel_source, T_source, xs_source, xg_source, &
         time_param, bottom_radial_source_flag, pore_pressure_flag, &
         pore_pres_fract, n_intervals, t_intervals, vel_intervals, &
         radial_source_flag, h_source, bottom_fissural_source_flag, &
         n_fissures, linear_vel_fissures, T_fissures, time_param_fissures, &
         xs_fissures, xg_fissures
    USE geometry_2d, ONLY : pi_g, bottom_source_area_ratio
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: time, lat_arc_perim_jk
    REAL(wp), INTENT(IN) :: lat_n_x_jk, lat_n_y_jk, cell_area_jk
    REAL(wp), INTENT(OUT) :: rate(n_eqns)
    REAL(wp) :: t_rem, t_coeff, h_dot
    REAL(wp) :: qp_source(n_vars+2)
    REAL(wp) :: r_Ri, r_rho_m, r_rho_c, r_red_grav
    REAL(wp) :: r_sp_heat_c, r_sp_heat_mix, exc_pore_pres

    rate=0.0_wp

    ! ----------- LATERAL RADIAL SOURCE (volume formulation) ------------------
    !
    ! Replaces the prior Dirichlet-BC implementation. For each in-arc ring
    ! cell, inject mass + momentum at a rate that integrates to the user-
    ! specified MFR. cell_arc_perim is rescaled in init_source so the
    ! per-cell rate rho_source * h_source * vel_source * cell_arc_perim_jk
    ! / cell_area_jk sums to MFR over all in-arc cells.

    IF ( radial_source_flag .AND. ( lat_arc_perim_jk .GT. 0.0_wp ) .AND.       &
         ( h_source .GT. 0.0_wp ) ) THEN

       ! Time gate: keep the current-main semantics used by eval_boundary_injection_explicit.
       ! TIME_PARAM(4) is an absolute cut-off for the lateral radial source.
       IF ( time .GE. time_param(4) ) THEN

          t_coeff = 0.0_wp

       ELSE

          t_rem = MOD( time , time_param(1) )
          t_coeff = 0.0_wp

          IF ( time_param(3) .EQ. 0.0_wp ) THEN
             IF ( t_rem .LE. time_param(2) ) t_coeff = 1.0_wp
          ELSE
             IF ( t_rem .LT. time_param(3) ) THEN
                t_coeff = 0.5_wp * ( 1.0_wp - COS( pi_g * t_rem /              &
                     time_param(3) ) )
             ELSEIF ( t_rem .LE. ( time_param(2) - time_param(3) ) ) THEN
                t_coeff = 1.0_wp
             ELSEIF ( t_rem .LE. time_param(2) ) THEN
                t_coeff = 0.5_wp * ( 1.0_wp + COS( pi_g * ( ( t_rem -          &
                     time_param(2) ) / time_param(3) + 1.0_wp ) ) )
             END IF
          END IF

       END IF

       IF ( lateral_coeff_step .GE. 0.0_wp ) t_coeff = lateral_coeff_step

       IF ( t_coeff .GT. 0.0_wp ) THEN

          qp_source(:) = 0.0_wp
          qp_source(idx_h) = 1.0_wp
          qp_source(idx_T) = T_source

          qp_source(idx_solid_first:idx_solid_last) = xs_source(1:n_solid)
          qp_source(idx_add_gas_first:idx_add_gas_last) = xg_source(1:n_add_gas)

          IF ( stoch_transport_flag ) qp_source(idx_stoch) = 0.0_wp
          IF ( pore_pressure_flag )   qp_source(idx_poreEqn) = 0.0_wp

          qp_source(idx_u) = 0.0_wp
          qp_source(idx_v) = 0.0_wp

          CALL mixt_var( qp_source, r_Ri, r_rho_m, r_rho_c, r_red_grav,        &
               r_sp_heat_c, r_sp_heat_mix )

          ! Per-cell injection rate (m/s) - integrates to MFR over all cells.
          h_dot = h_source * vel_source * lat_arc_perim_jk / cell_area_jk

          ! Mass equation.
          rate(1) = rate(1) + t_coeff * h_dot * r_rho_m

          ! Momentum equations - emit along the cell's outward unit normal.
          rate(2) = rate(2) + t_coeff * h_dot * r_rho_m * vel_source &
               * lat_n_x_jk
          rate(3) = rate(3) + t_coeff * h_dot * r_rho_m * vel_source &
               * lat_n_y_jk

       ! Thermal-energy equation.
          rate(4) = rate(4) + (t_coeff * &
               h_dot * r_rho_m * r_sp_heat_mix * T_source)

          ! Solid transport.
          rate(idx_solidEqn_first:idx_solidEqn_last) =                   &
               rate(idx_solidEqn_first:idx_solidEqn_last) +              &
               t_coeff * h_dot * r_rho_m * xs_source(1:n_solid)

          ! Additional gases.

          rate(idx_addGasEqn_first:idx_addGasEqn_last) =                  &
               rate(idx_addGasEqn_first:idx_addGasEqn_last) +             &
               t_coeff * h_dot * r_rho_m * xg_source(1:n_add_gas)


          IF ( pore_pressure_flag ) THEN
             exc_pore_pres = state%exc_pore_pres
             rate(idx_poreEqn) = rate(idx_poreEqn) +                &
                  t_coeff * ((state%h*state%rho_m) * pore_pres_fract * h_dot * &
                  r_rho_m * r_red_grav + exc_pore_pres * h_dot * r_rho_m)
          END IF

       END IF

    END IF
  END SUBROUTINE eval_radial_injection_explicit

  SUBROUTINE eval_radial_injection_semi_implicit(state, time, lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, &
       cell_area_jk, rate)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: time, lat_arc_perim_jk
    REAL(wp), INTENT(IN) :: lat_n_x_jk, lat_n_y_jk, cell_area_jk
    REAL(wp), INTENT(OUT) :: rate(n_eqns)
    CALL fatal_error('radial_injection: semi_implicit treatment is unsupported')
  END SUBROUTINE eval_radial_injection_semi_implicit

  SUBROUTINE eval_radial_injection_implicit(state, time, lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, &
       cell_area_jk, rate)
    TYPE(complex_physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: time, lat_arc_perim_jk
    REAL(wp), INTENT(IN) :: lat_n_x_jk, lat_n_y_jk, cell_area_jk
    COMPLEX(wp), INTENT(OUT) :: rate(n_eqns)
    CALL fatal_error('radial_injection: implicit treatment is unsupported')
  END SUBROUTINE eval_radial_injection_implicit

  SUBROUTINE eval_radial_injection_analytic(state, time, lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, &
       cell_area_jk, rate)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: time, lat_arc_perim_jk
    REAL(wp), INTENT(IN) :: lat_n_x_jk, lat_n_y_jk, cell_area_jk
    REAL(wp), INTENT(OUT) :: rate(n_eqns)
    CALL fatal_error('radial_injection: analytic treatment is unsupported')
  END SUBROUTINE eval_radial_injection_analytic
END MODULE mass_radial_injection_2d
