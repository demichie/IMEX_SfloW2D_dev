!> \brief Common physical source assembly entry points.
!> \details Owns the explicit STEP source-rate bridge independently of the
!> spatial flux and hydrostatic interfaces in equation_terms_2d.
!> The bounded mass-exchange transaction and all inventory updates remain
!> the responsibility of step_source_executor_2d.
MODULE source_terms_2d
  USE parameters_2d
  USE constitutive_parameters_2d
  USE state_conversion_2d, ONLY : mixt_var
  USE thermal_laws_2d, ONLY : eval_inlet_thermal_energy_source
  USE pore_pressure_laws_2d, ONLY : eval_pore_pressure_inlet_source
  USE geometry_2d, ONLY : pi_g, bottom_source_area_ratio
  USE mass_source_terms_2d, ONLY : eval_mass_exchange_terms
  USE rheology_2d, ONLY : curvature_acceleration
  USE diagnostics_2d, ONLY : fatal_error
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_step_source_explicit
  PUBLIC :: eval_geometry_curvature_explicit
  PUBLIC :: eval_geometry_curvature_semi_implicit
  PUBLIC :: eval_geometry_curvature_implicit
  PUBLIC :: solve_geometry_curvature_analytic
  PUBLIC :: bottom_coeff_step, bottom_vel_step, lateral_coeff_step, fissure_coeff_step
  PUBLIC :: eval_boundary_injection_explicit
  PUBLIC :: eval_boundary_injection_semi_implicit
  PUBLIC :: eval_boundary_injection_implicit
  PUBLIC :: solve_boundary_injection_analytic
  PUBLIC :: eval_bottom_circular_injection_explicit
  PUBLIC :: eval_bottom_circular_injection_semi_implicit
  PUBLIC :: eval_bottom_circular_injection_implicit
  PUBLIC :: solve_bottom_circular_injection_analytic
  PUBLIC :: eval_bottom_fissural_injection_explicit
  PUBLIC :: eval_bottom_fissural_injection_semi_implicit
  PUBLIC :: eval_bottom_fissural_injection_implicit
  PUBLIC :: solve_bottom_fissural_injection_analytic
  PUBLIC :: eval_radial_injection_explicit
  PUBLIC :: eval_radial_injection_semi_implicit
  PUBLIC :: eval_radial_injection_implicit
  PUBLIC :: solve_radial_injection_analytic

  ! Shared source-window coefficients updated once before the RK stages.
  ! They remain private to the physics source layer, but are imported by
  ! the legacy window initializer until its migration.
  REAL(wp) :: bottom_coeff_step = -1.0_wp
  REAL(wp) :: bottom_vel_step = -1.0_wp
  REAL(wp) :: lateral_coeff_step = -1.0_wp
  REAL(wp), ALLOCATABLE :: fissure_coeff_step(:)
CONTAINS

  SUBROUTINE eval_boundary_injection_explicit( time, vect_x , vect_y , source_bdry )

    USE parameters_2d, ONLY : h_source , vel_source , T_source , xs_source,    &
         xg_source, xl_source, time_param

    USE geometry_2d, ONLY : pi_g

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: time
    REAL(wp), INTENT(IN) :: vect_x
    REAL(wp), INTENT(IN) :: vect_y
    REAL(wp), INTENT(OUT) :: source_bdry(n_vars+2)

    REAL(wp) :: t_rem
    REAL(wp) :: t_coeff

    source_bdry = 0.0_wp
    source_bdry(4) = T_source

    IF ( time .GE. time_param(4) ) THEN

       RETURN

    END IF

    t_rem = MOD( time , time_param(1) )

    t_coeff = 0.0_wp

    IF ( time_param(3) .EQ. 0.0_wp ) THEN

       IF ( t_rem .LE. time_param(2) ) t_coeff = 1.0_wp

    ELSE

       IF ( t_rem .LT. time_param(3) ) THEN

          t_coeff = 0.5_wp * ( 1.0_wp - COS( pi_g * t_rem / time_param(3) ) )

       ELSEIF ( t_rem .LE. ( time_param(2) - time_param(3) ) ) THEN

          t_coeff = 1.0_wp

       ELSEIF ( t_rem .LE. time_param(2) ) THEN

          t_coeff = 0.5_wp * ( 1.0_wp + COS( pi_g * ( ( t_rem - time_param(2) ) &
               / time_param(3) + 1.0_wp ) ) )

       END IF

    END IF

    ! Preserve the common dry-state convention while the periodic source is
    ! inactive. In particular, no composition is attached to zero thickness.
    IF (t_coeff .LE. EPSILON(1.0_wp)) RETURN

    ! The exponents of t_coeff are such that Ri does not depend on t_coeff
    source_bdry(1) = t_coeff * h_source
    source_bdry(2) = t_coeff**1.5_wp * h_source * vel_source * vect_x
    source_bdry(3) = t_coeff**1.5_wp * h_source * vel_source * vect_y
    source_bdry(4) = T_source

    source_bdry(idx_solid_first:idx_solid_last) = xs_source(1:n_solid)
    source_bdry(idx_add_gas_first:idx_add_gas_last) = xg_source(1:n_add_gas)
    IF ( gas_flag .AND. liquid_flag ) source_bdry(n_vars) = xl_source

    source_bdry(idx_u) = t_coeff**0.5_wp * vel_source * vect_x
    source_bdry(idx_v) = t_coeff**0.5_wp * vel_source * vect_y

    RETURN

  END SUBROUTINE eval_boundary_injection_explicit

  SUBROUTINE eval_boundary_injection_semi_implicit(time, vect_x, vect_y, source_bdry)
    REAL(wp), INTENT(IN) :: time, vect_x, vect_y
    REAL(wp), INTENT(OUT) :: source_bdry(n_vars+2)
    CALL fatal_error('boundary_injection: semi_implicit treatment is unsupported')
  END SUBROUTINE eval_boundary_injection_semi_implicit

  SUBROUTINE eval_boundary_injection_implicit(time, vect_x, vect_y, source_bdry)
    REAL(wp), INTENT(IN) :: time, vect_x, vect_y
    REAL(wp), INTENT(OUT) :: source_bdry(n_vars+2)
    CALL fatal_error('boundary_injection: implicit treatment is unsupported')
  END SUBROUTINE eval_boundary_injection_implicit

  SUBROUTINE solve_boundary_injection_analytic(time, vect_x, vect_y, source_bdry)
    REAL(wp), INTENT(IN) :: time, vect_x, vect_y
    REAL(wp), INTENT(OUT) :: source_bdry(n_vars+2)
    CALL fatal_error('boundary_injection: analytic treatment is unsupported')
  END SUBROUTINE solve_boundary_injection_analytic

  SUBROUTINE eval_bottom_circular_injection_explicit(qpj, time, cell_fract_jk, q1, expl_term)
    USE parameters_2d, ONLY : vel_source, T_source, xs_source, xg_source, &
         xl_source, time_param, bottom_radial_source_flag, pore_pressure_flag, &
         pore_pres_fract, n_intervals, t_intervals, vel_intervals, &
         radial_source_flag, h_source, bottom_fissural_source_flag, &
         n_fissures, linear_vel_fissures, T_fissures, time_param_fissures, &
         xs_fissures, xg_fissures, xl_fissures
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
    qp_source(1) = 1.0_wp
    qp_source(2) = 0.0_wp
    qp_source(3) = 0.0_wp
    qp_source(4) = t_source

    qp_source(idx_solid_first:idx_solid_last) = xs_source(1:n_solid)
    qp_source(idx_add_gas_first:idx_add_gas_last) = xg_source(1:n_add_gas)
    IF ( gas_flag .AND. liquid_flag ) qp_source(n_vars) = xl_source

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

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! source term for the liquid phase
       expl_term(n_vars) = expl_term(n_vars) + t_coeff * h_dot * r_rho_m       &
            * xl_source

    END IF

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

  SUBROUTINE eval_bottom_fissural_injection_explicit(qpj, time, cell_fissure_fract, q1, expl_term)
    USE parameters_2d, ONLY : vel_source, T_source, xs_source, xg_source, &
         xl_source, time_param, bottom_radial_source_flag, pore_pressure_flag, &
         pore_pres_fract, n_intervals, t_intervals, vel_intervals, &
         radial_source_flag, h_source, bottom_fissural_source_flag, &
         n_fissures, linear_vel_fissures, T_fissures, time_param_fissures, &
         xs_fissures, xg_fissures, xl_fissures
    USE geometry_2d, ONLY : pi_g, bottom_source_area_ratio
    REAL(wp), INTENT(IN) :: qpj(n_vars+2), time, cell_fissure_fract(:), q1
    REAL(wp), INTENT(INOUT) :: expl_term(n_eqns)
    REAL(wp) :: t_rem, t_coeff, h_dot, vel_local
    REAL(wp) :: qp_source(n_vars+2)
    REAL(wp) :: r_Ri, r_rho_m, r_rho_c, r_red_grav
    REAL(wp) :: r_sp_heat_c, r_sp_heat_mix, exc_pore_pres
    INTEGER :: i_int, i_fissure

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
          qp_source(1) = 1.0_wp
          qp_source(4) = T_fissures(i_fissure)
          qp_source(idx_solid_first:idx_solid_last) =                          &
               xs_fissures(1:n_solid,i_fissure)
          qp_source(idx_add_gas_first:idx_add_gas_last) =                      &
               xg_fissures(1:n_add_gas,i_fissure)
          IF (gas_flag .AND. liquid_flag) qp_source(n_vars) =                  &
               xl_fissures(i_fissure)
          IF (stoch_transport_flag) qp_source(idx_stoch) = 0.0_wp
          IF (pore_pressure_flag) qp_source(idx_poreEqn) = 0.0_wp
          qp_source(idx_u) = 0.0_wp
          qp_source(idx_v) = 0.0_wp

          CALL mixt_var(qp_source, r_Ri, r_rho_m, r_rho_c, r_red_grav,       &
               r_sp_heat_c, r_sp_heat_mix)

          expl_term(1) = expl_term(1) + t_coeff * h_dot * r_rho_m
          expl_term(4) = expl_term(4) + eval_inlet_thermal_energy_source(      &
               t_coeff, h_dot, r_rho_m, r_sp_heat_mix, T_fissures(i_fissure))
          expl_term(idx_solidEqn_first:idx_solidEqn_last) =                   &
               expl_term(idx_solidEqn_first:idx_solidEqn_last) +              &
               t_coeff * h_dot * r_rho_m * xs_fissures(1:n_solid,i_fissure)
          expl_term(idx_addGasEqn_first:idx_addGasEqn_last) =                  &
               expl_term(idx_addGasEqn_first:idx_addGasEqn_last) +             &
               t_coeff * h_dot * r_rho_m * xg_fissures(1:n_add_gas,i_fissure)

          IF (gas_flag .AND. liquid_flag) THEN
             expl_term(n_vars) = expl_term(n_vars) + t_coeff * h_dot         &
                  * r_rho_m * xl_fissures(i_fissure)
          END IF

          IF (pore_pressure_flag) THEN
             exc_pore_pres = qpj(idx_pore)
             expl_term(idx_poreEqn) = expl_term(idx_poreEqn) +                &
                  eval_pore_pressure_inlet_source(t_coeff, q1,                 &
                  pore_pres_fract, h_dot, r_rho_m, r_red_grav, exc_pore_pres)
          END IF

       END DO

    END IF
  END SUBROUTINE eval_bottom_fissural_injection_explicit

  SUBROUTINE eval_bottom_fissural_injection_semi_implicit(qpj, time, cell_fissure_fract, q1, expl_term)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2), time, cell_fissure_fract(:), q1
    REAL(wp), INTENT(INOUT) :: expl_term(n_eqns)
    CALL fatal_error('bottom_fissural_injection: semi_implicit treatment is unsupported')
  END SUBROUTINE eval_bottom_fissural_injection_semi_implicit

  SUBROUTINE eval_bottom_fissural_injection_implicit(qpj, time, cell_fissure_fract, q1, expl_term)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2), time, cell_fissure_fract(:), q1
    REAL(wp), INTENT(INOUT) :: expl_term(n_eqns)
    CALL fatal_error('bottom_fissural_injection: implicit treatment is unsupported')
  END SUBROUTINE eval_bottom_fissural_injection_implicit

  SUBROUTINE solve_bottom_fissural_injection_analytic(qpj, time, cell_fissure_fract, q1, expl_term)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2), time, cell_fissure_fract(:), q1
    REAL(wp), INTENT(INOUT) :: expl_term(n_eqns)
    CALL fatal_error('bottom_fissural_injection: analytic treatment is unsupported')
  END SUBROUTINE solve_bottom_fissural_injection_analytic

  SUBROUTINE eval_radial_injection_explicit(qpj, time, lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, &
       cell_area_jk, q1, expl_term)
    USE parameters_2d, ONLY : vel_source, T_source, xs_source, xg_source, &
         xl_source, time_param, bottom_radial_source_flag, pore_pressure_flag, &
         pore_pres_fract, n_intervals, t_intervals, vel_intervals, &
         radial_source_flag, h_source, bottom_fissural_source_flag, &
         n_fissures, linear_vel_fissures, T_fissures, time_param_fissures, &
         xs_fissures, xg_fissures, xl_fissures
    USE geometry_2d, ONLY : pi_g, bottom_source_area_ratio
    REAL(wp), INTENT(IN) :: qpj(n_vars+2), time, lat_arc_perim_jk
    REAL(wp), INTENT(IN) :: lat_n_x_jk, lat_n_y_jk, cell_area_jk, q1
    REAL(wp), INTENT(INOUT) :: expl_term(n_eqns)
    REAL(wp) :: t_rem, t_coeff, h_dot, vel_local
    REAL(wp) :: qp_source(n_vars+2)
    REAL(wp) :: r_Ri, r_rho_m, r_rho_c, r_red_grav
    REAL(wp) :: r_sp_heat_c, r_sp_heat_mix, exc_pore_pres
    INTEGER :: i_int, i_fissure

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
          qp_source(1) = 1.0_wp
          qp_source(4) = T_source

          qp_source(idx_solid_first:idx_solid_last) = xs_source(1:n_solid)
          qp_source(idx_add_gas_first:idx_add_gas_last) = xg_source(1:n_add_gas)
          IF ( gas_flag .AND. liquid_flag ) qp_source(n_vars) = xl_source

          IF ( stoch_transport_flag ) qp_source(idx_stoch) = 0.0_wp
          IF ( pore_pressure_flag )   qp_source(idx_poreEqn) = 0.0_wp

          qp_source(idx_u) = 0.0_wp
          qp_source(idx_v) = 0.0_wp

          CALL mixt_var( qp_source, r_Ri, r_rho_m, r_rho_c, r_red_grav,        &
               r_sp_heat_c, r_sp_heat_mix )

          ! Per-cell injection rate (m/s) - integrates to MFR over all cells.
          h_dot = h_source * vel_source * lat_arc_perim_jk / cell_area_jk

          ! Mass equation.
          expl_term(1) = expl_term(1) + t_coeff * h_dot * r_rho_m

          ! Momentum equations - emit along the cell's outward unit normal.
          expl_term(2) = expl_term(2) + t_coeff * h_dot * r_rho_m * vel_source &
               * lat_n_x_jk
          expl_term(3) = expl_term(3) + t_coeff * h_dot * r_rho_m * vel_source &
               * lat_n_y_jk

       ! Thermal-energy equation.
          expl_term(4) = expl_term(4) + eval_inlet_thermal_energy_source(     &
               t_coeff, h_dot, r_rho_m, r_sp_heat_mix, T_source)

          ! Solid transport.
          expl_term(idx_solidEqn_first:idx_solidEqn_last) =                   &
               expl_term(idx_solidEqn_first:idx_solidEqn_last) +              &
               t_coeff * h_dot * r_rho_m * xs_source(1:n_solid)

          ! Additional gases.

          expl_term(idx_addGasEqn_first:idx_addGasEqn_last) =                  &
               expl_term(idx_addGasEqn_first:idx_addGasEqn_last) +             &
               t_coeff * h_dot * r_rho_m * xg_source(1:n_add_gas)

          IF ( gas_flag .AND. liquid_flag ) THEN
             expl_term(n_vars) = expl_term(n_vars) +                           &
                  t_coeff * h_dot * r_rho_m * xl_source
          END IF

          IF ( pore_pressure_flag ) THEN
             exc_pore_pres = qpj(idx_pore)
             expl_term(idx_poreEqn) = expl_term(idx_poreEqn) +                &
                  eval_pore_pressure_inlet_source(t_coeff, q1,                 &
                  pore_pres_fract, h_dot, r_rho_m, r_red_grav, exc_pore_pres)
          END IF

       END IF

    END IF
  END SUBROUTINE eval_radial_injection_explicit

  SUBROUTINE eval_radial_injection_semi_implicit(qpj, time, lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, &
       cell_area_jk, q1, expl_term)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2), time, lat_arc_perim_jk
    REAL(wp), INTENT(IN) :: lat_n_x_jk, lat_n_y_jk, cell_area_jk, q1
    REAL(wp), INTENT(INOUT) :: expl_term(n_eqns)
    CALL fatal_error('radial_injection: semi_implicit treatment is unsupported')
  END SUBROUTINE eval_radial_injection_semi_implicit

  SUBROUTINE eval_radial_injection_implicit(qpj, time, lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, &
       cell_area_jk, q1, expl_term)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2), time, lat_arc_perim_jk
    REAL(wp), INTENT(IN) :: lat_n_x_jk, lat_n_y_jk, cell_area_jk, q1
    REAL(wp), INTENT(INOUT) :: expl_term(n_eqns)
    CALL fatal_error('radial_injection: implicit treatment is unsupported')
  END SUBROUTINE eval_radial_injection_implicit

  SUBROUTINE solve_radial_injection_analytic(qpj, time, lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, &
       cell_area_jk, q1, expl_term)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2), time, lat_arc_perim_jk
    REAL(wp), INTENT(IN) :: lat_n_x_jk, lat_n_y_jk, cell_area_jk, q1
    REAL(wp), INTENT(INOUT) :: expl_term(n_eqns)
    CALL fatal_error('radial_injection: analytic treatment is unsupported')
  END SUBROUTINE solve_radial_injection_analytic


  PURE SUBROUTINE eval_geometry_curvature_explicit( Bprime_x, Bprime_y,         &
       Bsecond_xx, Bsecond_xy, Bsecond_yy, grav_coefficient,                 &
       mixture_mass, velocity_x, velocity_y, source_x, source_y )

    REAL(wp), INTENT(IN) :: Bprime_x, Bprime_y
    REAL(wp), INTENT(IN) :: Bsecond_xx, Bsecond_xy, Bsecond_yy
    REAL(wp), INTENT(IN) :: grav_coefficient
    REAL(wp), INTENT(IN) :: mixture_mass, velocity_x, velocity_y
    REAL(wp), INTENT(OUT) :: source_x, source_y
    REAL(wp) :: acceleration

    acceleration = curvature_acceleration( Bsecond_xx, Bsecond_xy,          &
         Bsecond_yy, velocity_x, velocity_y )
    source_x = -grav_coefficient * mixture_mass * acceleration * Bprime_x
    source_y = -grav_coefficient * mixture_mass * acceleration * Bprime_y

  END SUBROUTINE eval_geometry_curvature_explicit

  SUBROUTINE eval_geometry_curvature_semi_implicit(Bprime_x, Bprime_y, &
       Bsecond_xx, Bsecond_xy, Bsecond_yy, grav_coefficient, &
       mixture_mass, velocity_x, velocity_y, source_x, source_y)
    REAL(wp), INTENT(IN) :: Bprime_x, Bprime_y
    REAL(wp), INTENT(IN) :: Bsecond_xx, Bsecond_xy, Bsecond_yy
    REAL(wp), INTENT(IN) :: grav_coefficient
    REAL(wp), INTENT(IN) :: mixture_mass, velocity_x, velocity_y
    REAL(wp), INTENT(OUT) :: source_x, source_y
    CALL fatal_error('geometry_curvature: semi_implicit treatment is unsupported')
  END SUBROUTINE eval_geometry_curvature_semi_implicit

  SUBROUTINE eval_geometry_curvature_implicit(Bprime_x, Bprime_y, &
       Bsecond_xx, Bsecond_xy, Bsecond_yy, grav_coefficient, &
       mixture_mass, velocity_x, velocity_y, source_x, source_y)
    REAL(wp), INTENT(IN) :: Bprime_x, Bprime_y
    REAL(wp), INTENT(IN) :: Bsecond_xx, Bsecond_xy, Bsecond_yy
    REAL(wp), INTENT(IN) :: grav_coefficient
    REAL(wp), INTENT(IN) :: mixture_mass, velocity_x, velocity_y
    REAL(wp), INTENT(OUT) :: source_x, source_y
    CALL fatal_error('geometry_curvature: implicit treatment is unsupported')
  END SUBROUTINE eval_geometry_curvature_implicit

  SUBROUTINE solve_geometry_curvature_analytic(Bprime_x, Bprime_y, &
       Bsecond_xx, Bsecond_xy, Bsecond_yy, grav_coefficient, &
       mixture_mass, velocity_x, velocity_y, source_x, source_y)
    REAL(wp), INTENT(IN) :: Bprime_x, Bprime_y
    REAL(wp), INTENT(IN) :: Bsecond_xx, Bsecond_xy, Bsecond_yy
    REAL(wp), INTENT(IN) :: grav_coefficient
    REAL(wp), INTENT(IN) :: mixture_mass, velocity_x, velocity_y
    REAL(wp), INTENT(OUT) :: source_x, source_y
    CALL fatal_error('geometry_curvature: analytic treatment is unsupported')
  END SUBROUTINE solve_geometry_curvature_analytic

  !> \brief Evaluate the atomic explicit STEP mass-exchange rate.
  !> \details Returns rates only; it never mutates conservative state, bed or inventory.
  !> \param[in] qpj Cell-local primitive state.
  !> \param[in] B_zone Connected-water-region indicator.
  !> \param[in] B_prime_x Filtered bed x derivative.
  !> \param[in] B_prime_y Filtered bed y derivative.
  !> \param[in] erodible Per-class available substrate inventory.
  !> \param[in] dt Time-step length.
  !> \param[out] erosion_term Per-class substrate erosion rate.
  !> \param[out] deposition_term Per-class deposition rate.
  !> \param[out] continuous_phase_erosion_term Carrier entrainment rate.
  !> \param[out] continuous_phase_loss_term Limited carrier loss rate.
  !> \param[out] eqns_term Conservative equation source-rate vector.
  !> \param[out] topo_term Proposed bed-elevation rate.
  SUBROUTINE eval_step_source_explicit(qpj, B_zone, B_prime_x, B_prime_y, &
       erodible, dt, erosion_term, deposition_term, &
       continuous_phase_erosion_term, continuous_phase_loss_term, eqns_term, topo_term)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2), B_prime_x, B_prime_y
    REAL(wp), INTENT(IN) :: erodible(n_solid), dt
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(OUT) :: erosion_term(n_solid), deposition_term(n_solid)
    REAL(wp), INTENT(OUT) :: continuous_phase_erosion_term, continuous_phase_loss_term
    REAL(wp), INTENT(OUT) :: eqns_term(n_eqns), topo_term

    CALL eval_mass_exchange_terms(qpj, B_zone, B_prime_x, B_prime_y, &
         erodible, dt, erosion_term, deposition_term, &
         continuous_phase_erosion_term, continuous_phase_loss_term, eqns_term, topo_term)
  END SUBROUTINE eval_step_source_explicit
END MODULE source_terms_2d
