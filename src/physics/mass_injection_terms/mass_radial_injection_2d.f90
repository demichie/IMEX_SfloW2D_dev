!> @brief Implement the four radial lateral injection treatment interfaces.
!> @details Injection physics uses local state and geometric values; unsupported
!>          numerical treatments terminate with a diagnostic rather than silent rates.
MODULE mass_radial_injection_2d
    USE parameters_2d
    USE constitutive_parameters_2d
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
    USE state_new_mixture, ONLY : evaluate_mixture_density, evaluate_mixture_heat_capacity
    USE state_new_composition, ONLY : close_composition

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
    !> @brief Evaluate explicit radial lateral injection source.
    !> @details Computes cell-local rates using preconfigured source windows.
    !> @param[in] state Local physical cell state.
    !> @param[in] time Current source evaluation time [s].
    !> @param[in] lat_arc_perim_jk Local source arc perimeter [m].
    !> @param[in] lat_n_x_jk Local boundary normal x component.
    !> @param[in] lat_n_y_jk Local boundary normal y component.
    !> @param[in] cell_area_jk Local horizontal cell area [m2].
    !> @param[out] rate Conservative injection rate vector.
    SUBROUTINE eval_radial_injection_explicit(state, time, lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, cell_area_jk, rate)
        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: time, lat_arc_perim_jk
        REAL(wp), INTENT(IN) :: lat_n_x_jk, lat_n_y_jk, cell_area_jk
        REAL(wp), INTENT(OUT) :: rate(n_eqns)
        REAL(wp) :: t_rem, t_coeff, h_dot
        REAL(wp) :: r_rho_m, r_inv_rho_m, r_red_grav
        REAL(wp) :: source_xs(n_solid), source_xg(n_add_gas)
        REAL(wp) :: r_sp_heat_c, r_sp_heat_mix, exc_pore_pres

        rate=0.0_wp

        ! ----------- LATERAL RADIAL SOURCE (volume formulation) ------------------
        !
        ! Replaces the prior Dirichlet-BC implementation. For each in-arc ring
        ! cell, inject mass + momentum at a rate that integrates to the user-
        ! specified MFR. cell_arc_perim is rescaled in init_source so the
        ! per-cell rate rho_source * h_source * vel_source * cell_arc_perim_jk
        ! / cell_area_jk sums to MFR over all in-arc cells.

        IF ( radial_source_flag .AND. ( lat_arc_perim_jk .GT. 0.0_wp ) .AND.                                           &
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
                        t_coeff = 0.5_wp * ( 1.0_wp - COS( pi_g * t_rem /                                              &
                            time_param(3) ) )
                    ELSEIF ( t_rem .LE. ( time_param(2) - time_param(3) ) ) THEN
                        t_coeff = 1.0_wp
                    ELSEIF ( t_rem .LE. time_param(2) ) THEN
                        t_coeff = 0.5_wp * ( 1.0_wp + COS( pi_g * ( ( t_rem -                                          &
                            time_param(2) ) / time_param(3) + 1.0_wp ) ) )

                    END IF

                END IF

            END IF

            IF ( lateral_coeff_step .GE. 0.0_wp ) t_coeff = lateral_coeff_step

            IF ( t_coeff .GT. 0.0_wp ) THEN

                ! Source material is independent of the recipient cell.
                ! Close fractions before evaluating only the necessary properties.
                source_xs=xs_source(1:n_solid)
                source_xg=xg_source(1:n_add_gas)
                IF (.NOT.close_composition(source_xs,source_xg)) &
                    CALL fatal_error('Invalid injection material composition')
                CALL evaluate_mixture_density(T_source, &
                    source_xs,source_xg,r_rho_m,r_inv_rho_m)
                CALL evaluate_mixture_heat_capacity(source_xs,source_xg,r_sp_heat_mix,r_sp_heat_c)
                r_red_grav=(r_rho_m-rho_a_amb)*r_inv_rho_m*grav

                ! Per-cell injection rate (m/s) - integrates to MFR over all cells.
                h_dot = h_source * vel_source * lat_arc_perim_jk / cell_area_jk

                ! Mass equation.
                rate(1) = rate(1) + t_coeff * h_dot * r_rho_m

                ! Momentum equations - emit along the cell's outward unit normal.
                rate(2) = rate(2) + t_coeff * h_dot * r_rho_m * vel_source                                             &
                    * lat_n_x_jk
                rate(3) = rate(3) + t_coeff * h_dot * r_rho_m * vel_source                                             &
                    * lat_n_y_jk

            ! Thermal-energy equation.
                rate(4) = rate(4) + (t_coeff *                                                                         &
                    h_dot * r_rho_m * r_sp_heat_mix * T_source)

                ! Solid transport.
                rate(idx_solidEqn_first:idx_solidEqn_last) =                                                           &
                    rate(idx_solidEqn_first:idx_solidEqn_last) +                                                       &
                    t_coeff * h_dot * r_rho_m * xs_source(1:n_solid)

                ! Additional gases.

                rate(idx_addGasEqn_first:idx_addGasEqn_last) =                                                         &
                    rate(idx_addGasEqn_first:idx_addGasEqn_last) +                                                     &
                    t_coeff * h_dot * r_rho_m * xg_source(1:n_add_gas)


                IF ( pore_pressure_flag ) THEN
                    exc_pore_pres = state%exc_pore_pres
                    rate(idx_poreEqn) = rate(idx_poreEqn) +                                                            &
                        t_coeff * ((state%h*state%rho_m) * pore_pres_fract * h_dot *                                   &
                        r_rho_m * r_red_grav + exc_pore_pres * h_dot * r_rho_m)

                END IF

            END IF

        END IF
    END SUBROUTINE eval_radial_injection_explicit


    !> @brief Reject unsupported semi-implicit radial lateral injection source.
    !> @details This interface is preserved but fails explicitly when invoked.
    !> @param[in] state Local physical cell state.
    !> @param[in] time Current source evaluation time [s].
    !> @param[in] lat_arc_perim_jk Local source arc perimeter [m].
    !> @param[in] lat_n_x_jk Local boundary normal x component.
    !> @param[in] lat_n_y_jk Local boundary normal y component.
    !> @param[in] cell_area_jk Local horizontal cell area [m2].
    !> @param[out] rate Conservative injection rate vector.
    SUBROUTINE eval_radial_injection_semi_implicit(state, time, lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, cell_area_jk,&
        rate)
        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: time, lat_arc_perim_jk
        REAL(wp), INTENT(IN) :: lat_n_x_jk, lat_n_y_jk, cell_area_jk
        REAL(wp), INTENT(OUT) :: rate(n_eqns)
        CALL fatal_error('radial_injection: semi_implicit treatment is unsupported')
    END SUBROUTINE eval_radial_injection_semi_implicit


    !> @brief Reject unsupported implicit radial lateral injection source.
    !> @details This interface is preserved but fails explicitly when invoked.
    !> @param[in] state Local physical cell state.
    !> @param[in] time Current source evaluation time [s].
    !> @param[in] lat_arc_perim_jk Local source arc perimeter [m].
    !> @param[in] lat_n_x_jk Local boundary normal x component.
    !> @param[in] lat_n_y_jk Local boundary normal y component.
    !> @param[in] cell_area_jk Local horizontal cell area [m2].
    !> @param[out] rate Conservative injection rate vector.
    SUBROUTINE eval_radial_injection_implicit(state, time, lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, cell_area_jk, rate)
        TYPE(complex_physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: time, lat_arc_perim_jk
        REAL(wp), INTENT(IN) :: lat_n_x_jk, lat_n_y_jk, cell_area_jk
        COMPLEX(wp), INTENT(OUT) :: rate(n_eqns)
        CALL fatal_error('radial_injection: implicit treatment is unsupported')
    END SUBROUTINE eval_radial_injection_implicit


    !> @brief Reject unsupported analytic radial lateral injection source.
    !> @details This interface is preserved but fails explicitly when invoked.
    !> @param[in] state Local physical cell state.
    !> @param[in] time Current source evaluation time [s].
    !> @param[in] lat_arc_perim_jk Local source arc perimeter [m].
    !> @param[in] lat_n_x_jk Local boundary normal x component.
    !> @param[in] lat_n_y_jk Local boundary normal y component.
    !> @param[in] cell_area_jk Local horizontal cell area [m2].
    !> @param[out] rate Conservative injection rate vector.
    SUBROUTINE eval_radial_injection_analytic(state, time, lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, cell_area_jk, rate)
        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: time, lat_arc_perim_jk
        REAL(wp), INTENT(IN) :: lat_n_x_jk, lat_n_y_jk, cell_area_jk
        REAL(wp), INTENT(OUT) :: rate(n_eqns)
        CALL fatal_error('radial_injection: analytic treatment is unsupported')
    END SUBROUTINE eval_radial_injection_analytic
END MODULE mass_radial_injection_2d
