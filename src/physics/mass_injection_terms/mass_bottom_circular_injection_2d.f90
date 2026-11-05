!> @brief Implement the four bottom circular injection treatment interfaces.
!> @details Injection physics uses local state and geometric values; unsupported
!>          numerical treatments terminate with a diagnostic rather than silent rates.
MODULE mass_bottom_circular_injection_2d
    USE parameters_2d, ONLY : idx_h, idx_T
    USE parameters_2d
    USE constitutive_parameters_2d
    USE state_types, ONLY : physical_state_real, physical_state_complex
    USE mass_injection_coefficients_2d, ONLY : injection_material_type

    USE geometry_2d, ONLY : pi_g, bottom_source_area_ratio
    USE mass_injection_coefficients_2d
    USE diagnostics_2d, ONLY : fatal_error
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_bottom_circular_injection_explicit
    PUBLIC :: eval_bottom_circular_injection_semi_implicit
    PUBLIC :: eval_bottom_circular_injection_implicit
    PUBLIC :: eval_bottom_circular_injection_analytic
CONTAINS
    !> @brief Evaluate explicit bottom circular injection source.
    !> @details Computes cell-local rates using preconfigured source windows.
    !> @param[in] state Local physical cell state.
    !> @param[in] time Current source evaluation time [s].
    !> @param[in] cell_fract_jk Fraction of the cell covered by bottom circular injection.
    !> @param[out] rate Conservative injection rate vector.
    SUBROUTINE eval_bottom_circular_injection_explicit(state, time, cell_fract_jk, rate,material)
    USE parameters_2d, ONLY : vel_source, T_source, xs_source, xg_source,                                              &
            time_param, bottom_radial_source_flag, pore_pressure_flag,                                                 &
            pore_pres_fract, n_intervals, t_intervals, vel_intervals,                                                  &
            radial_source_flag, h_source, bottom_fissural_source_flag,                                                 &
            n_fissures, linear_vel_fissures, T_fissures, time_param_fissures,                                          &
            xs_fissures, xg_fissures
    USE geometry_2d, ONLY : pi_g, bottom_source_area_ratio
        TYPE(physical_state_real), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: time, cell_fract_jk
        REAL(wp), INTENT(OUT) :: rate(n_eqns)
        TYPE(injection_material_type), INTENT(IN) :: material
        REAL(wp) :: t_rem, t_coeff, h_dot, vel_local
        REAL(wp) :: r_rho_m, r_red_grav
        REAL(wp) :: source_xs(n_solid), source_xg(n_add_gas)
        REAL(wp) :: r_sp_heat_mix, exc_pore_pres
        INTEGER :: i_int

        rate=0.0_wp

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

                        t_coeff = 1.0_wp - ( t_rem - time_param(2) + time_param(3) ) /                                 &
                            time_param(3)

                    ELSE

                        t_coeff = 0.0_wp

                    END IF

                END IF

                IF ( bottom_coeff_step .GE. 0.0_wp ) t_coeff = bottom_coeff_step

                h_dot = cell_fract_jk * vel_source * bottom_source_area_ratio

            END IF

            ! Material properties were evaluated once before the cell sweep.
            source_xs=material%xs
            source_xg=material%xg
            r_rho_m=material%rho_m
            r_sp_heat_mix=material%cp_m
            r_red_grav=material%red_grav

            rate(1) = rate(1) + t_coeff * h_dot * r_rho_m
            rate(2) = rate(2) + 0.0_wp
            rate(3) = rate(3) + 0.0_wp

            rate(4) = rate(4) + (t_coeff *                                                                             &
                h_dot * r_rho_m * r_sp_heat_mix * t_source)

            ! source terms for the solid equations
            rate(idx_solidEqn_first:idx_solidEqn_last) =                                                               &
                rate(idx_solidEqn_first:idx_solidEqn_last) + t_coeff                                                   &
                * h_dot * r_rho_m * source_xs

            ! source terms for the additional gas equations
            rate(idx_addGasEqn_first:idx_addGasEqn_last) =                                                             &
                rate(idx_addGasEqn_first:idx_addGasEqn_last) + t_coeff                                                 &
                * h_dot * r_rho_m * source_xg


            IF ( pore_pressure_flag ) THEN

                exc_pore_pres = state%exc_pore_pres

                ! Pore-pressure contribution uses the existing cell mixture mass.
                ! the units of this source term are: kg^2 m^-3 s^-3
                rate(idx_poreEqn) = rate(idx_poreEqn) +                                                                &
                    t_coeff * ((state%h*state%rho_m) * pore_pres_fract * h_dot *                                       &
                            r_rho_m * r_red_grav + exc_pore_pres * h_dot * r_rho_m)

            END IF

        END IF   ! bottom_radial_source_flag
    END SUBROUTINE eval_bottom_circular_injection_explicit


    !> @brief Reject unsupported semi-implicit bottom circular injection source.
    !> @details This interface is preserved but fails explicitly when invoked.
    !> @param[in] state Local physical cell state.
    !> @param[in] time Current source evaluation time [s].
    !> @param[in] cell_fract_jk Fraction of the cell covered by bottom circular injection.
    !> @param[out] rate Conservative injection rate vector.
    SUBROUTINE eval_bottom_circular_injection_semi_implicit(state, time, cell_fract_jk, rate)
        TYPE(physical_state_real), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: time, cell_fract_jk
        REAL(wp), INTENT(OUT) :: rate(n_eqns)
        CALL fatal_error('bottom_circular_injection: semi_implicit treatment is unsupported')
    END SUBROUTINE eval_bottom_circular_injection_semi_implicit


    !> @brief Reject unsupported implicit bottom circular injection source.
    !> @details This interface is preserved but fails explicitly when invoked.
    !> @param[in] state Local physical cell state.
    !> @param[in] time Current source evaluation time [s].
    !> @param[in] cell_fract_jk Fraction of the cell covered by bottom circular injection.
    !> @param[out] rate Conservative injection rate vector.
    SUBROUTINE eval_bottom_circular_injection_implicit(state, time, cell_fract_jk, rate)
        TYPE(physical_state_complex), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: time, cell_fract_jk
        COMPLEX(wp), INTENT(OUT) :: rate(n_eqns)
        CALL fatal_error('bottom_circular_injection: implicit treatment is unsupported')
    END SUBROUTINE eval_bottom_circular_injection_implicit


    !> @brief Reject unsupported analytic bottom circular injection source.
    !> @details This interface is preserved but fails explicitly when invoked.
    !> @param[in] state Local physical cell state.
    !> @param[in] time Current source evaluation time [s].
    !> @param[in] cell_fract_jk Fraction of the cell covered by bottom circular injection.
    !> @param[out] rate Conservative injection rate vector.
    SUBROUTINE eval_bottom_circular_injection_analytic(state, time, cell_fract_jk, rate)
        TYPE(physical_state_real), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: time, cell_fract_jk
        REAL(wp), INTENT(OUT) :: rate(n_eqns)
        CALL fatal_error('bottom_circular_injection: analytic treatment is unsupported')
    END SUBROUTINE eval_bottom_circular_injection_analytic
END MODULE mass_bottom_circular_injection_2d
