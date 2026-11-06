!> @brief Implement the four bottom fissural injection treatment interfaces.
!> @details Injection physics uses local state and geometric values; unsupported
!>          numerical treatments terminate with a diagnostic rather than silent rates.
MODULE mass_bottom_fissural_injection_2d
    USE parameters_2d, ONLY : idx_h, idx_T
    USE parameters_2d
    USE constitutive_parameters_2d
    USE state_types, ONLY : physical_state_real, physical_state_complex
    USE mass_injection_coefficients_2d, ONLY : injection_material_type

    USE mass_injection_coefficients_2d
    USE diagnostics_2d, ONLY : fatal_error
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_bottom_fissural_injection_explicit
    PUBLIC :: eval_bottom_fissural_injection_semi_implicit
    PUBLIC :: eval_bottom_fissural_injection_implicit
    PUBLIC :: eval_bottom_fissural_injection_analytic
CONTAINS
    !> @brief Evaluate explicit bottom fissural injection source.
    !> @details Computes cell-local rates using preconfigured source windows.
    !> @param[in] state Local physical cell state.
    !> @param[in] time Current source evaluation time [s].
    !> @param[in] cell_fissure_fract Local coverage fraction for each fissure.
    !> @param[out] rate Conservative injection rate vector.
    SUBROUTINE eval_bottom_fissural_injection_explicit(state, time, cell_fissure_fract, rate,materials)
    USE parameters_2d, ONLY : pore_pressure_flag, pore_pres_fract, bottom_fissural_source_flag, n_fissures,            &
        linear_vel_fissures, T_fissures, time_param_fissures
        TYPE(physical_state_real), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: time, cell_fissure_fract(:)
        REAL(wp), INTENT(OUT) :: rate(n_eqns)
        TYPE(injection_material_type), INTENT(IN) :: materials(:)
        REAL(wp) :: t_rem, t_coeff, h_dot
        REAL(wp) :: r_rho_m, r_red_grav
        REAL(wp) :: source_xs(n_solid), source_xg(n_add_gas)
        REAL(wp) :: r_sp_heat_mix, exc_pore_pres
        INTEGER :: i_fissure

        rate=0.0_wp

        ! Each fissure uses its own clipped cell coverage, pulse and temperature,
        ! and the source composition converted at that temperature.
        IF ( bottom_fissural_source_flag ) THEN

            DO i_fissure = 1, n_fissures

                IF ( cell_fissure_fract(i_fissure) .LE. 0.0_wp ) CYCLE

                t_rem = MOD(time + time_param_fissures(4,i_fissure),                                                   &
                    time_param_fissures(1,i_fissure))

                IF (time_param_fissures(3,i_fissure) .EQ. 0.0_wp) THEN

                    IF (t_rem .LE. time_param_fissures(2,i_fissure)) THEN
                        t_coeff = 1.0_wp
                    ELSE
                        t_coeff = 0.0_wp

                    END IF
                ELSEIF (t_rem .LT. time_param_fissures(3,i_fissure)) THEN
                    t_coeff = t_rem / time_param_fissures(3,i_fissure)
                ELSEIF (t_rem .LE. time_param_fissures(2,i_fissure) -                                                  &
                    time_param_fissures(3,i_fissure)) THEN
                    t_coeff = 1.0_wp
                ELSEIF (t_rem .LE. time_param_fissures(2,i_fissure)) THEN
                    t_coeff = 1.0_wp - (t_rem - time_param_fissures(2,i_fissure) +                                     &
                        time_param_fissures(3,i_fissure)) /                                                            &
                        time_param_fissures(3,i_fissure)
                ELSE
                    t_coeff = 0.0_wp

                END IF

                IF ( ALLOCATED( fissure_coeff_step ) )                                                                 &
                    t_coeff = fissure_coeff_step(i_fissure)

                IF (t_coeff .LE. 0.0_wp) CYCLE

                h_dot = cell_fissure_fract(i_fissure) * linear_vel_fissures(i_fissure)
                ! Material properties were evaluated once before the cell sweep.
                source_xs=materials(i_fissure)%xs
                source_xg=materials(i_fissure)%xg
                r_rho_m=materials(i_fissure)%rho_m
                r_sp_heat_mix=materials(i_fissure)%cp_m
                r_red_grav=materials(i_fissure)%red_grav

                rate(1) = rate(1) + t_coeff * h_dot * r_rho_m
                rate(4) = rate(4) + (t_coeff *                                                                         &
                    h_dot * r_rho_m * r_sp_heat_mix * T_fissures(i_fissure))
                rate(idx_solidEqn_first:idx_solidEqn_last) =                                                           &
                    rate(idx_solidEqn_first:idx_solidEqn_last) +                                                       &
                    t_coeff * h_dot * r_rho_m * source_xs
                rate(idx_addGasEqn_first:idx_addGasEqn_last) =                                                         &
                    rate(idx_addGasEqn_first:idx_addGasEqn_last) +                                                     &
                    t_coeff * h_dot * r_rho_m * source_xg


                IF (pore_pressure_flag) THEN
                    exc_pore_pres = state%exc_pore_pres
                    rate(idx_poreEqn) = rate(idx_poreEqn) +                                                            &
                        t_coeff * ((state%h*state%rho_m) * pore_pres_fract * h_dot *                                   &
                        r_rho_m * r_red_grav + exc_pore_pres * h_dot * r_rho_m)

                END IF

            END DO

        END IF
    END SUBROUTINE eval_bottom_fissural_injection_explicit


    !> @brief Reject unsupported semi-implicit bottom fissural injection source.
    !> @details This interface is preserved but fails explicitly when invoked.
    !> @param[in] state Local physical cell state.
    !> @param[in] time Current source evaluation time [s].
    !> @param[in] cell_fissure_fract Local coverage fraction for each fissure.
    !> @param[out] rate Conservative injection rate vector.
    SUBROUTINE eval_bottom_fissural_injection_semi_implicit(state, time, cell_fissure_fract, rate)
        TYPE(physical_state_real), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: time, cell_fissure_fract(:)
        REAL(wp), INTENT(OUT) :: rate(n_eqns)
        CALL fatal_error('bottom_fissural_injection: semi_implicit treatment is unsupported')
    END SUBROUTINE eval_bottom_fissural_injection_semi_implicit


    !> @brief Reject unsupported implicit bottom fissural injection source.
    !> @details This interface is preserved but fails explicitly when invoked.
    !> @param[in] state Local physical cell state.
    !> @param[in] time Current source evaluation time [s].
    !> @param[in] cell_fissure_fract Local coverage fraction for each fissure.
    !> @param[out] rate Conservative injection rate vector.
    SUBROUTINE eval_bottom_fissural_injection_implicit(state, time, cell_fissure_fract, rate)
        TYPE(physical_state_complex), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: time, cell_fissure_fract(:)
        COMPLEX(wp), INTENT(OUT) :: rate(n_eqns)
        CALL fatal_error('bottom_fissural_injection: implicit treatment is unsupported')
    END SUBROUTINE eval_bottom_fissural_injection_implicit


    !> @brief Reject unsupported analytic bottom fissural injection source.
    !> @details This interface is preserved but fails explicitly when invoked.
    !> @param[in] state Local physical cell state.
    !> @param[in] time Current source evaluation time [s].
    !> @param[in] cell_fissure_fract Local coverage fraction for each fissure.
    !> @param[out] rate Conservative injection rate vector.
    SUBROUTINE eval_bottom_fissural_injection_analytic(state, time, cell_fissure_fract, rate)
        TYPE(physical_state_real), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: time, cell_fissure_fract(:)
        REAL(wp), INTENT(OUT) :: rate(n_eqns)
        CALL fatal_error('bottom_fissural_injection: analytic treatment is unsupported')
    END SUBROUTINE eval_bottom_fissural_injection_analytic
END MODULE mass_bottom_fissural_injection_2d
