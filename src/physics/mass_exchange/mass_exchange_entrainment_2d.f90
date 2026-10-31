!> @brief Evaluate entrainment process contributions to mass, energy and bed evolution.
!> @details Explicit treatment is supported; other treatment interfaces reject
!> unsupported configurations. Shared resource budgets belong to the coordinator.
MODULE mass_exchange_entrainment_2d
    USE parameters_2d, ONLY : wp, mass_exchange_rates_type, &
        reset_mass_exchange_rates
    USE diagnostics_2d, ONLY : fatal_error
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
    USE state_conversion_2d, ONLY : mixt_var_state
    USE constitutive_parameters_2d, ONLY : entrainment_flag, rho_a_amb, T_ambient, sp_heat_a
    IMPLICIT NONE
    PRIVATE

    PUBLIC :: eval_mass_exchange_entrainment_explicit
    PUBLIC :: eval_mass_exchange_entrainment_semi_implicit
    PUBLIC :: eval_mass_exchange_entrainment_implicit
    PUBLIC :: solve_mass_exchange_entrainment_analytic
    PUBLIC :: eval_air_entrainment_rate
    ! Four-entry physical treatment contract; only EXPLICIT is supported.

CONTAINS

    !> @brief Evaluate explicit ambient-air entrainment.
    !> @param[in] state Local real physical state.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_entrainment_explicit(state,r)

        TYPE(physical_state_type), INTENT(IN) :: state
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
        REAL(wp) :: speed_squared,Ri,rhom,rhoc,redgrav,cpc,cpm

        CALL reset_mass_exchange_rates(r)
        speed_squared=state%mod_vel2
        CALL mixt_var_state(state,Ri,rhom,rhoc,redgrav,cpc,cpm)

        ! Convert entrained ambient-air volume to mixture mass at ambient density.
        r%mixture_mass_rate=rho_a_amb*eval_entrainment_volume_kernel(state%h,Ri,speed_squared)
        r%energy_rate=T_ambient*sp_heat_a*r%mixture_mass_rate
    END SUBROUTINE

    !> @brief Reject unsupported semi-implicit treatment for entrainment.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] state Local real physical state.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_entrainment_semi_implicit(state,r)

        TYPE(physical_state_type), INTENT(IN) :: state
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('entrainment: unsupported physical treatment')
    END SUBROUTINE eval_mass_exchange_entrainment_semi_implicit

    !> @brief Reject unsupported implicit treatment for ambient-air entrainment.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] state Local real physical state.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_entrainment_implicit(state,r)

        TYPE(complex_physical_state_type), INTENT(IN) :: state
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('entrainment: unsupported physical treatment')
    END SUBROUTINE eval_mass_exchange_entrainment_implicit

    !> @brief Reject unsupported analytic treatment for ambient-air entrainment.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] state Local real physical state.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE solve_mass_exchange_entrainment_analytic(state,r)

        TYPE(physical_state_type), INTENT(IN) :: state
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('entrainment: unsupported physical treatment')
    END SUBROUTINE solve_mass_exchange_entrainment_analytic

    !> @brief Switch air entrainment off in dry or non-buoyant conditions.
    !> @param[in] r_h Mobile flow thickness [m].
    !> @param[in] r_Ri Mixture Richardson number [dimensionless].
    !> @param[in] mod_vel2 Squared slope-corrected speed [m2/s2].
    !> @return Physical rate from the indicated constitutive closure.
    REAL(wp) FUNCTION eval_entrainment_volume_kernel(r_h, r_Ri, mod_vel2)
        REAL(wp), INTENT(IN) :: r_h, r_Ri, mod_vel2

        IF ( entrainment_flag .AND.  ( r_h .GT. 0.0_wp ) .AND.                                     &
                  ( r_Ri .GT. 0.0_wp ) ) THEN

              eval_entrainment_volume_kernel = eval_air_entrainment_rate(r_Ri, mod_vel2)

        ELSE

              eval_entrainment_volume_kernel = 0.0_wp

        END IF

    END FUNCTION eval_entrainment_volume_kernel

    !> @brief Evaluate the Richardson-number-dependent ambient-air entrainment law.
    !> @param[in] richardson_number Richardson number controlling entrainment [dimensionless].
    !> @param[in] speed_squared Squared three-dimensional flow speed [m2/s2].
    !> @return Physical rate from the indicated constitutive closure.
    PURE FUNCTION eval_air_entrainment_rate(richardson_number, speed_squared)                      &
              RESULT(entrainment_rate)

        REAL(wp), INTENT(IN) :: richardson_number, speed_squared
        REAL(wp) :: entrainment_coefficient, entrainment_rate

        entrainment_coefficient = 0.075_wp / SQRT(1.0_wp + 718.0_wp *                              &
                  richardson_number**2.4_wp)
        entrainment_rate = entrainment_coefficient * SQRT(speed_squared)

    END FUNCTION eval_air_entrainment_rate

END MODULE mass_exchange_entrainment_2d
