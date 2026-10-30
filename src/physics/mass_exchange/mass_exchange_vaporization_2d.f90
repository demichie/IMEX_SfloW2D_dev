!> @brief Evaluate vaporization process contributions to mass, energy and bed evolution.
!> @details Explicit treatment is supported; other treatment interfaces reject
!> unsupported configurations. Shared resource budgets belong to the coordinator.
MODULE mass_exchange_vaporization_2d
    USE parameters_2d, ONLY : wp, n_solid, n_add_gas, gas_flag,                                    &
              liquid_vaporization_flag, mass_exchange_rates_type, reset_mass_exchange_rates
    USE diagnostics_2d, ONLY : fatal_error
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
    USE constitutive_parameters_2d, ONLY : rho_s, sp_heat_s, gamma_steam, sp_heat_g
    USE thermal_laws_2d, ONLY : eval_vaporization_mass_rate, eval_vaporization_energy_rate

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: eval_mass_exchange_vaporization_explicit
    PUBLIC :: eval_mass_exchange_vaporization_semi_implicit
    PUBLIC :: eval_mass_exchange_vaporization_implicit
    PUBLIC :: solve_mass_exchange_vaporization_analytic
    ! Four-entry physical treatment contract; only EXPLICIT is supported.

CONTAINS

    !> @brief Evaluate explicit steam generation from deposited hot solids.
    !> @param[in] B_zone Substrate zone identifier controlling vaporization eligibility.
    !> @param[in] T Mobile flow temperature [K].
    !> @param[in] deposition Actual post-limiting solid deposition rates by species [m/s].
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_vaporization_explicit(state,B_zone,deposition,r)

        TYPE(physical_state_type), INTENT(IN) :: state
        INTEGER, INTENT(IN) :: B_zone
        REAL(wp), INTENT(IN) :: deposition(n_solid)
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        IF (.NOT.(gas_flag .AND. liquid_vaporization_flag .AND. B_zone/=0)) RETURN
        IF (n_add_gas<1) CALL fatal_error('vaporization requires additional gas species')

        ! Use already-limited deposition: generated vapor must reflect deposited solids.
        r%mixture_mass_rate=eval_vaporization_mass_rate(state%T,rho_s,sp_heat_s,                         &
                  deposition,gamma_steam)
        r%add_gas_mass_rate(1)=r%mixture_mass_rate
        r%energy_rate=eval_vaporization_energy_rate(r%mixture_mass_rate,sp_heat_g(1))
    END SUBROUTINE

    !> @brief Reject unsupported semi-implicit treatment for vaporization.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] B_zone Substrate zone identifier controlling vaporization eligibility.
    !> @param[in] T Mobile flow temperature [K].
    !> @param[in] deposition Actual post-limiting solid deposition rates by species [m/s].
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_vaporization_semi_implicit(state,B_zone,deposition,r)

        TYPE(physical_state_type), INTENT(IN) :: state
        INTEGER, INTENT(IN) :: B_zone
        REAL(wp), INTENT(IN) :: deposition(n_solid)
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('vaporization: unsupported physical treatment')
    END SUBROUTINE eval_mass_exchange_vaporization_semi_implicit

    !> @brief Reject unsupported implicit treatment for steam generation from deposited hot solids.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] B_zone Substrate zone identifier controlling vaporization eligibility.
    !> @param[in] T Mobile flow temperature [K].
    !> @param[in] deposition Actual post-limiting solid deposition rates by species [m/s].
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_vaporization_implicit(state,B_zone,deposition,r)

        TYPE(complex_physical_state_type), INTENT(IN) :: state
        INTEGER, INTENT(IN) :: B_zone
        REAL(wp), INTENT(IN) :: deposition(n_solid)
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('vaporization: unsupported physical treatment')
    END SUBROUTINE eval_mass_exchange_vaporization_implicit

    !> @brief Reject unsupported analytic treatment for steam generation from deposited hot solids.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] B_zone Substrate zone identifier controlling vaporization eligibility.
    !> @param[in] T Mobile flow temperature [K].
    !> @param[in] deposition Actual post-limiting solid deposition rates by species [m/s].
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE solve_mass_exchange_vaporization_analytic(state,B_zone,deposition,r)

        TYPE(physical_state_type), INTENT(IN) :: state
        INTEGER, INTENT(IN) :: B_zone
        REAL(wp), INTENT(IN) :: deposition(n_solid)
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('vaporization: unsupported physical treatment')
    END SUBROUTINE solve_mass_exchange_vaporization_analytic

END MODULE mass_exchange_vaporization_2d
