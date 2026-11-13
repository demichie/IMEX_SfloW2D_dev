!> @brief Evaluate carrier liquid loss process contributions to mass, energy and bed evolution.
!> @details Explicit treatment is supported; other treatment interfaces reject
!> unsupported configurations. Shared resource budgets belong to the coordinator.
MODULE mass_exchange_carrier_liquid_loss_2d
    USE parameters_2d, ONLY : wp, liquid_flag, mass_exchange_rates_type,                                               &
            reset_mass_exchange_rates
    USE mass_exchange_parameters_2d, ONLY : mass_exchange_config_type
    USE diagnostics_2d, ONLY : fatal_error
    USE state_types, ONLY : physical_state_real, physical_state_complex

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: eval_mass_exchange_carrier_liquid_loss_explicit
    PUBLIC :: eval_mass_exchange_carrier_liquid_loss_semi_implicit
    PUBLIC :: eval_mass_exchange_carrier_liquid_loss_implicit
    PUBLIC :: solve_mass_exchange_carrier_liquid_loss_analytic

CONTAINS

        !> @brief Evaluate explicit prescribed liquid-carrier loss.
        !> @param[in] state Local real physical state.
        !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
        !> @param[in] mass_config Immutable layer carrier-liquid-loss configuration.
    SUBROUTINE eval_mass_exchange_carrier_liquid_loss_explicit(mass_config,state,r)

            TYPE(mass_exchange_config_type), INTENT(IN) :: mass_config
            TYPE(physical_state_real), INTENT(IN) :: state
            TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

            CALL reset_mass_exchange_rates(r)

            IF (liquid_flag) THEN
                r%carrier_withdrawal_volume_rate=MAX(0.0_wp,mass_config%carrier_liquid_loss%loss_rate)

                ! The liquid is the only carrier in this configuration.
                r%mixture_mass_rate=-state%rho_c*r%carrier_withdrawal_volume_rate
                r%energy_rate=state%T*state%cp_c*r%mixture_mass_rate

            END IF

    END SUBROUTINE eval_mass_exchange_carrier_liquid_loss_explicit

        !> @brief Reject unsupported semi-implicit treatment for carrier liquid loss.
        !> @details This entry point raises fatal_error; the rates are reset first.
        !> @param[in] state Local real physical state.
        !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_carrier_liquid_loss_semi_implicit(state,r)

            TYPE(physical_state_real), INTENT(IN) :: state
            TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

            CALL reset_mass_exchange_rates(r)
            CALL fatal_error('carrier_liquid_loss: unsupported semi_implicit method')
    END SUBROUTINE eval_mass_exchange_carrier_liquid_loss_semi_implicit

        !> @brief Reject unsupported implicit treatment for prescribed liquid-carrier loss.
        !> @details This entry point raises fatal_error; the rates are reset first.
        !> @param[in] state Local real physical state.
        !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_carrier_liquid_loss_implicit(state,r)

            TYPE(physical_state_complex), INTENT(IN) :: state
            TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

            CALL reset_mass_exchange_rates(r)
            CALL fatal_error('carrier_liquid_loss: unsupported implicit method')
    END SUBROUTINE eval_mass_exchange_carrier_liquid_loss_implicit

        !> @brief Reject unsupported analytic treatment for prescribed liquid-carrier loss.
        !> @details This entry point raises fatal_error; the rates are reset first.
        !> @param[in] state Local real physical state.
        !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE solve_mass_exchange_carrier_liquid_loss_analytic(state,r)

            TYPE(physical_state_real), INTENT(IN) :: state
            TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

            CALL reset_mass_exchange_rates(r)
            CALL fatal_error('carrier_liquid_loss: unsupported analytic method')
    END SUBROUTINE solve_mass_exchange_carrier_liquid_loss_analytic

END MODULE mass_exchange_carrier_liquid_loss_2d
