!> @brief Evaluate carrier liquid loss process contributions to mass, energy and bed evolution.
!> @details Explicit treatment is supported; other treatment interfaces reject
!> unsupported configurations. Shared resource budgets belong to the coordinator.
MODULE mass_exchange_carrier_liquid_loss_2d
    USE parameters_2d, ONLY : idx_T
    USE parameters_2d, ONLY : wp,n_vars,liquid_flag,mass_exchange_rates_type,                      &
        reset_mass_exchange_rates
    USE constitutive_parameters_2d, ONLY : loss_rate
    USE diagnostics_2d, ONLY : fatal_error
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
    USE state_conversion_2d, ONLY : mixt_var

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: eval_mass_exchange_carrier_liquid_loss_explicit
    PUBLIC :: eval_mass_exchange_carrier_liquid_loss_semi_implicit
    PUBLIC :: eval_mass_exchange_carrier_liquid_loss_implicit
    PUBLIC :: solve_mass_exchange_carrier_liquid_loss_analytic

CONTAINS

    !> @brief Evaluate explicit prescribed liquid-carrier loss.
    !> @param[in] state Local real physical state.
    !> @param[in] qp Primitive vector retained for mixture thermodynamic and species closures.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_carrier_liquid_loss_explicit(state,qp,r)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: qp(n_vars+2)
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
        REAL(wp) :: rho_c,T,cp_c,Ri,rhom,redgrav,cpm

        CALL reset_mass_exchange_rates(r)

        IF (liquid_flag .AND. ALLOCATED(loss_rate)) THEN
            T=state%T
            CALL mixt_var(qp,Ri,rhom,rho_c,redgrav,cp_c,cpm)
            r%carrier_withdrawal_volume_rate=MAX(0.0_wp,loss_rate)

            ! The liquid is the only carrier in this configuration.
            r%mixture_mass_rate=-rho_c*r%carrier_withdrawal_volume_rate
            r%energy_rate=T*cp_c*r%mixture_mass_rate
        END IF

    END SUBROUTINE eval_mass_exchange_carrier_liquid_loss_explicit

    !> @brief Reject unsupported semi-implicit treatment for carrier liquid loss.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] state Local real physical state.
    !> @param[in] qp Primitive vector retained for mixture thermodynamic and species closures.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_carrier_liquid_loss_semi_implicit(state,qp,r)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: qp(n_vars+2)
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('carrier_liquid_loss: unsupported semi_implicit method')
    END SUBROUTINE eval_mass_exchange_carrier_liquid_loss_semi_implicit

    !> @brief Reject unsupported implicit treatment for prescribed liquid-carrier loss.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] state Local real physical state.
    !> @param[in] qp Primitive vector retained for mixture thermodynamic and species closures.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_carrier_liquid_loss_implicit(state,qp,r)

        TYPE(complex_physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: qp(n_vars+2)
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('carrier_liquid_loss: unsupported implicit method')
    END SUBROUTINE eval_mass_exchange_carrier_liquid_loss_implicit

    !> @brief Reject unsupported analytic treatment for prescribed liquid-carrier loss.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] state Local real physical state.
    !> @param[in] qp Primitive vector retained for mixture thermodynamic and species closures.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE solve_mass_exchange_carrier_liquid_loss_analytic(state,qp,r)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: qp(n_vars+2)
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('carrier_liquid_loss: unsupported analytic method')
    END SUBROUTINE solve_mass_exchange_carrier_liquid_loss_analytic

END MODULE mass_exchange_carrier_liquid_loss_2d
