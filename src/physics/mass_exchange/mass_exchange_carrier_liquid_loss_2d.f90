!> @brief mass exchange carrier liquid loss 2d: physical mass-exchange components and rate evaluation.
!> @details Explicit physics is evaluated locally; unsupported treatment entry points
!> fail explicitly. The coordinator alone enforces shared carrier inventories.
MODULE mass_exchange_carrier_liquid_loss_2d
    USE parameters_2d, ONLY : wp,n_vars,liquid_flag,gas_flag,mass_exchange_rates_type,reset_mass_exchange_rates
    USE constitutive_parameters_2d, ONLY : loss_rate, rho_l, sp_heat_l
    USE diagnostics_2d, ONLY : fatal_error
    USE state_conversion_2d, ONLY : mixt_var
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_mass_exchange_carrier_liquid_loss_explicit
    PUBLIC :: eval_mass_exchange_carrier_liquid_loss_semi_implicit
    PUBLIC :: eval_mass_exchange_carrier_liquid_loss_implicit
    PUBLIC :: solve_mass_exchange_carrier_liquid_loss_analytic
CONTAINS
    !> @brief Evaluate explicit prescribed liquid-carrier loss.
    !> @param[in] qp Local primitive flow state; component order follows model configuration.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_carrier_liquid_loss_explicit(qp,r)
        REAL(wp), INTENT(IN) :: qp(n_vars+2)
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
        REAL(wp) :: rho_c,T,cp_c,Ri,rhom,redgrav,cpm
        CALL reset_mass_exchange_rates(r)
        IF (liquid_flag .AND. ALLOCATED(loss_rate)) THEN
            T=qp(4)
            CALL mixt_var(qp,Ri,rhom,rho_c,redgrav,cp_c,cpm)
            r%liquid_withdrawal_volume_rate=MAX(0.0_wp,loss_rate)
            ! In mixed-phase flow the removed liquid has its own density and enthalpy.
            IF (gas_flag) THEN
                r%mixture_mass_rate=-rho_l*MAX(0.0_wp,loss_rate)
                r%liquid_mass_rate=r%mixture_mass_rate
                r%energy_rate=T*sp_heat_l*r%mixture_mass_rate
            ELSE
                r%mixture_mass_rate=-rho_c*MAX(0.0_wp,loss_rate)
                r%energy_rate=T*cp_c*r%mixture_mass_rate
            END IF
        END IF
    END SUBROUTINE eval_mass_exchange_carrier_liquid_loss_explicit
    !> @brief Reject unsupported implicit treatment for carrier_liquid_loss_semi.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] qp Local primitive flow state; component order follows model configuration.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_carrier_liquid_loss_semi_implicit(qp,r)
        REAL(wp), INTENT(IN) :: qp(n_vars+2)
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('carrier_liquid_loss: unsupported semi_implicit method')
    END SUBROUTINE eval_mass_exchange_carrier_liquid_loss_semi_implicit
    !> @brief Reject unsupported implicit treatment for prescribed liquid-carrier loss.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] qp Local primitive flow state; component order follows model configuration.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_carrier_liquid_loss_implicit(qp,r)
        REAL(wp), INTENT(IN) :: qp(n_vars+2)
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('carrier_liquid_loss: unsupported implicit method')
    END SUBROUTINE eval_mass_exchange_carrier_liquid_loss_implicit
    !> @brief Reject unsupported analytic treatment for prescribed liquid-carrier loss.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] qp Local primitive flow state; component order follows model configuration.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE solve_mass_exchange_carrier_liquid_loss_analytic(qp,r)
        REAL(wp), INTENT(IN) :: qp(n_vars+2)
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('carrier_liquid_loss: unsupported analytic method')
    END SUBROUTINE solve_mass_exchange_carrier_liquid_loss_analytic

END MODULE mass_exchange_carrier_liquid_loss_2d
