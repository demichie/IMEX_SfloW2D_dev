!> @brief Evaluate carrier darcy loss process contributions to mass, energy and bed evolution.
!> @details Explicit treatment is supported; other treatment interfaces reject
!> unsupported configurations. Shared resource budgets belong to the coordinator.
MODULE mass_exchange_carrier_darcy_loss_2d
    USE parameters_2d, ONLY : wp, n_add_gas, pore_pressure_flag, gas_loss_flag, gas_flag, &
        mass_exchange_rates_type, reset_mass_exchange_rates
    USE constitutive_parameters_2d, ONLY : maximum_solid_packing
    USE pore_pressure_laws_2d, ONLY : eval_pore_pressure_gas_loss
    USE diagnostics_2d, ONLY : fatal_error
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
    USE state_conversion_2d, ONLY : mixt_var_state

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: eval_mass_exchange_carrier_darcy_loss_explicit
    PUBLIC :: eval_mass_exchange_carrier_darcy_loss_semi_implicit
    PUBLIC :: eval_mass_exchange_carrier_darcy_loss_implicit
    PUBLIC :: solve_mass_exchange_carrier_darcy_loss_analytic

CONTAINS

    !> @brief Evaluate explicit pressure-driven carrier gas loss.
    !> @param[in] state Local real physical state.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_carrier_darcy_loss_explicit(state,r)

        TYPE(physical_state_type), INTENT(IN) :: state
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
        REAL(wp) :: volume_loss,alpha,rho_c,cp_c
        REAL(wp) :: Ri,rhom,redgrav,cpm,gas_phase_mass_fraction

        CALL reset_mass_exchange_rates(r)
        IF (.NOT.(gas_flag .AND. pore_pressure_flag .AND. gas_loss_flag)) RETURN
        alpha=SUM(state%alphas)
        CALL mixt_var_state(state,Ri,rhom,rho_c,redgrav,cp_c,cpm)

        ! Recover carrier gas mass fraction before partitioning withdrawal among gas species.
        gas_phase_mass_fraction=MAX(0.0_wp,1.0_wp-                                                 &
                  SUM(state%xs))
        volume_loss=0.0_wp
        IF (gas_flag .AND. pore_pressure_flag .AND. gas_loss_flag .AND.                            &
                  alpha<maximum_solid_packing)                                                     &
                  CALL eval_pore_pressure_gas_loss(state%h,state%T,state%alphas,alpha,rho_c, &
                      state%exc_pore_pres,volume_loss)
        r%carrier_withdrawal_volume_rate=MAX(0.0_wp,volume_loss)
        r%mixture_mass_rate=-rho_c*r%carrier_withdrawal_volume_rate
        r%energy_rate=state%T*cp_c*r%mixture_mass_rate
        IF (n_add_gas>0 .AND. gas_phase_mass_fraction>0.0_wp)                                      &
                  r%add_gas_mass_rate=r%mixture_mass_rate*                                         &
                  state%xg/gas_phase_mass_fraction
    END SUBROUTINE eval_mass_exchange_carrier_darcy_loss_explicit

    !> @brief Reject unsupported semi-implicit treatment for carrier darcy loss.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] state Local real physical state.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_carrier_darcy_loss_semi_implicit(state,r)

        TYPE(physical_state_type), INTENT(IN) :: state
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('carrier_darcy_loss: unsupported semi_implicit method')
    END SUBROUTINE eval_mass_exchange_carrier_darcy_loss_semi_implicit

    !> @brief Reject unsupported implicit treatment for pressure-driven carrier gas loss.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] state Local real physical state.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_carrier_darcy_loss_implicit(state,r)

        TYPE(complex_physical_state_type), INTENT(IN) :: state
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('carrier_darcy_loss: unsupported implicit method')
    END SUBROUTINE eval_mass_exchange_carrier_darcy_loss_implicit

    !> @brief Reject unsupported analytic treatment for pressure-driven carrier gas loss.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] state Local real physical state.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE solve_mass_exchange_carrier_darcy_loss_analytic(state,r)

        TYPE(physical_state_type), INTENT(IN) :: state
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('carrier_darcy_loss: unsupported analytic method')
    END SUBROUTINE solve_mass_exchange_carrier_darcy_loss_analytic

END MODULE mass_exchange_carrier_darcy_loss_2d
