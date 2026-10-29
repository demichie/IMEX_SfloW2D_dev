!> @brief Evaluate lava deposition process contributions to mass, energy and bed evolution.
!> @details Explicit treatment is supported; other treatment interfaces reject
!> unsupported configurations. Shared resource budgets belong to the coordinator.
MODULE mass_exchange_lava_deposition_2d
    USE parameters_2d, ONLY : wp, mass_exchange_rates_type, reset_mass_exchange_rates
    USE diagnostics_2d, ONLY : fatal_error
    USE rheology_lava_2d, ONLY : lava_kinematic_viscosity
    USE constitutive_parameters_2d, ONLY : rho_l, sp_heat_l

    IMPLICIT NONE
    PRIVATE
    LOGICAL, PUBLIC :: lava_deposition_flag=.FALSE.
    CHARACTER(LEN=16), PUBLIC :: lava_deposition_model='CONSTANT'
    REAL(wp), PUBLIC :: lava_deposition_rate=0.0_wp
    REAL(wp), PUBLIC :: lava_deposition_visc_ref=1.0E-2_wp
    REAL(wp), PUBLIC :: lava_deposition_exponent=1.0_wp
    REAL(wp), PUBLIC :: lava_deposition_temp_ref=900.0_wp
    REAL(wp), PUBLIC :: lava_deposition_temp_width=100.0_wp
    PUBLIC :: eval_mass_exchange_lava_deposition_explicit
    PUBLIC :: eval_mass_exchange_lava_deposition_semi_implicit
    PUBLIC :: eval_mass_exchange_lava_deposition_implicit
    PUBLIC :: solve_mass_exchange_lava_deposition_analytic

CONTAINS

    !> @brief Evaluate explicit immobilization of mobile lava.
    !> @param[in] depth Mobile lava thickness [m].
    !> @param[in] T Mobile flow temperature [K].
    !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_lava_deposition_explicit(depth,T,dt,r)
        REAL(wp), INTENT(IN) :: depth,T,dt
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL eval_lava_deposition_volume_kernel(depth,T,dt, &
                  r%continuous_phase_deposition_rate)
        r%mixture_mass_rate=-rho_l*r%continuous_phase_deposition_rate
        r%energy_rate=sp_heat_l*T*r%mixture_mass_rate
        r%lava_deposit_rate=r%continuous_phase_deposition_rate
        r%dB_dt=r%continuous_phase_deposition_rate
    END SUBROUTINE

    !> @brief Reject unsupported semi-implicit treatment for lava deposition.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] depth Mobile lava thickness [m].
    !> @param[in] T Mobile flow temperature [K].
    !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_lava_deposition_semi_implicit(depth,T,dt,r)
        REAL(wp), INTENT(IN) :: depth,T,dt
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('lava_deposition: unsupported physical treatment')
    END SUBROUTINE eval_mass_exchange_lava_deposition_semi_implicit

    !> @brief Reject unsupported implicit treatment for immobilization of mobile lava.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] depth Mobile lava thickness [m].
    !> @param[in] T Mobile flow temperature [K].
    !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_lava_deposition_implicit(depth,T,dt,r)
        REAL(wp), INTENT(IN) :: depth,T,dt
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('lava_deposition: unsupported physical treatment')
    END SUBROUTINE eval_mass_exchange_lava_deposition_implicit

    !> @brief Reject unsupported analytic treatment for immobilization of mobile lava.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] depth Mobile lava thickness [m].
    !> @param[in] T Mobile flow temperature [K].
    !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE solve_mass_exchange_lava_deposition_analytic(depth,T,dt,r)
        REAL(wp), INTENT(IN) :: depth,T,dt
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('lava_deposition: unsupported physical treatment')
    END SUBROUTINE solve_mass_exchange_lava_deposition_analytic

    !> @brief Evaluate the configured temperature- or viscosity-dependent lava immobilization rate.
    !> @param[in] temperature Material temperature [K].
    !> @param[in] viscosity Lava kinematic viscosity [m2/s].
    !> @return Physical rate from the indicated constitutive closure.
    PURE REAL(wp) FUNCTION eval_lava_deposition_rate(temperature,viscosity) RESULT(rate)
        REAL(wp), INTENT(IN) :: temperature,viscosity
        REAL(wp) :: factor

        factor=1.0_wp

        ! Temperature and viscosity choices modulate the prescribed immobilization flux.
        SELECT CASE(TRIM(lava_deposition_model))
        CASE('TEMPERATURE')
              factor=MIN(1.0_wp,MAX(0.0_wp,(lava_deposition_temp_ref-temperature)/ &
                        lava_deposition_temp_width))
        CASE('VISCOSITY')
              factor=MIN(1.0_wp,(MAX(0.0_wp,viscosity)/lava_deposition_visc_ref)** &
                        lava_deposition_exponent)
        END SELECT
        rate=MAX(0.0_wp,lava_deposition_rate)*factor
    END FUNCTION eval_lava_deposition_rate

    !> @brief Cap lava immobilization by the mobile layer depth available over the time step.
    !> @param[in] depth Mobile lava thickness [m].
    !> @param[in] temperature Material temperature [K].
    !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
    !> @param[out] rate Requested immobilization thickness rate [m/s].
    SUBROUTINE eval_lava_deposition_volume_kernel(depth,temperature,dt,rate)
        REAL(wp), INTENT(IN) :: depth,temperature,dt
        REAL(wp), INTENT(OUT) :: rate

        rate=0.0_wp
        IF (depth<=0.0_wp) RETURN
        IF (dt<=0.0_wp) CALL fatal_error('lava deposition requires positive dt')
        rate=MIN(depth/dt,eval_lava_deposition_rate(temperature, &
                  lava_kinematic_viscosity(temperature)))
    END SUBROUTINE eval_lava_deposition_volume_kernel

END MODULE mass_exchange_lava_deposition_2d
