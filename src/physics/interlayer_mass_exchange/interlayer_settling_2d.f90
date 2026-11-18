!> @brief Per-species L2 to L1 settling rate from the shared deposition law.
!> @details The same settling velocity, rheology fluid properties and hindrance
!>          are used for basal deposition and interlayer settling. No bed
!>          deposition or carrier creation occurs in this kernel.
MODULE interlayer_settling_2d
    USE parameters_2d, ONLY : wp,n_solid
    USE diagnostics_2d, ONLY : fatal_error
    USE state_types, ONLY : physical_state_real,physical_state_complex
    USE rheology_parameters_2d, ONLY : rheology_config_type
    USE mass_exchange_parameters_2d, ONLY : mass_exchange_config_type
    USE interlayer_mass_exchange_parameters_2d, ONLY : interlayer_mass_exchange_config_type
    USE mass_exchange_deposition_2d, ONLY : eval_deposition_volume_kernel
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_interlayer_settling_explicit
    ! Four physical treatment entry points; unsupported methods fail explicitly.
    PUBLIC :: eval_interlayer_settling_semi_implicit
    PUBLIC :: eval_interlayer_settling_implicit
    PUBLIC :: solve_interlayer_settling_analytic

CONTAINS
    !> @brief Evaluate per-class positive solid volume fluxes from L2 to L1.
    !> @param[in] config Pair-wide interlayer activation and options.
    !> @param[in] rheology_l2 Donor L2 rheology for the existing settling law.
    !> @param[in] mass_l2 Donor L2 carrier/drag properties, not its activation.
    !> @param[in] state_l2 Local primitive-derived L2 physical state.
    !> @param[in] dt Common physical timestep [s].
    !> @param[out] settling Per-class downward solid volume flux [m/s].
    SUBROUTINE eval_interlayer_settling_explicit(config,rheology_l2,mass_l2,state_l2,dt,settling)
        TYPE(interlayer_mass_exchange_config_type), INTENT(IN) :: config
        TYPE(rheology_config_type), INTENT(IN) :: rheology_l2
        TYPE(mass_exchange_config_type), INTENT(IN) :: mass_l2
        TYPE(physical_state_real), INTENT(IN) :: state_l2
        REAL(wp), INTENT(IN) :: dt
        REAL(wp), INTENT(OUT) :: settling(n_solid)

        settling=0.0_wp
        IF (.NOT.config%settling_flag) RETURN
        IF (dt<=0.0_wp .OR. state_l2%h<=0.0_wp) RETURN
        CALL eval_deposition_volume_kernel(rheology_l2,mass_l2,state_l2%h,state_l2%T, &
            state_l2%alphas,SUM(state_l2%alphas),state_l2%rho_m,state_l2%rho_c, &
            dt,settling,enable_settling=.TRUE.)
    END SUBROUTINE eval_interlayer_settling_explicit
    !> @brief Reject unsupported semi-implicit treatment for interlayer settling.
    !> @details No candidate state is committed; report the unsupported method.
    SUBROUTINE eval_interlayer_settling_semi_implicit(config,rheology_l2,mass_l2,state_l2,dt,settling)
        TYPE(interlayer_mass_exchange_config_type), INTENT(IN) :: config
        TYPE(rheology_config_type), INTENT(IN) :: rheology_l2
        TYPE(mass_exchange_config_type), INTENT(IN) :: mass_l2
        TYPE(physical_state_real), INTENT(IN) :: state_l2
        REAL(wp), INTENT(IN) :: dt
        REAL(wp), INTENT(OUT) :: settling(n_solid)

        settling=0.0_wp
        CALL fatal_error('interlayer_settling: unsupported semi-implicit physical treatment')
    END SUBROUTINE eval_interlayer_settling_semi_implicit

    !> @brief Reject unsupported implicit treatment for interlayer settling.
    !> @details No candidate state is committed; report the unsupported method.
    SUBROUTINE eval_interlayer_settling_implicit(config,rheology_l2,mass_l2,state_l2,dt,settling)
        TYPE(interlayer_mass_exchange_config_type), INTENT(IN) :: config
        TYPE(rheology_config_type), INTENT(IN) :: rheology_l2
        TYPE(mass_exchange_config_type), INTENT(IN) :: mass_l2
        TYPE(physical_state_complex), INTENT(IN) :: state_l2
        REAL(wp), INTENT(IN) :: dt
        REAL(wp), INTENT(OUT) :: settling(n_solid)

        settling=0.0_wp
        CALL fatal_error('interlayer_settling: unsupported implicit physical treatment')
    END SUBROUTINE eval_interlayer_settling_implicit

    !> @brief Reject unsupported analytic treatment for interlayer settling.
    !> @details No candidate state is committed; report the unsupported method.
    SUBROUTINE solve_interlayer_settling_analytic(config,rheology_l2,mass_l2,state_l2,dt,settling)
        TYPE(interlayer_mass_exchange_config_type), INTENT(IN) :: config
        TYPE(rheology_config_type), INTENT(IN) :: rheology_l2
        TYPE(mass_exchange_config_type), INTENT(IN) :: mass_l2
        TYPE(physical_state_real), INTENT(IN) :: state_l2
        REAL(wp), INTENT(IN) :: dt
        REAL(wp), INTENT(OUT) :: settling(n_solid)

        settling=0.0_wp
        CALL fatal_error('interlayer_settling: unsupported analytic physical treatment')
    END SUBROUTINE solve_interlayer_settling_analytic

END MODULE interlayer_settling_2d
