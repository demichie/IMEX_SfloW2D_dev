!> @brief Dispatch thermal source models through a common four-treatment interface.
!> Add future models in the SELECT CASE statements below.
MODULE thermal_laws_2d
    USE parameters_2d, ONLY : wp
    USE diagnostics_2d, ONLY : fatal_error
    USE source_catalog_2d, ONLY : SOURCE_ID_THERMAL_COOLING
    USE state_types, ONLY : physical_state_real, physical_state_complex
    USE thermal_cooling_2d, ONLY : eval_thermal_cooling_explicit,                                                      &
            eval_thermal_cooling_semi_implicit, eval_thermal_cooling_implicit,                                         &
            eval_thermal_cooling_implicit_real, eval_thermal_cooling_analytic
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_thermal_explicit, eval_thermal_semi_implicit
    PUBLIC :: eval_thermal_implicit, eval_thermal_implicit_real, eval_thermal_analytic
CONTAINS

    !> @brief Dispatch and sum explicit thermal sources.
    !> @details The dispatcher retains each of the four numerical interfaces;
    !>          unsupported model treatments fail through the model implementation.
    !> @param[in] models Canonical thermal source names.
    !> @param[in] state Real physical cell state.
    !> @param[in] radiative Radiative heat-exchange coefficient.
    !> @param[in] convective Convective heat-exchange coefficient.
    !> @param[in] ambient Ambient temperature [K].
    !> @param[out] rate Assembled real thermal source rate.
    SUBROUTINE eval_thermal_explicit(model_id,state,radiative,convective,ambient,rate)
        INTEGER, INTENT(IN) :: model_id
        TYPE(physical_state_real), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: radiative,convective,ambient
        REAL(wp), INTENT(OUT) :: rate
        REAL(wp) :: process_rate
        rate=0.0_wp

        SELECT CASE (model_id)
        CASE (SOURCE_ID_THERMAL_COOLING)
            CALL eval_thermal_cooling_explicit(state,radiative,convective,ambient,process_rate)
        CASE DEFAULT
            CALL fatal_error('Unknown explicit thermal model ID')

        END SELECT
        rate=process_rate
    END SUBROUTINE eval_thermal_explicit


    !> @brief Dispatch and sum semi-implicit thermal contributions.
    !> @details The dispatcher retains each of the four numerical interfaces;
    !>          unsupported model treatments fail through the model implementation.
    !> @param[in] models Canonical thermal source names.
    !> @param[in] state Real physical cell state.
    !> @param[out] rate Assembled real thermal source rate.
    SUBROUTINE eval_thermal_semi_implicit(models,state,rate)
        CHARACTER(*), INTENT(IN) :: models(:)
        TYPE(physical_state_real), INTENT(IN) :: state
        REAL(wp), INTENT(OUT) :: rate
        REAL(wp) :: process_rate
        INTEGER :: i
        rate=0.0_wp

        DO i=1,SIZE(models)

        SELECT CASE (TRIM(models(i)))
        CASE ('thermal_cooling')
            CALL eval_thermal_cooling_semi_implicit(state,process_rate)
        CASE DEFAULT
            CALL fatal_error('Unknown semi-implicit thermal model: '//TRIM(models(i)))

        END SELECT

        rate=rate+process_rate

        END DO
    END SUBROUTINE eval_thermal_semi_implicit


    ! Complex-step source evaluation is side-effect free with respect to
    ! substrate, topography and other shared integration state.
    !> @brief Dispatch and sum complex thermal contributions for Newton.
    !> @details The dispatcher retains each of the four numerical interfaces;
    !>          unsupported model treatments fail through the model implementation.
    !> @param[in] models Canonical thermal source names.
    !> @param[in] state Complex physical cell state.
    !> @param[in] radiative Radiative heat-exchange coefficient.
    !> @param[in] convective Convective heat-exchange coefficient.
    !> @param[in] ambient Ambient temperature [K].
    !> @param[out] rate Assembled complex thermal source rate.
    SUBROUTINE eval_thermal_implicit(model_id,state,radiative,convective,ambient,rate)
        INTEGER, INTENT(IN) :: model_id
        TYPE(physical_state_complex), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: radiative,convective,ambient
        COMPLEX(wp), INTENT(OUT) :: rate
        COMPLEX(wp) :: process_rate
        rate=CMPLX(0.0_wp,0.0_wp,wp)

        SELECT CASE (model_id)
        CASE (SOURCE_ID_THERMAL_COOLING)
            CALL eval_thermal_cooling_implicit(state,radiative,convective,ambient,process_rate)
        CASE DEFAULT
            CALL fatal_error('Unknown implicit thermal model ID')

        END SELECT
        rate=process_rate
    END SUBROUTINE eval_thermal_implicit


    !> REAL implicit source for Newton residual evaluation.
    SUBROUTINE eval_thermal_implicit_real(model_id,state,radiative,convective,ambient,rate)
        INTEGER, INTENT(IN) :: model_id
        TYPE(physical_state_real), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: radiative,convective,ambient
        REAL(wp), INTENT(OUT) :: rate
        SELECT CASE(model_id)
        CASE(SOURCE_ID_THERMAL_COOLING)
            CALL eval_thermal_cooling_implicit_real(state,radiative,convective,ambient,rate)
        CASE DEFAULT
            CALL fatal_error('Unknown REAL implicit thermal model ID')
        END SELECT
    END SUBROUTINE eval_thermal_implicit_real

    !> @brief Advance thermal temperature analytically over the given interval.
    !> @details The dispatcher retains each of the four numerical interfaces;
    !>          unsupported model treatments fail through the model implementation.
    !> @param[in] models Canonical thermal source names.
    !> @param[in] initial_temperature Frozen mixture temperature; no full physical state is required.
    !> @param[in] interval Analytic substep duration [s].
    !> @param[in] heat_capacity Thermal capacity of the frozen mixture.
    !> @param[in] radiative Radiative heat-exchange coefficient.
    !> @param[in] convective Convective heat-exchange coefficient.
    !> @param[in] ambient Ambient temperature [K].
    !> @param[out] temperature_after Temperature after the thermal substep [K].
    SUBROUTINE eval_thermal_analytic(models,initial_temperature,interval,heat_capacity,                               &
            radiative,convective,ambient,temperature_after)
        CHARACTER(*), INTENT(IN) :: models(:)
        REAL(wp), INTENT(IN) :: initial_temperature
        REAL(wp), INTENT(IN) :: interval,heat_capacity,radiative,convective,ambient
        REAL(wp), INTENT(OUT) :: temperature_after
        REAL(wp) :: temperature_before
        INTEGER :: i
        temperature_after=initial_temperature

        DO i=1,SIZE(models)
            temperature_before=temperature_after

        SELECT CASE (TRIM(models(i)))
        CASE ('thermal_cooling')
            CALL eval_thermal_cooling_analytic(temperature_before,interval,heat_capacity,                              &
                radiative,convective,ambient,temperature_after)
        CASE DEFAULT
            CALL fatal_error('Unknown analytic thermal model: '//TRIM(models(i)))

        END SELECT

        END DO
    END SUBROUTINE eval_thermal_analytic

END MODULE thermal_laws_2d
