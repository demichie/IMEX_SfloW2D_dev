!********************************************************************************
!> \brief Rheology dispatch, model naming and friction-only momentum correction.
!> \details Dispatches common momentum-rate and analytic-velocity interfaces
!> to the selected constitutive law using canonical source names,
!> and provides a projection for momentum updates caused exclusively by
!> friction. The projection constrains the candidate update to
!> \f$\boldsymbol{m}_{after}=\alpha\boldsymbol{m}_{before}\f$,
!> \f$0\le\alpha\le1\f$, so it cannot rotate, reverse or increase the
!> momentum. Callers must pass the momentum immediately before and after the
!> isolated friction substep, not the full stage state.
!********************************************************************************

MODULE rheology_2d

    USE parameters_2d, ONLY : wp
    USE rheology_parameters_2d, ONLY : rheology_config_type
    USE state_types, ONLY : physical_state_real, physical_state_complex
    USE state_physical, ONLY : update_physical_state_kinematics
    USE diagnostics_2d, ONLY : fatal_error
    USE source_catalog_2d, ONLY : rheology_source_id,RHEOLOGY_ID_NONE,                                                 &
            RHEOLOGY_ID_VOELLMY,RHEOLOGY_ID_PLASTIC,RHEOLOGY_ID_LAVA,RHEOLOGY_ID_LAHAR,                                &
            RHEOLOGY_ID_EMPIRICAL,RHEOLOGY_ID_QUADRATIC,RHEOLOGY_ID_FROUDE,                                            &
            RHEOLOGY_ID_LUCAS,RHEOLOGY_ID_MUI,RHEOLOGY_ID_MUI_VOELLMY,RHEOLOGY_ID_NEWMODEL
    USE rheology_voellmy_2d, ONLY : eval_rheology_voellmy_explicit,                                                    &
                eval_rheology_voellmy_implicit_real,                                                                   &
                eval_rheology_voellmy_implicit,                                                                        &
                eval_rheology_voellmy_semi_implicit,                                                                   &
                eval_rheology_voellmy_analytic
    USE rheology_plastic_2d, ONLY : eval_rheology_plastic_explicit,                                                    &
                eval_rheology_plastic_implicit_real,                                                                   &
                eval_rheology_plastic_implicit,                                                                        &
                eval_rheology_plastic_analytic
    USE rheology_lava_2d, ONLY : eval_rheology_lava_explicit,                                                          &
                eval_rheology_lava_implicit_real,                                                                      &
                eval_rheology_lava_implicit,                                                                           &
                eval_rheology_lava_semi_implicit,                                                                      &
                eval_rheology_lava_analytic
    USE rheology_lahar_2d, ONLY : eval_rheology_lahar_explicit,                                                        &
                eval_rheology_lahar_implicit_real,                                                                     &
                eval_rheology_lahar_implicit,                                                                          &
                eval_rheology_lahar_semi_implicit,                                                                     &
                eval_rheology_lahar_analytic
    USE rheology_empirical_2d, ONLY : eval_rheology_empirical_explicit,                                                &
                eval_rheology_empirical_implicit_real,                                                                 &
                eval_rheology_empirical_implicit,                                                                      &
                eval_rheology_empirical_analytic
    USE rheology_quadratic_2d, ONLY : eval_rheology_quadratic_explicit,                                                &
                eval_rheology_quadratic_implicit_real,                                                                 &
                eval_rheology_quadratic_implicit,                                                                      &
                eval_rheology_quadratic_analytic
    USE rheology_newmodel_2d, ONLY : eval_rheology_newmodel_explicit,                                                  &
                eval_rheology_newmodel_implicit_real,                                                                  &
                eval_rheology_newmodel_implicit
    USE rheology_froude_2d, ONLY : eval_rheology_froude_explicit,                                                      &
                eval_rheology_froude_semi_implicit,                                                                    &
                eval_rheology_froude_analytic
    USE rheology_lucas_2d, ONLY : eval_rheology_lucas_explicit,                                                        &
                eval_rheology_lucas_semi_implicit,                                                                     &
                eval_rheology_lucas_analytic
    USE rheology_mui_2d, ONLY : eval_rheology_mui_explicit,                                                            &
                eval_rheology_mui_semi_implicit,                                                                       &
                eval_rheology_mui_analytic
    USE rheology_mui_voellmy_2d, ONLY : eval_rheology_mui_voellmy_explicit,                                            &
                eval_rheology_mui_voellmy_semi_implicit,                                                               &
                eval_rheology_mui_voellmy_implicit_real,                                                               &
                eval_rheology_mui_voellmy_implicit,                                                                    &
                eval_rheology_mui_voellmy_analytic

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: curvature_acceleration
    PUBLIC :: eval_rheology_explicit, eval_rheology_semi_implicit
    PUBLIC :: eval_rheology_implicit, eval_rheology_analytic
    PUBLIC :: eval_rheology_implicit_id, eval_rheology_implicit_real_id, eval_rheology_semi_implicit_id
    PUBLIC :: limit_rheology_momentum_reversal

        ! Uniform treatment surface for this source family.

CONTAINS

    !> Reject only motion reversal caused by an isolated rheology update.
    !> @brief Prevent a rheology-only momentum update from reversing the original flow.
    !> @param[in] momentum_before Two-dimensional momentum immediately before the rheology substep.
    !> @param[in,out] momentum_after Proposed momentum; its reversed projection is removed in place.
    PURE SUBROUTINE limit_rheology_momentum_reversal(momentum_before, momentum_after)
        REAL(wp), INTENT(IN) :: momentum_before(2)
        REAL(wp), INTENT(INOUT) :: momentum_after(2)
        REAL(wp) :: initial_norm2, reversed_component

        initial_norm2=DOT_PRODUCT(momentum_before,momentum_before)
        IF (initial_norm2<=0.0_wp) RETURN
        reversed_component=DOT_PRODUCT(momentum_after,momentum_before)

        IF (reversed_component<0.0_wp) THEN
        momentum_after=momentum_after-                                                                                 &
                (reversed_component/initial_norm2)*momentum_before

        END IF
    END SUBROUTINE limit_rheology_momentum_reversal


    !> @brief Evaluate the explicit two-component momentum rate for one resolved rheology source.
    !> @param[in] config Immutable constitutive coefficients of the active layer.
    !> @param[in] model_id Pre-resolved rheology source identifier.
    !> @param[in] state Real-valued physical cell state for source evaluation.
    !> @param[in] grav_coeff Topographic gravity correction in the momentum law.
    !> @param[in] use_pore_pressure Whether effective stress includes pore pressure.
    !> @param[in] use_stochastic Whether the stochastic rheology perturbation is active.
    !> @param[in] stochastic_value Cell-local stochastic rheology multiplier.
    !> @param[out] rate Two-component conservative momentum source [kg/(m s2)].
    SUBROUTINE eval_rheology_explicit(config,model_id, state, grav_coeff,                                              &
        use_pore_pressure, use_stochastic, stochastic_value,                                                           &
        rate)
        TYPE(rheology_config_type), INTENT(IN) :: config
        INTEGER, INTENT(IN) :: model_id
        TYPE(physical_state_real), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: grav_coeff
        REAL(wp), INTENT(IN) :: stochastic_value
        LOGICAL, INTENT(IN) :: use_pore_pressure,use_stochastic
        REAL(wp), INTENT(OUT) :: rate(2)
        REAL(wp) :: process_rate(2)
        rate=0.0_wp
        process_rate=0.0_wp
        SELECT CASE(model_id)
            CASE(RHEOLOGY_ID_VOELLMY)
                CALL eval_rheology_voellmy_explicit(config,state, grav_coeff, use_pore_pressure, use_stochastic,       &
                    stochastic_value, process_rate)
            CASE(RHEOLOGY_ID_PLASTIC)
                CALL eval_rheology_plastic_explicit(config,state, process_rate)
            CASE(RHEOLOGY_ID_LAVA)
                CALL eval_rheology_lava_explicit(config,state, process_rate)
            CASE(RHEOLOGY_ID_LAHAR)
                CALL eval_rheology_lahar_explicit(config,state, process_rate)
            CASE(RHEOLOGY_ID_EMPIRICAL)
                CALL eval_rheology_empirical_explicit(state, process_rate)
            CASE(RHEOLOGY_ID_QUADRATIC)
                CALL eval_rheology_quadratic_explicit(config,state, process_rate)
            CASE(RHEOLOGY_ID_NEWMODEL)
                CALL eval_rheology_newmodel_explicit(state, process_rate)
            CASE(RHEOLOGY_ID_FROUDE)
                CALL eval_rheology_froude_explicit(config,state, grav_coeff,                                           &
                    use_stochastic, stochastic_value, process_rate)
            CASE(RHEOLOGY_ID_LUCAS)
                CALL eval_rheology_lucas_explicit(config,state, grav_coeff,                                            &
                    use_stochastic, stochastic_value, process_rate)
            CASE(RHEOLOGY_ID_MUI)
                CALL eval_rheology_mui_explicit(config,state, grav_coeff, use_pore_pressure, process_rate)
            CASE(RHEOLOGY_ID_MUI_VOELLMY)
                CALL eval_rheology_mui_voellmy_explicit(config,state, grav_coeff, use_pore_pressure, use_stochastic,   &
                    stochastic_value, process_rate)
            CASE DEFAULT
                CALL fatal_error('Invalid rheology name in explicit dispatch')

        END SELECT
        rate=process_rate
    END SUBROUTINE eval_rheology_explicit


    !> @brief Sum frozen-state semi-implicit momentum rates for the selected rheologies.
    !> @param[in] config Canonical rheology coefficients for the current layer.
    !> @param[in] models Names of the rheology sources to evaluate.
    !> @param[in] state Frozen real physical state for this Newton evaluation.
    !> @param[in] grav_coeff Topographic gravity correction.
    !> @param[in] use_pore_pressure Whether the effective-stress correction is active.
    !> @param[in] use_stochastic Whether stochastic friction is enabled.
    !> @param[in] stochastic_value Cell-local stochastic rheology coefficient.
    !> @param[out] rate Summed real horizontal momentum-rate components.
    SUBROUTINE eval_rheology_semi_implicit(config,models, state, grav_coeff, use_pore_pressure, use_stochastic,        &
        stochastic_value, rate)
        TYPE(rheology_config_type), INTENT(IN) :: config
        CHARACTER(*), INTENT(IN) :: models(:)
        TYPE(physical_state_real), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: grav_coeff
        REAL(wp), INTENT(IN) :: stochastic_value
        LOGICAL, INTENT(IN) :: use_pore_pressure,use_stochastic
        REAL(wp), INTENT(OUT) :: rate(2)
        REAL(wp) :: process_rate(2)
        INTEGER :: i
        rate=0.0_wp

        DO i=1,SIZE(models)
            process_rate = 0.0_wp

            CALL eval_rheology_semi_implicit_id(config,rheology_source_id(models(i)),state,                            &
                grav_coeff,use_pore_pressure,use_stochastic,stochastic_value,process_rate)

            rate=rate+process_rate

        END DO
    END SUBROUTINE eval_rheology_semi_implicit

    !> Frozen Newton semi-implicit dispatch with pre-resolved source identity.
    !> @brief Evaluate a frozen real rheology rate using a pre-resolved source identifier.
    !> @param[in] config Current layer rheology coefficients.
    !> @param[in] model_id Rheology source ID resolved during configuration.
    !> @param[in] state Frozen physical state for the Newton residual.
    !> @param[in] grav_coeff Corrected gravitational acceleration factor.
    !> @param[in] use_pore_pressure Enable effective-stress pore-pressure correction.
    !> @param[in] use_stochastic Enable stochastic rheology coefficients.
    !> @param[in] stochastic_value Local stochastic coefficient for the source law.
    !> @param[out] rate Two-component semi-implicit momentum source.
    SUBROUTINE eval_rheology_semi_implicit_id(config,model_id,state,grav_coeff,use_pore_pressure,                      &
            use_stochastic,stochastic_value,rate)
        TYPE(rheology_config_type), INTENT(IN) :: config
        INTEGER, INTENT(IN) :: model_id
        TYPE(physical_state_real), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: grav_coeff,stochastic_value
        LOGICAL, INTENT(IN) :: use_pore_pressure,use_stochastic
        REAL(wp), INTENT(OUT) :: rate(2)
        REAL(wp) :: process_rate(2)
        process_rate=0.0_wp
        SELECT CASE(model_id)
            CASE(RHEOLOGY_ID_VOELLMY)
                CALL eval_rheology_voellmy_semi_implicit(config,state, grav_coeff, use_pore_pressure, process_rate(1), &
                    process_rate(2))
            CASE(RHEOLOGY_ID_LAVA)
                CALL eval_rheology_lava_semi_implicit(config,state, process_rate(1), process_rate(2))
            CASE(RHEOLOGY_ID_LAHAR)
                CALL eval_rheology_lahar_semi_implicit(config,state, process_rate(1), process_rate(2))
            CASE(RHEOLOGY_ID_FROUDE)
                CALL eval_rheology_froude_semi_implicit(config,state, grav_coeff, use_stochastic, stochastic_value,    &
                    process_rate(1), process_rate(2))
            CASE(RHEOLOGY_ID_LUCAS)
                CALL eval_rheology_lucas_semi_implicit(config,state, grav_coeff, use_stochastic, stochastic_value,     &
                    process_rate(1), process_rate(2))
            CASE(RHEOLOGY_ID_MUI)
                CALL eval_rheology_mui_semi_implicit(config,state, grav_coeff, use_pore_pressure, process_rate(1),     &
                    process_rate(2))
            CASE(RHEOLOGY_ID_MUI_VOELLMY)
                CALL eval_rheology_mui_voellmy_semi_implicit(config,state, grav_coeff,                                 &
                    use_pore_pressure, process_rate(1),                                                                &
                    process_rate(2))
            CASE(RHEOLOGY_ID_NONE,RHEOLOGY_ID_PLASTIC,RHEOLOGY_ID_EMPIRICAL,                                           &
                    RHEOLOGY_ID_QUADRATIC,RHEOLOGY_ID_NEWMODEL)
                CONTINUE
            CASE DEFAULT
                CALL fatal_error('Invalid rheology name in semi-implicit momentum-process_rate dispatch')

        END SELECT
        rate=process_rate
    END SUBROUTINE eval_rheology_semi_implicit_id


        ! Newton/complex-step contract: read-only state and coefficients; return rates only.
        ! Never advance shared inventories or geometry while evaluating this dispatch.
    !> @brief Sum complex rheology momentum rates for a Newton complex-step evaluation.
    !> @param[in] config Current layer rheology coefficients.
    !> @param[in] models Names of complex-differentiated rheology laws to evaluate.
    !> @param[in] state Unmodified complex physical cell state for differentiation.
    !> @param[in] use_stochastic Whether stochastic rheology coefficients are active.
    !> @param[in] stochastic_value Frozen cell-local stochastic coefficient.
    !> @param[out] rate Sum of complex two-component momentum rates.
    SUBROUTINE eval_rheology_implicit(config,models, state, use_stochastic, stochastic_value, rate)
        TYPE(rheology_config_type), INTENT(IN) :: config
        CHARACTER(*), INTENT(IN) :: models(:)
        TYPE(physical_state_complex), INTENT(IN) :: state
        LOGICAL, INTENT(IN) :: use_stochastic
        REAL(wp), INTENT(IN) :: stochastic_value
        COMPLEX(wp), INTENT(OUT) :: rate(2)
        COMPLEX(wp) :: process_rate(2)
        INTEGER :: i
        rate=CMPLX(0.0_wp,0.0_wp,wp)

        DO i=1,SIZE(models)
            process_rate = CMPLX(0.0_wp,0.0_wp,wp)

            CALL eval_rheology_implicit_id(config,rheology_source_id(models(i)),state,                                 &
                use_stochastic,stochastic_value,process_rate)

            rate=rate+process_rate

        END DO
    END SUBROUTINE eval_rheology_implicit

    !> Newton complex-rate dispatch with source identity resolved at initialization.
    !> @brief Evaluate a complex rheology momentum rate using a resolved source ID.
    !> @param[in] config Current layer rheology parameters.
    !> @param[in] model_id Rheology source identifier selected during setup.
    !> @param[in] state Complex physical cell state for Jacobian differentiation.
    !> @param[in] use_stochastic Whether the stochastic coefficient is applied.
    !> @param[in] stochastic_value Frozen cell-local stochastic coefficient.
    !> @param[out] rate Complex two-component conservative momentum rate.
    SUBROUTINE eval_rheology_implicit_id(config,model_id,state,use_stochastic,stochastic_value,rate)
        TYPE(rheology_config_type), INTENT(IN) :: config
        INTEGER, INTENT(IN) :: model_id
        TYPE(physical_state_complex), INTENT(IN) :: state
        LOGICAL, INTENT(IN) :: use_stochastic
        REAL(wp), INTENT(IN) :: stochastic_value
        COMPLEX(wp), INTENT(OUT) :: rate(2)
        COMPLEX(wp) :: process_rate(2)
        process_rate=CMPLX(0.0_wp,0.0_wp,wp)
        SELECT CASE(model_id)
            CASE(RHEOLOGY_ID_VOELLMY)
                CALL eval_rheology_voellmy_implicit(config,state, use_stochastic, stochastic_value, process_rate(1),   &
                    process_rate(2))
            CASE(RHEOLOGY_ID_MUI_VOELLMY)
                CALL eval_rheology_mui_voellmy_implicit(config,state, use_stochastic,                                  &
                    stochastic_value, process_rate(1),                                                                 &
                    process_rate(2))
            CASE(RHEOLOGY_ID_PLASTIC)
                CALL eval_rheology_plastic_implicit(config,state, process_rate(1), process_rate(2))
            CASE(RHEOLOGY_ID_LAVA)
                CALL eval_rheology_lava_implicit(config,state, process_rate(1), process_rate(2))
            CASE(RHEOLOGY_ID_LAHAR)
                CALL eval_rheology_lahar_implicit(config,state, process_rate(1), process_rate(2))
            CASE(RHEOLOGY_ID_EMPIRICAL)
                CALL eval_rheology_empirical_implicit(state, process_rate(1), process_rate(2))
            CASE(RHEOLOGY_ID_NEWMODEL)
                CALL eval_rheology_newmodel_implicit(state,process_rate(1),process_rate(2))
            CASE(RHEOLOGY_ID_QUADRATIC)
                CALL eval_rheology_quadratic_implicit(config,state, process_rate(1), process_rate(2))
            CASE(RHEOLOGY_ID_NONE, RHEOLOGY_ID_FROUDE, RHEOLOGY_ID_LUCAS, RHEOLOGY_ID_MUI)
                CONTINUE
            CASE DEFAULT
                CALL fatal_error('Invalid rheology name in implicit momentum-process_rate dispatch')

        END SELECT
        rate=process_rate
    END SUBROUTINE eval_rheology_implicit_id


    !> Evaluate the implicit rheology contribution in REAL for Newton residuals.
    !> @brief Evaluate the real counterpart of a Newton rheology rate without complex conversion.
    !> @param[in] config Layer-specific rheology coefficients.
    !> @param[in] model_id Selected rheology source identifier.
    !> @param[in] state Real physical cell state at the Newton evaluation point.
    !> @param[in] use_stochastic Whether stochastic rheology is active.
    !> @param[in] stochastic_value Frozen stochastic coefficient.
    !> @param[out] rate Real two-component rate matching the complex residual law.
    SUBROUTINE eval_rheology_implicit_real_id(config,model_id,state,use_stochastic,stochastic_value,rate)
        TYPE(rheology_config_type), INTENT(IN) :: config
        INTEGER, INTENT(IN) :: model_id
        TYPE(physical_state_real), INTENT(IN) :: state
        LOGICAL, INTENT(IN) :: use_stochastic
        REAL(wp), INTENT(IN) :: stochastic_value
        REAL(wp), INTENT(OUT) :: rate(2)
        rate=0.0_wp
        SELECT CASE(model_id)
        CASE(RHEOLOGY_ID_VOELLMY)
            CALL eval_rheology_voellmy_implicit_real(config,state,use_stochastic,stochastic_value,rate(1),rate(2))
        CASE(RHEOLOGY_ID_MUI_VOELLMY)
            CALL eval_rheology_mui_voellmy_implicit_real(config,state,use_stochastic,stochastic_value,rate(1),rate(2))
        CASE(RHEOLOGY_ID_PLASTIC)
            CALL eval_rheology_plastic_implicit_real(config,state,rate(1),rate(2))
        CASE(RHEOLOGY_ID_LAVA)
            CALL eval_rheology_lava_implicit_real(config,state,rate(1),rate(2))
        CASE(RHEOLOGY_ID_LAHAR)
            CALL eval_rheology_lahar_implicit_real(config,state,rate(1),rate(2))
        CASE(RHEOLOGY_ID_EMPIRICAL)
            CALL eval_rheology_empirical_implicit_real(state,rate(1),rate(2))
        CASE(RHEOLOGY_ID_NEWMODEL)
            CALL eval_rheology_newmodel_implicit_real(state,rate(1),rate(2))
        CASE(RHEOLOGY_ID_QUADRATIC)
            CALL eval_rheology_quadratic_implicit_real(config,state,rate(1),rate(2))
        CASE(RHEOLOGY_ID_NONE,RHEOLOGY_ID_FROUDE,RHEOLOGY_ID_LUCAS,RHEOLOGY_ID_MUI)
            CONTINUE
        CASE DEFAULT
            CALL fatal_error('Invalid rheology ID in REAL implicit dispatch')
        END SELECT
    END SUBROUTINE eval_rheology_implicit_real_id

    !> Integrate analytic rheology without copying state in the single-law Lie case.
    !> @brief Integrate the selected analytic rheologies over one Lie source substep.
    !> @param[in] config Current layer rheology coefficients.
    !> @param[in] models Ordered rheology law names to integrate sequentially.
    !> @param[in] dt Positive source integration interval [s].
    !> @param[in] state Physical state before analytic momentum integration.
    !> @param[in] speed_ratio Ratio of total to horizontal speed used by friction laws.
    !> @param[in] grav_coeff Topographic gravitational correction.
    !> @param[in] curvature_acceleration Normal acceleration induced by terrain curvature.
    !> @param[in] use_pore_pressure Enable effective-stress pore correction.
    !> @param[in] use_stochastic Enable stochastic rheological coefficients.
    !> @param[in] stochastic_value Frozen local stochastic coefficient.
    !> @param[out] velocity_after Horizontal velocity after all supported analytic substeps.
    !> @param[out] model_supported False if any requested analytic rheology is unavailable.
    !> @param[in] bx Terrain x-slope used to refresh kinematics between laws.
    !> @param[in] by Terrain y-slope used to refresh kinematics between laws.
    !> @param[in] bxx Terrain second x derivative for curvature kinematics.
    !> @param[in] bxy Terrain mixed second derivative for curvature kinematics.
    !> @param[in] byy Terrain second y derivative for curvature kinematics.
    SUBROUTINE eval_rheology_analytic(config,models, dt, state, speed_ratio, grav_coeff, curvature_acceleration,       &
        use_pore_pressure, use_stochastic, stochastic_value, velocity_after, model_supported, bx, by, bxx, bxy, byy)
        TYPE(rheology_config_type), INTENT(IN) :: config
        CHARACTER(*), INTENT(IN) :: models(:)
        REAL(wp), INTENT(IN) :: dt,speed_ratio,grav_coeff,curvature_acceleration,stochastic_value
        TYPE(physical_state_real), INTENT(IN) :: state
        LOGICAL, INTENT(IN) :: use_pore_pressure,use_stochastic
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        LOGICAL, INTENT(OUT) :: model_supported
        REAL(wp), INTENT(IN) :: bx,by,bxx,bxy,byy
        TYPE(physical_state_real) :: current_state
        INTEGER :: i

        velocity_after=[state%u,state%v]
        model_supported=.TRUE.
        IF (SIZE(models)==0) RETURN

        IF (SIZE(models)==1) THEN
            CALL dispatch_analytic_rheology(config,models(1),dt,state,speed_ratio,grav_coeff,                          &
                curvature_acceleration,use_pore_pressure,use_stochastic,stochastic_value,                              &
                velocity_after,model_supported)
            RETURN
        END IF

        ! Sequential laws require updated kinematics after every substep.
        ! Only this uncommon multi-law path copies the physical composition.
        current_state=state
        DO i=1,SIZE(models)
            CALL dispatch_analytic_rheology(config,models(i),dt,current_state,speed_ratio,grav_coeff,                  &
                curvature_acceleration,use_pore_pressure,use_stochastic,stochastic_value,                              &
                velocity_after,model_supported)
            IF (.NOT.model_supported) RETURN
            current_state%u=velocity_after(1)
            current_state%v=velocity_after(2)
            CALL update_physical_state_kinematics(current_state,bx,by,bxx,bxy,byy)
        END DO
    END SUBROUTINE eval_rheology_analytic

    !> Select the model-specific exact velocity flow for one frozen physical state.
    !> @brief Dispatch one analytic rheology source on the supplied frozen physical state.
    !> @param[in] config Layer rheology coefficients.
    !> @param[in] model Registered rheology source name to integrate.
    !> @param[in] dt Duration of this analytic Lie substep [s].
    !> @param[in] state Physical state at the start of this selected substep.
    !> @param[in] speed_ratio Ratio of total speed to horizontal speed.
    !> @param[in] grav_coeff Topographic gravitational factor.
    !> @param[in] curvature_acceleration Curvature-induced acceleration for the source law.
    !> @param[in] use_pore_pressure Whether the friction law accounts for pore pressure.
    !> @param[in] use_stochastic Whether stochastic rheology is active.
    !> @param[in] stochastic_value Local stochastic rheology coefficient.
    !> @param[out] velocity_after Updated horizontal velocity returned by the selected law.
    !> @param[out] model_supported False when this name has no supported analytic treatment.
    SUBROUTINE dispatch_analytic_rheology(config,model,dt,state,speed_ratio,grav_coeff,curvature_acceleration,         &
        use_pore_pressure,use_stochastic,stochastic_value,velocity_after,model_supported)
        TYPE(rheology_config_type), INTENT(IN) :: config
        CHARACTER(*), INTENT(IN) :: model
        REAL(wp), INTENT(IN) :: dt,speed_ratio,grav_coeff,curvature_acceleration,stochastic_value
        TYPE(physical_state_real), INTENT(IN) :: state
        LOGICAL, INTENT(IN) :: use_pore_pressure,use_stochastic
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        LOGICAL, INTENT(OUT) :: model_supported

        velocity_after=[state%u,state%v]
        model_supported=.TRUE.
        SELECT CASE(TRIM(model))
        CASE ('none')
            CONTINUE
        CASE ('rheology_voellmy')
            CALL eval_rheology_voellmy_analytic(config,dt, state, speed_ratio, grav_coeff, curvature_acceleration,     &
                use_pore_pressure, use_stochastic, stochastic_value, velocity_after)
        CASE ('rheology_plastic')
            CALL eval_rheology_plastic_analytic(config,dt, state, speed_ratio, velocity_after)
        CASE ('rheology_lava')
            CALL eval_rheology_lava_analytic(config,dt, state, speed_ratio, velocity_after)
        CASE ('rheology_lahar')
            CALL eval_rheology_lahar_analytic(config,dt, state, speed_ratio, velocity_after)
        CASE ('rheology_empirical')
            CALL eval_rheology_empirical_analytic(dt, state, speed_ratio, velocity_after)
        CASE ('rheology_quadratic')
            CALL eval_rheology_quadratic_analytic(config,dt, state, speed_ratio, velocity_after)
        CASE ('rheology_froude')
            CALL eval_rheology_froude_analytic(config,dt, state, speed_ratio, grav_coeff, curvature_acceleration,      &
                use_stochastic, stochastic_value, velocity_after)
        CASE ('rheology_lucas')
            CALL eval_rheology_lucas_analytic(config,dt, state, speed_ratio, grav_coeff, curvature_acceleration,       &
                use_stochastic, stochastic_value, velocity_after)
        CASE ('rheology_mui')
            CALL eval_rheology_mui_analytic(config,dt, state, speed_ratio, grav_coeff, curvature_acceleration,         &
                use_pore_pressure, velocity_after)
        CASE ('rheology_mui_voellmy')
            CALL eval_rheology_mui_voellmy_analytic(config,dt, state, speed_ratio, grav_coeff,                         &
                curvature_acceleration, use_pore_pressure, use_stochastic, stochastic_value, velocity_after)
        CASE DEFAULT
            model_supported = .FALSE.

        END SELECT
    END SUBROUTINE dispatch_analytic_rheology



        !> \brief Compute curvature-induced normal acceleration.
        !> \details Contract the topographic Hessian with the horizontal velocity twice: Bxx*u**2 +
        !> 2*Bxy*u*v + Byy*v**2.
        !> \param[in] Bsecond_xx Second x derivative of topography [m^-1].
        !> \param[in] Bsecond_xy Mixed derivative of topography [m^-1].
        !> \param[in] Bsecond_yy Second y derivative of topography [m^-1].
        !> \param[in] u x velocity [m s^-1].
        !> \param[in] v y velocity [m s^-1].
        !> \return Computed acceleration for the requested model and inputs.
    PURE FUNCTION curvature_acceleration(Bsecond_xx, Bsecond_xy, Bsecond_yy,                                           &
                u, v) RESULT(acceleration)

            REAL(wp), INTENT(IN) :: Bsecond_xx, Bsecond_xy, Bsecond_yy
            REAL(wp), INTENT(IN) :: u, v
            REAL(wp) :: acceleration

            acceleration = Bsecond_xx*u**2 +                                                                           &
                        2.0_wp*Bsecond_xy*u*v + Bsecond_yy*v**2

    END FUNCTION curvature_acceleration


        ! Additive rates from all active rheology processes at the same stage state.



        ! Exact flows are composed in order, never added as rates.

END MODULE rheology_2d
