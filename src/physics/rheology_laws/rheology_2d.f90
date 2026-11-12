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
                eval_rheology_voellmy_implicit_real, &
                eval_rheology_voellmy_implicit,                                                                        &
                eval_rheology_voellmy_semi_implicit,                                                                   &
                eval_rheology_voellmy_analytic
    USE rheology_plastic_2d, ONLY : eval_rheology_plastic_explicit,                                                    &
                eval_rheology_plastic_implicit_real, &
                eval_rheology_plastic_implicit,                                                                        &
                eval_rheology_plastic_analytic
    USE rheology_lava_2d, ONLY : eval_rheology_lava_explicit,                                                          &
                eval_rheology_lava_implicit_real, &
                eval_rheology_lava_implicit,                                                                           &
                eval_rheology_lava_semi_implicit,                                                                      &
                eval_rheology_lava_analytic
    USE rheology_lahar_2d, ONLY : eval_rheology_lahar_explicit,                                                        &
                eval_rheology_lahar_implicit_real, &
                eval_rheology_lahar_implicit,                                                                          &
                eval_rheology_lahar_semi_implicit,                                                                     &
                eval_rheology_lahar_analytic
    USE rheology_empirical_2d, ONLY : eval_rheology_empirical_explicit,                                                &
                eval_rheology_empirical_implicit_real, &
                eval_rheology_empirical_implicit,                                                                      &
                eval_rheology_empirical_analytic
    USE rheology_quadratic_2d, ONLY : eval_rheology_quadratic_explicit,                                                &
                eval_rheology_quadratic_implicit_real, &
                eval_rheology_quadratic_implicit,                                                                      &
                eval_rheology_quadratic_analytic
    USE rheology_newmodel_2d, ONLY : eval_rheology_newmodel_explicit,                                                  &
                eval_rheology_newmodel_implicit_real, &
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
                eval_rheology_mui_voellmy_implicit_real, &
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


    SUBROUTINE eval_rheology_explicit(config,model_id, state, grav_coeff, use_pore_pressure, use_stochastic, stochastic_value,  &
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
                CALL eval_rheology_voellmy_explicit(state, grav_coeff, use_pore_pressure, use_stochastic,              &
                    stochastic_value, process_rate)
            CASE(RHEOLOGY_ID_PLASTIC)
                CALL eval_rheology_plastic_explicit(state, process_rate)
            CASE(RHEOLOGY_ID_LAVA)
                CALL eval_rheology_lava_explicit(state, process_rate)
            CASE(RHEOLOGY_ID_LAHAR)
                CALL eval_rheology_lahar_explicit(state, process_rate)
            CASE(RHEOLOGY_ID_EMPIRICAL)
                CALL eval_rheology_empirical_explicit(state, process_rate)
            CASE(RHEOLOGY_ID_QUADRATIC)
                CALL eval_rheology_quadratic_explicit(config,state, process_rate)
            CASE(RHEOLOGY_ID_NEWMODEL)
                CALL eval_rheology_newmodel_explicit(state, process_rate)
            CASE(RHEOLOGY_ID_FROUDE)
                CALL eval_rheology_froude_explicit(state, grav_coeff, use_stochastic, stochastic_value, process_rate)
            CASE(RHEOLOGY_ID_LUCAS)
                CALL eval_rheology_lucas_explicit(state, grav_coeff, use_stochastic, stochastic_value, process_rate)
            CASE(RHEOLOGY_ID_MUI)
                CALL eval_rheology_mui_explicit(state, grav_coeff, use_pore_pressure, process_rate)
            CASE(RHEOLOGY_ID_MUI_VOELLMY)
                CALL eval_rheology_mui_voellmy_explicit(state, grav_coeff, use_pore_pressure, use_stochastic,          &
                    stochastic_value, process_rate)
            CASE DEFAULT
                CALL fatal_error('Invalid rheology name in explicit dispatch')

        END SELECT
        rate=process_rate
    END SUBROUTINE eval_rheology_explicit


    SUBROUTINE eval_rheology_semi_implicit(config,models, state, grav_coeff, use_pore_pressure, use_stochastic,               &
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

            CALL eval_rheology_semi_implicit_id(config,rheology_source_id(models(i)),state,                                   &
                grav_coeff,use_pore_pressure,use_stochastic,stochastic_value,process_rate)

            rate=rate+process_rate

        END DO
    END SUBROUTINE eval_rheology_semi_implicit

    !> Frozen Newton semi-implicit dispatch with pre-resolved source identity.
    SUBROUTINE eval_rheology_semi_implicit_id(config,model_id,state,grav_coeff,use_pore_pressure,                             &
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
                CALL eval_rheology_voellmy_semi_implicit(state, grav_coeff, use_pore_pressure, process_rate(1),        &
                    process_rate(2))
            CASE(RHEOLOGY_ID_LAVA)
                CALL eval_rheology_lava_semi_implicit(state, process_rate(1), process_rate(2))
            CASE(RHEOLOGY_ID_LAHAR)
                CALL eval_rheology_lahar_semi_implicit(state, process_rate(1), process_rate(2))
            CASE(RHEOLOGY_ID_FROUDE)
                CALL eval_rheology_froude_semi_implicit(state, grav_coeff, use_stochastic, stochastic_value,           &
                    process_rate(1), process_rate(2))
            CASE(RHEOLOGY_ID_LUCAS)
                CALL eval_rheology_lucas_semi_implicit(state, grav_coeff, use_stochastic, stochastic_value,            &
                    process_rate(1), process_rate(2))
            CASE(RHEOLOGY_ID_MUI)
                CALL eval_rheology_mui_semi_implicit(state, grav_coeff, use_pore_pressure, process_rate(1),            &
                    process_rate(2))
            CASE(RHEOLOGY_ID_MUI_VOELLMY)
                CALL eval_rheology_mui_voellmy_semi_implicit(state, grav_coeff, use_pore_pressure, process_rate(1),    &
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

            CALL eval_rheology_implicit_id(config,rheology_source_id(models(i)),state,                                        &
                use_stochastic,stochastic_value,process_rate)

            rate=rate+process_rate

        END DO
    END SUBROUTINE eval_rheology_implicit

    !> Newton complex-rate dispatch with source identity resolved at initialization.
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
                CALL eval_rheology_voellmy_implicit(state, use_stochastic, stochastic_value, process_rate(1),          &
                    process_rate(2))
            CASE(RHEOLOGY_ID_MUI_VOELLMY)
                CALL eval_rheology_mui_voellmy_implicit(state, use_stochastic, stochastic_value, process_rate(1),      &
                    process_rate(2))
            CASE(RHEOLOGY_ID_PLASTIC)
                CALL eval_rheology_plastic_implicit(state, process_rate(1), process_rate(2))
            CASE(RHEOLOGY_ID_LAVA)
                CALL eval_rheology_lava_implicit(state, process_rate(1), process_rate(2))
            CASE(RHEOLOGY_ID_LAHAR)
                CALL eval_rheology_lahar_implicit(state, process_rate(1), process_rate(2))
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
            CALL eval_rheology_voellmy_implicit_real(state,use_stochastic,stochastic_value,rate(1),rate(2))
        CASE(RHEOLOGY_ID_MUI_VOELLMY)
            CALL eval_rheology_mui_voellmy_implicit_real(state,use_stochastic,stochastic_value,rate(1),rate(2))
        CASE(RHEOLOGY_ID_PLASTIC)
            CALL eval_rheology_plastic_implicit_real(state,rate(1),rate(2))
        CASE(RHEOLOGY_ID_LAVA)
            CALL eval_rheology_lava_implicit_real(state,rate(1),rate(2))
        CASE(RHEOLOGY_ID_LAHAR)
            CALL eval_rheology_lahar_implicit_real(state,rate(1),rate(2))
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
    SUBROUTINE eval_rheology_analytic(config,models, dt, state, speed_ratio, grav_coeff, curvature_acceleration,              &
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
            CALL dispatch_analytic_rheology(config,models(1),dt,state,speed_ratio,grav_coeff,                                 &
                curvature_acceleration,use_pore_pressure,use_stochastic,stochastic_value,                              &
                velocity_after,model_supported)
            RETURN
        END IF

        ! Sequential laws require updated kinematics after every substep.
        ! Only this uncommon multi-law path copies the physical composition.
        current_state=state
        DO i=1,SIZE(models)
            CALL dispatch_analytic_rheology(config,models(i),dt,current_state,speed_ratio,grav_coeff,                         &
                curvature_acceleration,use_pore_pressure,use_stochastic,stochastic_value,                              &
                velocity_after,model_supported)
            IF (.NOT.model_supported) RETURN
            current_state%u=velocity_after(1)
            current_state%v=velocity_after(2)
            CALL update_physical_state_kinematics(current_state,bx,by,bxx,bxy,byy)
        END DO
    END SUBROUTINE eval_rheology_analytic

    !> Select the model-specific exact velocity flow for one frozen physical state.
    SUBROUTINE dispatch_analytic_rheology(config,model,dt,state,speed_ratio,grav_coeff,curvature_acceleration,                &
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
            CALL eval_rheology_voellmy_analytic(dt, state, speed_ratio, grav_coeff, curvature_acceleration,&
                use_pore_pressure, use_stochastic, stochastic_value, velocity_after)
        CASE ('rheology_plastic')
            CALL eval_rheology_plastic_analytic(dt, state, speed_ratio, velocity_after)
        CASE ('rheology_lava')
            CALL eval_rheology_lava_analytic(dt, state, speed_ratio, velocity_after)
        CASE ('rheology_lahar')
            CALL eval_rheology_lahar_analytic(dt, state, speed_ratio, velocity_after)
        CASE ('rheology_empirical')
            CALL eval_rheology_empirical_analytic(dt, state, speed_ratio, velocity_after)
        CASE ('rheology_quadratic')
            CALL eval_rheology_quadratic_analytic(config,dt, state, speed_ratio, velocity_after)
        CASE ('rheology_froude')
            CALL eval_rheology_froude_analytic(dt, state, speed_ratio, grav_coeff, curvature_acceleration,             &
                use_stochastic, stochastic_value, velocity_after)
        CASE ('rheology_lucas')
            CALL eval_rheology_lucas_analytic(dt, state, speed_ratio, grav_coeff, curvature_acceleration,              &
                use_stochastic, stochastic_value, velocity_after)
        CASE ('rheology_mui')
            CALL eval_rheology_mui_analytic(dt, state, speed_ratio, grav_coeff, curvature_acceleration,                &
                use_pore_pressure, velocity_after)
        CASE ('rheology_mui_voellmy')
            CALL eval_rheology_mui_voellmy_analytic(dt, state, speed_ratio, grav_coeff,                                &
                curvature_acceleration, use_pore_pressure, use_stochastic, stochastic_value, velocity_after)
        CASE DEFAULT
            model_supported = .FALSE.

        END SELECT
    END SUBROUTINE dispatch_analytic_rheology



        ! Both Newton components are evaluated on the same, unchanged physical state.

        !> \brief Dispatch the explicit real-valued rheology momentum source.
        !> \param[in] model_name Canonical rheology model name.
        !> \param[in] state Shared real-valued physical state.
        !> \param[in] mod_vel Total speed [m s^-1].
        !> \param[in] mod_hor_vel Horizontal speed [m s^-1].
        !> \param[in] mod_vel2 Squared total speed [m^2 s^-2].
        !> \param[in] grav_coeff Topographic gravity correction.
        !> \param[in] use_pore_pressure Enable effective-stress pore-pressure correction.
        !> \param[in] moving Whether the flow is moving.
        !> \param[in] use_stochastic Enable stochastic rheology coefficient.
        !> \param[in] stochastic_value Stochastic rheology coefficient.
        !> \param[out] rate Two-component momentum source.


        !> \brief Dispatch the semi-implicit real-valued rheology momentum source.
        !> \param[in] model_name Canonical rheology model name.
        !> \param[in] state Shared real-valued physical state.
        !> \param[in] mod_vel Total speed [m s^-1].
        !> \param[in] mod_hor_vel Horizontal speed [m s^-1].
        !> \param[in] grav_coeff Topographic gravity correction.
        !> \param[in] use_pore_pressure Enable pore-pressure correction.
        !> \param[in] use_stochastic Enable stochastic rheology coefficient.
        !> \param[in] stochastic_value Stochastic rheology coefficient.
        !> \param[out] rate Two-component momentum source.


        !> \brief Dispatch the complex-valued implicit rheology momentum source.
        !> \param[in] model_name Canonical rheology model name.
        !> \param[in] state Shared complex-valued physical state for complex-step differentiation.
        !> \param[in] mod_vel Total complex speed.
        !> \param[in] mod_hor_vel Horizontal complex speed.
        !> \param[in] mod_vel2 Squared total complex speed.
        !> \param[in] moving Whether the complex flow state is moving.
        !> \param[in] use_stochastic Enable stochastic rheology coefficient.
        !> \param[in] stochastic_value Stochastic rheology coefficient.
        !> \param[out] rate Two-component complex momentum source.


        !> \brief Dispatch exact velocity integration for the selected rheology.
        !> \param[in] model_name Canonical rheology model name.
        !> \param[in] dt Source integration interval [s].
        !> \param[in] state Shared real physical state before source integration.
        !> \param[in] speed_ratio Total-to-horizontal speed ratio.
        !> \param[in] grav_coeff Topographic gravity correction.
        !> \param[in] curvature_acceleration Curvature-induced normal acceleration [m s^-2].
        !> \param[in] use_pore_pressure Enable pore-pressure correction.
        !> \param[in] use_stochastic Enable stochastic rheology coefficient.
        !> \param[in] stochastic_value Stochastic rheology coefficient.
        !> \param[out] velocity_after Velocity after the analytic source substep.
        !> \param[out] model_supported Whether the selected model has an analytic treatment.


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
