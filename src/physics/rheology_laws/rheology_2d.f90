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
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
    USE diagnostics_2d, ONLY : fatal_error
    USE rheology_voellmy_2d, ONLY : eval_rheology_voellmy_explicit,                                &
            eval_rheology_voellmy_implicit,                                                        &
            eval_rheology_voellmy_semi_implicit,                                                   &
            solve_rheology_voellmy_analytic
    USE rheology_plastic_2d, ONLY : eval_rheology_plastic_explicit,                                &
            eval_rheology_plastic_implicit,                                                        &
            solve_rheology_plastic_analytic
    USE rheology_lava_2d, ONLY : eval_rheology_lava_explicit,                                      &
            eval_rheology_lava_implicit,                                                           &
            eval_rheology_lava_semi_implicit,                                                      &
            solve_rheology_lava_analytic
    USE rheology_lahar_2d, ONLY : eval_rheology_lahar_explicit,                                    &
            eval_rheology_lahar_implicit,                                                          &
            eval_rheology_lahar_semi_implicit,                                                     &
            solve_rheology_lahar_analytic
    USE rheology_empirical_2d, ONLY : eval_rheology_empirical_explicit,                            &
            eval_rheology_empirical_implicit,                                                      &
            solve_rheology_empirical_analytic
    USE rheology_quadratic_2d, ONLY : eval_rheology_quadratic_explicit,                            &
            eval_rheology_quadratic_implicit,                                                      &
            solve_rheology_quadratic_analytic
    USE rheology_newmodel_2d, ONLY : eval_rheology_newmodel_explicit,                              &
            eval_rheology_newmodel_implicit
    USE rheology_froude_2d, ONLY : eval_rheology_froude_explicit,                                  &
            eval_rheology_froude_semi_implicit,                                                    &
            solve_rheology_froude_analytic
    USE rheology_lucas_2d, ONLY : eval_rheology_lucas_explicit,                                    &
            eval_rheology_lucas_semi_implicit,                                                     &
            solve_rheology_lucas_analytic
    USE rheology_mui_2d, ONLY : eval_rheology_mui_explicit,                                        &
            eval_rheology_mui_semi_implicit,                                                       &
            solve_rheology_mui_analytic
    USE rheology_mui_voellmy_2d, ONLY : eval_rheology_mui_voellmy_explicit,                        &
            eval_rheology_mui_voellmy_semi_implicit,                                               &
            eval_rheology_mui_voellmy_implicit,                                                    &
            solve_rheology_mui_voellmy_analytic

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: eval_rheology_explicit
    PUBLIC :: eval_rheology_semi_implicit
    PUBLIC :: eval_rheology_implicit
    PUBLIC :: solve_rheology_analytic
    PUBLIC :: curvature_acceleration
    PUBLIC :: eval_rheology_explicit_processes, eval_rheology_semi_implicit_processes
    PUBLIC :: eval_rheology_implicit_processes, solve_rheology_analytic_processes

    ! Uniform treatment surface for this source family.

CONTAINS

    ! Both Newton components are evaluated on the same, unchanged physical state.

    !> \brief Dispatch the explicit real-valued rheology momentum source.
    !> \param[in] model_name Canonical rheology model name.
    !> \param[in] state Shared real-valued physical state.
    !> \param[in] mod_vel Total speed [m s^-1].
    !> \param[in] mod_hor_vel Horizontal speed [m s^-1].
    !> \param[in] mod_vel2 Squared total speed [m^2 s^-2].
    !> \param[in] grav_coeff Topographic gravity correction.
    !> \param[in] curvature_term Enable curvature-dependent normal stress.
    !> \param[in] B_second_xx, B_second_xy, B_second_yy Bed Hessian entries [m^-1].
    !> \param[in] use_pore_pressure Enable effective-stress pore-pressure correction.
    !> \param[in] moving Whether the flow is moving.
    !> \param[in] use_stochastic Enable stochastic rheology coefficient.
    !> \param[in] stochastic_value Stochastic rheology coefficient.
    !> \param[out] rate Two-component momentum source.
    SUBROUTINE eval_rheology_explicit(model_name, state, mod_vel, mod_hor_vel, mod_vel2, &
            grav_coeff, curvature_term, B_second_xx, B_second_xy, B_second_yy, &
            use_pore_pressure, moving, use_stochastic, stochastic_value, rate)

        CHARACTER(*), INTENT(IN) :: model_name
        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: mod_vel, mod_hor_vel, mod_vel2, grav_coeff, B_second_xx, B_second_xy, B_second_yy
        LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure, moving, use_stochastic
        REAL(wp), INTENT(IN) :: stochastic_value
        REAL(wp), INTENT(OUT) :: rate(2)

        SELECT CASE(TRIM(model_name))
        CASE('rheology_voellmy')
            CALL eval_rheology_voellmy_explicit(state, &
                    mod_vel, &
                    mod_vel2, &
                    grav_coeff, &
                    curvature_term, &
                    B_second_xx, &
                    B_second_xy, &
                    B_second_yy, &
                    use_pore_pressure, &
                    moving, &
                    use_stochastic, &
                    stochastic_value, &
                    rate)
        CASE('rheology_plastic')
            CALL eval_rheology_plastic_explicit(state, &
                    mod_vel, &
                    moving, &
                    rate)
        CASE('rheology_lava')
            CALL eval_rheology_lava_explicit(state, &
                    mod_vel, &
                    moving, &
                    rate)
        CASE('rheology_lahar')
            CALL eval_rheology_lahar_explicit(state, &
                    mod_vel, &
                    mod_hor_vel, &
                    moving, &
                    rate)
        CASE('rheology_empirical')
            CALL eval_rheology_empirical_explicit(state, &
                    mod_vel, &
                    mod_hor_vel, &
                    moving, &
                    rate)
        CASE('rheology_quadratic')
            CALL eval_rheology_quadratic_explicit(state, &
                    mod_hor_vel, &
                    mod_vel2, &
                    moving, &
                    rate)
        CASE('rheology_newmodel')
            CALL eval_rheology_newmodel_explicit(state, &
                    rate)
        CASE('rheology_froude')
            CALL eval_rheology_froude_explicit(state, &
                    mod_vel, &
                    mod_hor_vel, &
                    grav_coeff, &
                    use_stochastic, &
                    stochastic_value, &
                    rate)
        CASE('rheology_lucas')
            CALL eval_rheology_lucas_explicit(state, &
                    mod_vel, &
                    mod_hor_vel, &
                    grav_coeff, &
                    use_stochastic, &
                    stochastic_value, &
                    rate)
        CASE('rheology_mui')
            CALL eval_rheology_mui_explicit(state, &
                    mod_vel, &
                    grav_coeff, &
                    curvature_term, &
                    B_second_xx, &
                    B_second_xy, &
                    B_second_yy, &
                    use_pore_pressure, &
                    rate)
        CASE('rheology_mui_voellmy')
            CALL eval_rheology_mui_voellmy_explicit(state, &
                    mod_vel, &
                    mod_vel2, &
                    grav_coeff, &
                    curvature_term, &
                    B_second_xx, &
                    B_second_xy, &
                    B_second_yy, &
                    use_pore_pressure, &
                    moving, &
                    use_stochastic, &
                    stochastic_value, &
                    rate)
        CASE DEFAULT
            CALL fatal_error('Invalid rheology name in explicit dispatch')

        END SELECT
    END SUBROUTINE eval_rheology_explicit

    !> \brief Dispatch the semi-implicit real-valued rheology momentum source.
    !> \param[in] model_name Canonical rheology model name.
    !> \param[in] state Shared real-valued physical state.
    !> \param[in] mod_vel Total speed [m s^-1].
    !> \param[in] mod_hor_vel Horizontal speed [m s^-1].
    !> \param[in] grav_coeff Topographic gravity correction.
    !> \param[in] curvature_term Enable curvature-dependent normal stress.
    !> \param[in] B_second_xx, B_second_xy, B_second_yy Bed Hessian entries [m^-1].
    !> \param[in] use_pore_pressure Enable pore-pressure correction.
    !> \param[in] use_stochastic Enable stochastic rheology coefficient.
    !> \param[in] stochastic_value Stochastic rheology coefficient.
    !> \param[out] rate Two-component momentum source.
    SUBROUTINE eval_rheology_semi_implicit(model_name, state, mod_vel, mod_hor_vel, grav_coeff, curvature_term, &
            B_second_xx, B_second_xy, B_second_yy, use_pore_pressure, &
            use_stochastic, stochastic_value, rate)

        CHARACTER(*), INTENT(IN) :: model_name
        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: mod_vel, mod_hor_vel, grav_coeff, B_second_xx, B_second_xy, B_second_yy
        LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure, use_stochastic
        REAL(wp), INTENT(IN) :: stochastic_value
        REAL(wp), INTENT(OUT) :: rate(2)

        REAL(wp) :: curvature_term_value

        rate = 0.0_wp

        SELECT CASE(TRIM(model_name))
        CASE('rheology_voellmy')
            CALL eval_rheology_voellmy_semi_implicit(state, mod_vel, grav_coeff, curvature_term, &
                    B_second_xx, B_second_xy, B_second_yy, use_pore_pressure, rate(1), rate(2))
        CASE('rheology_lava')
            CALL eval_rheology_lava_semi_implicit(state, mod_vel, rate(1), rate(2))
        CASE('rheology_lahar')
            CALL eval_rheology_lahar_semi_implicit(state, mod_hor_vel, mod_hor_vel>0.0_wp, rate(1), rate(2))
        CASE('rheology_froude')
            CALL eval_rheology_froude_semi_implicit(state, mod_vel, mod_hor_vel, grav_coeff, &
                    0.0_wp, use_stochastic, stochastic_value, rate(1), rate(2))
        CASE('rheology_lucas')
            CALL eval_rheology_lucas_semi_implicit(state, mod_vel, mod_hor_vel, grav_coeff, &
                    0.0_wp, use_stochastic, stochastic_value, rate(1), rate(2))
        CASE('rheology_mui')
            curvature_term_value = 0.0_wp
            IF (curvature_term) curvature_term_value = curvature_acceleration(                     &
                    B_second_xx, B_second_xy, B_second_yy, state%u, state%v)
            CALL eval_rheology_mui_semi_implicit(state, mod_vel, grav_coeff, curvature_term_value, &
                    use_pore_pressure, rate(1), rate(2))
        CASE('rheology_mui_voellmy')
            curvature_term_value = 0.0_wp
            IF (curvature_term) curvature_term_value = curvature_acceleration(                     &
                    B_second_xx, B_second_xy, B_second_yy, state%u, state%v)
            CALL eval_rheology_mui_voellmy_semi_implicit(state, mod_vel, grav_coeff, curvature_term_value, &
                    use_pore_pressure, rate(1), rate(2))
        CASE('none', 'rheology_plastic', 'rheology_empirical',                                     &
                'rheology_quadratic', 'rheology_newmodel')
            CONTINUE
        CASE DEFAULT
            CALL fatal_error('Invalid rheology name in semi-implicit momentum-rate dispatch')

        END SELECT
    END SUBROUTINE eval_rheology_semi_implicit

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
    SUBROUTINE eval_rheology_implicit(model_name, state, mod_vel, mod_hor_vel, mod_vel2, &
            moving, use_stochastic, stochastic_value, rate)

        CHARACTER(*), INTENT(IN) :: model_name
        TYPE(complex_physical_state_type), INTENT(IN) :: state
        COMPLEX(wp), INTENT(IN) :: mod_vel, mod_hor_vel, mod_vel2
        LOGICAL, INTENT(IN) :: moving, use_stochastic
        REAL(wp), INTENT(IN) :: stochastic_value
        COMPLEX(wp), INTENT(OUT) :: rate(2)

        rate = CMPLX(0.0_wp, 0.0_wp, wp)

        SELECT CASE(TRIM(model_name))
        CASE('rheology_voellmy')
            CALL eval_rheology_voellmy_implicit(state, mod_vel, mod_vel2, moving, use_stochastic, stochastic_value, &
                    rate(1), rate(2))
        CASE('rheology_mui_voellmy')
            CALL eval_rheology_mui_voellmy_implicit(state, mod_vel, mod_vel2, moving, use_stochastic, &
                    stochastic_value, rate(1), rate(2))
        CASE('rheology_plastic')
            CALL eval_rheology_plastic_implicit(state, mod_vel, moving, rate(1), rate(2))
        CASE('rheology_lava')
            CALL eval_rheology_lava_implicit(state, moving, rate(1), rate(2))
        CASE('rheology_lahar')
            CALL eval_rheology_lahar_implicit(state, mod_vel, mod_hor_vel, moving, rate(1), rate(2))
        CASE('rheology_empirical')
            CALL eval_rheology_empirical_implicit(state, mod_hor_vel, mod_vel, moving, rate(1), rate(2))
        CASE('rheology_newmodel')
            CALL eval_rheology_newmodel_implicit(state,rate(1),rate(2))
        CASE('rheology_quadratic')
            CALL eval_rheology_quadratic_implicit(state, mod_hor_vel, mod_vel2, moving, rate(1), rate(2))
        CASE('none', 'rheology_froude', 'rheology_lucas', 'rheology_mui')
            CONTINUE
        CASE DEFAULT
            CALL fatal_error('Invalid rheology name in implicit momentum-rate dispatch')

        END SELECT
    END SUBROUTINE eval_rheology_implicit

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
    SUBROUTINE solve_rheology_analytic(model_name, dt, state, speed_ratio, grav_coeff, curvature_acceleration, &
            use_pore_pressure, use_stochastic, stochastic_value, velocity_after, model_supported)

        CHARACTER(*), INTENT(IN) :: model_name
        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: dt, speed_ratio, grav_coeff, curvature_acceleration, stochastic_value
        LOGICAL, INTENT(IN) :: use_pore_pressure, use_stochastic
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        LOGICAL, INTENT(OUT) :: model_supported
        REAL(wp) :: velocity_before(2)

        velocity_before = [state%u, state%v]
        velocity_after = velocity_before
        model_supported = .TRUE.

        SELECT CASE(TRIM(model_name))
        CASE ('none')
            CONTINUE
        CASE ('rheology_voellmy')
            CALL solve_rheology_voellmy_analytic(dt, state, speed_ratio, grav_coeff, curvature_acceleration, &
                    use_pore_pressure, use_stochastic, stochastic_value, velocity_after)
        CASE ('rheology_plastic')
            CALL solve_rheology_plastic_analytic(dt, state, speed_ratio, velocity_after)
        CASE ('rheology_lava')
            CALL solve_rheology_lava_analytic(dt, state, speed_ratio, velocity_after)
        CASE ('rheology_lahar')
            CALL solve_rheology_lahar_analytic(dt, state, speed_ratio, velocity_after)
        CASE ('rheology_empirical')
            CALL solve_rheology_empirical_analytic(dt, state, speed_ratio, velocity_after)
        CASE ('rheology_quadratic')
            CALL solve_rheology_quadratic_analytic(dt, state, speed_ratio, velocity_after)
        CASE ('rheology_froude')
            CALL solve_rheology_froude_analytic(dt, state, speed_ratio, grav_coeff, curvature_acceleration, &
                    use_stochastic, stochastic_value, velocity_after)
        CASE ('rheology_lucas')
            CALL solve_rheology_lucas_analytic(dt, state, speed_ratio, grav_coeff, curvature_acceleration, &
                    use_stochastic, stochastic_value, velocity_after)
        CASE ('rheology_mui')
            CALL solve_rheology_mui_analytic(dt, state, speed_ratio, grav_coeff, curvature_acceleration, &
                    use_pore_pressure, velocity_after)
        CASE ('rheology_mui_voellmy')
            CALL solve_rheology_mui_voellmy_analytic(dt, state, speed_ratio, grav_coeff, curvature_acceleration, &
                    use_pore_pressure, use_stochastic, stochastic_value, velocity_after)
        CASE DEFAULT
            model_supported = .FALSE.

        END SELECT
    END SUBROUTINE solve_rheology_analytic

    !> \brief Compute curvature-induced normal acceleration.
    !> \details Contract the topographic Hessian with the horizontal velocity twice: Bxx*u**2 +
    !> 2*Bxy*u*v + Byy*v**2.
    !> \param[in] Bsecond_xx Second x derivative of topography [m^-1].
    !> \param[in] Bsecond_xy Mixed derivative of topography [m^-1].
    !> \param[in] Bsecond_yy Second y derivative of topography [m^-1].
    !> \param[in] u x velocity [m s^-1].
    !> \param[in] v y velocity [m s^-1].
    !> \return Computed acceleration for the requested model and inputs.
    PURE FUNCTION curvature_acceleration(Bsecond_xx, Bsecond_xy, Bsecond_yy,                       &
            u, v) RESULT(acceleration)

        REAL(wp), INTENT(IN) :: Bsecond_xx, Bsecond_xy, Bsecond_yy
        REAL(wp), INTENT(IN) :: u, v
        REAL(wp) :: acceleration

        acceleration = Bsecond_xx*u**2 +                                                           &
                2.0_wp*Bsecond_xy*u*v + Bsecond_yy*v**2

    END FUNCTION curvature_acceleration


    ! Additive rates from all active rheology processes at the same stage state.
    SUBROUTINE eval_rheology_explicit_processes(models,state,mod_vel,mod_hor_vel,mod_vel2, &
         grav_coeff,curvature_term,B_second_xx,B_second_xy,B_second_yy, &
         use_pore_pressure,moving,use_stochastic,stochastic_value,rate)
        CHARACTER(*), INTENT(IN) :: models(:)
        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: mod_vel,mod_hor_vel,mod_vel2,grav_coeff
        REAL(wp), INTENT(IN) :: B_second_xx,B_second_xy,B_second_yy,stochastic_value
        LOGICAL, INTENT(IN) :: curvature_term,use_pore_pressure,moving,use_stochastic
        REAL(wp), INTENT(OUT) :: rate(2)
        REAL(wp) :: process_rate(2)
        INTEGER :: i
        rate=0.0_wp
        DO i=1,SIZE(models)
            CALL eval_rheology_explicit(TRIM(models(i)),state,mod_vel,mod_hor_vel,mod_vel2, &
                 grav_coeff,curvature_term,B_second_xx,B_second_xy,B_second_yy, &
                 use_pore_pressure,moving,use_stochastic,stochastic_value,process_rate)
            rate=rate+process_rate
        END DO
    END SUBROUTINE eval_rheology_explicit_processes

    SUBROUTINE eval_rheology_semi_implicit_processes(models,state,mod_vel,mod_hor_vel, &
         grav_coeff,curvature_term,B_second_xx,B_second_xy,B_second_yy, &
         use_pore_pressure,use_stochastic,stochastic_value,rate)
        CHARACTER(*), INTENT(IN) :: models(:)
        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: mod_vel,mod_hor_vel,grav_coeff
        REAL(wp), INTENT(IN) :: B_second_xx,B_second_xy,B_second_yy,stochastic_value
        LOGICAL, INTENT(IN) :: curvature_term,use_pore_pressure,use_stochastic
        REAL(wp), INTENT(OUT) :: rate(2)
        REAL(wp) :: process_rate(2)
        INTEGER :: i
        rate=0.0_wp
        DO i=1,SIZE(models)
            CALL eval_rheology_semi_implicit(TRIM(models(i)),state,mod_vel,mod_hor_vel, &
                 grav_coeff,curvature_term,B_second_xx,B_second_xy,B_second_yy, &
                 use_pore_pressure,use_stochastic,stochastic_value,process_rate)
            rate=rate+process_rate
        END DO
    END SUBROUTINE eval_rheology_semi_implicit_processes

    SUBROUTINE eval_rheology_implicit_processes(models,state,mod_vel,mod_hor_vel,mod_vel2, &
         moving,use_stochastic,stochastic_value,rate)
        CHARACTER(*), INTENT(IN) :: models(:)
        TYPE(complex_physical_state_type), INTENT(IN) :: state
        COMPLEX(wp), INTENT(IN) :: mod_vel,mod_hor_vel,mod_vel2
        LOGICAL, INTENT(IN) :: moving,use_stochastic
        REAL(wp), INTENT(IN) :: stochastic_value
        COMPLEX(wp), INTENT(OUT) :: rate(2)
        COMPLEX(wp) :: process_rate(2)
        INTEGER :: i
        rate=CMPLX(0.0_wp,0.0_wp,wp)
        DO i=1,SIZE(models)
            CALL eval_rheology_implicit(TRIM(models(i)),state,mod_vel,mod_hor_vel,mod_vel2, &
                 moving,use_stochastic,stochastic_value,process_rate)
            rate=rate+process_rate
        END DO
    END SUBROUTINE eval_rheology_implicit_processes

    ! Exact flows are composed in order, never added as rates.
    SUBROUTINE solve_rheology_analytic_processes(models,dt,state,speed_ratio,grav_coeff, &
         curvature_acceleration,use_pore_pressure,use_stochastic,stochastic_value, &
         velocity_after,model_supported)
        CHARACTER(*), INTENT(IN) :: models(:)
        REAL(wp), INTENT(IN) :: dt,speed_ratio,grav_coeff,curvature_acceleration,stochastic_value
        TYPE(physical_state_type), INTENT(IN) :: state
        LOGICAL, INTENT(IN) :: use_pore_pressure,use_stochastic
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        LOGICAL, INTENT(OUT) :: model_supported
        TYPE(physical_state_type) :: current_state
        INTEGER :: i
        current_state=state
        velocity_after=[state%u,state%v]
        model_supported=.TRUE.
        DO i=1,SIZE(models)
            CALL solve_rheology_analytic(TRIM(models(i)),dt,current_state,speed_ratio, &
                 grav_coeff,curvature_acceleration,use_pore_pressure,use_stochastic, &
                 stochastic_value,velocity_after,model_supported)
            IF (.NOT.model_supported) RETURN
            current_state%u=velocity_after(1)
            current_state%v=velocity_after(2)
        END DO
    END SUBROUTINE solve_rheology_analytic_processes

END MODULE rheology_2d
