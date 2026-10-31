!********************************************************************************
!> \brief Froude-dependent Coulomb basal-friction model (model 9).
!> \details The friction coefficient weakens exponentially with the Froude
!> number, \f$\mathrm{Fr}=|\boldsymbol{u}|/\sqrt{g' h}\f$:
!> \f[
!> \mu(\mathrm{Fr})=\mu_\infty+(\mu_0-\mu_\infty)
!>                 \exp(-\mathrm{Fr}/\mathrm{Fr}_0).
!> \f]
!> The horizontal source is
!> \f$\boldsymbol{S}=-\rho_m\mu h G(g'+a_c)
!> \boldsymbol{u}_h/|\boldsymbol{u}_h|\f$, where \f$G\f$ is the supplied
!> gravity coefficient and \f$a_c\f$ is the curvature acceleration. Optional
!> stochastic perturbations shift the Froude number for sufficiently deep
!> flow and clamp it at zero. The module provides
!> \ref eval_rheology_froude_semi_implicit to evaluate the two horizontal momentum rates.
!> \ref eval_rheology_froude_analytic integrates the frozen Froude-dependent
!> velocity source in time.
!********************************************************************************
MODULE rheology_froude_2d
    USE diagnostics_2d, ONLY : fatal_error
    USE parameters_2d, ONLY : wp
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
    USE rheology_parameters_2d, ONLY : mu_0, mu_inf, Fr_0
    USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs,                      &
            fatal_analytic_input_error, solve_exact_velocity_flow,                                 &
            EXACT_FLOW_FROUDE
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_rheology_froude_explicit
    PUBLIC :: eval_rheology_froude_semi_implicit
    PUBLIC :: eval_rheology_froude_implicit
    PUBLIC :: eval_rheology_froude_analytic
    ! Uniform source treatment entry points. Existing kernels remain unchanged.
CONTAINS
    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the explicit froude momentum contribution.
    !> \details Assemble the real-valued momentum source from the model's supported physical
    !> treatment contributions; unused model inputs are retained for the shared dispatch interface.
    !> \param[in] mod_vel Magnitude of three-dimensional velocity [m s^-1].
    !> \param[in] mod_hor_vel Magnitude of horizontal velocity [m s^-1].
    !> \param[in] grav_coeff Dimensionless gravity correction for local topography.
    !> \param[in] use_stochastic Whether stochastic rheology parameters are enabled.
    !> \param[in] stochastic_value Stochastic model parameter used in the friction law.
    !> \param[out] rate Conservative horizontal momentum source vector.
    SUBROUTINE eval_rheology_froude_explicit(state, grav_coeff, use_stochastic, stochastic_value, rate)
        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: grav_coeff
        REAL(wp), INTENT(IN) :: stochastic_value
        LOGICAL, INTENT(IN) :: use_stochastic
        REAL(wp), INTENT(OUT) :: rate(2)
        ! The implicit component is identically zero for Froude rheology.
        CALL eval_rheology_froude_semi_implicit(state, &
             grav_coeff, &
             use_stochastic, &
             stochastic_value, &
             rate(1), &
             rate(2))
    END SUBROUTINE eval_rheology_froude_explicit

    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the semi-implicit froude momentum contribution.
    !> \details Evaluate the contribution of a Froude-dependent basal friction law assigned to this
    !> treatment. The returned rates act on the conservative horizontal momenta; zero-valued
    !> implementations intentionally provide no contribution.
    !> \param[in] mod_vel Magnitude of three-dimensional velocity [m s^-1].
    !> \param[in] mod_hor_vel Magnitude of horizontal velocity [m s^-1].
    !> \param[in] grav_coeff Dimensionless gravity correction for local topography.
    !> \param[in] curvature_term Whether topographic curvature enters the source.
    !> \param[in] stochastic Model-specific input or output as defined by the constitutive law.
    !> \param[in] stochastic_value Stochastic model parameter used in the friction law.
    !> \param[out] rate_x x component of momentum source rate.
    !> \param[out] rate_y y component of momentum source rate.
    SUBROUTINE eval_rheology_froude_semi_implicit(state, grav_coeff, stochastic, stochastic_value, rate_x, rate_y)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: grav_coeff
        LOGICAL, INTENT(IN) :: stochastic
        REAL(wp), INTENT(IN) :: stochastic_value
        REAL(wp), INTENT(OUT) :: rate_x, rate_y
        REAL(wp) :: froude_number
        REAL(wp) :: friction_coefficient
        REAL(wp) :: friction_stress

        rate_x = 0.0_wp
        rate_y = 0.0_wp
        IF ( state%mod_vel .GT. 0.0_wp .AND. state%red_grav .GT. 0.0_wp ) THEN
            froude_number = state%mod_vel / SQRT(state%red_grav * state%h)
            IF ( stochastic ) THEN
                IF ( state%h .GT. 0.01_wp ) froude_number = froude_number + stochastic_value
                IF ( froude_number .LT. 0.0_wp ) froude_number = 0.0_wp

            END IF

            friction_coefficient = mu_inf + (mu_0 - mu_inf) * EXP(-froude_number/Fr_0)
            friction_stress = state%rho_m * ( friction_coefficient * state%h *                                  &
                    grav_coeff * ( state%red_grav + state%curvature_acceleration ) )
            rate_x = rate_x - friction_stress * state%u / state%mod_vel_hor
            rate_y = rate_y - friction_stress * state%v / state%mod_vel_hor

        END IF

    END SUBROUTINE eval_rheology_froude_semi_implicit

    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the implicit froude momentum contribution.
    !> \details Evaluate the part of a Froude-dependent basal friction law assigned to the implicit
    !> treatment. Complex-valued arguments retain complex-step sensitivities where implemented;
    !> zero-valued implementations intentionally contribute no term.
    !> \param[out] rate_x x component of momentum source rate.
    !> \param[out] rate_y y component of momentum source rate.
    SUBROUTINE eval_rheology_froude_implicit(state, rate_x, rate_y)

        TYPE(complex_physical_state_type), INTENT(IN) :: state
        COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y

        rate_x = CMPLX(0.0_wp,0.0_wp,wp)
        rate_y = CMPLX(0.0_wp,0.0_wp,wp)
    END SUBROUTINE eval_rheology_froude_implicit

    !> \param[in] state Shared real physical state before analytic integration.
    !> \brief Integrate the froude velocity source over a frozen-state substep.
    !> \details Evaluate the exact or time-integral-based decay of horizontal speed under a
    !> Froude-dependent basal friction law. The velocity direction is retained and the updated
    !> components are returned without evolving the frozen thermodynamic and geometric inputs.
    !> \param[in] dt Duration of the analytic source substep [s].
    !> \param[in] speed_ratio Ratio of three-dimensional to horizontal speed.
    !> s^-1].
    !> \param[in] grav_coeff Dimensionless gravity correction for local topography.
    !> \param[in] curvature_term Whether topographic curvature enters the source.
    !> \param[in] stochastic Model-specific input or output as defined by the constitutive law.
    !> \param[in] stochastic_value Stochastic model parameter used in the friction law.
    !> \param[out] velocity_after Two-component horizontal velocity after source integration [m
    !> s^-1].
    SUBROUTINE eval_rheology_froude_analytic(dt, state, speed_ratio, grav_coeff, curvature_term, &
            stochastic, stochastic_value, velocity_after)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: dt, speed_ratio, grav_coeff, curvature_term
        LOGICAL, INTENT(IN) :: stochastic
        REAL(wp), INTENT(IN) :: stochastic_value
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        REAL(wp) :: velocity_before(2)
        REAL(wp) :: speed_before, speed_after, shift, base_drag
        REAL(wp) :: froude_slope, froude_value_slope

        velocity_before = [state%u,state%v]
        velocity_after = velocity_before
        speed_before = SQRT(SUM(velocity_before**2))
        CALL validate_analytic_velocity_inputs(dt, state%h,                                              &
                speed_ratio, speed_before, 'Froude')
        IF (dt .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
        IF (state%h .LE. 0.0_wp .OR. state%red_grav .LE. 0.0_wp .OR.                                           &
                speed_ratio .LE. 0.0_wp) RETURN
        IF (Fr_0 .LE. 0.0_wp) THEN
            CALL fatal_analytic_input_error('Nonpositive Fr_0 in Froude analytic solver')

        END IF

        base_drag = grav_coeff * (state%red_grav + curvature_term)
        IF (base_drag .LE. 0.0_wp) RETURN
        shift = 0.0_wp
        IF (stochastic .AND. state%h .GT. 0.01_wp) shift = stochastic_value
        froude_value_slope = speed_ratio / SQRT(state%red_grav*state%h)
        froude_slope = froude_value_slope / Fr_0
        CALL solve_exact_velocity_flow(EXACT_FLOW_FROUDE,                                          &
                [base_drag, froude_slope, shift/Fr_0, mu_0, mu_inf],                               &
                dt, speed_before, speed_after)
        velocity_after = (speed_after / speed_before) * velocity_before
    END SUBROUTINE eval_rheology_froude_analytic

    ! Unsupported physical treatment, rejected by the source catalog.
    ! Valid zero component of the configured Newton treatment.
END MODULE rheology_froude_2d