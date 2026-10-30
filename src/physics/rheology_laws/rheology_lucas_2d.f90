!********************************************************************************
!> \brief Lucas velocity-weakening Coulomb basal-friction model (model 10).
!> \details For friction mod_vel \f$U_f\f$, the coefficient is constant below
!> the transition mod_vel \f$U_w\f$ (and for shallow flow), and weakens above
!> that threshold:
!> \f[
!> \mu(U_f)=
!> \begin{cases}
!> \mu_0, & U_f<U_w\ \mathrm{or}\ h<0.1,\\
!> \mu_\infty+(\mu_0-\mu_\infty)U_w/U_f, & \mathrm{otherwise}.
!> \end{cases}
!> \f]
!> The source is \f$\boldsymbol{S}=-\rho_m\mu hG(g'+a_c)
!> \boldsymbol{u}_h/|\boldsymbol{u}_h|\f$. When stochastic mode is enabled,
!> its increment is added to the mod_vel before evaluating the law. The module
!> exposes \ref eval_rheology_lucas_semi_implicit for the x/y momentum rates and
!> \ref solve_rheology_lucas_analytic for exact integration of the frozen
!> speed-decay law with non-velocity variables explicit.
!********************************************************************************
MODULE rheology_lucas_2d
    USE diagnostics_2d, ONLY : fatal_error
    USE parameters_2d, ONLY : wp
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
    USE rheology_parameters_2d, ONLY : mu_0, mu_inf, U_w
    USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs,                      &
            fatal_analytic_input_error, solve_exact_velocity_flow,                                 &
            EXACT_FLOW_LUCAS
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_rheology_lucas_explicit
    PUBLIC :: eval_rheology_lucas_semi_implicit
    PUBLIC :: eval_rheology_lucas_implicit
    PUBLIC :: solve_rheology_lucas_analytic
    ! Uniform source treatment entry points. Existing kernels remain unchanged.
CONTAINS
    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the explicit lucas momentum contribution.
    !> \details Assemble the real-valued momentum source from the model's supported physical
    !> treatment contributions; unused model inputs are retained for the shared dispatch interface.
    !> \param[in] mod_vel Magnitude of three-dimensional velocity [m s^-1].
    !> \param[in] mod_hor_vel Magnitude of horizontal velocity [m s^-1].
    !> \param[in] mod_vel2 Squared three-dimensional speed [m^2 s^-2].
    !> \param[in] grav_coeff Dimensionless gravity correction for local topography.
    !> \param[in] curvature_term Whether topographic curvature enters the source.
    !> \param[in] B_second_xx Second x derivative of topography [m^-1].
    !> \param[in] B_second_xy Mixed derivative of topography [m^-1].
    !> \param[in] B_second_yy Second y derivative of topography [m^-1].
    !> \param[in] use_pore_pressure Whether excess pore pressure modifies resistance.
    !> \param[in] moving Whether the local state is considered moving.
    !> \param[in] use_stochastic Whether stochastic rheology parameters are enabled.
    !> \param[in] stochastic_value Stochastic model parameter used in the friction law.
    !> \param[out] rate Conservative horizontal momentum source vector.
    SUBROUTINE eval_rheology_lucas_explicit(state, mod_vel, mod_hor_vel, mod_vel2, &
            grav_coeff, curvature_term, B_second_xx, B_second_xy, B_second_yy, &
            use_pore_pressure, moving, use_stochastic, stochastic_value, rate)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: mod_vel, mod_hor_vel, mod_vel2, grav_coeff
        REAL(wp), INTENT(IN) :: B_second_xx, B_second_xy, B_second_yy, stochastic_value
        LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure, moving, use_stochastic
        REAL(wp), INTENT(OUT) :: rate(2)
        REAL(wp) :: sx, sy, curvature_accel
        COMPLEX(wp) :: ix, iy
        TYPE(complex_physical_state_type) :: cs

        curvature_accel=0.0_wp

        IF (curvature_term) curvature_accel=B_second_xx*state%u**2 +                                     &
                2.0_wp*B_second_xy*state%u*state%v+B_second_yy*state%v**2
        CALL eval_rheology_lucas_semi_implicit(state, mod_vel, mod_hor_vel, grav_coeff, &
            0.0_wp, use_stochastic, stochastic_value, sx, sy)
        CALL eval_rheology_lucas_implicit(cs, ix, iy)
        rate=[sx,sy]+REAL([ix,iy],wp)
    END SUBROUTINE eval_rheology_lucas_explicit

    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the semi-implicit lucas momentum contribution.
    !> \details Evaluate the contribution of the Lucas basal-friction parameterization assigned to
    !> this treatment. The returned rates act on the conservative horizontal momenta; zero-valued
    !> implementations intentionally provide no contribution.
    !> \param[in] mod_vel Magnitude of three-dimensional velocity [m s^-1].
    !> \param[in] mod_hor_vel Magnitude of horizontal velocity [m s^-1].
    !> \param[in] grav_coeff Dimensionless gravity correction for local topography.
    !> \param[in] curvature_term Whether topographic curvature enters the source.
    !> \param[in] stochastic Model-specific input or output as defined by the constitutive law.
    !> \param[in] stochastic_value Stochastic model parameter used in the friction law.
    !> \param[out] rate_x x component of momentum source rate.
    !> \param[out] rate_y y component of momentum source rate.
    SUBROUTINE eval_rheology_lucas_semi_implicit(state, mod_vel, mod_hor_vel, grav_coeff, curvature_term, &
            stochastic, stochastic_value, rate_x, rate_y)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: mod_vel, mod_hor_vel, grav_coeff, curvature_term
        LOGICAL, INTENT(IN) :: stochastic
        REAL(wp), INTENT(IN) :: stochastic_value
        REAL(wp), INTENT(OUT) :: rate_x, rate_y
        REAL(wp) :: friction_speed
        REAL(wp) :: friction_coefficient
        REAL(wp) :: friction_stress

        rate_x = 0.0_wp
        rate_y = 0.0_wp
        IF ( mod_vel .GT. 0.0_wp ) THEN
            IF ( stochastic ) THEN
                friction_speed = mod_vel + stochastic_value

            ELSE
                friction_speed = mod_vel

            END IF

            IF ( friction_speed .LT. U_w ) THEN
                friction_coefficient = mu_0

            ELSE
                IF ( state%h .LT. 0.1_wp ) THEN
                        friction_coefficient = mu_0

                ELSE
                        friction_coefficient = mu_inf + (mu_0-mu_inf)/(friction_speed/U_w)

                END IF

            END IF

            friction_stress = state%rho_m * ( friction_coefficient * state%h *                                  &
                    grav_coeff * ( state%red_grav + curvature_term ) )
            rate_x = rate_x - friction_stress * state%u / mod_hor_vel
            rate_y = rate_y - friction_stress * state%v / mod_hor_vel

        END IF

    END SUBROUTINE eval_rheology_lucas_semi_implicit

    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the implicit lucas momentum contribution.
    !> \details Evaluate the part of the Lucas basal-friction parameterization assigned to the
    !> implicit treatment. Complex-valued arguments retain complex-step sensitivities where
    !> implemented; zero-valued implementations intentionally contribute no term.
    !> \param[out] rate_x x component of momentum source rate.
    !> \param[out] rate_y y component of momentum source rate.
    SUBROUTINE eval_rheology_lucas_implicit(state, rate_x, rate_y)

        TYPE(complex_physical_state_type), INTENT(IN) :: state
        COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y

        rate_x = CMPLX(0.0_wp,0.0_wp,wp)
        rate_y = CMPLX(0.0_wp,0.0_wp,wp)
    END SUBROUTINE eval_rheology_lucas_implicit

    !> \param[in] state Shared real physical state before analytic integration.
    !> \brief Integrate the lucas velocity source over a frozen-state substep.
    !> \details Evaluate the exact or time-integral-based decay of horizontal speed under the Lucas
    !> basal-friction parameterization. The velocity direction is retained and the updated
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
    SUBROUTINE solve_rheology_lucas_analytic(dt, state, speed_ratio, grav_coeff, curvature_term, &
            stochastic, stochastic_value, velocity_after)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: dt, speed_ratio, grav_coeff, curvature_term
        LOGICAL, INTENT(IN) :: stochastic
        REAL(wp), INTENT(IN) :: stochastic_value
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        REAL(wp) :: velocity_before(2)
        REAL(wp) :: speed_before, speed_after, shift, scale_factor, shallow

        velocity_before = [state%u,state%v]
        velocity_after = velocity_before
        speed_before = SQRT(SUM(velocity_before**2))
        CALL validate_analytic_velocity_inputs(dt, state%h,                                              &
                speed_ratio, speed_before, 'Lucas')
        IF (dt .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
        IF (U_w .LE. 0.0_wp) THEN
            CALL fatal_analytic_input_error('Nonpositive U_w in Lucas analytic solver')

        END IF

        scale_factor = grav_coeff * (state%red_grav + curvature_term)
        IF (scale_factor .LE. 0.0_wp) RETURN
        shift = 0.0_wp
        IF (stochastic) shift = stochastic_value
        shallow = 0.0_wp
        IF (state%h .LT. 0.1_wp) shallow = 1.0_wp
        CALL solve_exact_velocity_flow(EXACT_FLOW_LUCAS,                                           &
                [scale_factor, speed_ratio, shift, mu_0, mu_inf, U_w, shallow],                    &
                dt, speed_before, speed_after)
        velocity_after = (speed_after / speed_before) * velocity_before
    END SUBROUTINE solve_rheology_lucas_analytic

    ! Unsupported physical treatment, rejected by the source catalog.
    ! Valid zero component of the configured Newton treatment.
END MODULE rheology_lucas_2d