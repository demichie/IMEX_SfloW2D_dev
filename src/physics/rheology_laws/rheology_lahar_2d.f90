!********************************************************************************
!> \brief Lahar rheology: temperature-dependent viscosity, turbulent drag and yield stress.
!> \details With \f$T_C=T-273.15\f$, the carrier-fluid viscosity is
!> \f$\eta_f(T,C)=\alpha_1(T)\exp(\beta_1 C)\f$. Its piecewise temperature
!> factor is
!> \f[
!> \alpha_1(T)=\alpha_{1,\rm coeff}
!> \begin{cases}
!> 10^{-3}10^{1301/[998.333+8.1855(T_C-20)+0.00585(T_C-20)^2]-1.30223},
!>       & T_C<20,\\
!> 1.002\,10^{-3}10^{[1.3272(20-T_C)-0.001053(T_C-20)^2]/(T_C+105)},
!>       & T_C\ge20.
!> \end{cases}
!> \f]
!> The implicit momentum source
!> combines viscous and turbulent slope contributions,
!> \f[
!> S_v=\frac{\kappa\eta_f|\boldsymbol{u}|}
!>              {8\rho_m g h^2},\qquad
!> S_t=n_{td}^2|\boldsymbol{u}|^2h^{-4/3},\qquad
!> \boldsymbol{S}=-\frac{g\rho_mh(S_v+S_t)}
!>                       {|\boldsymbol{u}_h|}\boldsymbol{u}_h.
!> \f]
!> The semi-implicit yield term uses
!> \f$Y=\alpha_2[\exp(\beta_2 C)-1]\f$ and
!> \f$\boldsymbol{S}_y=-Y\boldsymbol{u}_h/|\boldsymbol{u}_h|\f$. The
!> \ref eval_rheology_lahar_implicit evaluates the viscous/turbulent drag and
!> \ref eval_rheology_lahar_semi_implicit evaluates the yield-stress term. The
!> real and complex overloads of \ref eval_lahar_fluid_viscosity implement
!> the same fluid-viscosity correlation.
!> \ref solve_rheology_lahar_analytic integrates the frozen polynomial
!> velocity source with depth, T and concentration explicit.
!********************************************************************************
MODULE rheology_lahar_2d
    USE diagnostics_2d, ONLY : fatal_error
    USE parameters_2d, ONLY : wp, sp, four_thirds, neg_four_thirds
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
    USE constitutive_parameters_2d, ONLY : grav, inv_grav
    USE rheology_parameters_2d, ONLY : T_ref
    USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs,                      &
            fatal_analytic_input_error, solve_exact_velocity_flow,                                 &
            EXACT_FLOW_QUADRATIC_POLYNOMIAL
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_rheology_lahar_explicit
    PUBLIC :: eval_rheology_lahar_semi_implicit
    PUBLIC :: eval_rheology_lahar_implicit
    PUBLIC :: solve_rheology_lahar_analytic
    PUBLIC :: alpha2, beta2, alpha1_coeff, beta1, Kappa, n_td, n_td2
    PUBLIC :: eval_lahar_fluid_viscosity
    REAL(wp) :: alpha2
    REAL(wp) :: beta2
    REAL(wp) :: alpha1_coeff
    REAL(wp) :: beta1
    REAL(wp) :: Kappa
    REAL(wp) :: n_td
    REAL(wp) :: n_td2
    INTERFACE eval_lahar_fluid_viscosity
        MODULE PROCEDURE eval_lahar_fluid_viscosity_real
        MODULE PROCEDURE eval_lahar_fluid_viscosity_complex
    END INTERFACE eval_lahar_fluid_viscosity

CONTAINS
    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the explicit lahar momentum contribution.
    !> \details Assemble the real-valued momentum source from the model's supported physical
    !> treatment contributions; unused model inputs are retained for the shared dispatch interface.
    !> \param[in] mod_vel Magnitude of three-dimensional velocity [m s^-1].
    !> \param[in] mod_hor_vel Magnitude of horizontal velocity [m s^-1].
    !> \param[in] moving Whether the local state is considered moving.
    !> \param[out] rate Conservative horizontal momentum source vector.
    SUBROUTINE eval_rheology_lahar_explicit(state, &
            mod_vel, &
            mod_hor_vel, &
            moving, &
            rate)
        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: mod_vel, mod_hor_vel
        LOGICAL, INTENT(IN) :: moving
        REAL(wp), INTENT(OUT) :: rate(2)
        REAL(wp) :: sx, sy
        COMPLEX(wp) :: ix, iy
        TYPE(complex_physical_state_type) :: cs
        CALL eval_rheology_lahar_semi_implicit(state, mod_hor_vel, moving, sx, sy)
        cs%h=CMPLX(state%h,0.0_wp,wp)
        cs%u=CMPLX(state%u,0.0_wp,wp)
        cs%v=CMPLX(state%v,0.0_wp,wp)
        cs%T=CMPLX(state%T,0.0_wp,wp)
        cs%rho_m=CMPLX(state%rho_m,0.0_wp,wp)
        cs%inv_rhom=CMPLX(state%inv_rhom,0.0_wp,wp)
        cs%alphas=CMPLX(state%alphas,0.0_wp,wp)
        CALL eval_rheology_lahar_implicit(cs, CMPLX(mod_vel,0.0_wp,wp), &
             CMPLX(mod_hor_vel,0.0_wp,wp), moving, ix, iy)
        rate=[sx,sy]+REAL([ix,iy],wp)
    END SUBROUTINE eval_rheology_lahar_explicit

    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the semi-implicit lahar momentum contribution.
    !> \details Evaluate the contribution of temperature-dependent carrier viscosity, turbulent
    !> resistance and a concentration-dependent yield stress assigned to this treatment. The
    !> returned rates act on the conservative horizontal momenta; zero-valued implementations
    !> intentionally provide no contribution.
    !> \param[in] mod_hor_vel Magnitude of horizontal velocity [m s^-1].
    !> \param[in] moving Whether the local state is considered moving.
    !> \param[out] rate_x x component of momentum source rate.
    !> \param[out] rate_y y component of momentum source rate.
    SUBROUTINE eval_rheology_lahar_semi_implicit(state, mod_hor_vel, moving, rate_x, rate_y)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: mod_hor_vel
        LOGICAL, INTENT(IN) :: moving
        REAL(wp), INTENT(OUT) :: rate_x, rate_y
        REAL(wp) :: yield_strength, yield_slope, friction_term
        REAL(wp), PARAMETER :: depth_threshold = 1.0E-20_wp

        REAL(wp) :: alphas_tot

        alphas_tot = SUM(state%alphas)
        rate_x = 0.0_wp
        rate_y = 0.0_wp
        IF (.NOT. moving) RETURN
        yield_strength = alpha2 * (EXP(beta2 * alphas_tot) - 1.0_wp)
        IF (state%h .GT. depth_threshold) THEN
            yield_slope = yield_strength / (grav * state%rho_m * state%h)

        ELSE
            yield_slope = yield_strength / (grav * state%rho_m * depth_threshold)

        END IF

        friction_term = grav * state%rho_m * state%h * yield_slope
        rate_x = -friction_term * state%u / mod_hor_vel
        rate_y = -friction_term * state%v / mod_hor_vel
    END SUBROUTINE eval_rheology_lahar_semi_implicit

    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the implicit lahar momentum contribution.
    !> \details Evaluate the part of temperature-dependent carrier viscosity, turbulent resistance
    !> and a concentration-dependent yield stress assigned to the implicit treatment. Complex-valued
    !> arguments retain complex-step sensitivities where implemented; zero-valued implementations
    !> intentionally contribute no term.
    !> \param[in] mod_vel Magnitude of three-dimensional velocity [m s^-1].
    !> \param[in] mod_hor_vel Magnitude of horizontal velocity [m s^-1].
    !> \param[in] moving Whether the local state is considered moving.
    !> \param[out] rate_x x component of momentum source rate.
    !> \param[out] rate_y y component of momentum source rate.
    SUBROUTINE eval_rheology_lahar_implicit(state, mod_vel, mod_hor_vel, moving, rate_x, rate_y)

        TYPE(complex_physical_state_type), INTENT(IN) :: state
        COMPLEX(wp), INTENT(IN) :: mod_vel, mod_hor_vel
        LOGICAL, INTENT(IN) :: moving
        COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y
        COMPLEX(wp) :: fluid_viscosity, inverse_depth, viscous_slope
        COMPLEX(wp) :: turbulent_slope, total_slope, friction_term
        REAL(wp) :: depth_threshold

        rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
        rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
        IF (wp .EQ. sp) THEN
            depth_threshold = 1.0E-10_wp

        ELSE
            depth_threshold = 1.0E-20_wp

        END IF

        fluid_viscosity = eval_lahar_fluid_viscosity(state%T, SUM(state%alphas))
        IF (REAL(state%h, wp) .GT. depth_threshold) THEN
            inverse_depth = 1.0_wp / state%h
            viscous_slope = Kappa * fluid_viscosity * mod_vel * 0.125_wp *                         &
                    state%inv_rhom * inv_grav * inverse_depth**2
            turbulent_slope = n_td2 * mod_vel**2 * inverse_depth**four_thirds

        ELSE
            viscous_slope = Kappa * fluid_viscosity * mod_vel /                                    &
                    (8.0_wp * state%rho_m * grav * depth_threshold**2)
            turbulent_slope = n_td2 * mod_vel**2 * depth_threshold**neg_four_thirds

        END IF

        total_slope = viscous_slope + turbulent_slope
        IF (.NOT. moving) RETURN
        friction_term = grav * state%rho_m * state%h * total_slope / mod_hor_vel
        rate_x = -state%u * friction_term
        rate_y = -state%v * friction_term
    END SUBROUTINE eval_rheology_lahar_implicit

    !> \param[in] state Shared real physical state before analytic integration.
    !> \brief Integrate the lahar velocity source over a frozen-state substep.
    !> \details Evaluate the exact or time-integral-based decay of horizontal speed under
    !> temperature-dependent carrier viscosity, turbulent resistance and a concentration-dependent
    !> yield stress. The velocity direction is retained and the updated components are returned
    !> without evolving the frozen thermodynamic and geometric inputs.
    !> \param[in] dt Duration of the analytic source substep [s].
    !> \param[in] speed_ratio Ratio of three-dimensional to horizontal speed.
    !> s^-1].
    !> \param[out] velocity_after Two-component horizontal velocity after source integration [m
    !> s^-1].
    SUBROUTINE solve_rheology_lahar_analytic(dt, state, speed_ratio, velocity_after)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: dt, speed_ratio
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        REAL(wp) :: velocity_before(2), alphas_tot
        REAL(wp) :: speed_before, speed_after, effective_depth, depth_threshold
        REAL(wp) :: fluid_viscosity, yield_strength
        REAL(wp) :: constant_drag, linear_drag, quadratic_drag

        velocity_before = [state%u,state%v]
        alphas_tot = SUM(state%alphas)
        velocity_after = velocity_before
        speed_before = SQRT(SUM(velocity_before**2))
        CALL validate_analytic_velocity_inputs(dt, state%h,                                              &
                speed_ratio, speed_before, 'lahar')
        IF (dt .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
        IF (state%rho_m .LE. 0.0_wp) THEN
            CALL fatal_analytic_input_error('Nonpositive mixture density in lahar analytic solver')

        END IF

        IF (wp .EQ. sp) THEN
            depth_threshold = 1.0E-10_wp

        ELSE
            depth_threshold = 1.0E-20_wp

        END IF

        effective_depth = MAX(state%h, depth_threshold)
        fluid_viscosity = eval_lahar_fluid_viscosity(state%T, alphas_tot)
        yield_strength = alpha2 * (EXP(beta2 * alphas_tot) - 1.0_wp)
        constant_drag = yield_strength / (state%rho_m *                                                   &
                MAX(state%h, 1.0E-20_wp))
        linear_drag = Kappa * fluid_viscosity * speed_ratio /                                      &
                (8.0_wp * state%rho_m * effective_depth**2)
        quadratic_drag = grav * n_td2 * speed_ratio**2 *                                           &
                effective_depth**neg_four_thirds
        CALL solve_exact_velocity_flow(EXACT_FLOW_QUADRATIC_POLYNOMIAL,                            &
                [constant_drag, linear_drag, quadratic_drag], dt,                                  &
                speed_before, speed_after)
        velocity_after = (speed_after / speed_before) * velocity_before
    END SUBROUTINE solve_rheology_lahar_analytic

    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the lahar carrier-fluid viscosity.
    !> \details Apply the piecewise temperature correlation and the exponential solid-concentration
    !> factor; the complex overload preserves complex-step derivatives.
    FUNCTION eval_lahar_fluid_viscosity_real(T, alphas_tot)                                        &
            RESULT(fluid_viscosity)

        REAL(wp), INTENT(IN) :: T, alphas_tot
        REAL(wp) :: fluid_viscosity
        COMPLEX(wp) :: complex_viscosity

        complex_viscosity = eval_lahar_fluid_viscosity_complex(                                    &
                CMPLX(T,0.0_wp,wp),CMPLX(alphas_tot,0.0_wp,wp))
        fluid_viscosity = REAL(complex_viscosity,wp)
    END FUNCTION eval_lahar_fluid_viscosity_real

    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the lahar carrier-fluid viscosity.
    !> \details Apply the piecewise temperature correlation and the exponential solid-concentration
    !> factor; the complex overload preserves complex-step derivatives.
    FUNCTION eval_lahar_fluid_viscosity_complex(T, alphas_tot)                                     &
            RESULT(fluid_viscosity)

        COMPLEX(wp), INTENT(IN) :: T, alphas_tot
        COMPLEX(wp) :: fluid_viscosity
        COMPLEX(wp) :: temperature_celsius, exponent_a, exponent_b, alpha1

        temperature_celsius = T - 273.15_wp
        IF (REAL(temperature_celsius, wp) .LT. 20.0_wp) THEN
            exponent_a = 1301.0_wp / (998.333_wp + 8.1855_wp *                                     &
                    (temperature_celsius - 20.0_wp) + 0.00585_wp *                                 &
                    (temperature_celsius - 20.0_wp)**2) - 1.30223_wp
            alpha1 = alpha1_coeff * 1.0E-3_wp * 10.0_wp**exponent_a

        ELSE
            exponent_b = (1.3272_wp * (20.0_wp - temperature_celsius)                              &
                    - 0.001053_wp * (temperature_celsius - 20.0_wp)**2)                            &
                    / (temperature_celsius + 105.0_wp)
            alpha1 = alpha1_coeff * 1.002E-3_wp * 10.0_wp**exponent_b

        END IF

        fluid_viscosity = alpha1 * EXP(beta1 * alphas_tot)
    END FUNCTION eval_lahar_fluid_viscosity_complex

END MODULE rheology_lahar_2d
