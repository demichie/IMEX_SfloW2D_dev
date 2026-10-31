!********************************************************************************
!> \brief Voellmy-Salm basal-friction model.
!> \details The model combines a Coulomb basal contribution and a
!> velocity-squared turbulent contribution. The implicit turbulent source is
!> \f[
!> \boldsymbol{S}_t=-\frac{\rho_m g'}{\xi}|\boldsymbol{u}|^2
!>                  \frac{\boldsymbol{u}_h}{|\boldsymbol{u}|}.
!> \f]
!> The semi-implicit Coulomb stress is
!> \f[
!> \tau_b=\max\!\left(0,\,
!> \mu\rho_m h(g'+a_c)\sqrt{G}
!> -\mathbf{1}_{p}\mu\sqrt{G}\,p_e\right),\qquad
!> \boldsymbol{S}_b=-\tau_b\boldsymbol{u}_h/|\boldsymbol{u}|,
!> \f]
!> where \f$a_c\f$ is the optional curvature acceleration and the pore
!> pressure subtraction is enabled only when requested. The module provides
!> \ref eval_rheology_voellmy_implicit for the turbulent term and
!> \ref eval_rheology_voellmy_semi_implicit for the Coulomb term; \ref mu and
!> \ref xi are their respective friction parameters. When enabled,
!> \f$a_c=c_{xx}u_x^2+2c_{xy}u_xu_y+c_{yy}u_y^2\f$.
!> \ref solve_rheology_voellmy_analytic exactly integrates the frozen
!> combined velocity source while holding non-velocity inputs fixed.
!********************************************************************************
MODULE rheology_voellmy_2d
    USE diagnostics_2d, ONLY : fatal_error
    USE parameters_2d, ONLY : wp
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type, complex_state_from_real
    USE rheology_parameters_2d, ONLY : xi
    USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs,                      &
            fatal_analytic_input_error, solve_exact_velocity_flow,                                 &
            EXACT_FLOW_QUADRATIC_POLYNOMIAL
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_rheology_voellmy_explicit
    PUBLIC :: eval_rheology_voellmy_semi_implicit
    PUBLIC :: eval_rheology_voellmy_implicit
    PUBLIC :: solve_rheology_voellmy_analytic
    PUBLIC :: mu
    REAL(wp) :: mu
    ! Uniform source treatment entry points. Existing kernels remain unchanged.
CONTAINS
    !> \param[in] state Shared physical state (complex for the implicit contribution).
    !> \brief Evaluate the explicit voellmy contribution to the horizontal momentum source.
    !> \details Compute the momentum source under the Voellmy Coulomb and turbulent basal resistance
    !> law, using the treatment-specific velocity and effective-stress inputs. Outputs are
    !> conservative x/y momentum source rates.
    !> \param[in] mod_vel Three-dimensional speed [m s^-1].
    !> \param[in] mod_vel2 Squared three-dimensional speed [m^2 s^-2].
    !> \param[in] grav_coeff Local dimensionless gravitational correction.
    !> \param[in] use_pore_pressure Enable pore-pressure correction.
    !> \param[in] moving True for a moving local state.
    !> \param[in] use_stochastic Enable model stochastic parameter.
    !> \param[in] stochastic_value Stochastic parameter value.
    !> \param[out] rate Two-component conservative momentum source.
    SUBROUTINE eval_rheology_voellmy_explicit(state, grav_coeff, use_pore_pressure, use_stochastic, stochastic_value, rate)
        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: grav_coeff
        REAL(wp), INTENT(IN) :: stochastic_value
        LOGICAL, INTENT(IN) :: use_pore_pressure, use_stochastic
        REAL(wp), INTENT(OUT) :: rate(2)
        REAL(wp) :: sx, sy
        COMPLEX(wp) :: ix, iy
        TYPE(complex_physical_state_type) :: cs

        CALL eval_rheology_voellmy_semi_implicit(state, &
             grav_coeff, &
             use_pore_pressure, &
             sx, &
             sy)
        CALL complex_state_from_real(state, cs)
        CALL eval_rheology_voellmy_implicit(cs, &
             use_stochastic, &
             &
                stochastic_value, &
             ix, &
             iy)
        rate=[sx,sy]+REAL([ix,iy],wp)
    END SUBROUTINE eval_rheology_voellmy_explicit

    !> \param[in] state Shared physical state (complex for the implicit contribution).
    !> \brief Evaluate the semi-implicit voellmy contribution to the horizontal momentum source.
    !> \details Return the contribution assigned to the semi-implicit treatment. This implementation
    !> supplies a zero source when the model has no term in this treatment.
    !> \param[in] mod_vel Three-dimensional speed [m s^-1].
    !> \param[in] grav_coeff Local dimensionless gravitational correction.
    !> \param[in] use_pore_pressure Enable pore-pressure correction.
    !> \param[out] rate_x x momentum source.
    !> \param[out] rate_y y momentum source.
    SUBROUTINE eval_rheology_voellmy_semi_implicit(state, grav_coeff, use_pore_pressure, rate_x, rate_y)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: grav_coeff
        LOGICAL, INTENT(IN) :: use_pore_pressure

        REAL(wp), INTENT(OUT) :: rate_x, rate_y
        REAL(wp) :: friction_stress

        ASSOCIATE(h => state%h, &
                  u => state%u, &
                  v => state%v, &
                  rhom => state%rho_m, &
                  red_grav => state%red_grav, &
                  exc_pore_pres => state%exc_pore_pres)
        rate_x = 0.0_wp
        rate_y = 0.0_wp
        IF (state%mod_vel .LE. 0.0_wp) RETURN

        friction_stress = mu * rhom *                                                              &
                (red_grav + state%curvature_acceleration) * h                                            &
                * SQRT(grav_coeff)
        friction_stress = MAX(0.0_wp, friction_stress)
        IF (use_pore_pressure) friction_stress = MAX(0.0_wp, friction_stress                       &
                - mu * SQRT(grav_coeff) * exc_pore_pres)
        rate_x = -friction_stress * u / state%mod_vel
        rate_y = -friction_stress * v / state%mod_vel
    
        END ASSOCIATE
    END SUBROUTINE eval_rheology_voellmy_semi_implicit

    !> \param[in] state Shared physical state (complex for the implicit contribution).
    !> \brief Evaluate the implicit voellmy contribution to the horizontal momentum source.
    !> \details Return the contribution assigned to the implicit treatment. This implementation
    !> supplies a zero source when the model has no term in this treatment.
    !> \param[in] mod_vel Three-dimensional speed [m s^-1].
    !> \param[in] mod_vel2 Squared three-dimensional speed [m^2 s^-2].
    !> \param[in] moving True for a moving local state.
    !> \param[in] use_stochastic_xi Enable stochastic turbulent resistance.
    !> \param[in] stochastic_value Stochastic parameter value.
    !> \param[out] rate_x x momentum source.
    !> \param[out] rate_y y momentum source.
    SUBROUTINE eval_rheology_voellmy_implicit(state, use_stochastic_xi, stochastic_value, rate_x, rate_y)

        TYPE(complex_physical_state_type), INTENT(IN) :: state

        LOGICAL, INTENT(IN) :: use_stochastic_xi
        REAL(wp), INTENT(IN) :: stochastic_value
        COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y
        REAL(wp) :: xi_effective
        COMPLEX(wp) :: turbulent_stress

        ASSOCIATE(u => state%u, &
                  v => state%v, &
                  rhom => state%rho_m, &
                  red_grav => state%red_grav)
        rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
        rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
        IF (.NOT. ((REAL(state%mod_vel,wp)/=0.0_wp).OR.(AIMAG(state%mod_vel)/=0.0_wp))) RETURN
        IF (use_stochastic_xi) THEN
            xi_effective = xi + stochastic_value
            IF (xi_effective .LT. 1.0_wp) xi_effective = 1.0_wp
            IF (xi_effective .GT. 1.0E5_wp) xi_effective = 1.0E5_wp

        ELSE
            xi_effective = xi

        END IF

        turbulent_stress = rhom * red_grav / xi_effective * state%mod_vel2
        rate_x = -turbulent_stress * (u / state%mod_vel)
        rate_y = -turbulent_stress * (v / state%mod_vel)
    
        END ASSOCIATE
    END SUBROUTINE eval_rheology_voellmy_implicit

    !> \param[in] state Shared real physical state before source integration.
    !> \brief Integrate the voellmy velocity source in time.
    !> \details Advance the scalar speed under the Voellmy Coulomb and turbulent basal resistance
    !> law with model coefficients fixed over dt; the horizontal direction is retained in the
    !> returned velocity vector.
    !> \param[in] dt Source integration interval [s].
    !> \param[in] speed_ratio Ratio of total to horizontal speed.
    !> \param[in] grav_coeff Local dimensionless gravitational correction.
    !> \param[in] curvature_acceleration Normal acceleration due to topographic curvature [m s^-2].
    !> \param[in] use_pore_pressure Enable pore-pressure correction.
    !> \param[in] use_stochastic_xi Enable stochastic turbulent resistance.
    !> \param[in] stochastic_xi_value Stochastic turbulent resistance parameter.
    !> \param[out] velocity_after Horizontal velocity vector after integration [m s^-1].
    SUBROUTINE solve_rheology_voellmy_analytic(dt, state, speed_ratio, grav_coeff, curvature_acceleration, &
            use_pore_pressure, use_stochastic_xi, stochastic_xi_value, velocity_after)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: dt, speed_ratio, grav_coeff, curvature_acceleration
        LOGICAL, INTENT(IN) :: use_pore_pressure, use_stochastic_xi
        REAL(wp), INTENT(IN) :: stochastic_xi_value
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        REAL(wp) :: velocity_before(2)
        REAL(wp) :: speed_before, speed_after, xi_effective
        REAL(wp) :: coulomb_stress, constant_drag, quadratic_drag

        velocity_before = [state%u, state%v]
        ASSOCIATE(h => state%h, &
                  rhom => state%rho_m, &
                  red_grav => state%red_grav, &
                  exc_pore_pres => state%exc_pore_pres)
        velocity_after = velocity_before
        speed_before = SQRT(SUM(velocity_before**2))
        CALL validate_analytic_velocity_inputs(dt, h,                                              &
                speed_ratio, speed_before, 'Voellmy')
        IF (dt .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
        IF (rhom .LE. 0.0_wp) THEN
            CALL fatal_analytic_input_error(                                                       &
                'Nonpositive mixture density in Voellmy analytic solver')

        END IF

        IF (grav_coeff .LE. 0.0_wp) THEN
            CALL fatal_analytic_input_error(                                                       &
                'Nonpositive gravity multiplier in Voellmy analytic solver')

        END IF

        xi_effective = xi
        IF (use_stochastic_xi) THEN
            xi_effective = xi + stochastic_xi_value
            xi_effective = MIN(1.0E5_wp, MAX(1.0_wp, xi_effective))

        ELSEIF (xi_effective .LE. 0.0_wp) THEN
            CALL fatal_analytic_input_error('Nonpositive xi in Voellmy analytic solver')

        END IF

        coulomb_stress = mu * rhom * (red_grav +                                                   &
                curvature_acceleration) * h * SQRT(grav_coeff)
        coulomb_stress = MAX(0.0_wp, coulomb_stress)
        IF (use_pore_pressure) coulomb_stress = MAX(0.0_wp, coulomb_stress -                       &
                mu * SQRT(grav_coeff) * exc_pore_pres)
        constant_drag = coulomb_stress / (rhom * h * speed_ratio)
        quadratic_drag = red_grav * speed_ratio / (xi_effective*h)
        CALL solve_exact_velocity_flow(EXACT_FLOW_QUADRATIC_POLYNOMIAL,                            &
                [constant_drag, 0.0_wp, quadratic_drag], dt,                                       &
                speed_before, speed_after)
        velocity_after = (speed_after / speed_before) * velocity_before
    
        END ASSOCIATE
    END SUBROUTINE solve_rheology_voellmy_analytic

    ! Unsupported physical treatment, rejected by the source catalog.
END MODULE rheology_voellmy_2d