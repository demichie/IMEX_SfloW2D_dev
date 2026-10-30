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
    USE rheology_parameters_2d, ONLY : xi
    USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs, &
                fatal_analytic_input_error, solve_exact_velocity_flow,                &
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

    SUBROUTINE eval_rheology_voellmy_explicit( &
                h, T, u, v, mod_vel, mod_hor_vel, &
                mod_vel2, rhom, inv_rhom, red_grav, &
                grav_coeff, alphas, curvature_term, B_second_xx, B_second_xy, &
                B_second_yy, use_pore_pressure, exc_pore_pres, moving, use_stochastic, &
                stochastic_value, rate )
        REAL(wp), INTENT(IN) :: h, T, u, v
        REAL(wp), INTENT(IN) :: mod_vel, mod_hor_vel, mod_vel2
        REAL(wp), INTENT(IN) :: rhom, inv_rhom, red_grav
        REAL(wp), INTENT(IN) :: grav_coeff, alphas(:)
        REAL(wp), INTENT(IN) :: B_second_xx, B_second_xy, B_second_yy
        REAL(wp), INTENT(IN) :: exc_pore_pres, stochastic_value
        LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure, moving, use_stochastic
        REAL(wp), INTENT(OUT) :: rate(2)
        REAL(wp) :: sx, sy
        COMPLEX(wp) :: ix, iy
        CALL eval_rheology_voellmy_semi_implicit( &
                h, u, v, mod_vel, rhom, red_grav, &
                grav_coeff, curvature_term, B_second_xx, B_second_xy, B_second_yy, &
                use_pore_pressure, exc_pore_pres, sx, sy )
        CALL eval_rheology_voellmy_implicit( &
                CMPLX(u,0.0_wp,wp), CMPLX(v,0.0_wp,wp), &
                CMPLX(rhom,0.0_wp,wp), CMPLX(red_grav,0.0_wp,wp), &
                CMPLX(mod_vel,0.0_wp,wp), CMPLX(mod_vel2,0.0_wp,wp), moving, use_stochastic, &
                stochastic_value, ix, iy )
        rate=[sx,sy]+REAL([ix,iy],wp)
    END SUBROUTINE eval_rheology_voellmy_explicit

    SUBROUTINE eval_rheology_voellmy_semi_implicit(h, u, v, &
                mod_vel, rhom, red_grav, grav_coeff,            &
                curvature_term, B_second_xx, B_second_xy, B_second_yy,                &
                use_pore_pressure, exc_pore_pres, rate_x, rate_y)
        REAL(wp), INTENT(IN) :: h, u, v, mod_vel
        REAL(wp), INTENT(IN) :: rhom, red_grav, grav_coeff
        LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure
        REAL(wp), INTENT(IN) :: B_second_xx, B_second_xy, B_second_yy
        REAL(wp), INTENT(IN) :: exc_pore_pres
        REAL(wp), INTENT(OUT) :: rate_x, rate_y
        REAL(wp) :: curvature_acceleration, friction_stress
        rate_x = 0.0_wp
        rate_y = 0.0_wp
        IF (mod_vel .LE. 0.0_wp) RETURN
        curvature_acceleration = 0.0_wp
        IF (curvature_term) THEN
                curvature_acceleration = B_second_xx * u**2                &
                        + 2.0_wp * B_second_xy * u * v                &
                        + B_second_yy * v**2
        END IF

        friction_stress = mu * rhom *                                  &
                    (red_grav + curvature_acceleration) * h              &
                    * SQRT(grav_coeff)
        friction_stress = MAX(0.0_wp, friction_stress)
        IF (use_pore_pressure) friction_stress = MAX(0.0_wp, friction_stress      &
                    - mu * SQRT(grav_coeff) * exc_pore_pres)
        rate_x = -friction_stress * u / mod_vel
        rate_y = -friction_stress * v / mod_vel
    END SUBROUTINE eval_rheology_voellmy_semi_implicit

    SUBROUTINE eval_rheology_voellmy_implicit(u, v, rhom, &
                red_grav, mod_vel, mod_vel2, moving, use_stochastic_xi,          &
                stochastic_value, rate_x, rate_y)
        COMPLEX(wp), INTENT(IN) :: u, v
        COMPLEX(wp), INTENT(IN) :: rhom, red_grav
        COMPLEX(wp), INTENT(IN) :: mod_vel, mod_vel2
        LOGICAL, INTENT(IN) :: moving, use_stochastic_xi
        REAL(wp), INTENT(IN) :: stochastic_value
        COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y
        REAL(wp) :: xi_effective
        COMPLEX(wp) :: turbulent_stress
        rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
        rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
        IF (.NOT. moving) RETURN
        IF (use_stochastic_xi) THEN
                xi_effective = xi + stochastic_value
                IF (xi_effective .LT. 1.0_wp) xi_effective = 1.0_wp
                IF (xi_effective .GT. 1.0E5_wp) xi_effective = 1.0E5_wp
        ELSE
                xi_effective = xi
        END IF

        turbulent_stress = rhom * red_grav / xi_effective * mod_vel2
        rate_x = -turbulent_stress * (u / mod_vel)
        rate_y = -turbulent_stress * (v / mod_vel)
    END SUBROUTINE eval_rheology_voellmy_implicit

    SUBROUTINE solve_rheology_voellmy_analytic(dt, h, &
                speed_ratio, velocity_before, rhom, red_grav,        &
                grav_coeff, curvature_acceleration, use_pore_pressure,       &
                exc_pore_pres, use_stochastic_xi, stochastic_xi_value,          &
                velocity_after)
        REAL(wp), INTENT(IN) :: dt, h, speed_ratio
        REAL(wp), INTENT(IN) :: velocity_before(2), rhom
        REAL(wp), INTENT(IN) :: red_grav, grav_coeff
        REAL(wp), INTENT(IN) :: curvature_acceleration, exc_pore_pres
        LOGICAL, INTENT(IN) :: use_pore_pressure, use_stochastic_xi
        REAL(wp), INTENT(IN) :: stochastic_xi_value
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        REAL(wp) :: speed_before, speed_after, xi_effective
        REAL(wp) :: coulomb_stress, constant_drag, quadratic_drag
        velocity_after = velocity_before
        speed_before = SQRT(SUM(velocity_before**2))
        CALL validate_analytic_velocity_inputs(dt, h,    &
                    speed_ratio, speed_before, 'Voellmy')
        IF (dt .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
        IF (rhom .LE. 0.0_wp) THEN
                CALL fatal_analytic_input_error('Nonpositive mixture density in Voellmy analytic solver')
        END IF

        IF (grav_coeff .LE. 0.0_wp) THEN
                CALL fatal_analytic_input_error('Nonpositive gravity multiplier in Voellmy analytic solver')
        END IF

        xi_effective = xi
        IF (use_stochastic_xi) THEN
                xi_effective = xi + stochastic_xi_value
                xi_effective = MIN(1.0E5_wp, MAX(1.0_wp, xi_effective))
        ELSEIF (xi_effective .LE. 0.0_wp) THEN
                CALL fatal_analytic_input_error('Nonpositive xi in Voellmy analytic solver')
        END IF

        coulomb_stress = mu * rhom * (red_grav +              &
                    curvature_acceleration) * h * SQRT(grav_coeff)
        coulomb_stress = MAX(0.0_wp, coulomb_stress)
        IF (use_pore_pressure) coulomb_stress = MAX(0.0_wp, coulomb_stress -    &
                    mu * SQRT(grav_coeff) * exc_pore_pres)
        constant_drag = coulomb_stress / (rhom * h * speed_ratio)
        quadratic_drag = red_grav * speed_ratio / (xi_effective*h)
        CALL solve_exact_velocity_flow(EXACT_FLOW_QUADRATIC_POLYNOMIAL,         &
                    [constant_drag, 0.0_wp, quadratic_drag], dt,          &
                    speed_before, speed_after)
        velocity_after = (speed_after / speed_before) * velocity_before
    END SUBROUTINE solve_rheology_voellmy_analytic

    !******************************************************************************
    !> \brief Integrate the Voellmy velocity source exactly.
    !> \details Combines the frozen Coulomb stress with quadratic turbulent
    !> drag. Depth, density, gravity, curvature, pore pressure, and the
    !> mod_vel-to-horizontal-mod_vel ratio remain fixed during integration.
    !> \param[in] dt Exact-flow integration interval [s].
    !> \param[in] h Explicit flow depth [m].
    !> \param[in] speed_ratio Explicit total-to-horizontal mod_vel ratio.
    !> \param[in] velocity_before Input horizontal velocity [m s^-1].
    !> \param[in] rhom Explicit mixture density [kg m^-3].
    !> \param[in] red_grav Explicit reduced gravity [m s^-2].
    !> \param[in] grav_coeff Explicit gravity multiplier [dimensionless].
    !> \param[in] curvature_acceleration Explicit curvature contribution [m s^-2].
    !> \param[in] use_pore_pressure Enable explicit pore-pressure reduction.
    !> \param[in] exc_pore_pres Explicit excess pore pressure [Pa].
    !> \param[in] use_stochastic_xi Enable stochastic \f$\xi\f$ adjustment.
    !> \param[in] stochastic_xi_value Explicit stochastic \f$\xi\f$ increment.
    !> \param[out] velocity_after Returned horizontal velocity [m s^-1].
    !******************************************************************************
    !******************************************************************************
    !> \brief Evaluate the implicit turbulent Voellmy source.
    !> \details Computes the velocity-squared stress using \f$\xi\f$ and returns
    !> its x/y components opposite to the velocity direction. Optional
    !> stochastic variation of \f$\xi\f$ is bounded to [1, 10^5].
    !> \param[in] u x velocity [m s^-1].
    !> \param[in] v y velocity [m s^-1].
    !> \param[in] rhom Mixture density \f$\rho_m\f$ [kg m^-3].
    !> \param[in] red_grav Reduced gravity \f$g'\f$ [m s^-2].
    !> \param[in] mod_vel mod_vel magnitude [m s^-1].
    !> \param[in] mod_vel2 Squared mod_vel used in the turbulent stress [m^2 s^-2].
    !> \param[in] moving Whether the state is moving.
    !> \param[in] use_stochastic_xi Enable the stochastic \f$\xi\f$ variation.
    !> \param[in] stochastic_value Stochastic increment to \f$\xi\f$.
    !> \param[out] rate_x x conservative momentum source rate [kg m^-1 s^-2].
    !> \param[out] rate_y y conservative momentum source rate [kg m^-1 s^-2].
    !> \note Complex-valued state and rates preserve complex-step derivatives.
    !******************************************************************************
    !******************************************************************************
    !> \brief Evaluate the semi-implicit Coulomb Voellmy momentum rate.
    !> \details Builds the basal stress from flow depth, reduced gravity and
    !> optional curvature acceleration, clips negative stress to zero, and
    !> optionally reduces it by excess pore pressure. The returned x/y sources
    !> oppose the horizontal flow direction.
    !> \param[in] h Flow depth [m].
    !> \param[in] u x velocity [m s^-1].
    !> \param[in] v y velocity [m s^-1].
    !> \param[in] mod_vel mod_vel magnitude [m s^-1].
    !> \param[in] rhom Mixture density [kg m^-3].
    !> \param[in] red_grav Reduced gravity [m s^-2].
    !> \param[in] grav_coeff Gravity multiplier [dimensionless].
    !> \param[in] curvature_term Enable curvature acceleration.
    !> \param[in] B_second_xx xx curvature coefficient [m^-1].
    !> \param[in] B_second_xy xy curvature coefficient [m^-1].
    !> \param[in] B_second_yy yy curvature coefficient [m^-1].
    !> \param[in] use_pore_pressure Enable pore-pressure stress reduction.
    !> \param[in] exc_pore_pres Excess pore pressure [Pa].
    !> \param[out] rate_x x momentum source rate.
    !> \param[out] rate_y y momentum source rate.
    !******************************************************************************
    ! Unsupported physical treatment, rejected by the source catalog.
END MODULE rheology_voellmy_2d