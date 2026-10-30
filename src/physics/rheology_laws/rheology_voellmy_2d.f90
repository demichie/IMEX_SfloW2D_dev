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
    !> \brief Evaluate the explicit voellmy contribution to the horizontal momentum source.
    !> \details Compute the momentum source under the Voellmy Coulomb and turbulent basal resistance law, using the treatment-specific velocity and effective-stress inputs. Outputs are conservative x/y momentum source rates.
    !> \param[in] h Flow depth [m].
    !> \param[in] T Absolute mixture temperature [K].
    !> \param[in] u Horizontal x velocity [m s^-1].
    !> \param[in] v Horizontal y velocity [m s^-1].
    !> \param[in] mod_vel Three-dimensional speed [m s^-1].
    !> \param[in] mod_hor_vel Horizontal speed [m s^-1].
    !> \param[in] mod_vel2 Squared three-dimensional speed [m^2 s^-2].
    !> \param[in] rhom Mixture density [kg m^-3].
    !> \param[in] inv_rhom Reciprocal mixture density [m^3 kg^-1].
    !> \param[in] red_grav Reduced gravitational acceleration [m s^-2].
    !> \param[in] grav_coeff Local dimensionless gravitational correction.
    !> \param[in] alphas Solid volume fractions by class.
    !> \param[in] curvature_term Enable curvature correction.
    !> \param[in] B_second_xx Topographic second derivative in x [m^-1].
    !> \param[in] B_second_xy Mixed topographic derivative [m^-1].
    !> \param[in] B_second_yy Topographic second derivative in y [m^-1].
    !> \param[in] use_pore_pressure Enable pore-pressure correction.
    !> \param[in] exc_pore_pres Excess pore pressure [Pa].
    !> \param[in] moving True for a moving local state.
    !> \param[in] use_stochastic Enable model stochastic parameter.
    !> \param[in] stochastic_value Stochastic parameter value.
    !> \param[out] rate Two-component conservative momentum source.
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

    !> \brief Evaluate the semi-implicit voellmy contribution to the horizontal momentum source.
    !> \details Return the contribution assigned to the semi-implicit treatment. This implementation supplies a zero source when the model has no term in this treatment.
    !> \param[in] h Flow depth [m].
    !> \param[in] u Horizontal x velocity [m s^-1].
    !> \param[in] v Horizontal y velocity [m s^-1].
    !> \param[in] mod_vel Three-dimensional speed [m s^-1].
    !> \param[in] rhom Mixture density [kg m^-3].
    !> \param[in] red_grav Reduced gravitational acceleration [m s^-2].
    !> \param[in] grav_coeff Local dimensionless gravitational correction.
    !> \param[in] curvature_term Enable curvature correction.
    !> \param[in] B_second_xx Topographic second derivative in x [m^-1].
    !> \param[in] B_second_xy Mixed topographic derivative [m^-1].
    !> \param[in] B_second_yy Topographic second derivative in y [m^-1].
    !> \param[in] use_pore_pressure Enable pore-pressure correction.
    !> \param[in] exc_pore_pres Excess pore pressure [Pa].
    !> \param[out] rate_x x momentum source.
    !> \param[out] rate_y y momentum source.
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

    !> \brief Evaluate the implicit voellmy contribution to the horizontal momentum source.
    !> \details Return the contribution assigned to the implicit treatment. This implementation supplies a zero source when the model has no term in this treatment.
    !> \param[in] u Horizontal x velocity [m s^-1].
    !> \param[in] v Horizontal y velocity [m s^-1].
    !> \param[in] rhom Mixture density [kg m^-3].
    !> \param[in] red_grav Reduced gravitational acceleration [m s^-2].
    !> \param[in] mod_vel Three-dimensional speed [m s^-1].
    !> \param[in] mod_vel2 Squared three-dimensional speed [m^2 s^-2].
    !> \param[in] moving True for a moving local state.
    !> \param[in] use_stochastic_xi Enable stochastic turbulent resistance.
    !> \param[in] stochastic_value Stochastic parameter value.
    !> \param[out] rate_x x momentum source.
    !> \param[out] rate_y y momentum source.
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

    !> \brief Integrate the voellmy velocity source in time.
    !> \details Advance the scalar speed under the Voellmy Coulomb and turbulent basal resistance law with model coefficients fixed over dt; the horizontal direction is retained in the returned velocity vector.
    !> \param[in] dt Source integration interval [s].
    !> \param[in] h Flow depth [m].
    !> \param[in] speed_ratio Ratio of total to horizontal speed.
    !> \param[in] velocity_before Horizontal velocity vector before source integration [m s^-1].
    !> \param[in] rhom Mixture density [kg m^-3].
    !> \param[in] red_grav Reduced gravitational acceleration [m s^-2].
    !> \param[in] grav_coeff Local dimensionless gravitational correction.
    !> \param[in] curvature_acceleration Normal acceleration due to topographic curvature [m s^-2].
    !> \param[in] use_pore_pressure Enable pore-pressure correction.
    !> \param[in] exc_pore_pres Excess pore pressure [Pa].
    !> \param[in] use_stochastic_xi Enable stochastic turbulent resistance.
    !> \param[in] stochastic_xi_value Stochastic turbulent resistance parameter.
    !> \param[out] velocity_after Horizontal velocity vector after integration [m s^-1].
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

    ! Unsupported physical treatment, rejected by the source catalog.
END MODULE rheology_voellmy_2d