!********************************************************************************
!> \brief Empirical depth- and mod_vel-dependent drag model.
!> \details The implemented friction coefficient is
!> \f$C(h,|\boldsymbol{u}|)=10^{-3}|\boldsymbol{u}|/(1+10h)\f$. The
!> horizontal source is
!> \f$\boldsymbol{S}=-\rho_m C\boldsymbol{u}_h/|\boldsymbol{u}_h|\f$.
!> \ref eval_rheology_empirical_implicit evaluates its complex-valued x/y
!> components for a moving state; the empirical constants are fixed in the
!> routine. \ref solve_rheology_empirical_analytic integrates its frozen
!> linear-drag velocity source exactly.
!********************************************************************************
MODULE rheology_empirical_2d
    USE diagnostics_2d, ONLY : fatal_error
    USE parameters_2d, ONLY : wp
    USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs,  &
            solve_exact_velocity_flow, EXACT_FLOW_LINEAR
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_rheology_empirical_explicit
    PUBLIC :: eval_rheology_empirical_semi_implicit
    PUBLIC :: eval_rheology_empirical_implicit
    PUBLIC :: solve_rheology_empirical_analytic
    ! Uniform source treatment entry points. Existing kernels remain unchanged.
CONTAINS
    !> \brief Evaluate the explicit empirical momentum contribution.
    !> \details Assemble the real-valued momentum source from the model's supported physical treatment contributions; unused model inputs are retained for the shared dispatch interface.
    !> \param[in] h Flow depth [m].
    !> \param[in] T Mixture temperature [K].
    !> \param[in] u Depth-averaged x velocity [m s^-1].
    !> \param[in] v Depth-averaged y velocity [m s^-1].
    !> \param[in] mod_vel Magnitude of three-dimensional velocity [m s^-1].
    !> \param[in] mod_hor_vel Magnitude of horizontal velocity [m s^-1].
    !> \param[in] mod_vel2 Squared three-dimensional speed [m^2 s^-2].
    !> \param[in] rhom Mixture density [kg m^-3].
    !> \param[in] inv_rhom Reciprocal mixture density [m^3 kg^-1].
    !> \param[in] red_grav Buoyancy-reduced gravity [m s^-2].
    !> \param[in] grav_coeff Dimensionless gravity correction for local topography.
    !> \param[in] alphas Solid volume fractions by class.
    !> \param[in] curvature_term Whether topographic curvature enters the source.
    !> \param[in] B_second_xx Second x derivative of topography [m^-1].
    !> \param[in] B_second_xy Mixed derivative of topography [m^-1].
    !> \param[in] B_second_yy Second y derivative of topography [m^-1].
    !> \param[in] use_pore_pressure Whether excess pore pressure modifies resistance.
    !> \param[in] exc_pore_pres Excess pore pressure [Pa].
    !> \param[in] moving Whether the local state is considered moving.
    !> \param[in] use_stochastic Whether stochastic rheology parameters are enabled.
    !> \param[in] stochastic_value Stochastic model parameter used in the friction law.
    !> \param[out] rate Conservative horizontal momentum source vector.
    SUBROUTINE eval_rheology_empirical_explicit( &
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

        CALL eval_rheology_empirical_semi_implicit( &
            sx, sy )
        CALL eval_rheology_empirical_implicit( &
            CMPLX(h,0.0_wp,wp), CMPLX(u,0.0_wp,wp), &
            CMPLX(v,0.0_wp,wp), CMPLX(mod_hor_vel,0.0_wp,wp), &
            CMPLX(mod_vel,0.0_wp,wp), CMPLX(rhom,0.0_wp,wp), moving, ix, iy )
        rate=[sx,sy]+REAL([ix,iy],wp)
    END SUBROUTINE eval_rheology_empirical_explicit
    !> \brief Evaluate the semi-implicit empirical momentum contribution.
    !> \details Evaluate the contribution of the empirical basal-friction parameterization assigned to this treatment. The returned rates act on the conservative horizontal momenta; zero-valued implementations intentionally provide no contribution.
    !> \param[out] rate_x x component of momentum source rate.
    !> \param[out] rate_y y component of momentum source rate.
    SUBROUTINE eval_rheology_empirical_semi_implicit(rate_x, rate_y)
        REAL(wp), INTENT(OUT) :: rate_x, rate_y

        rate_x = 0.0_wp
        rate_y = 0.0_wp
    END SUBROUTINE eval_rheology_empirical_semi_implicit
    !> \brief Evaluate the implicit empirical momentum contribution.
    !> \details Evaluate the part of the empirical basal-friction parameterization assigned to the implicit treatment. Complex-valued arguments retain complex-step sensitivities where implemented; zero-valued implementations intentionally contribute no term.
    !> \param[in] h Flow depth [m].
    !> \param[in] u Depth-averaged x velocity [m s^-1].
    !> \param[in] v Depth-averaged y velocity [m s^-1].
    !> \param[in] mod_hor_vel Magnitude of horizontal velocity [m s^-1].
    !> \param[in] mod_vel Magnitude of three-dimensional velocity [m s^-1].
    !> \param[in] rhom Mixture density [kg m^-3].
    !> \param[in] moving Whether the local state is considered moving.
    !> \param[out] rate_x x component of momentum source rate.
    !> \param[out] rate_y y component of momentum source rate.
    SUBROUTINE eval_rheology_empirical_implicit(h, u, v, &
            mod_hor_vel, mod_vel, rhom, moving, rate_x, rate_y)
        COMPLEX(wp), INTENT(IN) :: h, u, v
        COMPLEX(wp), INTENT(IN) :: mod_hor_vel, mod_vel, rhom
        LOGICAL, INTENT(IN) :: moving
        COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y
        COMPLEX(wp) :: friction_coefficient

        rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
        rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
        IF (.NOT. moving) RETURN
        friction_coefficient = 1.0E-3_wp / (1.0_wp + 10.0_wp * h) * mod_vel
        rate_x = -rhom * friction_coefficient * (u / mod_hor_vel)
        rate_y = -rhom * friction_coefficient * (v / mod_hor_vel)
    END SUBROUTINE eval_rheology_empirical_implicit
    !> \brief Integrate the empirical velocity source over a frozen-state substep.
    !> \details Evaluate the exact or time-integral-based decay of horizontal speed under the empirical basal-friction parameterization. The velocity direction is retained and the updated components are returned without evolving the frozen thermodynamic and geometric inputs.
    !> \param[in] dt Duration of the analytic source substep [s].
    !> \param[in] h Flow depth [m].
    !> \param[in] speed_ratio Ratio of three-dimensional to horizontal speed.
    !> \param[in] velocity_before Two-component horizontal velocity before source integration [m s^-1].
    !> \param[out] velocity_after Two-component horizontal velocity after source integration [m s^-1].
    SUBROUTINE solve_rheology_empirical_analytic(dt, h, &
            speed_ratio, velocity_before, velocity_after)
        REAL(wp), INTENT(IN) :: dt, h, speed_ratio
        REAL(wp), INTENT(IN) :: velocity_before(2)
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        REAL(wp) :: linear_drag, speed_before, speed_after

        velocity_after = velocity_before
        speed_before = SQRT(SUM(velocity_before**2))
        CALL validate_analytic_velocity_inputs(dt, h,    &
                speed_ratio, speed_before, 'empirical')
        IF (dt .EQ. 0.0_wp) RETURN
        IF (speed_before .EQ. 0.0_wp) RETURN
        linear_drag = 1.0E-3_wp * speed_ratio /                               &
                (h * (1.0_wp + 10.0_wp*h))
        CALL solve_exact_velocity_flow(EXACT_FLOW_LINEAR, [linear_drag],        &
                dt, speed_before, speed_after)
        velocity_after = (speed_after / speed_before) * velocity_before
    END SUBROUTINE solve_rheology_empirical_analytic

    ! Unsupported physical treatment, rejected by the source catalog.
    ! Valid zero component of the configured Newton treatment.
END MODULE rheology_empirical_2d