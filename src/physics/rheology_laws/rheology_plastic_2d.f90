!********************************************************************************
!> \brief Constant-strength plastic basal-friction model.
!> \details The model applies a mod_vel-independent resisting strength
!> \f$\tau\f$ opposite to the horizontal flow direction:
!> \f[
!> \boldsymbol{S}=-\rho_m\tau\,\boldsymbol{u}_h/|\boldsymbol{u}|.
!> \f]
!> The supplied mod_vel is the normalization used by the implementation. The
!> source is zero when the caller indicates that the flow is not moving.
!> \ref eval_rheology_plastic_implicit evaluates the complex-valued x/y rates
!> for the implicit solve, while \ref solve_rheology_plastic_analytic
!> integrates its frozen velocity source exactly.
!> \ref tau stores the model strength.
!********************************************************************************
MODULE rheology_plastic_2d
    USE diagnostics_2d, ONLY : fatal_error
    USE parameters_2d, ONLY : wp
    USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs,  &
            solve_exact_velocity_flow, EXACT_FLOW_CONSTANT
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_rheology_plastic_explicit
    PUBLIC :: eval_rheology_plastic_semi_implicit
    PUBLIC :: eval_rheology_plastic_implicit
    PUBLIC :: solve_rheology_plastic_analytic
    PUBLIC :: tau
    REAL(wp) :: tau
    ! Uniform source treatment entry points. Existing kernels remain unchanged.
CONTAINS

    SUBROUTINE eval_rheology_plastic_explicit( &
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
        CALL eval_rheology_plastic_semi_implicit( &
            sx, sy )
        CALL eval_rheology_plastic_implicit( &
            CMPLX(u,0.0_wp,wp), CMPLX(v,0.0_wp,wp), CMPLX(mod_vel,0.0_wp,wp), &
            CMPLX(rhom,0.0_wp,wp), moving, ix, iy )
        rate=[sx,sy]+REAL([ix,iy],wp)
    END SUBROUTINE eval_rheology_plastic_explicit

    SUBROUTINE eval_rheology_plastic_semi_implicit(rate_x, rate_y)
        REAL(wp), INTENT(OUT) :: rate_x, rate_y
        rate_x = 0.0_wp
        rate_y = 0.0_wp
    END SUBROUTINE eval_rheology_plastic_semi_implicit

    SUBROUTINE eval_rheology_plastic_implicit(u, v, mod_vel, &
            rhom, moving, rate_x, rate_y)
        COMPLEX(wp), INTENT(IN) :: u, v, mod_vel, rhom
        LOGICAL, INTENT(IN) :: moving
        COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y
        rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
        rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
        IF (.NOT. moving) RETURN
        rate_x = -rhom * tau * (u / mod_vel)
        rate_y = -rhom * tau * (v / mod_vel)
    END SUBROUTINE eval_rheology_plastic_implicit

    SUBROUTINE solve_rheology_plastic_analytic(dt, h,  &
            speed_ratio, velocity_before, velocity_after)
        REAL(wp), INTENT(IN) :: dt, h, speed_ratio
        REAL(wp), INTENT(IN) :: velocity_before(2)
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        REAL(wp) :: speed_before, speed_after
        velocity_after = velocity_before
        speed_before = SQRT(SUM(velocity_before**2))
        CALL validate_analytic_velocity_inputs(dt, h,    &
                speed_ratio, speed_before, 'plastic')
        IF (dt .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
        CALL solve_exact_velocity_flow(EXACT_FLOW_CONSTANT,                     &
                [tau/(h*speed_ratio)], dt, speed_before,    &
                speed_after)
        velocity_after = (speed_after / speed_before) * velocity_before
    END SUBROUTINE solve_rheology_plastic_analytic

    ! Unsupported physical treatment, rejected by the source catalog.
    ! Valid zero component of the configured Newton treatment.
END MODULE rheology_plastic_2d