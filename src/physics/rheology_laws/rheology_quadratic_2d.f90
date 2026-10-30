!********************************************************************************
!> \brief Quadratic basal-friction model.
!> \details The implicit drag is quadratic in the supplied mod_vel:
!> \f[
!> \boldsymbol{S}=-\rho_m f_q|\boldsymbol{u}|^2
!>                \frac{\boldsymbol{u}_h}{|\boldsymbol{u}_h|},
!> \f]
!> where \f$f_q\f$ is \ref friction_factor. The module provides
!> \ref eval_rheology_quadratic_implicit, which returns zero for a stationary
!> state and otherwise evaluates the complex-valued x/y sources.
!> \ref solve_rheology_quadratic_analytic integrates the frozen quadratic
!> velocity source and reconstructs both components along their input direction.
!********************************************************************************
MODULE rheology_quadratic_2d
    USE diagnostics_2d, ONLY : fatal_error
    USE parameters_2d, ONLY : wp
    USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs, &
            solve_exact_velocity_flow, EXACT_FLOW_QUADRATIC
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_rheology_quadratic_explicit
    PUBLIC :: eval_rheology_quadratic_semi_implicit
    PUBLIC :: eval_rheology_quadratic_implicit
    PUBLIC :: solve_rheology_quadratic_analytic
    PUBLIC :: friction_factor
    REAL(wp) :: friction_factor
    ! Uniform source treatment entry points. Existing kernels remain unchanged.
CONTAINS

    SUBROUTINE eval_rheology_quadratic_explicit( &
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
        REAL(wp) :: sx, sy, curvature_accel
        COMPLEX(wp) :: ix, iy
        curvature_accel=0.0_wp
        IF (curvature_term) curvature_accel=B_second_xx*u**2 + &
                2.0_wp*B_second_xy*u*v+B_second_yy*v**2
        CALL eval_rheology_quadratic_semi_implicit( &
            sx, sy )
        CALL eval_rheology_quadratic_implicit( &
            CMPLX(u,0.0_wp,wp), CMPLX(v,0.0_wp,wp), &
            CMPLX(mod_hor_vel,0.0_wp,wp), CMPLX(mod_vel2,0.0_wp,wp), &
            CMPLX(rhom,0.0_wp,wp), moving, ix, iy )
        rate=[sx,sy]+REAL([ix,iy],wp)
    END SUBROUTINE eval_rheology_quadratic_explicit

    SUBROUTINE eval_rheology_quadratic_semi_implicit(rate_x, rate_y)
        REAL(wp), INTENT(OUT) :: rate_x, rate_y
        rate_x = 0.0_wp
        rate_y = 0.0_wp
    END SUBROUTINE eval_rheology_quadratic_semi_implicit

    SUBROUTINE eval_rheology_quadratic_implicit(u, v,   &
            mod_hor_vel, mod_vel2, rhom, moving, rate_x, rate_y)
        COMPLEX(wp), INTENT(IN) :: u, v, mod_hor_vel
        COMPLEX(wp), INTENT(IN) :: mod_vel2, rhom
        LOGICAL, INTENT(IN) :: moving
        COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y
        rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
        rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
        IF (.NOT. moving) RETURN
        rate_x = -rhom * (u / mod_hor_vel)            &
                * friction_factor * mod_vel2
        rate_y = -rhom * (v / mod_hor_vel)            &
                * friction_factor * mod_vel2
    END SUBROUTINE eval_rheology_quadratic_implicit

    SUBROUTINE solve_rheology_quadratic_analytic(dt, h, &
            speed_ratio, velocity_before, velocity_after)
        REAL(wp), INTENT(IN) :: dt, h, speed_ratio
        REAL(wp), INTENT(IN) :: velocity_before(2)
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        REAL(wp) :: speed_before, speed_after, quadratic_drag
        velocity_after = velocity_before
        speed_before = SQRT(SUM(velocity_before**2))
        CALL validate_analytic_velocity_inputs(dt, h,    &
                speed_ratio, speed_before, 'quadratic')
        IF (dt .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
        quadratic_drag = friction_factor * speed_ratio**2 / h
        CALL solve_exact_velocity_flow(EXACT_FLOW_QUADRATIC, [quadratic_drag],  &
                dt, speed_before, speed_after)
        velocity_after = (speed_after / speed_before) * velocity_before
    END SUBROUTINE solve_rheology_quadratic_analytic

    ! Unsupported physical treatment, rejected by the source catalog.
    ! Valid zero component of the configured Newton treatment.
END MODULE rheology_quadratic_2d