!********************************************************************************
!> \brief temperature-dependent viscous and yield-stress lava rheology.
!> \details The implicit viscous drag uses a temperature-dependent kinematic
!> viscosity and a depth-averaged shear factor:
!> \f[
!> \gamma=\frac{3\nu_{\rm ref}}{h_*}
!>   \exp[-b(T-T_{\rm ref})],\qquad
!> \boldsymbol{S}_{v}=-\rho_m\gamma\boldsymbol{u}_h,
!> \f]
!> where \f$h_*=\max(h,10^{-10})\f$ in the depth branch selected by the
!> implementation. The semi-implicit yield-stress source is
!> \f$\boldsymbol{S}_{y}=-\tau_0\boldsymbol{u}_h/|\boldsymbol{u}|\f$.
!> \ref eval_rheology_lava_implicit evaluates the complex-valued viscous rate;
!> \ref eval_rheology_lava_semi_implicit evaluates the real-valued yield rate.
!> \ref nu_ref, \ref visc_par and \ref tau0 are the model parameters. Thermal
!> cooling is owned by the thermal-source module, not by this rheology.
!> \ref solve_rheology_lava_analytic integrates the frozen affine velocity
!> source exactly with T and other non-velocity inputs explicit.
!********************************************************************************
MODULE rheology_lava_2d
    USE diagnostics_2d, ONLY : fatal_error
    USE parameters_2d, ONLY : wp
    USE rheology_parameters_2d, ONLY : T_ref
    USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs,  &
            fatal_analytic_input_error, solve_exact_velocity_flow, EXACT_FLOW_AFFINE
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_rheology_lava_explicit
    PUBLIC :: eval_rheology_lava_semi_implicit
    PUBLIC :: eval_rheology_lava_implicit
    PUBLIC :: solve_rheology_lava_analytic
    PUBLIC :: nu_ref, visc_par, tau0, lava_kinematic_viscosity
    REAL(wp) :: nu_ref
    REAL(wp) :: visc_par
    REAL(wp) :: tau0
    ! Uniform source treatment entry points. Existing kernels remain unchanged.
CONTAINS
    !> \brief Kinematic viscosity shared by lava rheology and deposition.
    !> \param[in] T lava temperature [K].
    !> \return Kinematic viscosity [m2/s].
    PURE REAL(wp) FUNCTION lava_kinematic_viscosity(T) RESULT(nu)
        REAL(wp), INTENT(IN) :: T
        nu=nu_ref*EXP(-visc_par*(T-T_ref))
    END FUNCTION lava_kinematic_viscosity

    SUBROUTINE eval_rheology_lava_explicit( &
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
        CALL eval_rheology_lava_semi_implicit( &
            u, v, mod_vel, sx, sy )
        CALL eval_rheology_lava_implicit( &
            CMPLX(h,0.0_wp,wp), CMPLX(T,0.0_wp,wp), &
            CMPLX(rhom,0.0_wp,wp), CMPLX(u,0.0_wp,wp), &
            CMPLX(v,0.0_wp,wp), moving, ix, iy )
        rate=[sx,sy]+REAL([ix,iy],wp)
    END SUBROUTINE eval_rheology_lava_explicit

    SUBROUTINE eval_rheology_lava_semi_implicit(u, v, mod_vel, rate_x, rate_y)
        REAL(wp), INTENT(IN) :: u, v, mod_vel
        REAL(wp), INTENT(OUT) :: rate_x, rate_y
        rate_x = 0.0_wp
        rate_y = 0.0_wp
        IF (mod_vel .LE. 0.0_wp) RETURN
        rate_x = -tau0 * u / mod_vel
        rate_y = -tau0 * v / mod_vel
    END SUBROUTINE eval_rheology_lava_semi_implicit

    SUBROUTINE eval_rheology_lava_implicit(h, T, rhom, &
            u, v, moving, rate_x, rate_y)
        COMPLEX(wp), INTENT(IN) :: h, T, rhom
        COMPLEX(wp), INTENT(IN) :: u, v
        LOGICAL, INTENT(IN) :: moving
        COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y
        COMPLEX(wp) :: gamma
        REAL(wp), PARAMETER :: h_threshold = 1.0E-10_wp
        rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
        rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
        IF (REAL(h, wp) .GT. h_threshold) THEN
            gamma = 3.0_wp * nu_ref / h * EXP(-visc_par * (T - T_ref))
        ELSE
            gamma = 3.0_wp * nu_ref / h_threshold * EXP(-visc_par * (T - T_ref))
        END IF

        IF (moving) THEN
            rate_x = -rhom * gamma * u
            rate_y = -rhom * gamma * v
        END IF
    END SUBROUTINE eval_rheology_lava_implicit

    SUBROUTINE solve_rheology_lava_analytic(dt, h,     &
            T, rhom, speed_ratio, velocity_before,           &
            velocity_after)
        REAL(wp), INTENT(IN) :: dt, h, T
        REAL(wp), INTENT(IN) :: rhom, speed_ratio
        REAL(wp), INTENT(IN) :: velocity_before(2)
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        REAL(wp), PARAMETER :: h_threshold = 1.0E-10_wp
        REAL(wp) :: speed_before, speed_after, effective_depth, gamma
        REAL(wp) :: linear_drag, yield_drag
        velocity_after = velocity_before
        speed_before = SQRT(SUM(velocity_before**2))
        CALL validate_analytic_velocity_inputs(dt, h,    &
                speed_ratio, speed_before, 'lava')
        IF (dt .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
        IF (rhom .LE. 0.0_wp) THEN
            CALL fatal_analytic_input_error('Nonpositive mixture density in lava analytic solver')
        END IF

        effective_depth = MAX(h, h_threshold)
        gamma = 3.0_wp * nu_ref / effective_depth *                             &
                EXP(-visc_par * (T - T_ref))
        linear_drag = gamma / h
        yield_drag = tau0 / (rhom * h * speed_ratio)
        CALL solve_exact_velocity_flow(EXACT_FLOW_AFFINE,                       &
                [yield_drag, linear_drag], dt, speed_before,         &
                speed_after)
        velocity_after = (speed_after / speed_before) * velocity_before
    END SUBROUTINE solve_rheology_lava_analytic

    ! Unsupported physical treatment, rejected by the source catalog.
END MODULE rheology_lava_2d