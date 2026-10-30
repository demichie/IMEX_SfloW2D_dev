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
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
    USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs,                      &
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
    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the explicit quadratic contribution to the horizontal momentum source.
    !> \details Compute the momentum source under the quadratic velocity-dependent resistance law,
    !> using the treatment-specific velocity and effective-stress inputs. Outputs are conservative
    !> x/y momentum source rates.
    !> \param[in] mod_vel Three-dimensional speed [m s^-1].
    !> \param[in] mod_hor_vel Horizontal speed [m s^-1].
    !> \param[in] mod_vel2 Squared three-dimensional speed [m^2 s^-2].
    !> \param[in] grav_coeff Local dimensionless gravitational correction.
    !> \param[in] curvature_term Enable curvature correction.
    !> \param[in] B_second_xx Topographic second derivative in x [m^-1].
    !> \param[in] B_second_xy Mixed topographic derivative [m^-1].
    !> \param[in] B_second_yy Topographic second derivative in y [m^-1].
    !> \param[in] use_pore_pressure Enable pore-pressure correction.
    !> \param[in] moving True for a moving local state.
    !> \param[in] use_stochastic Enable model stochastic parameter.
    !> \param[in] stochastic_value Stochastic parameter value.
    !> \param[out] rate Two-component conservative momentum source.
    SUBROUTINE eval_rheology_quadratic_explicit(state, mod_vel, mod_hor_vel, mod_vel2, &
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
        CALL eval_rheology_quadratic_semi_implicit(state,sx,sy)
        cs%u = CMPLX(state%u,0.0_wp,wp)
        cs%v = CMPLX(state%v,0.0_wp,wp)
        cs%rho_m = CMPLX(state%rho_m,0.0_wp,wp)
        CALL eval_rheology_quadratic_implicit(cs, CMPLX(mod_hor_vel,0.0_wp,wp), &
             CMPLX(mod_vel2,0.0_wp,wp), moving, ix, iy)
        rate=[sx,sy]+REAL([ix,iy],wp)
    END SUBROUTINE eval_rheology_quadratic_explicit

    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the semi-implicit quadratic contribution to the horizontal momentum source.
    !> \details Return the contribution assigned to the semi-implicit treatment. This implementation
    !> supplies a zero source when the model has no term in this treatment.
    !> \param[out] rate_x x momentum source.
    !> \param[out] rate_y y momentum source.
    SUBROUTINE eval_rheology_quadratic_semi_implicit(state, rate_x, rate_y)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(OUT) :: rate_x, rate_y

        rate_x = 0.0_wp
        rate_y = 0.0_wp
    END SUBROUTINE eval_rheology_quadratic_semi_implicit

    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the implicit quadratic contribution to the horizontal momentum source.
    !> \details Return the contribution assigned to the implicit treatment. This implementation
    !> supplies a zero source when the model has no term in this treatment.
    !> \param[in] mod_hor_vel Horizontal speed [m s^-1].
    !> \param[in] mod_vel2 Squared three-dimensional speed [m^2 s^-2].
    !> \param[in] moving True for a moving local state.
    !> \param[out] rate_x x momentum source.
    !> \param[out] rate_y y momentum source.
    SUBROUTINE eval_rheology_quadratic_implicit(state, mod_hor_vel, mod_vel2, &
            moving, rate_x, rate_y)

        TYPE(complex_physical_state_type), INTENT(IN) :: state
        COMPLEX(wp), INTENT(IN) :: mod_hor_vel, mod_vel2
        LOGICAL, INTENT(IN) :: moving
        COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y

        rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
        rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
        IF (.NOT. moving) RETURN
        rate_x = -state%rho_m * (state%u / mod_hor_vel)                                                         &
                * friction_factor * mod_vel2
        rate_y = -state%rho_m * (state%v / mod_hor_vel)                                                         &
                * friction_factor * mod_vel2
    END SUBROUTINE eval_rheology_quadratic_implicit

    !> \param[in] state Shared real state at the start of the analytic source step.
    !> \brief Integrate the quadratic velocity source in time.
    !> \details Advance the scalar speed under the quadratic velocity-dependent resistance law with
    !> model coefficients fixed over dt; the horizontal direction is retained in the returned
    !> velocity vector.
    !> \param[in] dt Source integration interval [s].
    !> \param[in] speed_ratio Ratio of total to horizontal speed.
    !> \param[out] velocity_after Horizontal velocity vector after integration [m s^-1].
    SUBROUTINE solve_rheology_quadratic_analytic(dt, state, speed_ratio, velocity_after)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: dt, speed_ratio
        REAL(wp) :: velocity_before(2)
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        REAL(wp) :: speed_before, speed_after, quadratic_drag

        velocity_before = [state%u,state%v]
        velocity_after = velocity_before
        speed_before = SQRT(SUM(velocity_before**2))
        CALL validate_analytic_velocity_inputs(dt, state%h,                                              &
                speed_ratio, speed_before, 'quadratic')
        IF (dt .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
        quadratic_drag = friction_factor * speed_ratio**2 / state%h
        CALL solve_exact_velocity_flow(EXACT_FLOW_QUADRATIC, [quadratic_drag],                     &
                dt, speed_before, speed_after)
        velocity_after = (speed_after / speed_before) * velocity_before
    END SUBROUTINE solve_rheology_quadratic_analytic

    ! Unsupported physical treatment, rejected by the source catalog.
    ! Valid zero component of the configured Newton treatment.
END MODULE rheology_quadratic_2d