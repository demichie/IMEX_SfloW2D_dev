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
!> \ref eval_rheology_quadratic_analytic integrates the frozen quadratic
!> velocity source and reconstructs both components along their input direction.
!********************************************************************************
MODULE rheology_quadratic_2d
    USE diagnostics_2d, ONLY : fatal_error
    USE parameters_2d, ONLY : wp
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type, complex_state_from_real
    USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs,                      &
            solve_exact_velocity_flow, EXACT_FLOW_QUADRATIC
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_rheology_quadratic_explicit
    PUBLIC :: eval_rheology_quadratic_semi_implicit
    PUBLIC :: eval_rheology_quadratic_implicit
    PUBLIC :: eval_rheology_quadratic_analytic
    PUBLIC :: friction_factor
    REAL(wp) :: friction_factor
    ! Uniform source treatment entry points. Existing kernels remain unchanged.
CONTAINS
    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the explicit quadratic contribution to the horizontal momentum source.
    !> \details Compute the momentum source under the quadratic velocity-dependent resistance law,
    !> using the treatment-specific velocity and effective-stress inputs. Outputs are conservative
    !> x/y momentum source rates.
    !> \param[in] mod_hor_vel Horizontal speed [m s^-1].
    !> \param[in] mod_vel2 Squared three-dimensional speed [m^2 s^-2].
    !> \param[in] moving True for a moving local state.
    !> \param[out] rate Two-component conservative momentum source.
    SUBROUTINE eval_rheology_quadratic_explicit(state, rate)
        TYPE(physical_state_type), INTENT(IN) :: state


        REAL(wp), INTENT(OUT) :: rate(2)
        REAL(wp) :: sx, sy
        COMPLEX(wp) :: ix, iy
        TYPE(complex_physical_state_type) :: cs
        CALL eval_rheology_quadratic_semi_implicit(state,sx,sy)
        CALL complex_state_from_real(state, cs)
        CALL eval_rheology_quadratic_implicit(cs, &
             ix, &
             iy)
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
    SUBROUTINE eval_rheology_quadratic_implicit(state, rate_x, rate_y)

        TYPE(complex_physical_state_type), INTENT(IN) :: state


        COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y

        rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
        rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
        IF (.NOT. ((REAL(state%mod_vel,wp)/=0.0_wp).OR.(AIMAG(state%mod_vel)/=0.0_wp))) RETURN
        rate_x = -state%rho_m * (state%u / state%mod_vel_hor)                                                         &
                * friction_factor * state%mod_vel2
        rate_y = -state%rho_m * (state%v / state%mod_vel_hor)                                                         &
                * friction_factor * state%mod_vel2
    END SUBROUTINE eval_rheology_quadratic_implicit

    !> \param[in] state Shared real state at the start of the analytic source step.
    !> \brief Integrate the quadratic velocity source in time.
    !> \details Advance the scalar speed under the quadratic velocity-dependent resistance law with
    !> model coefficients fixed over dt; the horizontal direction is retained in the returned
    !> velocity vector.
    !> \param[in] dt Source integration interval [s].
    !> \param[in] speed_ratio Ratio of total to horizontal speed.
    !> \param[out] velocity_after Horizontal velocity vector after integration [m s^-1].
    SUBROUTINE eval_rheology_quadratic_analytic(dt, state, speed_ratio, velocity_after)

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
    END SUBROUTINE eval_rheology_quadratic_analytic

    ! Unsupported physical treatment, rejected by the source catalog.
    ! Valid zero component of the configured Newton treatment.
END MODULE rheology_quadratic_2d