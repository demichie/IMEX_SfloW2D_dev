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
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type, complex_state_from_real, &
         real_part_physical_state
    USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs,                      &
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
    !> \param[in] state Shared physical state (real or complex according to treatment).
    !> \brief Evaluate the explicit plastic contribution to the horizontal momentum source.
    !> \details Compute the momentum source under the basal plastic yield-resistance model, using
    !> the treatment-specific velocity and effective-stress inputs. Outputs are conservative x/y
    !> momentum source rates.
    !> \param[in] mod_vel Three-dimensional speed [m s^-1].
    !> \param[in] moving True for a moving local state.
    !> \param[out] rate Two-component conservative momentum source.
    SUBROUTINE eval_rheology_plastic_explicit(state, rate)
        TYPE(physical_state_type), INTENT(IN) :: state


        REAL(wp), INTENT(OUT) :: rate(2)
        REAL(wp) :: sx, sy
        COMPLEX(wp) :: ix, iy
        TYPE(complex_physical_state_type) :: cs
        CALL eval_rheology_plastic_semi_implicit(state, sx, sy)
        CALL complex_state_from_real(state, cs)
        CALL eval_rheology_plastic_implicit(cs, &
             ix, &
             iy)
        rate=[sx,sy]+REAL([ix,iy],wp)
    END SUBROUTINE eval_rheology_plastic_explicit

    !> \param[in] state Shared physical state (real or complex according to treatment).
    !> \brief Evaluate the semi-implicit plastic contribution to the horizontal momentum source.
    !> \details Return the contribution assigned to the semi-implicit treatment. This implementation
    !> supplies a zero source when the model has no term in this treatment.
    !> \param[out] rate_x x momentum source.
    !> \param[out] rate_y y momentum source.
    SUBROUTINE eval_rheology_plastic_semi_implicit(state, rate_x, rate_y)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(OUT) :: rate_x, rate_y

        rate_x = 0.0_wp
        rate_y = 0.0_wp
    END SUBROUTINE eval_rheology_plastic_semi_implicit

    !> \param[in] state Shared physical state (real or complex according to treatment).
    !> \brief Evaluate the implicit plastic contribution to the horizontal momentum source.
    !> \details Return the contribution assigned to the implicit treatment. This implementation
    !> supplies a zero source when the model has no term in this treatment.
    !> \param[in] mod_vel Three-dimensional speed [m s^-1].
    !> \param[in] moving True for a moving local state.
    !> \param[out] rate_x x momentum source.
    !> \param[out] rate_y y momentum source.
    SUBROUTINE eval_rheology_plastic_implicit(state, rate_x, rate_y)

        TYPE(complex_physical_state_type), INTENT(IN) :: state


        COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y

        rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
        rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
        IF (.NOT. ((REAL(state%mod_vel,wp)/=0.0_wp).OR.(AIMAG(state%mod_vel)/=0.0_wp))) RETURN
        rate_x = -state%rho_m * tau * (state%u / state%mod_vel)
        rate_y = -state%rho_m * tau * (state%v / state%mod_vel)
    END SUBROUTINE eval_rheology_plastic_implicit

    !> \param[in] state Shared real physical state at the start of the analytic substep.
    !> \brief Integrate the plastic velocity source in time.
    !> \details Advance the scalar speed under the basal plastic yield-resistance model with model
    !> coefficients fixed over dt; the horizontal direction is retained in the returned velocity
    !> vector.
    !> \param[in] dt Source integration interval [s].
    !> \param[in] speed_ratio Ratio of total to horizontal speed.
    !> \param[out] velocity_after Horizontal velocity vector after integration [m s^-1].
    SUBROUTINE solve_rheology_plastic_analytic(dt, state, speed_ratio, velocity_after)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: dt, speed_ratio
        REAL(wp) :: velocity_before(2)
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        REAL(wp) :: speed_before, speed_after

        velocity_before = [state%u,state%v]
        velocity_after = velocity_before
        speed_before = SQRT(SUM(velocity_before**2))
        CALL validate_analytic_velocity_inputs(dt, state%h,                                              &
                speed_ratio, speed_before, 'plastic')
        IF (dt .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
        CALL solve_exact_velocity_flow(EXACT_FLOW_CONSTANT,                                        &
                [tau/(state%h*speed_ratio)], dt, speed_before,                                           &
                speed_after)
        velocity_after = (speed_after / speed_before) * velocity_before
    END SUBROUTINE solve_rheology_plastic_analytic

    ! Unsupported physical treatment, rejected by the source catalog.
    ! Valid zero component of the configured Newton treatment.
END MODULE rheology_plastic_2d