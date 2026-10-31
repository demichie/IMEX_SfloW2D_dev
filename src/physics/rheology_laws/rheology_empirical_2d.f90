!********************************************************************************
!> \brief Empirical depth- and mod_vel-dependent drag model.
!> \details The implemented friction coefficient is
!> \f$C(h,|\boldsymbol{u}|)=10^{-3}|\boldsymbol{u}|/(1+10h)\f$. The
!> horizontal source is
!> \f$\boldsymbol{S}=-\rho_m C\boldsymbol{u}_h/|\boldsymbol{u}_h|\f$.
!> \ref eval_rheology_empirical_implicit evaluates its complex-valued x/y
!> components for a moving state; the empirical constants are fixed in the
!> routine. \ref eval_rheology_empirical_analytic integrates its frozen
!> linear-drag velocity source exactly.
!********************************************************************************
MODULE rheology_empirical_2d
    USE diagnostics_2d, ONLY : fatal_error
    USE parameters_2d, ONLY : wp
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type, complex_state_from_real
    USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs,                      &
            solve_exact_velocity_flow, EXACT_FLOW_LINEAR
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_rheology_empirical_explicit
    PUBLIC :: eval_rheology_empirical_semi_implicit
    PUBLIC :: eval_rheology_empirical_implicit
    PUBLIC :: eval_rheology_empirical_analytic
    ! Uniform source treatment entry points. Existing kernels remain unchanged.
CONTAINS
    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the explicit empirical momentum contribution.
    !> \details Assemble the real-valued momentum source from the model's supported physical
    !> treatment contributions; unused model inputs are retained for the shared dispatch interface.
    !> \param[in] mod_vel Magnitude of three-dimensional velocity [m s^-1].
    !> \param[in] mod_hor_vel Magnitude of horizontal velocity [m s^-1].
    !> \param[in] moving Whether the local state is considered moving.
    !> \param[out] rate Conservative horizontal momentum source vector.
    SUBROUTINE eval_rheology_empirical_explicit(state, rate)
        TYPE(physical_state_type), INTENT(IN) :: state


        REAL(wp), INTENT(OUT) :: rate(2)
        REAL(wp) :: sx, sy
        COMPLEX(wp) :: ix, iy
        TYPE(complex_physical_state_type) :: cs
        CALL eval_rheology_empirical_semi_implicit(state, sx, sy)
        CALL complex_state_from_real(state, cs)
        CALL eval_rheology_empirical_implicit(cs, &
             ix, &
             iy)
        rate=[sx,sy]+REAL([ix,iy],wp)
    END SUBROUTINE eval_rheology_empirical_explicit

    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the semi-implicit empirical momentum contribution.
    !> \details Evaluate the contribution of the empirical basal-friction parameterization assigned
    !> to this treatment. The returned rates act on the conservative horizontal momenta; zero-valued
    !> implementations intentionally provide no contribution.
    !> \param[out] rate_x x component of momentum source rate.
    !> \param[out] rate_y y component of momentum source rate.
    SUBROUTINE eval_rheology_empirical_semi_implicit(state, rate_x, rate_y)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(OUT) :: rate_x, rate_y

        rate_x = 0.0_wp
        rate_y = 0.0_wp
    END SUBROUTINE eval_rheology_empirical_semi_implicit

    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the implicit empirical momentum contribution.
    !> \details Evaluate the part of the empirical basal-friction parameterization assigned to the
    !> implicit treatment. Complex-valued arguments retain complex-step sensitivities where
    !> implemented; zero-valued implementations intentionally contribute no term.
    !> \param[in] mod_hor_vel Magnitude of horizontal velocity [m s^-1].
    !> \param[in] mod_vel Magnitude of three-dimensional velocity [m s^-1].
    !> \param[in] moving Whether the local state is considered moving.
    !> \param[out] rate_x x component of momentum source rate.
    !> \param[out] rate_y y component of momentum source rate.
    SUBROUTINE eval_rheology_empirical_implicit(state, rate_x, rate_y)

        TYPE(complex_physical_state_type), INTENT(IN) :: state


        COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y
        COMPLEX(wp) :: friction_coefficient

        rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
        rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
        IF (.NOT. ((REAL(state%mod_vel,wp)/=0.0_wp).OR.(AIMAG(state%mod_vel)/=0.0_wp))) RETURN
        friction_coefficient = 1.0E-3_wp / (1.0_wp + 10.0_wp * state%h) * state%mod_vel
        rate_x = -state%rho_m * friction_coefficient * (state%u / state%mod_vel_hor)
        rate_y = -state%rho_m * friction_coefficient * (state%v / state%mod_vel_hor)
    END SUBROUTINE eval_rheology_empirical_implicit

    !> \param[in] state Shared real physical state before the analytic substep.
    !> \brief Integrate the empirical velocity source over a frozen-state substep.
    !> \details Evaluate the exact or time-integral-based decay of horizontal speed under the
    !> empirical basal-friction parameterization. The velocity direction is retained and the updated
    !> components are returned without evolving the frozen thermodynamic and geometric inputs.
    !> \param[in] dt Duration of the analytic source substep [s].
    !> \param[in] speed_ratio Ratio of three-dimensional to horizontal speed.
    !> s^-1].
    !> \param[out] velocity_after Two-component horizontal velocity after source integration [m
    !> s^-1].
    SUBROUTINE eval_rheology_empirical_analytic(dt, state, speed_ratio, velocity_after)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: dt, speed_ratio
        REAL(wp) :: velocity_before(2)
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        REAL(wp) :: linear_drag, speed_before, speed_after

        velocity_before = [state%u,state%v]
        velocity_after = velocity_before
        speed_before = SQRT(SUM(velocity_before**2))
        CALL validate_analytic_velocity_inputs(dt, state%h,                                              &
                speed_ratio, speed_before, 'empirical')
        IF (dt .EQ. 0.0_wp) RETURN
        IF (speed_before .EQ. 0.0_wp) RETURN
        linear_drag = 1.0E-3_wp * speed_ratio /                                                    &
                (state%h * (1.0_wp + 10.0_wp*state%h))
        CALL solve_exact_velocity_flow(EXACT_FLOW_LINEAR, [linear_drag],                           &
                dt, speed_before, speed_after)
        velocity_after = (speed_after / speed_before) * velocity_before
    END SUBROUTINE eval_rheology_empirical_analytic

    ! Unsupported physical treatment, rejected by the source catalog.
    ! Valid zero component of the configured Newton treatment.
END MODULE rheology_empirical_2d