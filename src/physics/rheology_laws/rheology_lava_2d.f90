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
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type, complex_state_from_real
    USE rheology_parameters_2d, ONLY : T_ref
    USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs,                      &
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

    SUBROUTINE eval_rheology_lava_explicit(state, rate)
        TYPE(physical_state_type), INTENT(IN) :: state


        REAL(wp), INTENT(OUT) :: rate(2)
        REAL(wp) :: sx, sy
        COMPLEX(wp) :: ix, iy
        TYPE(complex_physical_state_type) :: cs
        CALL eval_rheology_lava_semi_implicit(state, &
             sx, &
             sy)
        CALL complex_state_from_real(state, cs)
        CALL eval_rheology_lava_implicit(cs, &
             ix, &
             iy)
        rate=[sx,sy]+REAL([ix,iy],wp)
    END SUBROUTINE eval_rheology_lava_explicit


    SUBROUTINE eval_rheology_lava_semi_implicit(state, rate_x, rate_y)

        TYPE(physical_state_type), INTENT(IN) :: state

        REAL(wp), INTENT(OUT) :: rate_x, rate_y

        rate_x = 0.0_wp
        rate_y = 0.0_wp
        IF (state%mod_vel .LE. 0.0_wp) RETURN
        rate_x = -tau0 * state%u / state%mod_vel
        rate_y = -tau0 * state%v / state%mod_vel
    END SUBROUTINE eval_rheology_lava_semi_implicit


    SUBROUTINE eval_rheology_lava_implicit(state, rate_x, rate_y)

        TYPE(complex_physical_state_type), INTENT(IN) :: state

        COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y
        COMPLEX(wp) :: gamma
        REAL(wp), PARAMETER :: h_threshold = 1.0E-10_wp

        rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
        rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
        IF (REAL(state%h, wp) .GT. h_threshold) THEN
            gamma = 3.0_wp * nu_ref / state%h * EXP(-visc_par * (state%T - T_ref))

        ELSE
            gamma = 3.0_wp * nu_ref / h_threshold * EXP(-visc_par * (state%T - T_ref))

        END IF

        IF (((REAL(state%mod_vel,wp)/=0.0_wp).OR.(AIMAG(state%mod_vel)/=0.0_wp))) THEN
            rate_x = -state%rho_m * gamma * state%u
            rate_y = -state%rho_m * gamma * state%v

        END IF

    END SUBROUTINE eval_rheology_lava_implicit


    SUBROUTINE solve_rheology_lava_analytic(dt, state, speed_ratio, velocity_after)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: dt, speed_ratio
        REAL(wp) :: velocity_before(2)
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        REAL(wp), PARAMETER :: h_threshold = 1.0E-10_wp
        REAL(wp) :: speed_before, speed_after, effective_depth, gamma
        REAL(wp) :: linear_drag, yield_drag

        velocity_before = [state%u,state%v]
        velocity_after = velocity_before
        speed_before = SQRT(SUM(velocity_before**2))
        CALL validate_analytic_velocity_inputs(dt, state%h,                                              &
                speed_ratio, speed_before, 'lava')
        IF (dt .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
        IF (state%rho_m .LE. 0.0_wp) THEN
            CALL fatal_analytic_input_error('Nonpositive mixture density in lava analytic solver')

        END IF

        effective_depth = MAX(state%h, h_threshold)
        gamma = 3.0_wp * nu_ref / effective_depth *                                                &
                EXP(-visc_par * (state%T - T_ref))
        linear_drag = gamma / state%h
        yield_drag = tau0 / (state%rho_m * state%h * speed_ratio)
        CALL solve_exact_velocity_flow(EXACT_FLOW_AFFINE,                                          &
                [yield_drag, linear_drag], dt, speed_before,                                       &
                speed_after)
        velocity_after = (speed_after / speed_before) * velocity_before
    END SUBROUTINE solve_rheology_lava_analytic

    !> \brief Kinematic viscosity shared by lava rheology and deposition.
    !> \return Kinematic viscosity [m2/s].
    PURE REAL(wp) FUNCTION lava_kinematic_viscosity(T) RESULT(nu)

        REAL(wp), INTENT(IN) :: T

        nu=nu_ref*EXP(-visc_par*(T-T_ref))
    END FUNCTION lava_kinematic_viscosity

    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the explicit lava momentum contribution.
    !> \details Assemble the real-valued momentum source from the model's supported physical
    !> treatment contributions; unused model inputs are retained for the shared dispatch interface.
    !> \param[in] mod_vel Magnitude of three-dimensional velocity [m s^-1].
    !> \param[in] moving Whether the local state is considered moving.
    !> \param[out] rate Conservative horizontal momentum source vector.

    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the semi-implicit lava momentum contribution.
    !> \details Evaluate the contribution of temperature-dependent viscous drag and a yield-stress
    !> contribution assigned to this treatment. The returned rates act on the conservative
    !> horizontal momenta; zero-valued implementations intentionally provide no contribution.
    !> \param[in] mod_vel Magnitude of three-dimensional velocity [m s^-1].
    !> \param[out] rate_x x component of momentum source rate.
    !> \param[out] rate_y y component of momentum source rate.

    !> \param[in] state Shared physical state (complex for implicit treatment).
    !> \brief Evaluate the implicit lava momentum contribution.
    !> \details Evaluate the part of temperature-dependent viscous drag and a yield-stress
    !> contribution assigned to the implicit treatment. Complex-valued arguments retain complex-step
    !> sensitivities where implemented; zero-valued implementations intentionally contribute no
    !> term.
    !> \param[in] moving Whether the local state is considered moving.
    !> \param[out] rate_x x component of momentum source rate.
    !> \param[out] rate_y y component of momentum source rate.

    !> \param[in] state Shared real physical state before the analytic substep.
    !> \brief Integrate the lava velocity source over a frozen-state substep.
    !> \details Evaluate the exact or time-integral-based decay of horizontal speed under
    !> temperature-dependent viscous drag and a yield-stress contribution. The velocity direction is
    !> retained and the updated components are returned without evolving the frozen thermodynamic
    !> and geometric inputs.
    !> \param[in] dt Duration of the analytic source substep [s].
    !> \param[in] speed_ratio Ratio of three-dimensional to horizontal speed.
    !> s^-1].
    !> \param[out] velocity_after Two-component horizontal velocity after source integration [m
    !> s^-1].

    ! Unsupported physical treatment, rejected by the source catalog.
END MODULE rheology_lava_2d