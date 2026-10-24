!********************************************************************************
!> \brief Constant-strength plastic basal-friction model.
!> \details The model applies a speed-independent resisting strength
!> \f$\tau\f$ opposite to the horizontal flow direction:
!> \f[
!> \boldsymbol{S}=-\rho_m\tau\,\boldsymbol{u}_h/|\boldsymbol{u}|.
!> \f]
!> The supplied speed is the normalization used by the implementation. The
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

  SUBROUTINE eval_rheology_plastic_explicit()
    CALL fatal_error('rheology_plastic: explicit treatment is unsupported')
  END SUBROUTINE eval_rheology_plastic_explicit

  SUBROUTINE eval_rheology_plastic_semi_implicit(rate_x, rate_y)
    REAL(wp), INTENT(OUT) :: rate_x, rate_y
    rate_x = 0.0_wp
    rate_y = 0.0_wp
  END SUBROUTINE eval_rheology_plastic_semi_implicit

  SUBROUTINE eval_rheology_plastic_implicit(velocity_x, velocity_y, speed, &
       mixture_density, moving, rate_x, rate_y)
    COMPLEX(wp), INTENT(IN) :: velocity_x, velocity_y, speed, mixture_density
    LOGICAL, INTENT(IN) :: moving
    COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y
    rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
    rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
    IF (.NOT. moving) RETURN
    rate_x = -mixture_density * tau * (velocity_x / speed)
    rate_y = -mixture_density * tau * (velocity_y / speed)
  END SUBROUTINE eval_rheology_plastic_implicit

  SUBROUTINE solve_rheology_plastic_analytic(integration_step, flow_depth,  &
       speed_ratio, velocity_before, velocity_after)
    REAL(wp), INTENT(IN) :: integration_step, flow_depth, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2)
    REAL(wp), INTENT(OUT) :: velocity_after(2)
    REAL(wp) :: speed_before, speed_after
    velocity_after = velocity_before
    speed_before = SQRT(SUM(velocity_before**2))
    CALL validate_analytic_velocity_inputs(integration_step, flow_depth,    &
         speed_ratio, speed_before, 'plastic')
    IF (integration_step .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
    CALL solve_exact_velocity_flow(EXACT_FLOW_CONSTANT,                     &
         [tau/(flow_depth*speed_ratio)], integration_step, speed_before,    &
         speed_after)
    velocity_after = (speed_after / speed_before) * velocity_before
  END SUBROUTINE solve_rheology_plastic_analytic

  !******************************************************************************
  !> \brief Integrate the constant plastic-friction velocity source exactly.
  !> \details Holds depth and the speed-to-horizontal-speed ratio fixed. The
  !> speed is reduced by the constant friction acceleration until it reaches
  !> zero; the horizontal velocity direction is preserved:
  !> \f$r^+=\max(0,r_0-\lambda\tau/(h\beta))\f$.
  !> \param[in] integration_step Exact-flow interval, normally dt * a_ii.
  !> \param[in] flow_depth Explicit flow depth.
  !> \param[in] speed_ratio Explicit ratio of total to horizontal speed.
  !> \param[in] velocity_before Horizontal velocity before integration.
  !> \param[out] velocity_after Horizontal velocity after integration.
  !******************************************************************************
  !******************************************************************************
  !> \brief Evaluate the implicit plastic-friction source.
  !> \details Returns \f$-\rho_m\tau(u_x,u_y)/|\boldsymbol{u}|\f$ for a
  !> moving state and zero otherwise. Complex arithmetic preserves derivatives
  !> with respect to the state in the implicit Jacobian.
  !> \param[in] velocity_x x velocity [m s^-1].
  !> \param[in] velocity_y y velocity [m s^-1].
  !> \param[in] speed Speed used to normalize the direction [m s^-1].
  !> \param[in] mixture_density Mixture density \f$\rho_m\f$ [kg m^-3].
  !> \param[in] moving Whether the state is moving.
  !> \param[out] rate_x x conservative momentum source rate [kg m^-1 s^-2].
  !> \param[out] rate_y y conservative momentum source rate [kg m^-1 s^-2].
  !> \note Complex-valued state and rates preserve complex-step derivatives.
  !******************************************************************************
  ! Unsupported physical treatment, rejected by the source catalog.
  ! Valid zero component of the configured Newton treatment.
END MODULE rheology_plastic_2d