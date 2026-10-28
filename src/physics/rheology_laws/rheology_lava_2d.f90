!********************************************************************************
!> \brief Temperature-dependent viscous and yield-stress lava rheology.
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
!> source exactly with temperature and other non-velocity inputs explicit.
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
  !> \param[in] temperature Lava temperature [K].
  !> \return Kinematic viscosity [m2/s].
  PURE REAL(wp) FUNCTION lava_kinematic_viscosity(temperature) RESULT(nu)
    REAL(wp), INTENT(IN) :: temperature
    nu=nu_ref*EXP(-visc_par*(temperature-T_ref))
  END FUNCTION lava_kinematic_viscosity


  SUBROUTINE eval_rheology_lava_explicit( &
       flow_depth, temperature, velocity_x, velocity_y, speed, horizontal_speed, &
       speed_squared, mixture_density, inverse_density, reduced_gravity, &
       gravity_coefficient, solid_fractions, curvature_term, curvature_xx, curvature_xy, &
       curvature_yy, use_pore_pressure, excess_pore_pressure, moving, use_stochastic, &
       stochastic_value, rate )
    REAL(wp), INTENT(IN) :: flow_depth, temperature, velocity_x, velocity_y
    REAL(wp), INTENT(IN) :: speed, horizontal_speed, speed_squared
    REAL(wp), INTENT(IN) :: mixture_density, inverse_density, reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient, solid_fractions(:)
    REAL(wp), INTENT(IN) :: curvature_xx, curvature_xy, curvature_yy
    REAL(wp), INTENT(IN) :: excess_pore_pressure, stochastic_value
    LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure, moving, use_stochastic
    REAL(wp), INTENT(OUT) :: rate(2)
    REAL(wp) :: sx, sy
    COMPLEX(wp) :: ix, iy
    CALL eval_rheology_lava_semi_implicit( &
       velocity_x, velocity_y, speed, sx, sy )
    CALL eval_rheology_lava_implicit( &
       CMPLX(flow_depth,0.0_wp,wp), CMPLX(temperature,0.0_wp,wp), &
       CMPLX(mixture_density,0.0_wp,wp), CMPLX(velocity_x,0.0_wp,wp), &
       CMPLX(velocity_y,0.0_wp,wp), moving, ix, iy )
    rate=[sx,sy]+REAL([ix,iy],wp)
  END SUBROUTINE eval_rheology_lava_explicit

  SUBROUTINE eval_rheology_lava_semi_implicit(velocity_x, velocity_y, speed, rate_x, rate_y)
    REAL(wp), INTENT(IN) :: velocity_x, velocity_y, speed
    REAL(wp), INTENT(OUT) :: rate_x, rate_y
    rate_x = 0.0_wp
    rate_y = 0.0_wp
    IF (speed .LE. 0.0_wp) RETURN
    rate_x = -tau0 * velocity_x / speed
    rate_y = -tau0 * velocity_y / speed
  END SUBROUTINE eval_rheology_lava_semi_implicit

  SUBROUTINE eval_rheology_lava_implicit(flow_depth, temperature, mixture_density, &
       velocity_x, velocity_y, moving, rate_x, rate_y)
    COMPLEX(wp), INTENT(IN) :: flow_depth, temperature, mixture_density
    COMPLEX(wp), INTENT(IN) :: velocity_x, velocity_y
    LOGICAL, INTENT(IN) :: moving
    COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y
    COMPLEX(wp) :: gamma
    REAL(wp), PARAMETER :: h_threshold = 1.0E-10_wp
    rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
    rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
    IF (REAL(flow_depth, wp) .GT. h_threshold) THEN
       gamma = 3.0_wp * nu_ref / flow_depth * EXP(-visc_par * (temperature - T_ref))
    ELSE
       gamma = 3.0_wp * nu_ref / h_threshold * EXP(-visc_par * (temperature - T_ref))
    END IF
    IF (moving) THEN
       rate_x = -mixture_density * gamma * velocity_x
       rate_y = -mixture_density * gamma * velocity_y
    END IF
  END SUBROUTINE eval_rheology_lava_implicit

  SUBROUTINE solve_rheology_lava_analytic(integration_step, flow_depth,     &
       temperature, mixture_density, speed_ratio, velocity_before,           &
       velocity_after)
    REAL(wp), INTENT(IN) :: integration_step, flow_depth, temperature
    REAL(wp), INTENT(IN) :: mixture_density, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2)
    REAL(wp), INTENT(OUT) :: velocity_after(2)
    REAL(wp), PARAMETER :: h_threshold = 1.0E-10_wp
    REAL(wp) :: speed_before, speed_after, effective_depth, gamma
    REAL(wp) :: linear_drag, yield_drag
    velocity_after = velocity_before
    speed_before = SQRT(SUM(velocity_before**2))
    CALL validate_analytic_velocity_inputs(integration_step, flow_depth,    &
         speed_ratio, speed_before, 'lava')
    IF (integration_step .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
    IF (mixture_density .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive mixture density in lava analytic solver')
    END IF
    effective_depth = MAX(flow_depth, h_threshold)
    gamma = 3.0_wp * nu_ref / effective_depth *                             &
         EXP(-visc_par * (temperature - T_ref))
    linear_drag = gamma / flow_depth
    yield_drag = tau0 / (mixture_density * flow_depth * speed_ratio)
    CALL solve_exact_velocity_flow(EXACT_FLOW_AFFINE,                       &
         [yield_drag, linear_drag], integration_step, speed_before,         &
         speed_after)
    velocity_after = (speed_after / speed_before) * velocity_before
  END SUBROUTINE solve_rheology_lava_analytic

  !******************************************************************************
  !> \brief Integrate the lava-rheology velocity source exactly.
  !> \details Holds depth, temperature, density, and speed ratio fixed. The
  !> affine speed-decay law is integrated exactly, with stopping at zero.
  !> \param[in] integration_step Exact-flow interval [s], normally dt * a_ii.
  !> \param[in] flow_depth Explicit flow depth [m].
  !> \param[in] temperature Explicit temperature [K].
  !> \param[in] mixture_density Explicit mixture density [kg m^-3].
  !> \param[in] speed_ratio Explicit total-to-horizontal speed ratio.
  !> \param[in] velocity_before Input horizontal velocity [m s^-1].
  !> \param[out] velocity_after Returned horizontal velocity [m s^-1].
  !******************************************************************************
  !******************************************************************************
  !> \brief Evaluate the implicit temperature-dependent viscous drag.
  !> \details Computes the depth-regularized drag coefficient and returns the
  !> x/y momentum sources for a moving state. Complex inputs preserve
  !> temperature, depth and velocity sensitivities in the implicit Jacobian.
  !> \param[in] flow_depth Flow depth [m].
  !> \param[in] temperature Lava temperature [K].
  !> \param[in] mixture_density Mixture density [kg m^-3].
  !> \param[in] velocity_x x velocity [m s^-1].
  !> \param[in] velocity_y y velocity [m s^-1].
  !> \param[in] moving Whether the state is moving.
  !> \param[out] rate_x x conservative momentum source rate [kg m^-1 s^-2].
  !> \param[out] rate_y y conservative momentum source rate [kg m^-1 s^-2].
  !> \note Complex-valued state and rates preserve complex-step derivatives.
  !******************************************************************************
  !******************************************************************************
  !> \brief Evaluate the semi-implicit lava yield-stress momentum rate.
  !> \details Returns \f$-\tau_0(u_x,u_y)/|\boldsymbol{u}|\f$ for positive
  !> speed and zero otherwise.
  !> \param[in] velocity_x x velocity [m s^-1].
  !> \param[in] velocity_y y velocity [m s^-1].
  !> \param[in] speed Speed magnitude [m s^-1].
  !> \param[out] rate_x x momentum source rate.
  !> \param[out] rate_y y momentum source rate.
  !******************************************************************************
  ! Unsupported physical treatment, rejected by the source catalog.
END MODULE rheology_lava_2d