!********************************************************************************
!> \brief Lahar rheology: temperature-dependent viscosity, turbulent drag and yield stress.
!> \details With \f$T_C=T-273.15\f$, the carrier-fluid viscosity is
!> \f$\eta_f(T,C)=\alpha_1(T)\exp(\beta_1 C)\f$. Its piecewise temperature
!> factor is
!> \f[
!> \alpha_1(T)=\alpha_{1,\rm coeff}
!> \begin{cases}
!> 10^{-3}10^{1301/[998.333+8.1855(T_C-20)+0.00585(T_C-20)^2]-1.30223},
!>       & T_C<20,\\
!> 1.002\,10^{-3}10^{[1.3272(20-T_C)-0.001053(T_C-20)^2]/(T_C+105)},
!>       & T_C\ge20.
!> \end{cases}
!> \f]
!> The implicit momentum source
!> combines viscous and turbulent slope contributions,
!> \f[
!> S_v=\frac{\kappa\eta_f|\boldsymbol{u}|}
!>              {8\rho_m g h^2},\qquad
!> S_t=n_{td}^2|\boldsymbol{u}|^2h^{-4/3},\qquad
!> \boldsymbol{S}=-\frac{g\rho_mh(S_v+S_t)}
!>                       {|\boldsymbol{u}_h|}\boldsymbol{u}_h.
!> \f]
!> The semi-implicit yield term uses
!> \f$Y=\alpha_2[\exp(\beta_2 C)-1]\f$ and
!> \f$\boldsymbol{S}_y=-Y\boldsymbol{u}_h/|\boldsymbol{u}_h|\f$. The
!> \ref eval_rheology_lahar_implicit evaluates the viscous/turbulent drag and
!> \ref eval_rheology_lahar_semi_implicit evaluates the yield-stress term. The
!> real and complex overloads of \ref eval_lahar_fluid_viscosity implement
!> the same fluid-viscosity correlation.
!> \ref solve_rheology_lahar_analytic integrates the frozen polynomial
!> velocity source with depth, temperature and concentration explicit.
!********************************************************************************
MODULE rheology_lahar_2d
  USE diagnostics_2d, ONLY : fatal_error
  USE parameters_2d, ONLY : wp, sp, four_thirds, neg_four_thirds
  USE constitutive_parameters_2d, ONLY : grav, inv_grav
  USE rheology_parameters_2d, ONLY : T_ref
  USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs, &
       fatal_analytic_input_error, solve_exact_velocity_flow,                &
       EXACT_FLOW_QUADRATIC_POLYNOMIAL
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_rheology_lahar_explicit
  PUBLIC :: eval_rheology_lahar_semi_implicit
  PUBLIC :: eval_rheology_lahar_implicit
  PUBLIC :: solve_rheology_lahar_analytic
  PUBLIC :: alpha2, beta2, alpha1_coeff, beta1, Kappa, n_td, n_td2
  PUBLIC :: eval_lahar_fluid_viscosity
  REAL(wp) :: alpha2
  REAL(wp) :: beta2
  REAL(wp) :: alpha1_coeff
  REAL(wp) :: beta1
  REAL(wp) :: Kappa
  REAL(wp) :: n_td
  REAL(wp) :: n_td2
  INTERFACE eval_lahar_fluid_viscosity
     MODULE PROCEDURE eval_lahar_fluid_viscosity_real
     MODULE PROCEDURE eval_lahar_fluid_viscosity_complex
  END INTERFACE eval_lahar_fluid_viscosity
  ! Uniform source treatment entry points. Existing kernels remain unchanged.
CONTAINS

  SUBROUTINE eval_rheology_lahar_explicit( &
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
    CALL eval_rheology_lahar_semi_implicit( &
       flow_depth, velocity_x, velocity_y, horizontal_speed, mixture_density, &
       SUM(solid_fractions), moving, sx, sy )
    CALL eval_rheology_lahar_implicit( &
       CMPLX(flow_depth,0.0_wp,wp), CMPLX(velocity_x,0.0_wp,wp), &
       CMPLX(velocity_y,0.0_wp,wp), CMPLX(speed,0.0_wp,wp), &
       CMPLX(horizontal_speed,0.0_wp,wp), CMPLX(temperature,0.0_wp,wp), &
       CMPLX(mixture_density,0.0_wp,wp), CMPLX(inverse_density,0.0_wp,wp), &
       CMPLX(solid_fractions,0.0_wp,wp), moving, ix, iy )
    rate=[sx,sy]+REAL([ix,iy],wp)
  END SUBROUTINE eval_rheology_lahar_explicit

  SUBROUTINE eval_rheology_lahar_semi_implicit(flow_depth, velocity_x, velocity_y, &
       horizontal_speed, mixture_density, solid_fraction_sum, moving,          &
       rate_x, rate_y)
    REAL(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y
    REAL(wp), INTENT(IN) :: horizontal_speed, mixture_density, solid_fraction_sum
    LOGICAL, INTENT(IN) :: moving
    REAL(wp), INTENT(OUT) :: rate_x, rate_y
    REAL(wp) :: yield_strength, yield_slope, friction_term
    REAL(wp), PARAMETER :: depth_threshold = 1.0E-20_wp
    rate_x = 0.0_wp
    rate_y = 0.0_wp
    IF (.NOT. moving) RETURN
    yield_strength = alpha2 * (EXP(beta2 * solid_fraction_sum) - 1.0_wp)
    IF (flow_depth .GT. depth_threshold) THEN
       yield_slope = yield_strength / (grav * mixture_density * flow_depth)
    ELSE
       yield_slope = yield_strength / (grav * mixture_density * depth_threshold)
    END IF
    friction_term = grav * mixture_density * flow_depth * yield_slope
    rate_x = -friction_term * velocity_x / horizontal_speed
    rate_y = -friction_term * velocity_y / horizontal_speed
  END SUBROUTINE eval_rheology_lahar_semi_implicit

  SUBROUTINE eval_rheology_lahar_implicit(flow_depth, velocity_x, velocity_y,    &
       speed, horizontal_speed, temperature, mixture_density, inverse_density, &
       solid_fractions, moving, rate_x, rate_y)
    COMPLEX(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y
    COMPLEX(wp), INTENT(IN) :: speed, horizontal_speed, temperature
    COMPLEX(wp), INTENT(IN) :: mixture_density, inverse_density
    COMPLEX(wp), INTENT(IN) :: solid_fractions(:)
    LOGICAL, INTENT(IN) :: moving
    COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y
    COMPLEX(wp) :: fluid_viscosity, inverse_depth, viscous_slope
    COMPLEX(wp) :: turbulent_slope, total_slope, friction_term
    REAL(wp) :: depth_threshold
    rate_x = CMPLX(0.0_wp, 0.0_wp, wp)
    rate_y = CMPLX(0.0_wp, 0.0_wp, wp)
    IF (wp .EQ. sp) THEN
       depth_threshold = 1.0E-10_wp
    ELSE
       depth_threshold = 1.0E-20_wp
    END IF
    fluid_viscosity = eval_lahar_fluid_viscosity(temperature, SUM(solid_fractions))
    IF (REAL(flow_depth, wp) .GT. depth_threshold) THEN
       inverse_depth = 1.0_wp / flow_depth
       viscous_slope = Kappa * fluid_viscosity * speed * 0.125_wp *          &
            inverse_density * inv_grav * inverse_depth**2
       turbulent_slope = n_td2 * speed**2 * inverse_depth**four_thirds
    ELSE
       viscous_slope = Kappa * fluid_viscosity * speed /                    &
            (8.0_wp * mixture_density * grav * depth_threshold**2)
       turbulent_slope = n_td2 * speed**2 * depth_threshold**neg_four_thirds
    END IF
    total_slope = viscous_slope + turbulent_slope
    IF (.NOT. moving) RETURN
    friction_term = grav * mixture_density * flow_depth * total_slope / horizontal_speed
    rate_x = -velocity_x * friction_term
    rate_y = -velocity_y * friction_term
  END SUBROUTINE eval_rheology_lahar_implicit

  SUBROUTINE solve_rheology_lahar_analytic(integration_step, flow_depth,   &
       temperature, solid_fraction_sum, mixture_density, speed_ratio,       &
       velocity_before, velocity_after)
    REAL(wp), INTENT(IN) :: integration_step, flow_depth, temperature
    REAL(wp), INTENT(IN) :: solid_fraction_sum, mixture_density, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2)
    REAL(wp), INTENT(OUT) :: velocity_after(2)
    REAL(wp) :: speed_before, speed_after, effective_depth, depth_threshold
    REAL(wp) :: fluid_viscosity, yield_strength
    REAL(wp) :: constant_drag, linear_drag, quadratic_drag
    velocity_after = velocity_before
    speed_before = SQRT(SUM(velocity_before**2))
    CALL validate_analytic_velocity_inputs(integration_step, flow_depth,    &
         speed_ratio, speed_before, 'lahar')
    IF (integration_step .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
    IF (mixture_density .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive mixture density in lahar analytic solver')
    END IF
    IF (wp .EQ. sp) THEN
       depth_threshold = 1.0E-10_wp
    ELSE
       depth_threshold = 1.0E-20_wp
    END IF
    effective_depth = MAX(flow_depth, depth_threshold)
    fluid_viscosity = eval_lahar_fluid_viscosity(temperature, solid_fraction_sum)
    yield_strength = alpha2 * (EXP(beta2 * solid_fraction_sum) - 1.0_wp)
    constant_drag = yield_strength / (mixture_density *                    &
         MAX(flow_depth, 1.0E-20_wp))
    linear_drag = Kappa * fluid_viscosity * speed_ratio /                   &
         (8.0_wp * mixture_density * effective_depth**2)
    quadratic_drag = grav * n_td2 * speed_ratio**2 *                       &
         effective_depth**neg_four_thirds
    CALL solve_exact_velocity_flow(EXACT_FLOW_QUADRATIC_POLYNOMIAL,         &
         [constant_drag, linear_drag, quadratic_drag], integration_step,    &
         speed_before, speed_after)
    velocity_after = (speed_after / speed_before) * velocity_before
  END SUBROUTINE solve_rheology_lahar_analytic

  !******************************************************************************
  !> \brief Integrate the lahar velocity source exactly.
  !> \details Holds temperature, concentration, depth, density, and speed
  !> ratio fixed while integrating the polynomial speed-decay law.
  !> \param[in] integration_step Exact-flow integration interval [s].
  !> \param[in] flow_depth Explicit flow depth [m].
  !> \param[in] temperature Explicit absolute temperature [K].
  !> \param[in] solid_fraction_sum Explicit total solid volume fraction.
  !> \param[in] mixture_density Explicit mixture density [kg m^-3].
  !> \param[in] speed_ratio Explicit total-to-horizontal speed ratio.
  !> \param[in] velocity_before Input horizontal velocity [m s^-1].
  !> \param[out] velocity_after Returned horizontal velocity [m s^-1].
  !******************************************************************************
  !******************************************************************************
  !> \brief Evaluate the implicit viscous and turbulent lahar source.
  !> \details Computes the temperature- and concentration-dependent fluid
  !> viscosity, adds the viscous and turbulent contributions to the basal
  !> slope, and returns their x/y momentum sources. A depth floor regularizes
  !> the source in shallow states; complex arithmetic preserves state
  !> sensitivities in the implicit Jacobian.
  !> \param[in] flow_depth Flow depth [m].
  !> \param[in] velocity_x x velocity [m s^-1].
  !> \param[in] velocity_y y velocity [m s^-1].
  !> \param[in] speed Speed magnitude [m s^-1].
  !> \param[in] horizontal_speed Horizontal speed magnitude [m s^-1].
  !> \param[in] temperature Absolute temperature [K].
  !> \param[in] mixture_density Mixture density [kg m^-3].
  !> \param[in] inverse_density Reciprocal mixture density [m^3 kg^-1].
  !> \param[in] solid_fractions Solid volume fractions by class.
  !> \param[in] moving Whether the state is moving.
  !> \param[out] rate_x x conservative momentum source rate [kg m^-1 s^-2].
  !> \param[out] rate_y y conservative momentum source rate [kg m^-1 s^-2].
  !******************************************************************************
  !******************************************************************************
  !> \brief Evaluate the semi-implicit lahar yield-stress momentum rate.
  !> \details Computes \f$Y=\alpha_2[\exp(\beta_2 C)-1]\f$ and returns
  !> \f$-Y\boldsymbol{u}_h/|\boldsymbol{u}_h|\f$ for moving flow.
  !> \param[in] flow_depth Flow depth [m].
  !> \param[in] velocity_x x velocity [m s^-1].
  !> \param[in] velocity_y y velocity [m s^-1].
  !> \param[in] horizontal_speed Horizontal speed magnitude [m s^-1].
  !> \param[in] mixture_density Mixture density [kg m^-3].
  !> \param[in] solid_fraction_sum Sum of solid volume fractions.
  !> \param[in] moving Whether the state is moving.
  !> \param[out] rate_x x conservative momentum source rate [kg m^-1 s^-2].
  !> \param[out] rate_y y conservative momentum source rate [kg m^-1 s^-2].
  !******************************************************************************
  !******************************************************************************
  !> \brief Evaluate carrier-fluid viscosity for real state variables.
  !> \details Wraps the complex correlation for real state variables.
  !> \param[in] temperature Absolute temperature [K].
  !> \param[in] solid_fraction_sum Sum of solid volume fractions.
  !> \return Carrier-fluid dynamic viscosity.
  !******************************************************************************
  FUNCTION eval_lahar_fluid_viscosity_real(temperature, solid_fraction_sum)  &
       RESULT(fluid_viscosity)
    REAL(wp), INTENT(IN) :: temperature, solid_fraction_sum
    REAL(wp) :: fluid_viscosity
    COMPLEX(wp) :: complex_viscosity
    complex_viscosity = eval_lahar_fluid_viscosity_complex(                  &
         CMPLX(temperature,0.0_wp,wp),CMPLX(solid_fraction_sum,0.0_wp,wp))
    fluid_viscosity = REAL(complex_viscosity,wp)
  END FUNCTION eval_lahar_fluid_viscosity_real
  !******************************************************************************
  !> \brief Evaluate carrier-fluid viscosity for complex state variables.
  !> \details Complex overload of the viscosity correlation. The branch in
  !> the piecewise temperature law is selected from the real temperature;
  !> algebra within that branch remains complex for derivative evaluation.
  !> \param[in] temperature Absolute temperature [K].
  !> \param[in] solid_fraction_sum Sum of solid volume fractions.
  !> \return Complex carrier-fluid dynamic viscosity.
  !******************************************************************************
  FUNCTION eval_lahar_fluid_viscosity_complex(temperature, solid_fraction_sum) &
       RESULT(fluid_viscosity)
    COMPLEX(wp), INTENT(IN) :: temperature, solid_fraction_sum
    COMPLEX(wp) :: fluid_viscosity
    COMPLEX(wp) :: temperature_celsius, exponent_a, exponent_b, alpha1
    temperature_celsius = temperature - 273.15_wp
    IF (REAL(temperature_celsius, wp) .LT. 20.0_wp) THEN
       exponent_a = 1301.0_wp / (998.333_wp + 8.1855_wp *                    &
            (temperature_celsius - 20.0_wp) + 0.00585_wp *                   &
            (temperature_celsius - 20.0_wp)**2) - 1.30223_wp
       alpha1 = alpha1_coeff * 1.0E-3_wp * 10.0_wp**exponent_a
    ELSE
       exponent_b = (1.3272_wp * (20.0_wp - temperature_celsius)             &
            - 0.001053_wp * (temperature_celsius - 20.0_wp)**2)              &
            / (temperature_celsius + 105.0_wp)
       alpha1 = alpha1_coeff * 1.002E-3_wp * 10.0_wp**exponent_b
    END IF
    fluid_viscosity = alpha1 * EXP(beta1 * solid_fraction_sum)
  END FUNCTION eval_lahar_fluid_viscosity_complex
  ! Unsupported physical treatment, rejected by the source catalog.
END MODULE rheology_lahar_2d