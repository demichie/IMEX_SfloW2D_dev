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
!> \ref eval_lahar_implicit_source evaluates the viscous/turbulent drag and
!> \ref eval_lahar_semi_impl_source evaluates the yield-stress term. The
!> real and complex overloads of \ref eval_lahar_fluid_viscosity implement
!> the same fluid-viscosity correlation, while
!> \ref eval_lahar_deposition_rate supplies the capped settling rate.
!> \ref solve_lahar_analytic_velocity solves the algebraic DIRK velocity
!> stage as a quadratic with depth, temperature and concentration explicit.
!********************************************************************************

MODULE rheology_lahar_2d

  USE parameters_2d, ONLY : wp, sp, four_thirds, neg_four_thirds
  USE constitutive_parameters_2d, ONLY : grav, inv_grav
  USE rheology_parameters_2d, ONLY : T_ref
  USE rheology_analytic_utils_2d, ONLY : solve_quadratic_in_interval,       &
       validate_analytic_velocity_inputs, fatal_analytic_input_error

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: alpha2, beta2, alpha1_coeff, beta1, Kappa, n_td, n_td2
  PUBLIC :: eval_lahar_fluid_viscosity, eval_lahar_implicit_source
  PUBLIC :: eval_lahar_semi_impl_source, eval_lahar_deposition_rate
  PUBLIC :: solve_lahar_analytic_velocity

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

CONTAINS

  !******************************************************************************
  !> \brief Solve a complete lahar DIRK stage analytically for velocity.
  !> \details Holds temperature, concentration, depth, density, and speed
  !> ratio fixed. The implicit viscous and turbulent terms and semi-implicit
  !> yield stress give a scalar quadratic in the stage speed:
  !> \f$\lambda a_2r^2+(1+\lambda a_1)r+\lambda c-r_0=0\f$.
  !> \param integration_step Effective DIRK integration step (dt * a_ii).
  !> \param flow_depth Explicit flow depth.
  !> \param temperature Explicit temperature.
  !> \param solid_fraction_sum Explicit total solid volume fraction.
  !> \param mixture_density Explicit mixture density.
  !> \param speed_ratio Explicit total-to-horizontal speed ratio.
  !> \param velocity_before Input horizontal velocity.
  !> \param velocity_after Returned horizontal velocity.
  !******************************************************************************
  SUBROUTINE solve_lahar_analytic_velocity(integration_step, flow_depth,   &
       temperature, solid_fraction_sum, mixture_density, speed_ratio,       &
       velocity_before, velocity_after)

    REAL(wp), INTENT(IN) :: integration_step, flow_depth, temperature
    REAL(wp), INTENT(IN) :: solid_fraction_sum, mixture_density, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2)
    REAL(wp), INTENT(OUT) :: velocity_after(2)

    REAL(wp) :: speed_before, speed_after, effective_depth, depth_threshold
    REAL(wp) :: fluid_viscosity, yield_strength
    REAL(wp) :: constant_drag, linear_drag, quadratic_drag
    LOGICAL :: root_found

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

    CALL solve_quadratic_in_interval(integration_step*quadratic_drag,       &
         1.0_wp + integration_step*linear_drag,                             &
         integration_step*constant_drag - speed_before, 0.0_wp, speed_before, &
         speed_after, root_found)
    IF (.NOT. root_found) speed_after = 0.0_wp
    velocity_after = (speed_after / speed_before) * velocity_before

  END SUBROUTINE solve_lahar_analytic_velocity

  !******************************************************************************
  !> \brief Evaluate carrier-fluid viscosity for real state variables.
  !> \details Converts temperature to degrees Celsius, evaluates the
  !> piecewise \f$\alpha_1(T)\f$ correlation, and multiplies by
  !> \f$\exp(\beta_1 C)\f$ for the total solid fraction.
  !> \param temperature Absolute temperature in kelvin.
  !> \param solid_fraction_sum Sum of solid volume fractions.
  !> \return Carrier-fluid dynamic viscosity.
  !******************************************************************************
  FUNCTION eval_lahar_fluid_viscosity_real(temperature, solid_fraction_sum)  &
       RESULT(fluid_viscosity)

    REAL(wp), INTENT(IN) :: temperature, solid_fraction_sum
    REAL(wp) :: fluid_viscosity
    REAL(wp) :: temperature_celsius, exponent_a, exponent_b, alpha1

    temperature_celsius = temperature - 273.15_wp
    IF (temperature_celsius .LT. 20.0_wp) THEN
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

  END FUNCTION eval_lahar_fluid_viscosity_real

  !******************************************************************************
  !> \brief Evaluate carrier-fluid viscosity for complex state variables.
  !> \details Complex overload of the viscosity correlation. The branch in
  !> the piecewise temperature law is selected from the real temperature;
  !> algebra within that branch remains complex for derivative evaluation.
  !> \param temperature Absolute temperature in kelvin.
  !> \param solid_fraction_sum Sum of solid volume fractions.
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

  !******************************************************************************
  !> \brief Evaluate the deposition-limited settling rate.
  !> \details Caps the settling flux by the amount of solid available in the
  !> flow during the current time step:
  !> \f$D=\min(Cw_s,hC/\Delta t)\f$.
  !> \param solid_fraction Solid volume fraction.
  !> \param settling_velocity Settling velocity.
  !> \param flow_depth Flow depth.
  !> \param time_step Time-step size.
  !> \return Deposition rate per unit area.
  !******************************************************************************
  PURE FUNCTION eval_lahar_deposition_rate(solid_fraction, settling_velocity, &
       flow_depth, time_step) RESULT(deposition_rate)

    REAL(wp), INTENT(IN) :: solid_fraction, settling_velocity
    REAL(wp), INTENT(IN) :: flow_depth, time_step
    REAL(wp) :: deposition_rate

    deposition_rate = MIN(solid_fraction * settling_velocity,                 &
         flow_depth * solid_fraction / time_step)

  END FUNCTION eval_lahar_deposition_rate

  !******************************************************************************
  !> \brief Evaluate the implicit viscous and turbulent lahar source.
  !> \details Computes the temperature- and concentration-dependent fluid
  !> viscosity, adds the viscous and turbulent contributions to the basal
  !> slope, and returns their x/y momentum sources. A depth floor regularizes
  !> the source in shallow states; complex arithmetic preserves state
  !> sensitivities in the implicit Jacobian.
  !> \param flow_depth Flow depth.
  !> \param velocity_x x component of horizontal velocity.
  !> \param velocity_y y component of horizontal velocity.
  !> \param speed Speed magnitude.
  !> \param horizontal_speed Horizontal speed magnitude.
  !> \param temperature Absolute temperature.
  !> \param mixture_density Mixture density.
  !> \param inverse_density Reciprocal mixture density.
  !> \param solid_fractions Solid volume fractions.
  !> \param moving Whether the state is moving.
  !> \param source_x Returned x momentum source.
  !> \param source_y Returned y momentum source.
  !******************************************************************************
  SUBROUTINE eval_lahar_implicit_source(flow_depth, velocity_x, velocity_y,    &
       speed, horizontal_speed, temperature, mixture_density, inverse_density, &
       solid_fractions, moving, source_x, source_y)

    COMPLEX(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y
    COMPLEX(wp), INTENT(IN) :: speed, horizontal_speed, temperature
    COMPLEX(wp), INTENT(IN) :: mixture_density, inverse_density
    COMPLEX(wp), INTENT(IN) :: solid_fractions(:)
    LOGICAL, INTENT(IN) :: moving
    COMPLEX(wp), INTENT(OUT) :: source_x, source_y

    COMPLEX(wp) :: fluid_viscosity, inverse_depth, viscous_slope
    COMPLEX(wp) :: turbulent_slope, total_slope, friction_term
    REAL(wp) :: depth_threshold

    source_x = CMPLX(0.0_wp, 0.0_wp, wp)
    source_y = CMPLX(0.0_wp, 0.0_wp, wp)
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
    source_x = -velocity_x * friction_term
    source_y = -velocity_y * friction_term

  END SUBROUTINE eval_lahar_implicit_source

  !******************************************************************************
  !> \brief Evaluate the semi-implicit lahar yield-stress source.
  !> \details Computes \f$Y=\alpha_2[\exp(\beta_2 C)-1]\f$ and returns
  !> \f$-Y\boldsymbol{u}_h/|\boldsymbol{u}_h|\f$ for moving flow.
  !> \param flow_depth Flow depth.
  !> \param velocity_x x component of horizontal velocity.
  !> \param velocity_y y component of horizontal velocity.
  !> \param horizontal_speed Horizontal speed magnitude.
  !> \param mixture_density Mixture density.
  !> \param solid_fraction_sum Sum of solid volume fractions.
  !> \param moving Whether the state is moving.
  !> \param source_x Returned x momentum source.
  !> \param source_y Returned y momentum source.
  !******************************************************************************
  SUBROUTINE eval_lahar_semi_impl_source(flow_depth, velocity_x, velocity_y, &
       horizontal_speed, mixture_density, solid_fraction_sum, moving,          &
       source_x, source_y)

    REAL(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y
    REAL(wp), INTENT(IN) :: horizontal_speed, mixture_density, solid_fraction_sum
    LOGICAL, INTENT(IN) :: moving
    REAL(wp), INTENT(OUT) :: source_x, source_y

    REAL(wp) :: yield_strength, yield_slope, friction_term
    REAL(wp), PARAMETER :: depth_threshold = 1.0E-20_wp

    source_x = 0.0_wp
    source_y = 0.0_wp
    IF (.NOT. moving) RETURN

    yield_strength = alpha2 * (EXP(beta2 * solid_fraction_sum) - 1.0_wp)
    IF (flow_depth .GT. depth_threshold) THEN
       yield_slope = yield_strength / (grav * mixture_density * flow_depth)
    ELSE
       yield_slope = yield_strength / (grav * mixture_density * depth_threshold)
    END IF
    friction_term = grav * mixture_density * flow_depth * yield_slope
    source_x = -friction_term * velocity_x / horizontal_speed
    source_y = -friction_term * velocity_y / horizontal_speed

  END SUBROUTINE eval_lahar_semi_impl_source

END MODULE rheology_lahar_2d
