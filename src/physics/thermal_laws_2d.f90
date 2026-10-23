!********************************************************************************
!> \brief Thermal laws and exact lava cooling integration.
!> \details Provides thermal-energy rates for inlet injection and vaporization,
!>          the lava radiative/convective cooling law and its exact-flow
!>          integrator.
!********************************************************************************

MODULE thermal_laws_2d

  USE parameters_2d, ONLY : wp, rheology_id_lava
  USE diagnostics_2d, ONLY : fatal_error
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_inlet_thermal_energy_source
  PUBLIC :: eval_lava_cooling_energy_rate
  PUBLIC :: eval_rheology_cooling_energy_rate
  PUBLIC :: solve_lava_analytic_temperature
  PUBLIC :: advance_lava_cooling_temperature
  PUBLIC :: eval_vaporization_mass_rate
  PUBLIC :: eval_vaporization_energy_rate

  REAL(wp), PARAMETER :: vaporization_liquid_temperature = 290.0_wp
  REAL(wp), PARAMETER :: vaporization_boiling_temperature = 373.15_wp
  REAL(wp), PARAMETER :: vaporization_latent_heat = 2264705.0_wp
  REAL(wp), PARAMETER :: vaporization_liquid_heat_capacity = 4184.0_wp

CONTAINS

  !******************************************************************************
  !> \brief Compute the sensible thermal-energy source from injected material.
  !> \details Returns
  !>          \f$S_T=c_t\dot{H}\rho_m c_{p,m}T_s\f$,
  !>          the injected material's enthalpy rate per unit horizontal area.
  !> \param[in] time_coefficient Dimensionless temporal activation or step-averaged source weight.
  !> \param[in] injection_rate Equivalent injected layer-thickness rate [m s^-1].
  !> \param[in] mixture_density Density of the injected mixture [kg m^-3].
  !> \param[in] mixture_heat_capacity Specific heat capacity of the injected mixture
  !>                                  [J kg^-1 K^-1].
  !> \param[in] source_temperature Temperature of the injected material [K].
  !> \return Thermal-energy source per unit area [W m^-2].
  !******************************************************************************
  PURE FUNCTION eval_inlet_thermal_energy_source(time_coefficient, injection_rate, &
       mixture_density, mixture_heat_capacity, source_temperature) RESULT(source)

    REAL(wp), INTENT(IN) :: time_coefficient
    REAL(wp), INTENT(IN) :: injection_rate
    REAL(wp), INTENT(IN) :: mixture_density
    REAL(wp), INTENT(IN) :: mixture_heat_capacity
    REAL(wp), INTENT(IN) :: source_temperature
    REAL(wp) :: source

    source = time_coefficient * injection_rate * mixture_density *             &
         mixture_heat_capacity * source_temperature

  END FUNCTION eval_inlet_thermal_energy_source

  !******************************************************************************
  !> \brief Evaluate the lava radiative and convective cooling energy rate.
  !> \details Uses the net Stefan-type radiation and linear convection law
  !>          \f$S_T=-k_r(T^4-T_a^4)-k_c(T-T_a)\f$.
  !>          The complex temperature argument preserves complex-step
  !>          derivatives used by the implicit solver.
  !> \param[in] temperature Current lava temperature, real or complex [K].
  !> \param[in] radiative_coefficient Radiative coefficient k_r [W m^-2 K^-4].
  !> \param[in] convective_coefficient Convective coefficient k_c [W m^-2 K^-1].
  !> \param[in] ambient_temperature Ambient temperature T_a [K].
  !> \return Net cooling energy rate per unit area [W m^-2].
  !******************************************************************************
  PURE FUNCTION eval_lava_cooling_energy_rate(temperature, radiative_coefficient, &
       convective_coefficient, ambient_temperature) RESULT(source)

    COMPLEX(wp), INTENT(IN) :: temperature
    REAL(wp), INTENT(IN) :: radiative_coefficient
    REAL(wp), INTENT(IN) :: convective_coefficient
    REAL(wp), INTENT(IN) :: ambient_temperature
    COMPLEX(wp) :: source

    source = -radiative_coefficient * (temperature**4 - ambient_temperature**4) &
         - convective_coefficient * (temperature - ambient_temperature)

  END FUNCTION eval_lava_cooling_energy_rate

  !******************************************************************************
  !> \brief Evaluate the cooling energy rate associated with the selected rheology.
  !> \details Lava uses the radiative/convective cooling law; other models
  !>          contribute no rheology-specific cooling source.
  !> \param[in] model_id Selected rheology model ID.
  !> \param[in] temperature Current temperature [K].
  !> \param[in] radiative_coefficient Radiative coefficient [W m^-2 K^-4].
  !> \param[in] convective_coefficient Convective coefficient [W m^-2 K^-1].
  !> \param[in] ambient_temperature Ambient temperature [K].
  !> \return Rheology-specific thermal-energy rate per unit area [W m^-2].
  !******************************************************************************
  PURE FUNCTION eval_rheology_cooling_energy_rate(model_id, temperature,     &
       radiative_coefficient, convective_coefficient, ambient_temperature)   &
       RESULT(rate)

    INTEGER, INTENT(IN) :: model_id
    COMPLEX(wp), INTENT(IN) :: temperature
    REAL(wp), INTENT(IN) :: radiative_coefficient, convective_coefficient
    REAL(wp), INTENT(IN) :: ambient_temperature
    COMPLEX(wp) :: rate

    SELECT CASE (model_id)
    CASE (rheology_id_lava)
       rate = eval_lava_cooling_energy_rate(temperature, radiative_coefficient, &
            convective_coefficient, ambient_temperature)
    CASE DEFAULT
       rate = CMPLX(0.0_wp, 0.0_wp, wp)
    END SELECT

  END FUNCTION eval_rheology_cooling_energy_rate

  !******************************************************************************
  !> \brief Solve the exact-flow lava cooling temperature over a time interval.
  !> \details Separates the autonomous temperature equation and evaluates its
  !>          analytic time-temperature primitive. A bisection inversion finds
  !>          the temperature after the requested integration interval. The
  !>          areal thermal capacity and cooling coefficients are held fixed:
  !>          \f$C_A\frac{dT}{dt}=-k_r(T^4-T_a^4)-k_c(T-T_a)\f$.
  !>          The result is constrained to the closed interval between the
  !>          initial and ambient temperatures.
  !> \param[in] temperature Initial lava temperature T(0) [K].
  !> \param[in] integration_step Integration interval, normally dt*a_ii [s].
  !> \param[in] areal_heat_capacity Areal heat capacity C_A [J m^-2 K^-1].
  !> \param[in] radiative_coefficient Radiative cooling coefficient k_r [W m^-2 K^-4].
  !> \param[in] convective_coefficient Convective cooling coefficient k_c [W m^-2 K^-1].
  !> \param[in] ambient_temperature Ambient temperature T_a [K].
  !> \param[out] temperature_after Exact-flow temperature after integration_step [K].
  !> \note Requires finite inputs, positive initial and ambient temperatures, positive heat
  !>       capacity, nonnegative integration step, and nonnegative cooling coefficients.
  !******************************************************************************
  SUBROUTINE solve_lava_analytic_temperature(temperature, integration_step,   &
       areal_heat_capacity, radiative_coefficient, convective_coefficient,     &
       ambient_temperature, temperature_after)

    REAL(wp), INTENT(IN) :: temperature, integration_step, areal_heat_capacity
    REAL(wp), INTENT(IN) :: radiative_coefficient, convective_coefficient
    REAL(wp), INTENT(IN) :: ambient_temperature
    REAL(wp), INTENT(OUT) :: temperature_after

    REAL(wp) :: beta, cubic_root, lower, upper, midpoint
    REAL(wp) :: initial_primitive, target_primitive, midpoint_primitive
    REAL(wp) :: elapsed_scale, tolerance
    INTEGER :: iteration

    IF (.NOT. ieee_is_finite(temperature) .OR.                              &
         .NOT. ieee_is_finite(integration_step) .OR.                        &
         .NOT. ieee_is_finite(areal_heat_capacity) .OR.                    &
         .NOT. ieee_is_finite(radiative_coefficient) .OR.                  &
         .NOT. ieee_is_finite(convective_coefficient) .OR.                 &
         .NOT. ieee_is_finite(ambient_temperature)) THEN
       CALL fatal_error('Nonfinite input to exact lava cooling integrator')
    END IF
    IF (temperature .LE. 0.0_wp .OR. ambient_temperature .LE. 0.0_wp) THEN
       CALL fatal_error('Nonpositive temperature in exact lava cooling integrator')
    END IF
    IF (integration_step .LT. 0.0_wp .OR. areal_heat_capacity .LE. 0.0_wp) THEN
       CALL fatal_error('Invalid time or heat capacity in exact lava cooling integrator')
    END IF
    IF (radiative_coefficient .LT. 0.0_wp .OR.                            &
         convective_coefficient .LT. 0.0_wp) THEN
       CALL fatal_error('Negative cooling coefficient in exact lava cooling integrator')
    END IF

    temperature_after = temperature
    IF (integration_step .EQ. 0.0_wp .OR.                                   &
         (radiative_coefficient .EQ. 0.0_wp .AND.                          &
         convective_coefficient .EQ. 0.0_wp) .OR.                          &
         temperature .EQ. ambient_temperature) RETURN

    IF (radiative_coefficient .EQ. 0.0_wp) THEN
       temperature_after = ambient_temperature + (temperature -             &
            ambient_temperature) * EXP(-convective_coefficient *          &
            integration_step / areal_heat_capacity)
       RETURN
    END IF

    beta = convective_coefficient / (radiative_coefficient *               &
         ambient_temperature**3)
    IF (.NOT. ieee_is_finite(beta)) THEN
       CALL fatal_error('Cooling coefficient ratio overflow in lava integrator')
    END IF
    CALL lava_cooling_cubic_root(beta, cubic_root)

    initial_primitive = lava_cooling_time_primitive(                        &
         temperature / ambient_temperature, beta, cubic_root)
    elapsed_scale = radiative_coefficient * ambient_temperature**3 *        &
         integration_step / areal_heat_capacity
    IF (.NOT. ieee_is_finite(elapsed_scale)) THEN
       temperature_after = ambient_temperature
       RETURN
    END IF
    target_primitive = initial_primitive - elapsed_scale

    IF (temperature .GT. ambient_temperature) THEN
       lower = 1.0_wp
       upper = temperature / ambient_temperature
    ELSE
       lower = temperature / ambient_temperature
       upper = 1.0_wp
    END IF

    tolerance = 8.0_wp * EPSILON(1.0_wp) *                                  &
         MAX(1.0_wp, ABS(lower), ABS(upper))
    DO iteration = 1, 256
       midpoint = lower + 0.5_wp * (upper - lower)
       IF (midpoint .EQ. lower .OR. midpoint .EQ. upper .OR.                &
            ABS(upper - lower) .LE. tolerance) EXIT
       IF (midpoint .EQ. 1.0_wp) THEN
          IF (temperature .GT. ambient_temperature) THEN
             lower = midpoint
          ELSE
             upper = midpoint
          END IF
          CYCLE
       END IF

       midpoint_primitive = lava_cooling_time_primitive(midpoint, beta,       &
            cubic_root)
       IF (temperature .GT. ambient_temperature) THEN
          IF (midpoint_primitive .LT. target_primitive) THEN
             lower = midpoint
          ELSE
             upper = midpoint
          END IF
       ELSE
          IF (midpoint_primitive .GT. target_primitive) THEN
             lower = midpoint
          ELSE
             upper = midpoint
          END IF
       END IF
    END DO
    IF (iteration .GT. 256) THEN
       CALL fatal_error('Exact lava cooling inversion failed to converge')
    END IF

    midpoint = lower + 0.5_wp * (upper - lower)
    temperature_after = ambient_temperature * midpoint
    IF (temperature .GT. ambient_temperature) THEN
       temperature_after = MAX(ambient_temperature,                        &
            MIN(temperature, temperature_after))
    ELSE
       temperature_after = MIN(ambient_temperature,                        &
            MAX(temperature, temperature_after))
    END IF

  END SUBROUTINE solve_lava_analytic_temperature

  !******************************************************************************
  !> \brief Backward-compatible name for the exact-flow lava cooling solver.
  !> \param[in] temperature Initial lava temperature [K].
  !> \param[in] integration_step Integration interval [s].
  !> \param[in] areal_heat_capacity Lava heat capacity per unit area [J m^-2 K^-1].
  !> \param[in] radiative_coefficient Radiative coefficient [W m^-2 K^-4].
  !> \param[in] convective_coefficient Convective coefficient [W m^-2 K^-1].
  !> \param[in] ambient_temperature Ambient temperature [K].
  !> \param[out] temperature_after Exact-flow temperature after the interval [K].
  !******************************************************************************
  SUBROUTINE advance_lava_cooling_temperature(temperature, integration_step,  &
       areal_heat_capacity, radiative_coefficient, convective_coefficient,    &
       ambient_temperature, temperature_after)

    REAL(wp), INTENT(IN) :: temperature, integration_step, areal_heat_capacity
    REAL(wp), INTENT(IN) :: radiative_coefficient, convective_coefficient
    REAL(wp), INTENT(IN) :: ambient_temperature
    REAL(wp), INTENT(OUT) :: temperature_after

    CALL solve_lava_analytic_temperature(temperature, integration_step,      &
         areal_heat_capacity, radiative_coefficient, convective_coefficient,  &
         ambient_temperature, temperature_after)

  END SUBROUTINE advance_lava_cooling_temperature

  !******************************************************************************
  !> \brief Compute the vaporization mass rate from hot-solid deposition.
  !> \details Deposition heat above the boiling point supplies liquid warming
  !>          and latent vaporization energy. The generated steam rate is
  !>          multiplied by the configured steam fraction.
  !>          \f$\dot{m}_v=f_v\max(0,T_f-T_b)\sum_i\rho_i c_iD_i/
  !>          [c_l(T_b-T_l)+L_v]\f$.
  !> \param[in] flow_temperature Flow temperature T_f [K].
  !> \param[in] solid_density Per-class solid densities rho_i [kg m^-3].
  !> \param[in] solid_heat_capacity Per-class solid specific heat capacities c_i [J kg^-1 K^-1].
  !> \param[in] deposition_rate Per-class deposition thickness rates D_i [m s^-1].
  !> \param[in] steam_fraction Fraction f_v of the available deposition heat assigned to vapor.
  !> \return Generated vapor mass rate per unit area [kg m^-2 s^-1].
  !******************************************************************************
  PURE FUNCTION eval_vaporization_mass_rate(flow_temperature, solid_density,   &
       solid_heat_capacity, deposition_rate, steam_fraction) RESULT(rate)

    REAL(wp), INTENT(IN) :: flow_temperature
    REAL(wp), INTENT(IN) :: solid_density(:)
    REAL(wp), INTENT(IN) :: solid_heat_capacity(:)
    REAL(wp), INTENT(IN) :: deposition_rate(:)
    REAL(wp), INTENT(IN) :: steam_fraction
    REAL(wp) :: rate

    rate = steam_fraction * MAX(0.0_wp, flow_temperature -                    &
         vaporization_boiling_temperature) *                                 &
         SUM(solid_density * solid_heat_capacity * deposition_rate) /         &
         (vaporization_liquid_heat_capacity *                                 &
         (vaporization_boiling_temperature -                                  &
         vaporization_liquid_temperature) + vaporization_latent_heat)

  END FUNCTION eval_vaporization_mass_rate

  !******************************************************************************
  !> \brief Compute the sensible-energy rate carried by generated vapor.
  !> \details Evaluates \f$\dot{E}_v=\dot{m}_vc_gT_b\f$ at the boiling
  !>          temperature.
  !> \param[in] vapor_mass_rate Generated vapor mass rate per unit area [kg m^-2 s^-1].
  !> \param[in] gas_heat_capacity Steam specific heat capacity c_g [J kg^-1 K^-1].
  !> \return Vapor sensible-energy rate per unit area [W m^-2].
  !******************************************************************************
  PURE FUNCTION eval_vaporization_energy_rate(vapor_mass_rate,               &
       gas_heat_capacity) RESULT(rate)

    REAL(wp), INTENT(IN) :: vapor_mass_rate, gas_heat_capacity
    REAL(wp) :: rate

    rate = vapor_mass_rate * gas_heat_capacity *                             &
         vaporization_boiling_temperature

  END FUNCTION eval_vaporization_energy_rate

  !******************************************************************************
  !> \brief Find the unique real root below -1 of the normalized cooling cubic.
  !> \details For beta >= 0, solves
  !>          \f$r^3+r^2+r+(1+\beta)=0\f$ by bisection on the unique real root
  !>          \f$r<-1\f$.
  !>          This root factors the normalized quartic cooling polynomial.
  !> \param[in] beta Nonnegative ratio k_c/(k_r*T_a^3).
  !> \param[out] root Unique real cubic root below -1.
  !> \note For beta=0, the root is exactly -1.
  !******************************************************************************
  SUBROUTINE lava_cooling_cubic_root(beta, root)

    REAL(wp), INTENT(IN) :: beta
    REAL(wp), INTENT(OUT) :: root

    REAL(wp) :: lower, upper, midpoint, scale, polynomial
    INTEGER :: iteration

    IF (beta .EQ. 0.0_wp) THEN
       root = -1.0_wp
       RETURN
    END IF

    scale = MAX(1.0_wp, (1.0_wp + beta)**(1.0_wp / 3.0_wp))
    lower = -2.0_wp
    upper = -1.0_wp / scale
    DO iteration = 1, 256
       midpoint = lower + 0.5_wp * (upper - lower)
       IF (midpoint .EQ. lower .OR. midpoint .EQ. upper) EXIT
       polynomial = midpoint**3 + midpoint**2/scale + midpoint/scale**2 +  &
            1.0_wp
       IF (polynomial .LT. 0.0_wp) THEN
          lower = midpoint
       ELSE
          upper = midpoint
       END IF
    END DO
    IF (iteration .GT. 256) THEN
       CALL fatal_error('Cooling cubic root calculation failed to converge')
    END IF
    root = scale * (lower + 0.5_wp * (upper - lower))

  END SUBROUTINE lava_cooling_cubic_root

  !******************************************************************************
  !> \brief Evaluate the analytic primitive for the normalized cooling law.
  !> \details Returns an antiderivative of
  !>          \f$1/[(1-x)(x^3+x^2+x+1+\beta)]\f$, used to convert temperature
  !>          change into elapsed dimensionless time in the lava cooling ODE.
  !> \param[in] x Positive temperature ratio T/T_a, with x != 1.
  !> \param[in] beta Nonnegative ratio k_c/(k_r*T_a^3).
  !> \param[in] cubic_root Unique real root below -1 of the normalized cubic.
  !> \return Analytic primitive value, up to an irrelevant additive constant.
  !******************************************************************************
  PURE FUNCTION lava_cooling_time_primitive(x, beta, cubic_root) RESULT(value)

    REAL(wp), INTENT(IN) :: x, beta, cubic_root
    REAL(wp) :: value

    REAL(wp) :: quadratic_linear, quadratic_constant, discriminant
    REAL(wp) :: partial_a, partial_b, partial_c, partial_d, residual

    quadratic_linear = 1.0_wp + cubic_root
    quadratic_constant = 1.0_wp + cubic_root + cubic_root**2
    discriminant = 4.0_wp * quadratic_constant - quadratic_linear**2

    partial_a = 1.0_wp / (4.0_wp + beta)
    partial_b = 1.0_wp / ((cubic_root - 1.0_wp) *                        &
         (3.0_wp * cubic_root**2 + 2.0_wp * cubic_root + 1.0_wp))
    partial_c = -partial_a - partial_b
    partial_d = -partial_a * (cubic_root + 2.0_wp) -                     &
         partial_b * (2.0_wp * cubic_root + 1.0_wp)
    residual = partial_d - 0.5_wp * partial_c * quadratic_linear

    value = partial_a * LOG(ABS(x - 1.0_wp)) +                            &
         partial_b * LOG(ABS(x - cubic_root)) +                           &
         0.5_wp * partial_c * LOG(x**2 + quadratic_linear*x +             &
         quadratic_constant) + 2.0_wp * residual / SQRT(discriminant) *   &
         ATAN((2.0_wp*x + quadratic_linear) / SQRT(discriminant))

  END FUNCTION lava_cooling_time_primitive

END MODULE thermal_laws_2d
