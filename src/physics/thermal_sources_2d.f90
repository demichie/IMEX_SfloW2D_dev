!********************************************************************************
!> \brief Thermal source laws and lava cooling stage solvers.
!> \details Provides thermal-energy sources for inlet injection, mass exchange
!>          and vaporization, as well as the lava radiative/convective cooling
!>          law. It provides both an exact-flow integrator for the autonomous
!>          cooling ODE and a closed-form solution of its scalar DIRK stage
!>          equation; callers select which contract applies.
!********************************************************************************

MODULE thermal_sources_2d

  USE parameters_2d, ONLY : wp
  USE diagnostics_2d, ONLY : fatal_error
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_inlet_thermal_energy_source
  PUBLIC :: eval_lava_cooling_source
  PUBLIC :: advance_lava_cooling_temperature
  PUBLIC :: solve_lava_cooling_stage_temperature
  PUBLIC :: eval_mass_exchange_thermal_energy
  PUBLIC :: eval_vaporization_sources

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
  !> \brief Evaluate the lava radiative and convective cooling source.
  !> \details Uses the net Stefan-type radiation and linear convection law
  !>          \f$S_T=-k_r(T^4-T_a^4)-k_c(T-T_a)\f$.
  !>          The complex temperature argument preserves complex-step
  !>          derivatives used by the implicit solver.
  !> \param[in] temperature Current lava temperature, real or complex [K].
  !> \param[in] radiative_coefficient Radiative coefficient k_r [W m^-2 K^-4].
  !> \param[in] convective_coefficient Convective coefficient k_c [W m^-2 K^-1].
  !> \param[in] ambient_temperature Ambient temperature T_a [K].
  !> \return Net cooling energy source per unit area [W m^-2].
  !******************************************************************************
  PURE FUNCTION eval_lava_cooling_source(temperature, radiative_coefficient,   &
       convective_coefficient, ambient_temperature) RESULT(source)

    COMPLEX(wp), INTENT(IN) :: temperature
    REAL(wp), INTENT(IN) :: radiative_coefficient
    REAL(wp), INTENT(IN) :: convective_coefficient
    REAL(wp), INTENT(IN) :: ambient_temperature
    COMPLEX(wp) :: source

    source = -radiative_coefficient * (temperature**4 - ambient_temperature**4) &
         - convective_coefficient * (temperature - ambient_temperature)

  END FUNCTION eval_lava_cooling_source

  !******************************************************************************
  !> \brief Advance the lava cooling ODE exactly over a time interval.
  !> \details Separates the autonomous temperature equation and evaluates its
  !>          analytic time-temperature primitive. A bisection inversion finds
  !>          the temperature after the requested DIRK integration step. The
  !>          areal thermal capacity and cooling coefficients are held fixed:
  !>          \f$C_A\frac{dT}{dt}=-k_r(T^4-T_a^4)-k_c(T-T_a)\f$.
  !>          The result is constrained to the closed interval between the
  !>          initial and ambient temperatures.
  !> \param[in] temperature Initial lava temperature T(0) [K].
  !> \param[in] integration_step Effective DIRK integration interval, normally dt*a_ii [s].
  !> \param[in] areal_heat_capacity Areal heat capacity C_A [J m^-2 K^-1].
  !> \param[in] radiative_coefficient Radiative cooling coefficient k_r [W m^-2 K^-4].
  !> \param[in] convective_coefficient Convective cooling coefficient k_c [W m^-2 K^-1].
  !> \param[in] ambient_temperature Ambient temperature T_a [K].
  !> \param[out] temperature_after Exact-flow temperature after integration_step [K].
  !> \note Requires finite inputs, positive initial and ambient temperatures, positive heat
  !>       capacity, nonnegative integration step, and nonnegative cooling coefficients.
  !******************************************************************************
  SUBROUTINE advance_lava_cooling_temperature(temperature, integration_step,    &
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

  END SUBROUTINE advance_lava_cooling_temperature

  !******************************************************************************
  !> \brief Solve the scalar DIRK stage equation for lava cooling.
  !> \details Solves
  !>          \f$T=T_k+\tau S_T(T)/C_A\f$, where
  !>          \f$\tau=dt\,a_{ii}\f$ and
  !>          \f$S_T=-k_r(T^4-T_a^4)-k_c(T-T_a)\f$. The equation reduces to
  !>          a monotone quartic of the form \f$aT^4+bT-c=0\f$ and is evaluated
  !>          with Ferrari/Cardano radicals, not by integrating the continuous
  !>          cooling ODE or iterating Newton.
  !> \param[in] temperature_known Stage-known temperature after prior RK terms [K].
  !> \param[in] integration_step Effective DIRK coefficient \f$dt\,a_{ii}\f$ [s].
  !> \param[in] areal_heat_capacity Areal heat capacity \f$C_A\f$ [J m^-2 K^-1].
  !> \param[in] radiative_coefficient Radiative coefficient \f$k_r\f$ [W m^-2 K^-4].
  !> \param[in] convective_coefficient Convective coefficient \f$k_c\f$ [W m^-2 K^-1].
  !> \param[in] ambient_temperature Ambient temperature \f$T_a\f$ [K].
  !> \param[out] temperature_stage DIRK stage temperature satisfying the stage equation [K].
  !> \note Requires finite nonnegative temperature_known and integration_step, positive
  !>       heat capacity and ambient temperature, and nonnegative cooling coefficients.
  !******************************************************************************
  SUBROUTINE solve_lava_cooling_stage_temperature(temperature_known,          &
       integration_step, areal_heat_capacity, radiative_coefficient,          &
       convective_coefficient, ambient_temperature, temperature_stage)

    REAL(wp), INTENT(IN) :: temperature_known, integration_step
    REAL(wp), INTENT(IN) :: areal_heat_capacity
    REAL(wp), INTENT(IN) :: radiative_coefficient, convective_coefficient
    REAL(wp), INTENT(IN) :: ambient_temperature
    REAL(wp), INTENT(OUT) :: temperature_stage

    REAL(wp) :: scale, ratio, scaled_radiation, a, b, c, w
    REAL(wp) :: log_p, p, u, x

    IF (.NOT. ieee_is_finite(temperature_known) .OR.                       &
         .NOT. ieee_is_finite(integration_step) .OR.                       &
         .NOT. ieee_is_finite(areal_heat_capacity) .OR.                   &
         .NOT. ieee_is_finite(radiative_coefficient) .OR.                 &
         .NOT. ieee_is_finite(convective_coefficient) .OR.                &
         .NOT. ieee_is_finite(ambient_temperature)) THEN
       CALL fatal_error('Nonfinite input to lava cooling DIRK stage solver')
    END IF
    IF (temperature_known .LT. 0.0_wp .OR. integration_step .LT. 0.0_wp .OR. &
         areal_heat_capacity .LE. 0.0_wp .OR. ambient_temperature .LE. 0.0_wp) THEN
       CALL fatal_error('Invalid state or integration data in lava cooling DIRK stage solver')
    END IF
    IF (radiative_coefficient .LT. 0.0_wp .OR.                               &
         convective_coefficient .LT. 0.0_wp) THEN
       CALL fatal_error('Negative coefficient in lava cooling DIRK stage solver')
    END IF

    temperature_stage = temperature_known
    IF (integration_step .EQ. 0.0_wp .OR.                                   &
         (radiative_coefficient .EQ. 0.0_wp .AND.                           &
         convective_coefficient .EQ. 0.0_wp)) RETURN

    scale = MAX(temperature_known, ambient_temperature)
    ratio = integration_step / areal_heat_capacity
    scaled_radiation = radiative_coefficient * scale**3
    a = ratio * scaled_radiation
    b = 1.0_wp + ratio * convective_coefficient
    c = temperature_known / scale + ratio *                                  &
         (scaled_radiation * (ambient_temperature / scale)**4 +              &
         convective_coefficient * ambient_temperature / scale)
    IF (.NOT. ieee_is_finite(a) .OR. .NOT. ieee_is_finite(b) .OR.             &
         .NOT. ieee_is_finite(c) .OR. b .LE. 0.0_wp .OR. c .LE. 0.0_wp) THEN
       CALL fatal_error('Invalid quartic coefficients in lava cooling DIRK stage solver')
    END IF

    IF (a .EQ. 0.0_wp) THEN
       x = c / b
    ELSE
       w = c / b
       log_p = LOG(b) - 0.25_wp * LOG(a) - 0.75_wp * LOG(c)
       IF (log_p .GE. 0.0_wp) THEN
          u = EXP(-4.0_wp * log_p)
          x = w * solve_u_y4_plus_y_minus_one(u)
       ELSE
          p = EXP(log_p)
          x = EXP(0.25_wp * (LOG(c) - LOG(a))) *                           &
               solve_x4_plus_p_x_minus_one(p)
       END IF
    END IF

    IF (.NOT. ieee_is_finite(x) .OR. x .LT. 0.0_wp .OR. x .GT. 1.0_wp +      &
         128.0_wp * EPSILON(1.0_wp)) THEN
       CALL fatal_error('Quartic root outside the physical interval in lava cooling DIRK stage')
    END IF
    temperature_stage = scale * MIN(1.0_wp, MAX(0.0_wp, x))

  END SUBROUTINE solve_lava_cooling_stage_temperature

  !******************************************************************************
  !> \brief Solve \f$u y^4+y-1=0\f$ for its unique positive root.
  !******************************************************************************
  PURE FUNCTION solve_u_y4_plus_y_minus_one(u) RESULT(y)

    REAL(wp), INTENT(IN) :: u
    REAL(wp) :: y
    REAL(wp) :: v, delta, s, c

    IF (u .EQ. 0.0_wp) THEN
       y = 1.0_wp
       RETURN
    END IF

    v = u**(1.0_wp/3.0_wp)
    delta = 1.0_wp/256.0_wp + v**3/27.0_wp
    s = SIGN(ABS(1.0_wp/16.0_wp + SQRT(delta))**(1.0_wp/3.0_wp),       &
         1.0_wp/16.0_wp + SQRT(delta)) +                               &
         SIGN(ABS(1.0_wp/16.0_wp - SQRT(delta))**(1.0_wp/3.0_wp),       &
         1.0_wp/16.0_wp - SQRT(delta))
    c = 1.0_wp / (2.0_wp * SQRT(2.0_wp*s))
    y = 2.0_wp / (s + c) / (SQRT(2.0_wp*s) +                            &
         SQRT(2.0_wp*s + 4.0_wp*v/(s+c)))

  END FUNCTION solve_u_y4_plus_y_minus_one

  !******************************************************************************
  !> \brief Solve \f$y^4+p y-1=0\f$ for its unique positive root when 0<=p<=1.
  !******************************************************************************
  PURE FUNCTION solve_x4_plus_p_x_minus_one(p) RESULT(y)

    REAL(wp), INTENT(IN) :: p
    REAL(wp) :: y
    REAL(wp) :: z, m, n, discriminant

    IF (p .EQ. 0.0_wp) THEN
       y = 1.0_wp
       RETURN
    END IF

    z = 2.0_wp / SQRT(3.0_wp) * SINH(ASINH(3.0_wp*p**2*SQRT(3.0_wp) / 16.0_wp) / &
         3.0_wp)
    m = SQRT(2.0_wp*z)
    n = -p / (2.0_wp*m)
    discriminant = MAX(0.0_wp, m**2 - 4.0_wp*(z+n))
    y = 0.5_wp * (SQRT(discriminant) - m)

  END FUNCTION solve_x4_plus_p_x_minus_one

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

  !******************************************************************************
  !> \brief Compute the thermal-energy source associated with mass exchange.
  !> \details Accounts for energy removed with deposited solids and lost carrier
  !>          material, energy supplied by eroded solids and carrier material,
  !>          and ambient-air entrainment:
  !>          \f[
  !>          S_T=-T_f\left(\sum_i\rho_i c_iD_i+\rho_c c_cL_c\right)
  !>              +T_s\left(\sum_i\rho_i c_iE_i+\rho_{sc}c_cE_c\right)
  !>              +T_a c_a\rho_aE_a.
  !>          \f]
  !> \param[in] flow_temperature Temperature of the flowing mixture T_f [K].
  !> \param[in] solid_density Per-class solid densities rho_i [kg m^-3].
  !> \param[in] solid_heat_capacity Per-class solid specific heat capacities c_i [J kg^-1 K^-1].
  !> \param[in] deposition_rate Per-class deposition thickness rates D_i [m s^-1].
  !> \param[in] erosion_rate Per-class solid erosion thickness rates E_i [m s^-1].
  !> \param[in] carrier_density Flow carrier density rho_c [kg m^-3].
  !> \param[in] carrier_heat_capacity Carrier specific heat capacity c_c [J kg^-1 K^-1].
  !> \param[in] carrier_loss_rate Carrier-layer loss rate L_c [m s^-1].
  !> \param[in] substrate_temperature Substrate temperature T_s for incoming eroded material [K].
  !> \param[in] substrate_carrier_density Carrier density of the substrate material rho_sc
  !>                                      [kg m^-3].
  !> \param[in] carrier_erosion_rate Carrier-layer erosion/entrainment rate E_c [m s^-1].
  !> \param[in] ambient_temperature Ambient-air temperature T_a [K].
  !> \param[in] air_heat_capacity Ambient-air specific heat capacity c_a [J kg^-1 K^-1].
  !> \param[in] ambient_air_density Ambient-air density rho_a [kg m^-3].
  !> \param[in] air_entrainment_rate Ambient-air entrainment rate E_a [m s^-1].
  !> \return Net mass-exchange thermal-energy source per unit area [W m^-2].
  !******************************************************************************
  PURE FUNCTION eval_mass_exchange_thermal_energy(flow_temperature,            &
       solid_density, solid_heat_capacity, deposition_rate, erosion_rate,      &
       carrier_density, carrier_heat_capacity, carrier_loss_rate,              &
       substrate_temperature, substrate_carrier_density,                       &
       carrier_erosion_rate, ambient_temperature, air_heat_capacity,            &
       ambient_air_density, air_entrainment_rate) RESULT(source)

    REAL(wp), INTENT(IN) :: flow_temperature
    REAL(wp), INTENT(IN) :: solid_density(:)
    REAL(wp), INTENT(IN) :: solid_heat_capacity(:)
    REAL(wp), INTENT(IN) :: deposition_rate(:)
    REAL(wp), INTENT(IN) :: erosion_rate(:)
    REAL(wp), INTENT(IN) :: carrier_density
    REAL(wp), INTENT(IN) :: carrier_heat_capacity
    REAL(wp), INTENT(IN) :: carrier_loss_rate
    REAL(wp), INTENT(IN) :: substrate_temperature
    REAL(wp), INTENT(IN) :: substrate_carrier_density
    REAL(wp), INTENT(IN) :: carrier_erosion_rate
    REAL(wp), INTENT(IN) :: ambient_temperature
    REAL(wp), INTENT(IN) :: air_heat_capacity
    REAL(wp), INTENT(IN) :: ambient_air_density
    REAL(wp), INTENT(IN) :: air_entrainment_rate
    REAL(wp) :: source

    source = -flow_temperature * ( SUM( solid_density * solid_heat_capacity *  &
         deposition_rate ) + carrier_density * carrier_heat_capacity *         &
         carrier_loss_rate ) + substrate_temperature * ( SUM( solid_density *  &
         solid_heat_capacity * erosion_rate ) + substrate_carrier_density *    &
         carrier_heat_capacity * carrier_erosion_rate ) + ambient_temperature &
         * air_heat_capacity * ambient_air_density * air_entrainment_rate

  END FUNCTION eval_mass_exchange_thermal_energy

  !******************************************************************************
  !> \brief Compute vapor mass and sensible-energy sources from hot-solid deposition.
  !> \details Deposition heat above the boiling point supplies liquid warming
  !>          and latent vaporization energy. The generated steam rate is
  !>          multiplied by the configured steam fraction; its energy source is
  !>          evaluated at the boiling temperature:
  !>          \f[
  !>          \dot{m}_v=\frac{f_v\max(0,T_f-T_b)\sum_i\rho_i c_iD_i}
  !>                         {c_l(T_b-T_l)+L_v},\qquad
  !>          S_v=\dot{m}_v c_gT_b.
  !>          \f]
  !> \param[in] flow_temperature Flow temperature T_f [K].
  !> \param[in] solid_density Per-class solid densities rho_i [kg m^-3].
  !> \param[in] solid_heat_capacity Per-class solid specific heat capacities c_i [J kg^-1 K^-1].
  !> \param[in] deposition_rate Per-class deposition thickness rates D_i [m s^-1].
  !> \param[in] steam_fraction Fraction f_v of the available deposition heat assigned to vapor.
  !> \param[in] gas_heat_capacity Steam specific heat capacity c_g [J kg^-1 K^-1].
  !> \param[out] vapor_mass_rate Generated vapor mass rate per unit area [kg m^-2 s^-1].
  !> \param[out] vapor_energy_source Vapor sensible-energy source per unit area [W m^-2].
  !******************************************************************************
  PURE SUBROUTINE eval_vaporization_sources(flow_temperature, solid_density,    &
       solid_heat_capacity, deposition_rate, steam_fraction, gas_heat_capacity,&
       vapor_mass_rate, vapor_energy_source)

    REAL(wp), INTENT(IN) :: flow_temperature
    REAL(wp), INTENT(IN) :: solid_density(:)
    REAL(wp), INTENT(IN) :: solid_heat_capacity(:)
    REAL(wp), INTENT(IN) :: deposition_rate(:)
    REAL(wp), INTENT(IN) :: steam_fraction
    REAL(wp), INTENT(IN) :: gas_heat_capacity
    REAL(wp), INTENT(OUT) :: vapor_mass_rate
    REAL(wp), INTENT(OUT) :: vapor_energy_source

    REAL(wp), PARAMETER :: liquid_temperature = 290.0_wp
    REAL(wp), PARAMETER :: boiling_temperature = 373.15_wp
    REAL(wp), PARAMETER :: latent_heat = 2264705.0_wp
    REAL(wp), PARAMETER :: liquid_heat_capacity = 4184.0_wp

    vapor_mass_rate = steam_fraction * MAX(0.0_wp, flow_temperature -         &
         boiling_temperature) * SUM(solid_density * solid_heat_capacity *      &
         deposition_rate) / (liquid_heat_capacity * (boiling_temperature -    &
         liquid_temperature) + latent_heat)
    vapor_energy_source = vapor_mass_rate * gas_heat_capacity *                &
         boiling_temperature

  END SUBROUTINE eval_vaporization_sources

END MODULE thermal_sources_2d
