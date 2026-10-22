!********************************************************************************
!> \brief Constitutive laws used by mass-exchange source assembly.
!> \details Provides cell-local erosion, deposition, carrier-loss and
!>          entrainment rates independently of source assembly and bookkeeping.
!********************************************************************************

MODULE mass_exchange_laws_2d

  USE parameters_2d, ONLY : wp, rheology_id_lahar
  USE rheology_lahar_2d, ONLY : eval_lahar_fluid_viscosity

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_lahar_deposition_rate
  PUBLIC :: eval_hindered_settling_deposition_rate
  PUBLIC :: eval_settling_fluid_properties
  PUBLIC :: eval_solid_deposition_rate
  PUBLIC :: eval_erosion_rates
  PUBLIC :: eval_continuous_phase_loss_rate
  PUBLIC :: eval_air_entrainment_rate

CONTAINS

  !******************************************************************************
  !> \brief Select the fluid properties used by the settling calculation.
  !> \details Lahar settling uses the lahar rheology's temperature- and
  !>          concentration-dependent viscosity and mixture density. Other
  !>          models retain the carrier-fluid properties supplied by the caller.
  !> \param[in] model_id Selected rheology model ID.
  !> \param[in] temperature Mixture temperature [K].
  !> \param[in] total_solid_fraction Sum of solid volume fractions.
  !> \param[in] mixture_density Mixture density [kg m^-3].
  !> \param[in] carrier_density_fallback Carrier density for non-lahar models.
  !> \param[in] kinematic_viscosity_fallback Kinematic viscosity for
  !>        non-lahar models [m^2 s^-1].
  !> \param[out] carrier_density Density used for settling [kg m^-3].
  !> \param[out] inverse_kinematic_viscosity Reciprocal kinematic viscosity
  !>        used for settling [s m^-2].
  !******************************************************************************
  SUBROUTINE eval_settling_fluid_properties(model_id, temperature,           &
       total_solid_fraction, mixture_density, carrier_density_fallback,     &
       kinematic_viscosity_fallback, carrier_density,                       &
       inverse_kinematic_viscosity)

    INTEGER, INTENT(IN) :: model_id
    REAL(wp), INTENT(IN) :: temperature, total_solid_fraction, mixture_density
    REAL(wp), INTENT(IN) :: carrier_density_fallback
    REAL(wp), INTENT(IN) :: kinematic_viscosity_fallback
    REAL(wp), INTENT(OUT) :: carrier_density, inverse_kinematic_viscosity
    REAL(wp) :: fluid_viscosity

    carrier_density = carrier_density_fallback

    IF (model_id .EQ. rheology_id_lahar) THEN
       fluid_viscosity = eval_lahar_fluid_viscosity(temperature,            &
            total_solid_fraction)
       carrier_density = mixture_density
       inverse_kinematic_viscosity = mixture_density / fluid_viscosity
    ELSE
       inverse_kinematic_viscosity = 1.0_wp / kinematic_viscosity_fallback
    END IF

  END SUBROUTINE eval_settling_fluid_properties

  !******************************************************************************
  !> \brief Dispatch the solid deposition law selected for settling.
  !> \param[in] model_id Selected rheology model ID.
  !> \param[in] solid_fraction Solid volume fraction for this class.
  !> \param[in] settling_velocity Settling velocity [m s^-1].
  !> \param[in] flow_depth Flow depth [m].
  !> \param[in] time_step Time-step size [s].
  !> \param[in] total_solid_fraction Sum of solid volume fractions.
  !> \param[in] maximum_solid_fraction Hindered-settling packing threshold.
  !> \param[in] hindered_settling_exponent Michaels-Bolger exponent.
  !> \return Deposition rate [m s^-1].
  !******************************************************************************
  PURE FUNCTION eval_solid_deposition_rate(model_id, solid_fraction,         &
       settling_velocity, flow_depth, time_step, total_solid_fraction,       &
       maximum_solid_fraction, hindered_settling_exponent) RESULT(rate)

    INTEGER, INTENT(IN) :: model_id
    REAL(wp), INTENT(IN) :: solid_fraction, settling_velocity
    REAL(wp), INTENT(IN) :: flow_depth, time_step, total_solid_fraction
    REAL(wp), INTENT(IN) :: maximum_solid_fraction, hindered_settling_exponent
    REAL(wp) :: rate

    IF (model_id .EQ. rheology_id_lahar) THEN
       rate = eval_lahar_deposition_rate(solid_fraction, settling_velocity,   &
            flow_depth, time_step)
    ELSE
       rate = eval_hindered_settling_deposition_rate(solid_fraction,         &
            settling_velocity, flow_depth, time_step, total_solid_fraction,  &
            maximum_solid_fraction, hindered_settling_exponent)
    END IF

  END FUNCTION eval_solid_deposition_rate

  !******************************************************************************
  !> \brief Evaluate the deposition-limited lahar settling rate.
  !> \details Caps the settling flux by the amount of solid available in the
  !> flow during the current time step:
  !> \f$D=\min(Cw_s,hC/\Delta t)\f$.
  !> \param[in] solid_fraction Solid volume fraction.
  !> \param[in] settling_velocity Settling velocity.
  !> \param[in] flow_depth Flow depth.
  !> \param[in] time_step Time-step size.
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
  !> \brief Evaluate deposition from hindered settling with an availability cap.
  !> \details Applies the Michaels-Bolger settling correction, then limits
  !> deposition to the solid mass available during the current time step.
  !> \param[in] solid_fraction Solid volume fraction for this class.
  !> \param[in] settling_velocity Settling velocity [m s^-1].
  !> \param[in] flow_depth Flow depth [m].
  !> \param[in] time_step Time-step size [s].
  !> \param[in] total_solid_fraction Sum of solid volume fractions.
  !> \param[in] maximum_solid_fraction Hindered-settling packing threshold.
  !> \param[in] hindered_settling_exponent Michaels-Bolger exponent.
  !> \return Deposition rate [m s^-1].
  !******************************************************************************
  PURE FUNCTION eval_hindered_settling_deposition_rate(solid_fraction,       &
       settling_velocity, flow_depth, time_step, total_solid_fraction,       &
       maximum_solid_fraction, hindered_settling_exponent) RESULT(rate)

    REAL(wp), INTENT(IN) :: solid_fraction, settling_velocity
    REAL(wp), INTENT(IN) :: flow_depth, time_step, total_solid_fraction
    REAL(wp), INTENT(IN) :: maximum_solid_fraction, hindered_settling_exponent
    REAL(wp) :: rate

    rate = solid_fraction * settling_velocity
    rate = rate * (1.0_wp - MIN(1.0_wp, total_solid_fraction /               &
         maximum_solid_fraction))**hindered_settling_exponent
    rate = MIN(rate, flow_depth * solid_fraction / time_step)

  END FUNCTION eval_hindered_settling_deposition_rate

  !******************************************************************************
  !> \brief Evaluate solid erosion rates limited by substrate inventory.
  !> \details Applies the empirical erosion law, partitions erosion across
  !> erodible components, limits each rate to available inventory and returns
  !> the associated continuous-phase erosion rate.
  !> \param[in] flow_depth Flow depth [m].
  !> \param[in] speed Three-dimensional flow-speed magnitude [m s^-1].
  !> \param[in] total_solid_fraction Sum of solid volume fractions.
  !> \param[in] erosion_coefficient Empirical erosion coefficient.
  !> \param[in] erodible_porosity Substrate porosity.
  !> \param[in] erodible_fraction Per-class fraction of erodible material.
  !> \param[in] erodible_inventory Available per-class substrate thickness [m].
  !> \param[in] time_step Time-step size [s].
  !> \param[out] solid_erosion_rate Per-class solid erosion rates [m s^-1].
  !> \param[out] continuous_phase_erosion_rate Continuous-phase erosion rate [m s^-1].
  !******************************************************************************
  PURE SUBROUTINE eval_erosion_rates(flow_depth, speed, total_solid_fraction, &
       erosion_coefficient, erodible_porosity, erodible_fraction,            &
       erodible_inventory, time_step, solid_erosion_rate,                    &
       continuous_phase_erosion_rate)

    REAL(wp), INTENT(IN) :: flow_depth, speed, total_solid_fraction
    REAL(wp), INTENT(IN) :: erosion_coefficient, erodible_porosity
    REAL(wp), INTENT(IN) :: erodible_fraction(:), erodible_inventory(:)
    REAL(wp), INTENT(IN) :: time_step
    REAL(wp), INTENT(OUT) :: solid_erosion_rate(:)
    REAL(wp), INTENT(OUT) :: continuous_phase_erosion_rate

    REAL(wp) :: total_erosion, total_solid_erosion

    IF (erosion_coefficient .GT. 0.0_wp) THEN
       total_erosion = erosion_coefficient * speed * flow_depth *            &
            (1.0_wp - total_solid_fraction)
       total_solid_erosion = total_erosion * (1.0_wp - erodible_porosity)
       solid_erosion_rate = erodible_fraction * total_solid_erosion
    ELSE
       solid_erosion_rate = 0.0_wp
    END IF

    solid_erosion_rate = MAX(0.0_wp, MIN(solid_erosion_rate,                 &
         erodible_inventory / time_step))
    total_solid_erosion = SUM(solid_erosion_rate)
    total_erosion = total_solid_erosion / (1.0_wp - erodible_porosity)
    continuous_phase_erosion_rate = total_erosion * erodible_porosity

  END SUBROUTINE eval_erosion_rates

  !******************************************************************************
  !> \brief Limit the requested continuous-phase loss to available capacity.
  !> \details Applies the solid-packing cap and the carrier-volume reserve after
  !>          accounting for the associated solid deposition.
  !> \param[in] requested_loss_rate Requested carrier loss rate [m s^-1].
  !> \param[in] flow_depth Flow depth [m].
  !> \param[in] total_solid_fraction Sum of solid volume fractions.
  !> \param[in] maximum_solid_packing Maximum solid packing fraction.
  !> \param[in] porosity_coefficient Carrier volume associated with deposited solids.
  !> \param[in] deposition_rate Per-class solid deposition rates [m s^-1].
  !> \param[in] time_step Time-step size [s].
  !> \return Limited continuous-phase loss rate [m s^-1].
  !******************************************************************************
  PURE FUNCTION eval_continuous_phase_loss_rate(requested_loss_rate,         &
       flow_depth, total_solid_fraction, maximum_solid_packing,               &
       porosity_coefficient, deposition_rate, time_step) RESULT(loss_rate)

    REAL(wp), INTENT(IN) :: requested_loss_rate, flow_depth
    REAL(wp), INTENT(IN) :: total_solid_fraction, maximum_solid_packing
    REAL(wp), INTENT(IN) :: porosity_coefficient, deposition_rate(:)
    REAL(wp), INTENT(IN) :: time_step
    REAL(wp) :: loss_rate

    loss_rate = MIN(requested_loss_rate, flow_depth * MAX(0.0_wp,            &
         maximum_solid_packing - total_solid_fraction) / time_step)
    loss_rate = MIN(loss_rate, MAX(0.0_wp, (flow_depth *                     &
         (1.0_wp - total_solid_fraction) - porosity_coefficient *            &
         (flow_depth * total_solid_fraction - time_step *                    &
         SUM(deposition_rate))) / time_step))

  END FUNCTION eval_continuous_phase_loss_rate

  !******************************************************************************
  !> \brief Evaluate the empirical air-entrainment rate.
  !> \param[in] richardson_number Richardson number.
  !> \param[in] speed_squared Squared flow-speed magnitude [m2 s-2].
  !> \return Air-entrainment rate [m s-1].
  !******************************************************************************
  PURE FUNCTION eval_air_entrainment_rate(richardson_number, speed_squared)  &
       RESULT(entrainment_rate)

    REAL(wp), INTENT(IN) :: richardson_number, speed_squared
    REAL(wp) :: entrainment_coefficient, entrainment_rate

    entrainment_coefficient = 0.075_wp / SQRT(1.0_wp + 718.0_wp *            &
         richardson_number**2.4_wp)
    entrainment_rate = entrainment_coefficient * SQRT(speed_squared)

  END FUNCTION eval_air_entrainment_rate

END MODULE mass_exchange_laws_2d
