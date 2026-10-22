!********************************************************************************
!> \brief Constitutive laws used by mass-exchange source assembly.
!> \details Provides cell-local deposition rates independently of source
!>          assembly, state updates and inventory bookkeeping.
!********************************************************************************

MODULE mass_exchange_laws_2d

  USE parameters_2d, ONLY : wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_lahar_deposition_rate
  PUBLIC :: eval_hindered_settling_deposition_rate

CONTAINS

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

END MODULE mass_exchange_laws_2d
