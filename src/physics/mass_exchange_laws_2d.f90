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

END MODULE mass_exchange_laws_2d
