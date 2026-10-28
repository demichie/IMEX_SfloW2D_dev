!> \brief Prescribed immobilization of liquid lava into a bed deposit.
!> \details Only EXPLICIT is supported; all other treatments fail explicitly.
MODULE mass_exchange_lava_deposition_2d
  USE parameters_2d, ONLY : wp
  USE diagnostics_2d, ONLY : fatal_error
  USE source_input_2d, ONLY : eval_lava_deposition_rate
  USE rheology_lava_2d, ONLY : lava_kinematic_viscosity
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_mass_exchange_lava_deposition_explicit
  PUBLIC :: eval_mass_exchange_lava_deposition_semi_implicit
  PUBLIC :: eval_mass_exchange_lava_deposition_implicit
  PUBLIC :: solve_mass_exchange_lava_deposition_analytic
CONTAINS
  !> \brief Evaluate an available-depth-limited lava immobilization rate.
  !> \param[in] depth Mobile lava thickness [m].
  !> \param[in] temperature Mobile lava temperature [K].
  !> \param[in] dt Physical step duration [s].
  !> \param[out] rate Immobile lava thickness gained per unit time [m/s].
  SUBROUTINE eval_mass_exchange_lava_deposition_explicit(depth,temperature,dt,rate)
    REAL(wp), INTENT(IN) :: depth,temperature,dt
    REAL(wp), INTENT(OUT) :: rate
    rate=0.0_wp
    IF (depth<=0.0_wp) RETURN
    IF (dt<=0.0_wp) CALL fatal_error('lava deposition requires positive dt')
    rate=MIN(depth/dt,eval_lava_deposition_rate(temperature, &
         lava_kinematic_viscosity(temperature)))
  END SUBROUTINE eval_mass_exchange_lava_deposition_explicit

  !> \brief Unsupported semi-implicit lava deposition.
  SUBROUTINE eval_mass_exchange_lava_deposition_semi_implicit()
    CALL fatal_error('unsupported semi_implicit lava deposition')
  END SUBROUTINE eval_mass_exchange_lava_deposition_semi_implicit

  !> \brief Unsupported implicit lava deposition.
  SUBROUTINE eval_mass_exchange_lava_deposition_implicit()
    CALL fatal_error('unsupported implicit lava deposition')
  END SUBROUTINE eval_mass_exchange_lava_deposition_implicit

  !> \brief Unsupported analytic lava deposition.
  SUBROUTINE solve_mass_exchange_lava_deposition_analytic()
    CALL fatal_error('unsupported analytic lava deposition')
  END SUBROUTINE solve_mass_exchange_lava_deposition_analytic
END MODULE mass_exchange_lava_deposition_2d
