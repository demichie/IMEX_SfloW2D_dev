!> \brief Prescribed immobilization of liquid lava into a bed deposit.
!> \details Only EXPLICIT is supported; all other treatments fail explicitly.
MODULE mass_exchange_lava_deposition_2d
  USE parameters_2d, ONLY : wp
  USE diagnostics_2d, ONLY : fatal_error
  USE rheology_lava_2d, ONLY : lava_kinematic_viscosity
  USE mass_exchange_rates_2d, ONLY: mass_exchange_rates_type, reset_mass_exchange_rates
  USE constitutive_parameters_2d, ONLY: rho_l, sp_heat_l
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_mass_exchange_lava_deposition_explicit
  LOGICAL, PUBLIC :: lava_deposition_flag=.FALSE.
  CHARACTER(LEN=16), PUBLIC :: lava_deposition_model='CONSTANT'
  REAL(wp), PUBLIC :: lava_deposition_rate=0.0_wp
  REAL(wp), PUBLIC :: lava_deposition_visc_ref=1.0E-2_wp
  REAL(wp), PUBLIC :: lava_deposition_exponent=1.0_wp
  REAL(wp), PUBLIC :: lava_deposition_temp_ref=900.0_wp
  REAL(wp), PUBLIC :: lava_deposition_temp_width=100.0_wp
  PUBLIC :: eval_mass_exchange_lava_deposition_explicit
  PUBLIC :: eval_mass_exchange_lava_deposition_semi_implicit
  PUBLIC :: eval_mass_exchange_lava_deposition_implicit
  PUBLIC :: solve_mass_exchange_lava_deposition_analytic
CONTAINS
  SUBROUTINE eval_mass_exchange_lava_deposition_explicit(depth,T,dt,r)
    REAL(wp), INTENT(IN) :: depth,T,dt
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL eval_lava_deposition_volume_kernel(depth,T,dt, &
         r%continuous_phase_deposition_rate)
    r%mixture_mass_rate=-rho_l*r%continuous_phase_deposition_rate
    r%energy_rate=sp_heat_l*T*r%mixture_mass_rate
    r%dB_dt=r%continuous_phase_deposition_rate
  END SUBROUTINE


  !> \brief Unsupported semi-implicit lava deposition.
  SUBROUTINE eval_mass_exchange_lava_deposition_semi_implicit(depth,T,dt,r)
    REAL(wp), INTENT(IN) :: depth,T,dt
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('lava_deposition: unsupported physical treatment')
  END SUBROUTINE eval_mass_exchange_lava_deposition_semi_implicit


  !> \brief Unsupported implicit lava deposition.
  SUBROUTINE eval_mass_exchange_lava_deposition_implicit(depth,T,dt,r)
    REAL(wp), INTENT(IN) :: depth,T,dt
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('lava_deposition: unsupported physical treatment')
  END SUBROUTINE eval_mass_exchange_lava_deposition_implicit


  !> \brief Unsupported analytic lava deposition.
  SUBROUTINE solve_mass_exchange_lava_deposition_analytic(depth,T,dt,r)
    REAL(wp), INTENT(IN) :: depth,T,dt
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('lava_deposition: unsupported physical treatment')
  END SUBROUTINE solve_mass_exchange_lava_deposition_analytic

  !> \brief Requested immobile lava deposition thickness rate.
  !> \details Caller enforces available thickness over dt.
  !> \param[in] temperature Current mobile temperature [K].
  !> \param[in] viscosity Lava kinematic viscosity [m2/s].
  !> \return Prescribed deposition rate [m/s].
  PURE REAL(wp) FUNCTION eval_lava_deposition_rate(temperature,viscosity) RESULT(rate)
    REAL(wp), INTENT(IN) :: temperature,viscosity
    REAL(wp) :: factor
    factor=1.0_wp
    SELECT CASE(TRIM(lava_deposition_model))
    CASE('TEMPERATURE')
       factor=MIN(1.0_wp,MAX(0.0_wp,(lava_deposition_temp_ref-temperature)/ &
            lava_deposition_temp_width))
    CASE('VISCOSITY')
       factor=MIN(1.0_wp,(MAX(0.0_wp,viscosity)/lava_deposition_visc_ref)** &
            lava_deposition_exponent)
    END SELECT
    rate=MAX(0.0_wp,lava_deposition_rate)*factor
  END FUNCTION eval_lava_deposition_rate


  !> \brief Evaluate an available-depth-limited lava immobilization rate.
  !> \param[in] depth Mobile lava thickness [m].
  !> \param[in] temperature Mobile lava temperature [K].
  !> \param[in] dt Physical step duration [s].
  !> \param[out] rate Immobile lava thickness gained per unit time [m/s].
  SUBROUTINE eval_lava_deposition_volume_kernel(depth,temperature,dt,rate)
    REAL(wp), INTENT(IN) :: depth,temperature,dt
    REAL(wp), INTENT(OUT) :: rate
    rate=0.0_wp
    IF (depth<=0.0_wp) RETURN
    IF (dt<=0.0_wp) CALL fatal_error('lava deposition requires positive dt')
    rate=MIN(depth/dt,eval_lava_deposition_rate(temperature, &
         lava_kinematic_viscosity(temperature)))
  END SUBROUTINE eval_lava_deposition_volume_kernel
END MODULE mass_exchange_lava_deposition_2d
