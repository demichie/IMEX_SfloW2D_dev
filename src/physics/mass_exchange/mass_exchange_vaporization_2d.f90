!> Vaporization source contribution to conserved mass, energy and gas.
MODULE mass_exchange_vaporization_2d
  USE parameters_2d, ONLY : wp, n_solid, n_add_gas, gas_flag, &
       liquid_vaporization_flag, mass_exchange_rates_type, reset_mass_exchange_rates
  USE diagnostics_2d, ONLY : fatal_error
  USE constitutive_parameters_2d, ONLY : rho_s, sp_heat_s, gamma_steam, sp_heat_g
  USE thermal_laws_2d, ONLY : eval_vaporization_mass_rate, eval_vaporization_energy_rate

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_mass_exchange_vaporization_explicit
  PUBLIC :: eval_mass_exchange_vaporization_semi_implicit
  PUBLIC :: eval_mass_exchange_vaporization_implicit
  PUBLIC :: solve_mass_exchange_vaporization_analytic
  ! Four-entry physical treatment contract; only EXPLICIT is supported.
CONTAINS




  SUBROUTINE eval_mass_exchange_vaporization_explicit(B_zone,T,deposition,r)
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(IN) :: T,deposition(n_solid)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    IF (.NOT.(gas_flag .AND. liquid_vaporization_flag .AND. B_zone/=0)) RETURN
    IF (n_add_gas<1) CALL fatal_error('vaporization requires additional gas species')
    r%mixture_mass_rate=eval_vaporization_mass_rate(T,rho_s,sp_heat_s, &
         deposition,gamma_steam)
    r%add_gas_mass_rate(1)=r%mixture_mass_rate
    r%energy_rate=eval_vaporization_energy_rate(r%mixture_mass_rate,sp_heat_g(1))
  END SUBROUTINE


  SUBROUTINE eval_mass_exchange_vaporization_semi_implicit(B_zone,T,deposition,r)
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(IN) :: T,deposition(n_solid)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('vaporization: unsupported physical treatment')
  END SUBROUTINE eval_mass_exchange_vaporization_semi_implicit

  SUBROUTINE eval_mass_exchange_vaporization_implicit(B_zone,T,deposition,r)
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(IN) :: T,deposition(n_solid)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('vaporization: unsupported physical treatment')
  END SUBROUTINE eval_mass_exchange_vaporization_implicit

  SUBROUTINE solve_mass_exchange_vaporization_analytic(B_zone,T,deposition,r)
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(IN) :: T,deposition(n_solid)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('vaporization: unsupported physical treatment')
  END SUBROUTINE solve_mass_exchange_vaporization_analytic


END MODULE mass_exchange_vaporization_2d
