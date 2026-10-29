!> Vaporization source contribution to conserved mass, energy and gas.
MODULE mass_exchange_vaporization_2d
  USE parameters_2d, ONLY : wp, n_eqns, n_solid, gas_flag, liquid_vaporization_flag
  USE diagnostics_2d, ONLY : fatal_error
  USE constitutive_parameters_2d, ONLY : rho_s, sp_heat_s, gamma_steam, sp_heat_g
  USE thermal_laws_2d, ONLY : eval_vaporization_mass_rate, eval_vaporization_energy_rate
  USE mass_exchange_rates_2d, ONLY: mass_exchange_rates_type, reset_mass_exchange_rates

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_mass_exchange_vaporization_process_explicit
  PUBLIC :: eval_mass_exchange_vaporization_explicit
  PUBLIC :: eval_mass_exchange_vaporization_rates_explicit
  PUBLIC :: eval_mass_exchange_vaporization_semi_implicit
  PUBLIC :: eval_mass_exchange_vaporization_implicit
  PUBLIC :: solve_mass_exchange_vaporization_analytic
  ! Four-entry physical treatment contract; only EXPLICIT is supported.
CONTAINS

  !> Physical outputs: identical mass added to mixture and added gas #1,
  !> plus the energy of the vapor. No conservative equation indices here.
  SUBROUTINE eval_mass_exchange_vaporization_rates_explicit(B_zone, r_T, &
       deposition_term, mixture_mass_rate, add_gas_mass_rate, energy_rate)
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(IN) :: r_T, deposition_term(n_solid)
    REAL(wp), INTENT(OUT) :: mixture_mass_rate, add_gas_mass_rate, energy_rate
    mixture_mass_rate = 0.0_wp
    add_gas_mass_rate = 0.0_wp
    energy_rate = 0.0_wp
    IF (gas_flag .AND. liquid_vaporization_flag .AND. B_zone /= 0) THEN
       mixture_mass_rate = eval_vaporization_mass_rate(r_T, rho_s, sp_heat_s, &
            deposition_term, gamma_steam)
       add_gas_mass_rate = mixture_mass_rate
       energy_rate = eval_vaporization_energy_rate(mixture_mass_rate, sp_heat_g(1))
    END IF
  END SUBROUTINE eval_mass_exchange_vaporization_rates_explicit

  !> Compatibility entry point until all process interfaces use physical rates.
  SUBROUTINE eval_mass_exchange_vaporization_explicit(B_zone, r_T, deposition_term, eqns_term)
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(IN) :: r_T, deposition_term(n_solid)
    REAL(wp), INTENT(INOUT) :: eqns_term(n_eqns)
    REAL(wp) :: mixture_mass_rate, add_gas_mass_rate, energy_rate
    CALL eval_mass_exchange_vaporization_rates_explicit(B_zone, r_T, deposition_term, &
         mixture_mass_rate, add_gas_mass_rate, energy_rate)
    eqns_term(1) = eqns_term(1) + mixture_mass_rate
    eqns_term(4) = eqns_term(4) + energy_rate
    IF (gas_flag .AND. liquid_vaporization_flag .AND. B_zone /= 0) &
         eqns_term(5+n_solid) = eqns_term(5+n_solid) + add_gas_mass_rate
  END SUBROUTINE eval_mass_exchange_vaporization_explicit

  SUBROUTINE eval_mass_exchange_vaporization_semi_implicit(B_zone, r_T, deposition_term, eqns_term)
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(IN) :: r_T, deposition_term(n_solid)
    REAL(wp), INTENT(INOUT) :: eqns_term(n_eqns)
    CALL fatal_error('vaporization: semi_implicit treatment is unsupported')
  END SUBROUTINE eval_mass_exchange_vaporization_semi_implicit

  SUBROUTINE eval_mass_exchange_vaporization_implicit(B_zone, r_T, deposition_term, eqns_term)
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(IN) :: r_T, deposition_term(n_solid)
    REAL(wp), INTENT(INOUT) :: eqns_term(n_eqns)
    CALL fatal_error('vaporization: implicit treatment is unsupported')
  END SUBROUTINE eval_mass_exchange_vaporization_implicit

  SUBROUTINE solve_mass_exchange_vaporization_analytic(B_zone, r_T, deposition_term, eqns_term)
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(IN) :: r_T, deposition_term(n_solid)
    REAL(wp), INTENT(INOUT) :: eqns_term(n_eqns)
    CALL fatal_error('vaporization: analytic treatment is unsupported')
  END SUBROUTINE solve_mass_exchange_vaporization_analytic




  SUBROUTINE eval_mass_exchange_vaporization_process_explicit(B_zone,T,deposition,r)
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(IN) :: T,deposition(n_solid)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    REAL(wp) :: gas_mass
    CALL reset_mass_exchange_rates(r)
    CALL eval_mass_exchange_vaporization_rates_explicit(B_zone,T,deposition, &
         r%mixture_mass_rate,gas_mass,r%energy_rate)
    IF (SIZE(r%add_gas_mass_rate)>0) r%add_gas_mass_rate(1)=gas_mass
  END SUBROUTINE
END MODULE mass_exchange_vaporization_2d
