!> Vaporization source contribution to conserved mass, energy and gas.
MODULE mass_exchange_vaporization_2d
  USE parameters_2d, ONLY : wp, n_eqns, n_solid, gas_flag, liquid_vaporization_flag
  USE diagnostics_2d, ONLY : fatal_error
  USE constitutive_parameters_2d, ONLY : rho_s, sp_heat_s, gamma_steam, sp_heat_g
  USE thermal_laws_2d, ONLY : eval_vaporization_mass_rate, eval_vaporization_energy_rate
  IMPLICIT NONE
  PRIVATE
  ! Four-entry physical treatment contract; only EXPLICIT is supported.
  PUBLIC :: eval_mass_exchange_vaporization_explicit
  PUBLIC :: eval_mass_exchange_vaporization_semi_implicit
  PUBLIC :: eval_mass_exchange_vaporization_implicit
  PUBLIC :: solve_mass_exchange_vaporization_analytic
CONTAINS
  SUBROUTINE eval_mass_exchange_vaporization_explicit(B_zone, r_T, deposition_term, eqns_term)
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(IN) :: r_T, deposition_term(n_solid)
    REAL(wp), INTENT(INOUT) :: eqns_term(n_eqns)
    REAL(wp) :: vapor_mass_rate, vapor_energy_rate
    IF ( gas_flag .AND. liquid_vaporization_flag .AND. ( B_zone .NE. 0 ) ) THEN

       vapor_mass_rate = eval_vaporization_mass_rate(r_T, rho_s, sp_heat_s,   &
            deposition_term, gamma_steam)
       vapor_energy_rate = eval_vaporization_energy_rate(vapor_mass_rate,     &
            sp_heat_g(1))

       eqns_term(1) = eqns_term(1) + vapor_mass_rate
       eqns_term(4) = eqns_term(4) + vapor_energy_rate
       eqns_term(5+n_solid) = eqns_term(5+n_solid) + vapor_mass_rate

    END IF
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
END MODULE mass_exchange_vaporization_2d
