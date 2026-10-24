!> Vaporization source contribution to conserved mass, energy and gas.
MODULE mass_exchange_vaporization_2d
  USE parameters_2d, ONLY : wp, n_eqns, n_solid, gas_flag, liquid_vaporization_flag
  USE constitutive_parameters_2d, ONLY : rho_s, sp_heat_s, gamma_steam, sp_heat_g
  USE thermal_laws_2d, ONLY : eval_vaporization_mass_rate, eval_vaporization_energy_rate
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: add_mass_exchange_vaporization
CONTAINS
  SUBROUTINE add_mass_exchange_vaporization(B_zone, r_T, deposition_term, eqns_term)
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
  END SUBROUTINE add_mass_exchange_vaporization
END MODULE mass_exchange_vaporization_2d
