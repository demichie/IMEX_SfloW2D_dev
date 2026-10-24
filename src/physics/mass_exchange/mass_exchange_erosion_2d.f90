!> Erosion process: limit solid and carrier entrainment by available bed inventory.
MODULE mass_exchange_erosion_2d
  USE parameters_2d, ONLY : wp, n_solid
  USE constitutive_parameters_2d, ONLY : erosion_coeff, erodible_porosity, erodible_fract
  USE mass_exchange_laws_2d, ONLY : eval_erosion_rates
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_mass_exchange_erosion
CONTAINS
  SUBROUTINE eval_mass_exchange_erosion(r_h, mod_vel, alphas_tot, erodible, dt, &
       erosion_term, continuous_phase_erosion_term)
    REAL(wp), INTENT(IN) :: r_h, mod_vel, alphas_tot, erodible(n_solid), dt
    REAL(wp), INTENT(OUT) :: erosion_term(n_solid), continuous_phase_erosion_term
    IF (erosion_coeff .GT. 0.0_wp) THEN
       CALL eval_erosion_rates(r_h, mod_vel, alphas_tot, erosion_coeff,       &
            erodible_porosity, erodible_fract(1:n_solid),                    &
            erodible(1:n_solid), dt, erosion_term(1:n_solid),                &
            continuous_phase_erosion_term)
    ELSE
       erosion_term = 0.0_wp
       continuous_phase_erosion_term = 0.0_wp
    END IF
  END SUBROUTINE eval_mass_exchange_erosion
END MODULE mass_exchange_erosion_2d
