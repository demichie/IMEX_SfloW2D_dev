!> Air entrainment contribution to mass exchange.
MODULE mass_exchange_entrainment_2d
  USE parameters_2d, ONLY : wp
  USE constitutive_parameters_2d, ONLY : entrainment_flag
  USE mass_exchange_laws_2d, ONLY : eval_air_entrainment_rate
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_mass_exchange_air_entrainment
CONTAINS
  REAL(wp) FUNCTION eval_mass_exchange_air_entrainment(r_h, r_Ri, mod_vel2)
    REAL(wp), INTENT(IN) :: r_h, r_Ri, mod_vel2
    IF ( entrainment_flag .AND.  ( r_h .GT. 0.0_wp ) .AND.                      &
         ( r_Ri .GT. 0.0_wp ) ) THEN

       eval_mass_exchange_air_entrainment = eval_air_entrainment_rate(r_Ri, mod_vel2)

    ELSE

       eval_mass_exchange_air_entrainment = 0.0_wp

    END IF
  END FUNCTION eval_mass_exchange_air_entrainment
END MODULE mass_exchange_entrainment_2d
