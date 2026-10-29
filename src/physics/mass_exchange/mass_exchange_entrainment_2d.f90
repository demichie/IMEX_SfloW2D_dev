!> Air entrainment contribution to mass exchange.
MODULE mass_exchange_entrainment_2d
  USE parameters_2d, ONLY : wp
  USE diagnostics_2d, ONLY : fatal_error
  USE constitutive_parameters_2d, ONLY : entrainment_flag

  USE mass_exchange_rates_2d, ONLY: mass_exchange_rates_type, reset_mass_exchange_rates
  USE constitutive_parameters_2d, ONLY: rho_a_amb, T_ambient, sp_heat_a
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_mass_exchange_entrainment_rates_explicit
  PUBLIC :: eval_mass_exchange_entrainment_explicit
  PUBLIC :: eval_mass_exchange_entrainment_semi_implicit
  PUBLIC :: eval_mass_exchange_entrainment_implicit
  PUBLIC :: solve_mass_exchange_entrainment_analytic
  PUBLIC :: eval_air_entrainment_rate
  ! Four-entry physical treatment contract; only EXPLICIT is supported.
CONTAINS

  REAL(wp) FUNCTION eval_mass_exchange_entrainment_explicit(r_h, r_Ri, mod_vel2)
    REAL(wp), INTENT(IN) :: r_h, r_Ri, mod_vel2
    IF ( entrainment_flag .AND.  ( r_h .GT. 0.0_wp ) .AND.                      &
         ( r_Ri .GT. 0.0_wp ) ) THEN

       eval_mass_exchange_entrainment_explicit = eval_air_entrainment_rate(r_Ri, mod_vel2)

    ELSE

       eval_mass_exchange_entrainment_explicit = 0.0_wp

    END IF
  END FUNCTION eval_mass_exchange_entrainment_explicit

  REAL(wp) FUNCTION eval_mass_exchange_entrainment_semi_implicit(r_h, r_Ri, mod_vel2)
    REAL(wp), INTENT(IN) :: r_h, r_Ri, mod_vel2
    CALL fatal_error('entrainment: semi_implicit treatment is unsupported')
    eval_mass_exchange_entrainment_semi_implicit = 0.0_wp
  END FUNCTION eval_mass_exchange_entrainment_semi_implicit

  REAL(wp) FUNCTION eval_mass_exchange_entrainment_implicit(r_h, r_Ri, mod_vel2)
    REAL(wp), INTENT(IN) :: r_h, r_Ri, mod_vel2
    CALL fatal_error('entrainment: implicit treatment is unsupported')
    eval_mass_exchange_entrainment_implicit = 0.0_wp
  END FUNCTION eval_mass_exchange_entrainment_implicit

  REAL(wp) FUNCTION solve_mass_exchange_entrainment_analytic(r_h, r_Ri, mod_vel2)
    REAL(wp), INTENT(IN) :: r_h, r_Ri, mod_vel2
    CALL fatal_error('entrainment: analytic treatment is unsupported')
    solve_mass_exchange_entrainment_analytic = 0.0_wp
  END FUNCTION solve_mass_exchange_entrainment_analytic

  PURE FUNCTION eval_air_entrainment_rate(richardson_number, speed_squared)  &
       RESULT(entrainment_rate)

    REAL(wp), INTENT(IN) :: richardson_number, speed_squared
    REAL(wp) :: entrainment_coefficient, entrainment_rate

    entrainment_coefficient = 0.075_wp / SQRT(1.0_wp + 718.0_wp *            &
         richardson_number**2.4_wp)
    entrainment_rate = entrainment_coefficient * SQRT(speed_squared)

  END FUNCTION eval_air_entrainment_rate



  SUBROUTINE eval_mass_exchange_entrainment_rates_explicit(h,Ri,speed_squared,r)
    REAL(wp), INTENT(IN) :: h,Ri,speed_squared
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    r%mixture_mass_rate=rho_a_amb*eval_mass_exchange_entrainment_explicit(h,Ri,speed_squared)
    r%energy_rate=T_ambient*sp_heat_a*r%mixture_mass_rate
  END SUBROUTINE
END MODULE mass_exchange_entrainment_2d
