!> Air entrainment contribution to mass exchange.
MODULE mass_exchange_entrainment_2d
  USE parameters_2d, ONLY : wp
  USE diagnostics_2d, ONLY : fatal_error
  USE constitutive_parameters_2d, ONLY : entrainment_flag

  PRIVATE
  PUBLIC :: eval_air_entrainment_rate
  ! Four-entry physical treatment contract; only EXPLICIT is supported.
  PUBLIC :: eval_mass_exchange_entrainment_explicit
  PUBLIC :: eval_mass_exchange_entrainment_semi_implicit
  PUBLIC :: eval_mass_exchange_entrainment_implicit
  PUBLIC :: solve_mass_exchange_entrainment_analytic
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
  PURE FUNCTION eval_air_entrainment_rate(richardson_number, speed_squared)  &
       RESULT(entrainment_rate)

    REAL(wp), INTENT(IN) :: richardson_number, speed_squared
    REAL(wp) :: entrainment_coefficient, entrainment_rate

    entrainment_coefficient = 0.075_wp / SQRT(1.0_wp + 718.0_wp *            &
         richardson_number**2.4_wp)
    entrainment_rate = entrainment_coefficient * SQRT(speed_squared)

  END FUNCTION eval_air_entrainment_rate

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
END MODULE mass_exchange_entrainment_2d
