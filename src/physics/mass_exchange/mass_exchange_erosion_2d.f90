!> Erosion process: limit solid and carrier entrainment by available bed inventory.
MODULE mass_exchange_erosion_2d
  USE parameters_2d, ONLY : wp, n_solid
  USE diagnostics_2d, ONLY : fatal_error
  USE constitutive_parameters_2d, ONLY : erosion_coeff, erodible_porosity, erodible_fract

  PRIVATE
  PUBLIC :: eval_erosion_rates
  PUBLIC :: eval_mass_exchange_erosion
  ! Four-entry physical treatment contract; only EXPLICIT is supported.
  PUBLIC :: eval_erosion_explicit
  PUBLIC :: eval_erosion_semi_implicit
  PUBLIC :: eval_erosion_implicit
  PUBLIC :: eval_erosion_analytic
  INTERFACE eval_erosion_explicit
    MODULE PROCEDURE eval_mass_exchange_erosion
  END INTERFACE eval_erosion_explicit
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
  PURE SUBROUTINE eval_erosion_rates(flow_depth, speed, total_solid_fraction, &
       erosion_coefficient, erodible_porosity, erodible_fraction,            &
       erodible_inventory, time_step, solid_erosion_rate,                    &
       continuous_phase_erosion_rate)

    REAL(wp), INTENT(IN) :: flow_depth, speed, total_solid_fraction
    REAL(wp), INTENT(IN) :: erosion_coefficient, erodible_porosity
    REAL(wp), INTENT(IN) :: erodible_fraction(:), erodible_inventory(:)
    REAL(wp), INTENT(IN) :: time_step
    REAL(wp), INTENT(OUT) :: solid_erosion_rate(:)
    REAL(wp), INTENT(OUT) :: continuous_phase_erosion_rate

    REAL(wp) :: total_erosion, total_solid_erosion

    IF (erosion_coefficient .GT. 0.0_wp) THEN
       total_erosion = erosion_coefficient * speed * flow_depth *            &
            (1.0_wp - total_solid_fraction)
       total_solid_erosion = total_erosion * (1.0_wp - erodible_porosity)
       solid_erosion_rate = erodible_fraction * total_solid_erosion
    ELSE
       solid_erosion_rate = 0.0_wp
    END IF

    solid_erosion_rate = MAX(0.0_wp, MIN(solid_erosion_rate,                 &
         erodible_inventory / time_step))
    total_solid_erosion = SUM(solid_erosion_rate)
    total_erosion = total_solid_erosion / (1.0_wp - erodible_porosity)
    continuous_phase_erosion_rate = total_erosion * erodible_porosity

  END SUBROUTINE eval_erosion_rates

  SUBROUTINE eval_erosion_semi_implicit(r_h, mod_vel, alphas_tot, erodible, dt, erosion_term, continuous_phase_erosion_term)
    REAL(wp), INTENT(IN) :: r_h, mod_vel, alphas_tot, erodible(n_solid), dt
    REAL(wp), INTENT(OUT) :: erosion_term(n_solid), continuous_phase_erosion_term
    CALL fatal_error('erosion: semi_implicit treatment is unsupported')
  END SUBROUTINE eval_erosion_semi_implicit

  SUBROUTINE eval_erosion_implicit(r_h, mod_vel, alphas_tot, erodible, dt, erosion_term, continuous_phase_erosion_term)
    REAL(wp), INTENT(IN) :: r_h, mod_vel, alphas_tot, erodible(n_solid), dt
    REAL(wp), INTENT(OUT) :: erosion_term(n_solid), continuous_phase_erosion_term
    CALL fatal_error('erosion: implicit treatment is unsupported')
  END SUBROUTINE eval_erosion_implicit

  SUBROUTINE eval_erosion_analytic(r_h, mod_vel, alphas_tot, erodible, dt, erosion_term, continuous_phase_erosion_term)
    REAL(wp), INTENT(IN) :: r_h, mod_vel, alphas_tot, erodible(n_solid), dt
    REAL(wp), INTENT(OUT) :: erosion_term(n_solid), continuous_phase_erosion_term
    CALL fatal_error('erosion: analytic treatment is unsupported')
  END SUBROUTINE eval_erosion_analytic
END MODULE mass_exchange_erosion_2d
