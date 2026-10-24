!> Coupled carrier-phase loss: prescribed liquid loss, deposition and pore gas.
!> Apply the shared physical limiter only once, after all contributions.
MODULE mass_exchange_carrier_loss_2d
  USE parameters_2d, ONLY : wp, n_solid, liquid_flag, pore_pressure_flag, gas_loss_flag
  USE diagnostics_2d, ONLY : fatal_error
  USE constitutive_parameters_2d, ONLY : loss_rate, coeff_porosity, maximum_solid_packing
  USE pore_pressure_laws_2d, ONLY : eval_pore_pressure_gas_loss

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_continuous_phase_loss_rate
  PUBLIC :: eval_total_carrier_loss_rate
  ! Four-entry physical treatment contract; only EXPLICIT is supported.
  PUBLIC :: eval_carrier_loss_explicit
  PUBLIC :: eval_carrier_loss_semi_implicit
  PUBLIC :: eval_carrier_loss_implicit
  PUBLIC :: eval_carrier_loss_analytic
  INTERFACE eval_carrier_loss_explicit
    MODULE PROCEDURE eval_total_carrier_loss_rate
  END INTERFACE eval_carrier_loss_explicit
CONTAINS
  SUBROUTINE eval_total_carrier_loss_rate(r_h, r_T, r_alphas, alphas_tot, &
       r_rho_c, r_exc_pore_pres, &
       deposition_term, dt, continuous_phase_loss_term)
    REAL(wp), INTENT(IN) :: r_h, r_T, r_alphas(n_solid), alphas_tot
    REAL(wp), INTENT(IN) :: r_rho_c, r_exc_pore_pres, deposition_term(n_solid), dt
    REAL(wp), INTENT(OUT) :: continuous_phase_loss_term
    REAL(wp) :: pore_pressure_term

    IF ( liquid_flag ) THEN

       ! set the rate of loss of continuous phase
       continuous_phase_loss_term = loss_rate

    ELSE

       continuous_phase_loss_term = 0.0_wp

    END IF

    ! add the loss associated with solid deposition
    continuous_phase_loss_term =  continuous_phase_loss_term +                  &
         coeff_porosity * SUM( deposition_term(1:n_solid) )

    IF ( pore_pressure_flag .AND. gas_loss_flag ) THEN

       IF ( alphas_tot .LT. maximum_solid_packing ) THEN

         CALL eval_pore_pressure_gas_loss(r_h, r_T, r_alphas, alphas_tot,      &
              r_rho_c, r_exc_pore_pres, pore_pressure_term)

         continuous_phase_loss_term = continuous_phase_loss_term +             &
               pore_pressure_term

       END IF

    END IF

    ! Limit carrier loss by solid packing and the available carrier reserve.
    continuous_phase_loss_term = eval_continuous_phase_loss_rate(             &
         continuous_phase_loss_term, r_h, alphas_tot, maximum_solid_packing,  &
         coeff_porosity, deposition_term(1:n_solid), dt)
  END SUBROUTINE eval_total_carrier_loss_rate
  PURE FUNCTION eval_continuous_phase_loss_rate(requested_loss_rate,         &
       flow_depth, total_solid_fraction, maximum_solid_packing,               &
       porosity_coefficient, deposition_rate, time_step) RESULT(loss_rate)

    REAL(wp), INTENT(IN) :: requested_loss_rate, flow_depth
    REAL(wp), INTENT(IN) :: total_solid_fraction, maximum_solid_packing
    REAL(wp), INTENT(IN) :: porosity_coefficient, deposition_rate(:)
    REAL(wp), INTENT(IN) :: time_step
    REAL(wp) :: loss_rate

    loss_rate = MIN(requested_loss_rate, flow_depth * MAX(0.0_wp,            &
         maximum_solid_packing - total_solid_fraction) / time_step)
    loss_rate = MIN(loss_rate, MAX(0.0_wp, (flow_depth *                     &
         (1.0_wp - total_solid_fraction) - porosity_coefficient *            &
         (flow_depth * total_solid_fraction - time_step *                    &
         SUM(deposition_rate))) / time_step))

  END FUNCTION eval_continuous_phase_loss_rate

  SUBROUTINE eval_carrier_loss_semi_implicit(r_h, r_T, r_alphas, alphas_tot, r_rho_c, r_exc_pore_pres, &
       deposition_term, dt, continuous_phase_loss_term)
    REAL(wp), INTENT(IN) :: r_h, r_T, r_alphas(n_solid), alphas_tot
    REAL(wp), INTENT(IN) :: r_rho_c, r_exc_pore_pres, deposition_term(n_solid), dt
    REAL(wp), INTENT(OUT) :: continuous_phase_loss_term
    CALL fatal_error('carrier_loss: semi_implicit treatment is unsupported')
  END SUBROUTINE eval_carrier_loss_semi_implicit

  SUBROUTINE eval_carrier_loss_implicit(r_h, r_T, r_alphas, alphas_tot, r_rho_c, r_exc_pore_pres, &
       deposition_term, dt, continuous_phase_loss_term)
    REAL(wp), INTENT(IN) :: r_h, r_T, r_alphas(n_solid), alphas_tot
    REAL(wp), INTENT(IN) :: r_rho_c, r_exc_pore_pres, deposition_term(n_solid), dt
    REAL(wp), INTENT(OUT) :: continuous_phase_loss_term
    CALL fatal_error('carrier_loss: implicit treatment is unsupported')
  END SUBROUTINE eval_carrier_loss_implicit

  SUBROUTINE eval_carrier_loss_analytic(r_h, r_T, r_alphas, alphas_tot, r_rho_c, r_exc_pore_pres, &
       deposition_term, dt, continuous_phase_loss_term)
    REAL(wp), INTENT(IN) :: r_h, r_T, r_alphas(n_solid), alphas_tot
    REAL(wp), INTENT(IN) :: r_rho_c, r_exc_pore_pres, deposition_term(n_solid), dt
    REAL(wp), INTENT(OUT) :: continuous_phase_loss_term
    CALL fatal_error('carrier_loss: analytic treatment is unsupported')
  END SUBROUTINE eval_carrier_loss_analytic
END MODULE mass_exchange_carrier_loss_2d
