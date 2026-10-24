!> \brief Evaluate carrier-phase loss contributions as one physical batch.
!> \details Combine prescribed liquid, deposition-related and pore-gas
!> losses before applying the shared limiter exactly once.
MODULE mass_exchange_carrier_loss_2d
  USE parameters_2d, ONLY : wp, n_solid, liquid_flag, pore_pressure_flag, gas_loss_flag
  USE diagnostics_2d, ONLY : fatal_error
  USE constitutive_parameters_2d, ONLY : loss_rate, coeff_porosity, maximum_solid_packing
  USE pore_pressure_laws_2d, ONLY : eval_pore_pressure_gas_loss

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_mass_exchange_carrier_loss_explicit
  PUBLIC :: eval_mass_exchange_carrier_loss_semi_implicit
  PUBLIC :: eval_mass_exchange_carrier_loss_implicit
  PUBLIC :: solve_mass_exchange_carrier_loss_analytic
  PUBLIC :: eval_mass_exchange_carrier_liquid_loss_explicit
  PUBLIC :: eval_mass_exchange_carrier_liquid_loss_semi_implicit
  PUBLIC :: eval_mass_exchange_carrier_liquid_loss_implicit
  PUBLIC :: solve_mass_exchange_carrier_liquid_loss_analytic
  PUBLIC :: eval_mass_exchange_carrier_darcy_loss_explicit
  PUBLIC :: eval_mass_exchange_carrier_darcy_loss_semi_implicit
  PUBLIC :: eval_mass_exchange_carrier_darcy_loss_implicit
  PUBLIC :: solve_mass_exchange_carrier_darcy_loss_analytic
  PUBLIC :: eval_continuous_phase_loss_rate
  ! Four-entry physical treatment contract; only EXPLICIT is supported.
  ! Catalog-addressable carrier-loss contributors; only explicit is supported.
CONTAINS

  SUBROUTINE eval_mass_exchange_carrier_loss_explicit(r_h, r_T, r_alphas, alphas_tot, &
       r_rho_c, r_exc_pore_pres, &
       deposition_term, dt, continuous_phase_loss_term)
    REAL(wp), INTENT(IN) :: r_h, r_T, r_alphas(n_solid), alphas_tot
    REAL(wp), INTENT(IN) :: r_rho_c, r_exc_pore_pres, deposition_term(n_solid), dt
    REAL(wp), INTENT(OUT) :: continuous_phase_loss_term
    ! Compose independent requested rates, then apply the shared limiter once.
    continuous_phase_loss_term = eval_mass_exchange_carrier_liquid_loss_explicit() &
         + deposition_induced_carrier_loss(deposition_term)
    continuous_phase_loss_term = continuous_phase_loss_term + &
         eval_mass_exchange_carrier_darcy_loss_explicit(r_h, r_T, r_alphas, &
         alphas_tot, r_rho_c, r_exc_pore_pres)

    ! Limit carrier loss by solid packing and the available carrier reserve.
    continuous_phase_loss_term = eval_continuous_phase_loss_rate(             &
         continuous_phase_loss_term, r_h, alphas_tot, maximum_solid_packing,  &
         coeff_porosity, deposition_term(1:n_solid), dt)
  END SUBROUTINE eval_mass_exchange_carrier_loss_explicit

  SUBROUTINE eval_mass_exchange_carrier_loss_semi_implicit(r_h, r_T, r_alphas, alphas_tot, r_rho_c, r_exc_pore_pres, &
       deposition_term, dt, continuous_phase_loss_term)
    REAL(wp), INTENT(IN) :: r_h, r_T, r_alphas(n_solid), alphas_tot
    REAL(wp), INTENT(IN) :: r_rho_c, r_exc_pore_pres, deposition_term(n_solid), dt
    REAL(wp), INTENT(OUT) :: continuous_phase_loss_term
    CALL fatal_error('carrier_loss: semi_implicit treatment is unsupported')
  END SUBROUTINE eval_mass_exchange_carrier_loss_semi_implicit

  SUBROUTINE eval_mass_exchange_carrier_loss_implicit(r_h, r_T, r_alphas, alphas_tot, r_rho_c, r_exc_pore_pres, &
       deposition_term, dt, continuous_phase_loss_term)
    REAL(wp), INTENT(IN) :: r_h, r_T, r_alphas(n_solid), alphas_tot
    REAL(wp), INTENT(IN) :: r_rho_c, r_exc_pore_pres, deposition_term(n_solid), dt
    REAL(wp), INTENT(OUT) :: continuous_phase_loss_term
    CALL fatal_error('carrier_loss: implicit treatment is unsupported')
  END SUBROUTINE eval_mass_exchange_carrier_loss_implicit

  SUBROUTINE solve_mass_exchange_carrier_loss_analytic(r_h, r_T, r_alphas, alphas_tot, r_rho_c, r_exc_pore_pres, &
       deposition_term, dt, continuous_phase_loss_term)
    REAL(wp), INTENT(IN) :: r_h, r_T, r_alphas(n_solid), alphas_tot
    REAL(wp), INTENT(IN) :: r_rho_c, r_exc_pore_pres, deposition_term(n_solid), dt
    REAL(wp), INTENT(OUT) :: continuous_phase_loss_term
    CALL fatal_error('carrier_loss: analytic treatment is unsupported')
  END SUBROUTINE solve_mass_exchange_carrier_loss_analytic

  REAL(wp) FUNCTION eval_mass_exchange_carrier_liquid_loss_explicit() RESULT(rate)
    ! Requested constituent carrier loss; limiting is applied by the batch.
    rate = 0.0_wp
    IF (liquid_flag) rate = loss_rate
  END FUNCTION eval_mass_exchange_carrier_liquid_loss_explicit

  REAL(wp) FUNCTION eval_mass_exchange_carrier_liquid_loss_semi_implicit() RESULT(rate)
    ! Unsupported physical treatment.
    CALL fatal_error('mass_exchange_carrier_liquid_loss: semi_implicit is unsupported')
    rate = 0.0_wp
  END FUNCTION eval_mass_exchange_carrier_liquid_loss_semi_implicit

  REAL(wp) FUNCTION eval_mass_exchange_carrier_liquid_loss_implicit() RESULT(rate)
    ! Unsupported physical treatment.
    CALL fatal_error('mass_exchange_carrier_liquid_loss: implicit is unsupported')
    rate = 0.0_wp
  END FUNCTION eval_mass_exchange_carrier_liquid_loss_implicit

  REAL(wp) FUNCTION solve_mass_exchange_carrier_liquid_loss_analytic() RESULT(rate)
    ! Unsupported physical treatment.
    CALL fatal_error('mass_exchange_carrier_liquid_loss: analytic is unsupported')
    rate = 0.0_wp
  END FUNCTION solve_mass_exchange_carrier_liquid_loss_analytic

  REAL(wp) FUNCTION deposition_induced_carrier_loss(deposition_term) RESULT(rate)
    REAL(wp), INTENT(IN) :: deposition_term(n_solid)
    ! Requested constituent carrier loss; limiting is applied by the batch.
    rate = coeff_porosity * SUM(deposition_term(1:n_solid))
  END FUNCTION deposition_induced_carrier_loss


  REAL(wp) FUNCTION eval_mass_exchange_carrier_darcy_loss_explicit(r_h, r_T, r_alphas, alphas_tot, &
       r_rho_c, r_exc_pore_pres) RESULT(rate)
    REAL(wp), INTENT(IN) :: r_h, r_T, r_alphas(n_solid), alphas_tot
    REAL(wp), INTENT(IN) :: r_rho_c, r_exc_pore_pres
    ! Requested constituent carrier loss; limiting is applied by the batch.
    rate = 0.0_wp
    IF (pore_pressure_flag .AND. gas_loss_flag) THEN
       IF (alphas_tot .LT. maximum_solid_packing) THEN
          CALL eval_pore_pressure_gas_loss(r_h, r_T, r_alphas, alphas_tot, &
               r_rho_c, r_exc_pore_pres, rate)
       END IF
    END IF
  END FUNCTION eval_mass_exchange_carrier_darcy_loss_explicit

  REAL(wp) FUNCTION eval_mass_exchange_carrier_darcy_loss_semi_implicit(r_h, r_T, r_alphas, alphas_tot, &
       r_rho_c, r_exc_pore_pres) RESULT(rate)
    REAL(wp), INTENT(IN) :: r_h, r_T, r_alphas(n_solid), alphas_tot
    REAL(wp), INTENT(IN) :: r_rho_c, r_exc_pore_pres
    ! Unsupported physical treatment.
    CALL fatal_error('mass_exchange_carrier_darcy_loss: semi_implicit is unsupported')
    rate = 0.0_wp
  END FUNCTION eval_mass_exchange_carrier_darcy_loss_semi_implicit

  REAL(wp) FUNCTION eval_mass_exchange_carrier_darcy_loss_implicit(r_h, r_T, r_alphas, alphas_tot, &
       r_rho_c, r_exc_pore_pres) RESULT(rate)
    REAL(wp), INTENT(IN) :: r_h, r_T, r_alphas(n_solid), alphas_tot
    REAL(wp), INTENT(IN) :: r_rho_c, r_exc_pore_pres
    ! Unsupported physical treatment.
    CALL fatal_error('mass_exchange_carrier_darcy_loss: implicit is unsupported')
    rate = 0.0_wp
  END FUNCTION eval_mass_exchange_carrier_darcy_loss_implicit

  REAL(wp) FUNCTION solve_mass_exchange_carrier_darcy_loss_analytic(r_h, r_T, r_alphas, alphas_tot, &
       r_rho_c, r_exc_pore_pres) RESULT(rate)
    REAL(wp), INTENT(IN) :: r_h, r_T, r_alphas(n_solid), alphas_tot
    REAL(wp), INTENT(IN) :: r_rho_c, r_exc_pore_pres
    ! Unsupported physical treatment.
    CALL fatal_error('mass_exchange_carrier_darcy_loss: analytic is unsupported')
    rate = 0.0_wp
  END FUNCTION solve_mass_exchange_carrier_darcy_loss_analytic

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
















END MODULE mass_exchange_carrier_loss_2d
