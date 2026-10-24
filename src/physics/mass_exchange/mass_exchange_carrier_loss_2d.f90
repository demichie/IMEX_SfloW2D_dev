!> Coupled carrier-phase loss: prescribed liquid loss, deposition and pore gas.
!> Apply the shared physical limiter only once, after all contributions.
MODULE mass_exchange_carrier_loss_2d
  USE parameters_2d, ONLY : wp, n_solid, liquid_flag, pore_pressure_flag, gas_loss_flag
  USE constitutive_parameters_2d, ONLY : loss_rate, coeff_porosity, maximum_solid_packing
  USE mass_exchange_laws_2d, ONLY : eval_continuous_phase_loss_rate
  USE pore_pressure_laws_2d, ONLY : eval_pore_pressure_gas_loss
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_total_carrier_loss_rate
CONTAINS
  SUBROUTINE eval_total_carrier_loss_rate(r_h, r_T, r_alphas, alphas_tot, &
       r_rho_c, r_exc_pore_pres, deposition_term, dt, continuous_phase_loss_term)
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
END MODULE mass_exchange_carrier_loss_2d
