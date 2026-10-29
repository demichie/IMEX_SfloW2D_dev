MODULE mass_exchange_carrier_darcy_loss_2d
  USE parameters_2d, ONLY : wp,n_solid,n_add_gas,n_vars,idx_pore, &
       idx_solid_first,idx_solid_last,idx_add_gas_first,idx_add_gas_last, &
       pore_pressure_flag,gas_loss_flag,gas_flag,liquid_flag, &
       mass_exchange_rates_type,reset_mass_exchange_rates
  USE constitutive_parameters_2d, ONLY : maximum_solid_packing
  USE pore_pressure_laws_2d, ONLY : eval_pore_pressure_gas_loss
  USE diagnostics_2d, ONLY : fatal_error
  USE state_conversion_2d, ONLY : primitive_to_volume_fractions, mixt_var
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_mass_exchange_carrier_darcy_loss_explicit
  PUBLIC :: eval_mass_exchange_carrier_darcy_loss_semi_implicit
  PUBLIC :: eval_mass_exchange_carrier_darcy_loss_implicit
  PUBLIC :: solve_mass_exchange_carrier_darcy_loss_analytic
CONTAINS
  SUBROUTINE eval_mass_exchange_carrier_darcy_loss_explicit(qp,r)
    REAL(wp), INTENT(IN) :: qp(n_vars+2)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    REAL(wp) :: volume_loss,h,T,alpha,rho_c,cp_c,pore
    REAL(wp) :: alphas(n_solid),alphag(n_add_gas),alphal
    REAL(wp) :: Ri,rhom,redgrav,cpm,gas_phase_mass_fraction
    REAL(wp) :: gas_species_mass_fraction(n_add_gas)
    CALL reset_mass_exchange_rates(r)
    IF (.NOT.(gas_flag .AND. pore_pressure_flag .AND. gas_loss_flag)) RETURN
    h=qp(1)
    T=qp(4)
    pore=qp(idx_pore)
    CALL primitive_to_volume_fractions(qp,alphas,alphag,alphal)
    alpha=SUM(alphas)
    CALL mixt_var(qp,Ri,rhom,rho_c,redgrav,cp_c,cpm)
    gas_species_mass_fraction=qp(idx_add_gas_first:idx_add_gas_last)
    gas_phase_mass_fraction=MAX(0.0_wp,1.0_wp- &
         SUM(qp(idx_solid_first:idx_solid_last)))
    IF (liquid_flag) gas_phase_mass_fraction= &
         MAX(0.0_wp,gas_phase_mass_fraction-qp(n_vars))
    volume_loss=0.0_wp
    IF (gas_flag .AND. pore_pressure_flag .AND. gas_loss_flag .AND. &
         alpha<maximum_solid_packing) &
         CALL eval_pore_pressure_gas_loss(h,T,alphas,alpha,rho_c,pore,volume_loss)
    r%gas_withdrawal_volume_rate=MAX(0.0_wp,volume_loss)
    r%mixture_mass_rate=-rho_c*r%gas_withdrawal_volume_rate
    r%energy_rate=T*cp_c*r%mixture_mass_rate
    IF (n_add_gas>0 .AND. gas_phase_mass_fraction>0.0_wp) &
         r%add_gas_mass_rate=r%mixture_mass_rate* &
         gas_species_mass_fraction/gas_phase_mass_fraction
  END SUBROUTINE eval_mass_exchange_carrier_darcy_loss_explicit

  SUBROUTINE eval_mass_exchange_carrier_darcy_loss_semi_implicit(qp,r)
    REAL(wp), INTENT(IN) :: qp(n_vars+2)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('carrier_darcy_loss: unsupported semi_implicit method')
  END SUBROUTINE eval_mass_exchange_carrier_darcy_loss_semi_implicit

  SUBROUTINE eval_mass_exchange_carrier_darcy_loss_implicit(qp,r)
    REAL(wp), INTENT(IN) :: qp(n_vars+2)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('carrier_darcy_loss: unsupported implicit method')
  END SUBROUTINE eval_mass_exchange_carrier_darcy_loss_implicit

  SUBROUTINE solve_mass_exchange_carrier_darcy_loss_analytic(qp,r)
    REAL(wp), INTENT(IN) :: qp(n_vars+2)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('carrier_darcy_loss: unsupported analytic method')
  END SUBROUTINE solve_mass_exchange_carrier_darcy_loss_analytic

END MODULE mass_exchange_carrier_darcy_loss_2d
