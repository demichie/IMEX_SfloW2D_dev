MODULE mass_exchange_carrier_darcy_loss_2d
  USE parameters_2d, ONLY: wp,n_solid,n_add_gas,pore_pressure_flag,gas_loss_flag, &
       mass_exchange_rates_type,reset_mass_exchange_rates
  USE constitutive_parameters_2d, ONLY: maximum_solid_packing
  USE pore_pressure_laws_2d, ONLY: eval_pore_pressure_gas_loss
  USE diagnostics_2d, ONLY: fatal_error
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_mass_exchange_carrier_darcy_loss_explicit
  PUBLIC :: eval_mass_exchange_carrier_darcy_loss_semi_implicit
  PUBLIC :: eval_mass_exchange_carrier_darcy_loss_implicit
  PUBLIC :: solve_mass_exchange_carrier_darcy_loss_analytic
CONTAINS
  SUBROUTINE eval_mass_exchange_carrier_darcy_loss_explicit(h,T,alphas,alpha,rho_c,pore,cp_c,gas_species_fraction,r)
    REAL(wp), INTENT(IN) :: h,T,alphas(n_solid),alpha,rho_c,pore,cp_c
    REAL(wp), INTENT(IN) :: gas_species_fraction(n_add_gas)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    REAL(wp) :: volume_loss
    CALL reset_mass_exchange_rates(r)
    volume_loss=0.0_wp
    IF (pore_pressure_flag .AND. gas_loss_flag .AND. alpha<maximum_solid_packing) &
         CALL eval_pore_pressure_gas_loss(h,T,alphas,alpha,rho_c,pore,volume_loss)
    r%mixture_mass_rate=-rho_c*MAX(0.0_wp,volume_loss)
    r%energy_rate=T*cp_c*r%mixture_mass_rate
    IF (n_add_gas>0) r%add_gas_mass_rate= &
         r%mixture_mass_rate*gas_species_fraction
  END SUBROUTINE eval_mass_exchange_carrier_darcy_loss_explicit

  SUBROUTINE eval_mass_exchange_carrier_darcy_loss_semi_implicit(h,T,alphas,alpha,rho_c,pore,cp_c,gas_species_fraction,r)
    REAL(wp), INTENT(IN) :: h,T,alphas(n_solid),alpha,rho_c,pore,cp_c
    REAL(wp), INTENT(IN) :: gas_species_fraction(n_add_gas)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('carrier_darcy_loss: unsupported semi_implicit method')
  END SUBROUTINE eval_mass_exchange_carrier_darcy_loss_semi_implicit

  SUBROUTINE eval_mass_exchange_carrier_darcy_loss_implicit(h,T,alphas,alpha,rho_c,pore,cp_c,gas_species_fraction,r)
    REAL(wp), INTENT(IN) :: h,T,alphas(n_solid),alpha,rho_c,pore,cp_c
    REAL(wp), INTENT(IN) :: gas_species_fraction(n_add_gas)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('carrier_darcy_loss: unsupported implicit method')
  END SUBROUTINE eval_mass_exchange_carrier_darcy_loss_implicit

  SUBROUTINE solve_mass_exchange_carrier_darcy_loss_analytic(h,T,alphas,alpha,rho_c,pore,cp_c,gas_species_fraction,r)
    REAL(wp), INTENT(IN) :: h,T,alphas(n_solid),alpha,rho_c,pore,cp_c
    REAL(wp), INTENT(IN) :: gas_species_fraction(n_add_gas)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('carrier_darcy_loss: unsupported analytic method')
  END SUBROUTINE solve_mass_exchange_carrier_darcy_loss_analytic

END MODULE mass_exchange_carrier_darcy_loss_2d
