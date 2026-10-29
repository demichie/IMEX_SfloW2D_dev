!> Cell-local mass exchange: process dispatch, shared-resource limiting, conservative assembly.
MODULE mass_exchange_2d
  USE parameters_2d
  USE constitutive_parameters_2d
  USE state_conversion_2d, ONLY: primitive_to_volume_fractions, mixt_var
  USE pore_pressure_laws_2d, ONLY: eval_pore_pressure_mass_exchange_source
  USE mass_exchange_erosion_2d, ONLY: eval_mass_exchange_erosion_explicit
  USE mass_exchange_deposition_2d, ONLY: eval_mass_exchange_deposition_explicit
  USE mass_exchange_carrier_liquid_loss_2d, ONLY: eval_mass_exchange_carrier_liquid_loss_explicit
  USE mass_exchange_carrier_darcy_loss_2d, ONLY: eval_mass_exchange_carrier_darcy_loss_explicit
  USE mass_exchange_entrainment_2d, ONLY: eval_mass_exchange_entrainment_explicit
  USE mass_exchange_vaporization_2d, ONLY: eval_mass_exchange_vaporization_explicit
  USE mass_exchange_lava_deposition_2d, ONLY: lava_deposition_flag, eval_mass_exchange_lava_deposition_explicit
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_mass_exchange_terms
CONTAINS
  SUBROUTINE eval_mass_exchange_terms(qp,B_zone,Bx,By,erodible,dt,rates,eq)
    REAL(wp), INTENT(IN) :: qp(n_vars+2),Bx,By,erodible(n_solid),dt
    INTEGER, INTENT(IN) :: B_zone
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: rates
    REAL(wp), INTENT(OUT) :: eq(n_eqns)
    TYPE(mass_exchange_rates_type) :: process,deposit_rates,liquid_rates,darcy_rates
    REAL(wp) :: alphas(n_solid),alphag(n_add_gas),alphal,alpha
    REAL(wp) :: u,v,w,speed2,h,T,pore,stoch
    REAL(wp) :: Ri,rhom,rhoc,redgrav,cpc,cpm,substrate_T
    REAL(wp) :: species(n_add_gas),gas_fraction
    REAL(wp) :: gas_available,liquid_available,gas_factor,liquid_factor
    REAL(wp) :: gas_withdrawals(2),liquid_withdrawals(2),carrier_volume
    REAL(wp) :: liquid_volume_rate,gas_volume_rate
    LOGICAL :: solids_enabled
    CALL reset_mass_exchange_rates(rates)
    eq=0.0_wp
    h=qp(1)
    IF (h<=EPSILON(1.0_wp)) RETURN
    u=qp(idx_u)
    v=qp(idx_v)
    T=qp(4)
    ! Processes requiring a gas/liquid carrier are evaluated only when that
    ! physical phase exists. Lava deposition is dispatched independently below.
    IF (gas_flag .OR. liquid_flag) THEN
    CALL primitive_to_volume_fractions(qp,alphas,alphag,alphal)
    alpha=SUM(alphas)
    solids_enabled=alpha>alphastot_min
    pore=0.0_wp
    stoch=0.0_wp
    IF (pore_pressure_flag) pore=qp(idx_pore)
    IF (stoch_transport_flag) stoch=qp(idx_stoch)
    species=qp(idx_add_gas_first:idx_add_gas_last)
    gas_fraction=0.0_wp
    IF (gas_flag) THEN
      gas_fraction=1.0_wp-SUM(qp(idx_solid_first:idx_solid_last))
      IF (liquid_flag) gas_fraction=gas_fraction-qp(n_vars)
      gas_fraction=MAX(0.0_wp,gas_fraction)
    END IF
    w=0.0_wp
    IF (slope_correction_flag) w=u*Bx+v*By
    speed2=u*u+v*v+w*w
    CALL mixt_var(qp,Ri,rhom,rhoc,redgrav,cpc,cpm)
    substrate_T=T_erodible
    IF (erodible_deposit_flag) substrate_T=T
    IF (solids_enabled) THEN
      CALL eval_mass_exchange_erosion_explicit(h,SQRT(speed2),alpha,erodible,dt,substrate_T,cpc,process)
      CALL accumulate_process(process,u,v,rates,eq)
      CALL eval_mass_exchange_deposition_explicit(h,T,alphas,alpha,rhom,rhoc,cpc,dt, &
           species,gas_fraction,deposit_rates)
    ELSE
      CALL reset_mass_exchange_rates(deposit_rates)
    END IF
    CALL eval_mass_exchange_carrier_liquid_loss_explicit(rhoc,T,cpc,liquid_rates)
    CALL eval_mass_exchange_carrier_darcy_loss_explicit(h,T,alphas,alpha,rhoc,pore,cpc, &
         species,gas_fraction,darcy_rates)
    ! Resource budgets use gross withdrawals, never signed net mass exchange.
    ! Gas and liquid are different inventories when both phases are present.
    ! The solid deposition rate is limited only when its required pore phase
    ! cannot be supplied; a later Darcy/liquid sink cannot consume it twice.
    carrier_volume=deposit_rates%continuous_phase_deposition_rate
    gas_volume_rate=MAX(0.0_wp,-darcy_rates%mixture_mass_rate/MAX(rhoc,TINY(1.0_wp)))
    liquid_volume_rate=MAX(0.0_wp,-liquid_rates%mixture_mass_rate/ &
         MAX(MERGE(rho_l,rhoc,gas_flag),TINY(1.0_wp)))
    gas_factor=1.0_wp
    liquid_factor=1.0_wp
    IF (gas_flag) THEN
      ! Subtract dispersed solids and liquid from the available gas volume.
      gas_available=MIN(h*MAX(0.0_wp,1.0_wp-alpha-alphal)/dt, &
           h*MAX(0.0_wp,maximum_solid_packing-alpha)/dt)
      gas_withdrawals=[-carrier_volume,-gas_volume_rate]
      gas_factor=shared_withdrawal_factor(gas_withdrawals,gas_available)
      CALL scale_process(deposit_rates,gas_factor)
      CALL scale_process(darcy_rates,gas_factor)
    ELSE
      ! In a liquid-only run, pore filling and prescribed loss share one pool.
      liquid_available=MIN(h*MAX(0.0_wp,1.0_wp-alpha)/dt, &
           h*MAX(0.0_wp,maximum_solid_packing-alpha)/dt)
      liquid_withdrawals=[-carrier_volume,-liquid_volume_rate]
      liquid_factor=shared_withdrawal_factor(liquid_withdrawals,liquid_available)
      CALL scale_process(deposit_rates,liquid_factor)
      CALL scale_process(liquid_rates,liquid_factor)
    END IF
    IF (gas_flag .AND. liquid_flag) THEN
      ! Prescribed liquid loss consumes dispersed liquid, NOT carrier gas.
      liquid_available=h*MAX(0.0_wp,alphal)/dt
      liquid_withdrawals=[-liquid_volume_rate,0.0_wp]
      liquid_factor=shared_withdrawal_factor(liquid_withdrawals,liquid_available)
      CALL scale_process(liquid_rates,liquid_factor)
    END IF
    CALL accumulate_process(deposit_rates,u,v,rates,eq)
    CALL accumulate_process(liquid_rates,u,v,rates,eq)
    CALL accumulate_process(darcy_rates,u,v,rates,eq)
    CALL eval_mass_exchange_entrainment_explicit(h,Ri,speed2,process)
    CALL accumulate_process(process,u,v,rates,eq)
    CALL eval_mass_exchange_vaporization_explicit(B_zone,T,deposit_rates%solid_deposition_rate,process)
    CALL accumulate_process(process,u,v,rates,eq)
    IF (stoch_transport_flag) eq(idx_stochEqn)=eq(1)*stoch
    IF (pore_pressure_flag) eq(idx_poreEqn)= &
         eval_pore_pressure_mass_exchange_source(eq(1),pore)
    END IF
    ! Lava deposition follows the same process-rate accumulation contract.
    IF (lava_deposition_flag) THEN
      CALL eval_mass_exchange_lava_deposition_explicit(h,T,dt,process)
      CALL accumulate_process(process,u,v,rates,eq)
    END IF
  END SUBROUTINE eval_mass_exchange_terms

  ! Limit simultaneous withdrawal demands from one physical inventory.
  ! Positive source rates are excluded from the budget, not netted against sinks.
  PURE REAL(wp) FUNCTION shared_withdrawal_factor(signed_rates,available) RESULT(f)
    REAL(wp), INTENT(IN) :: signed_rates(:),available
    REAL(wp) :: total_withdrawal
    total_withdrawal=SUM(MAX(0.0_wp,-signed_rates))
    f=1.0_wp
    IF (total_withdrawal>0.0_wp) &
         f=MIN(1.0_wp,MAX(0.0_wp,available)/total_withdrawal)
  END FUNCTION shared_withdrawal_factor

  SUBROUTINE accumulate_process(p,u,v,total,eq)
    TYPE(mass_exchange_rates_type), INTENT(IN) :: p
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: total
    REAL(wp), INTENT(IN) :: u,v
    REAL(wp), INTENT(INOUT) :: eq(n_eqns)
    eq(1)=eq(1)+p%mixture_mass_rate
    eq(2)=eq(2)+u*MIN(0.0_wp,p%mixture_mass_rate)
    eq(3)=eq(3)+v*MIN(0.0_wp,p%mixture_mass_rate)
    eq(4)=eq(4)+p%energy_rate
    eq(idx_solidEqn_first:idx_solidEqn_last)= &
         eq(idx_solidEqn_first:idx_solidEqn_last)+p%solid_mass_rate
    IF (gas_flag .AND. liquid_flag) eq(n_eqns)=eq(n_eqns)+p%liquid_mass_rate
    IF (n_add_gas>0) eq(idx_addGasEqn_first:idx_addGasEqn_last)= &
         eq(idx_addGasEqn_first:idx_addGasEqn_last)+p%add_gas_mass_rate
    total%mixture_mass_rate=total%mixture_mass_rate+p%mixture_mass_rate
    total%liquid_mass_rate=total%liquid_mass_rate+p%liquid_mass_rate
    total%energy_rate=total%energy_rate+p%energy_rate
    total%solid_mass_rate=total%solid_mass_rate+p%solid_mass_rate
    total%add_gas_mass_rate=total%add_gas_mass_rate+p%add_gas_mass_rate
    total%solid_erosion_rate=total%solid_erosion_rate+p%solid_erosion_rate
    total%solid_deposition_rate=total%solid_deposition_rate+p%solid_deposition_rate
    total%continuous_phase_erosion_rate= &
         total%continuous_phase_erosion_rate+p%continuous_phase_erosion_rate
    total%continuous_phase_deposition_rate= &
         total%continuous_phase_deposition_rate+p%continuous_phase_deposition_rate
    total%dB_dt=total%dB_dt+p%dB_dt
    total%lava_deposit_rate=total%lava_deposit_rate+p%lava_deposit_rate
  END SUBROUTINE accumulate_process

  SUBROUTINE scale_process(p,factor)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: p
    REAL(wp), INTENT(IN) :: factor
    p%mixture_mass_rate=p%mixture_mass_rate*factor
    p%liquid_mass_rate=p%liquid_mass_rate*factor
    p%energy_rate=p%energy_rate*factor
    p%solid_mass_rate=p%solid_mass_rate*factor
    p%add_gas_mass_rate=p%add_gas_mass_rate*factor
    p%solid_erosion_rate=p%solid_erosion_rate*factor
    p%solid_deposition_rate=p%solid_deposition_rate*factor
    p%continuous_phase_erosion_rate=p%continuous_phase_erosion_rate*factor
    p%continuous_phase_deposition_rate=p%continuous_phase_deposition_rate*factor
    p%dB_dt=p%dB_dt*factor
    p%lava_deposit_rate=p%lava_deposit_rate*factor
  END SUBROUTINE scale_process
END MODULE mass_exchange_2d
