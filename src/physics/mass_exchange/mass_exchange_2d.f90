!> @brief mass exchange 2d: physical mass-exchange components and rate evaluation.
!> @details Explicit physics is evaluated locally; unsupported treatment entry points
!> fail explicitly. The coordinator alone enforces shared carrier inventories.
MODULE mass_exchange_2d
    USE parameters_2d, ONLY : wp, n_vars, n_eqns, n_solid, n_add_gas, &
              idx_u, idx_v, idx_pore, idx_stoch, idx_stochEqn, idx_poreEqn, &
              idx_solidEqn_first, idx_solidEqn_last, &
              idx_addGasEqn_first, idx_addGasEqn_last, gas_flag, liquid_flag, &
              pore_pressure_flag, stoch_transport_flag, mass_exchange_rates_type, &
              reset_mass_exchange_rates
    USE constitutive_parameters_2d, ONLY : maximum_solid_packing
    USE state_conversion_2d, ONLY : primitive_to_volume_fractions
    USE pore_pressure_laws_2d, ONLY : eval_pore_pressure_mass_exchange_source
    USE mass_exchange_erosion_2d, ONLY : eval_mass_exchange_erosion_explicit
    USE mass_exchange_deposition_2d, ONLY : eval_mass_exchange_deposition_explicit
    USE mass_exchange_carrier_liquid_loss_2d, ONLY : eval_mass_exchange_carrier_liquid_loss_explicit
    USE mass_exchange_carrier_darcy_loss_2d, ONLY : eval_mass_exchange_carrier_darcy_loss_explicit
    USE mass_exchange_entrainment_2d, ONLY : eval_mass_exchange_entrainment_explicit
    USE mass_exchange_vaporization_2d, ONLY : eval_mass_exchange_vaporization_explicit
    USE mass_exchange_lava_deposition_2d, ONLY : lava_deposition_flag, eval_mass_exchange_lava_deposition_explicit

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: eval_mass_exchange_terms

CONTAINS

    !> @brief Coordinate cell-local mass exchange, shared phase budgets and conservative source assembly.
    !> @details Applies available-phase limits before conservative assembly.
    !> Positive sources never cancel gross withdrawals in a shared inventory.
    !> @param[in] qp Local primitive flow state; component order follows model configuration.
    !> @param[in] B_zone Substrate zone identifier controlling vaporization eligibility.
    !> @param[in] Bx Bed elevation gradient in the x direction [dimensionless].
    !> @param[in] By Bed elevation gradient in the y direction [dimensionless].
    !> @param[in] erodible Available erodible solid inventory by species [m].
    !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
    !> @param[in,out] total_rates Accumulated process rates and bed changes for this cell.
    !> @param[out] eqns_term Conservative mass, momentum, energy and species source vector.
    SUBROUTINE eval_mass_exchange_terms(qp,B_zone,Bx,By,erodible,dt,total_rates,eqns_term)
        REAL(wp), INTENT(IN) :: qp(n_vars+2),Bx,By,erodible(n_solid),dt
        INTEGER, INTENT(IN) :: B_zone
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: total_rates
        REAL(wp), INTENT(OUT) :: eqns_term(n_eqns)
        TYPE(mass_exchange_rates_type) :: process_rates,deposition_rates,liquid_loss_rates,darcy_loss_rates
        REAL(wp) :: alphas(n_solid), alphag(n_add_gas), alphal, alpha
        REAL(wp) :: u,v,h,T,pore,stoch
        REAL(wp) :: gas_available,liquid_available,gas_factor,liquid_factor
        REAL(wp) :: gas_withdrawals(2),liquid_withdrawals(2),carrier_volume
        REAL(wp) :: solid_erosion,carrier_erosion,solid_deposition
        REAL(wp) :: liquid_loss,gas_loss,packing_factor

        CALL reset_mass_exchange_rates(total_rates)
        eqns_term=0.0_wp
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
            pore=0.0_wp
            stoch=0.0_wp
            IF (pore_pressure_flag) pore=qp(idx_pore)
            IF (stoch_transport_flag) stoch=qp(idx_stoch)
            CALL eval_mass_exchange_erosion_explicit(qp,Bx,By,alpha,erodible, &
                      dt,process_rates)
            CALL accumulate_process(process_rates,u,v,total_rates,eqns_term)
            CALL eval_mass_exchange_deposition_explicit(qp,dt,deposition_rates)
            CALL eval_mass_exchange_carrier_liquid_loss_explicit(qp,liquid_loss_rates)
            CALL eval_mass_exchange_carrier_darcy_loss_explicit(qp,darcy_loss_rates)

            ! Resource budgets use gross withdrawals, never signed net mass exchange.
            ! Gas and liquid are different inventories when both phases are present.
            ! The solid deposition rate is limited only when its required pore phase
            ! cannot be supplied; a later Darcy/liquid sink cannot consume it twice.
            carrier_volume=MAX(0.0_wp,deposition_rates%continuous_phase_deposition_rate)
            solid_deposition=SUM(MAX(0.0_wp,deposition_rates%solid_deposition_rate))
            solid_erosion=SUM(MAX(0.0_wp,process_rates%solid_erosion_rate))
            carrier_erosion=MAX(0.0_wp,process_rates%continuous_phase_erosion_rate)
            gas_loss=MAX(0.0_wp,darcy_loss_rates%gas_withdrawal_volume_rate)
            liquid_loss=MAX(0.0_wp,liquid_loss_rates%liquid_withdrawal_volume_rate)
            gas_factor=1.0_wp
            liquid_factor=1.0_wp

            ! Limit the dispersed-liquid sink first in two-carrier flows. It
            ! shares neither the gas inventory nor the deposition pore gas.
            IF (gas_flag .AND. liquid_flag) THEN

                liquid_available=h*MAX(0.0_wp,alphal)/dt
                liquid_factor=shared_withdrawal_factor([-liquid_loss],liquid_available)
                CALL scale_process(liquid_loss_rates,liquid_factor)
                liquid_loss=liquid_loss*liquid_factor
            END IF

            ! Gas withdrawal and pore filling draw on the same gas inventory;
            ! do not offset physical sinks with any positive source.
            IF (gas_flag) THEN

                ! Subtract dispersed solids and liquid from the available gas volume.
                gas_available=h*MAX(0.0_wp,1.0_wp-alpha-alphal)/dt
                gas_withdrawals=[-carrier_volume,-gas_loss]
                gas_factor=shared_withdrawal_factor(gas_withdrawals,gas_available)

                ! Erosion increases the mobile solid and carrier volumes.
                ! Solid deposition reduces both solid volume and total depth:
                ! the packing constraint must use these coupled increments.
                packing_factor=packing_withdrawal_factor(h,alpha,solid_erosion, &
                          carrier_erosion,solid_deposition,carrier_volume, &
                          gas_loss,liquid_loss,maximum_solid_packing,dt)
                gas_factor=MIN(gas_factor,packing_factor)
                CALL scale_process(deposition_rates,gas_factor)
                CALL scale_process(darcy_loss_rates,gas_factor)

            ELSE

                ! In a liquid-only run, pore filling and prescribed loss share one pool.
                liquid_available=h*MAX(0.0_wp,1.0_wp-alpha)/dt
                liquid_withdrawals=[-carrier_volume,-liquid_loss]
                liquid_factor=shared_withdrawal_factor(liquid_withdrawals,liquid_available)
                packing_factor=packing_withdrawal_factor(h,alpha,solid_erosion, &
                          carrier_erosion,solid_deposition,carrier_volume, &
                          liquid_loss,0.0_wp,maximum_solid_packing,dt)
                liquid_factor=MIN(liquid_factor,packing_factor)
                CALL scale_process(deposition_rates,liquid_factor)
                CALL scale_process(liquid_loss_rates,liquid_factor)
            END IF

            ! Assemble only the rates after all common-inventory limiters have acted.
            CALL accumulate_process(deposition_rates,u,v,total_rates,eqns_term)
            CALL accumulate_process(liquid_loss_rates,u,v,total_rates,eqns_term)
            CALL accumulate_process(darcy_loss_rates,u,v,total_rates,eqns_term)
            CALL eval_mass_exchange_entrainment_explicit(qp,u,v,Bx,By,process_rates)
            CALL accumulate_process(process_rates,u,v,total_rates,eqns_term)
            CALL eval_mass_exchange_vaporization_explicit(B_zone,T,deposition_rates%solid_deposition_rate,process_rates)
            CALL accumulate_process(process_rates,u,v,total_rates,eqns_term)

            ! Scalar transport and pore-pressure sources inherit the net mass exchange.
            IF (stoch_transport_flag) eqns_term(idx_stochEqn)=eqns_term(1)*stoch
            IF (pore_pressure_flag) eqns_term(idx_poreEqn)= &
                      eval_pore_pressure_mass_exchange_source(eqns_term(1),pore)
        END IF  ! gas_flag .OR. liquid_flag

        ! Lava deposition follows the same process_rates-rate accumulation contract.
        IF (lava_deposition_flag) THEN
            CALL eval_mass_exchange_lava_deposition_explicit(h,T,dt,process_rates)
            CALL accumulate_process(process_rates,u,v,total_rates,eqns_term)
        END IF

    END SUBROUTINE eval_mass_exchange_terms

    ! Limit simultaneous withdrawal demands from one physical inventory.
    ! Positive source total_rates are excluded from the budget, not netted against sinks.

    !> @brief Limit simultaneous positive withdrawal demands against a single available inventory.
    !> @details Positive signed rates are excluded from gross withdrawal budgets.
    !> @param[in] signed_rates Signed process flux demands; negative values mean withdrawals.
    !> @param[in] available Maximum withdrawal flux available from the shared phase inventory.
    !> @return Dimensionless limiter in [0, 1].
    PURE REAL(wp) FUNCTION shared_withdrawal_factor(signed_rates,available) RESULT(f)
        REAL(wp), INTENT(IN) :: signed_rates(:),available
        REAL(wp) :: total_withdrawal

        total_withdrawal=SUM(MAX(0.0_wp,-signed_rates))
        f=1.0_wp
        IF (total_withdrawal>0.0_wp) &
                  f=MIN(1.0_wp,MAX(0.0_wp,available)/total_withdrawal)
    END FUNCTION shared_withdrawal_factor

    !> @brief Limit a coupled deposition/carrier withdrawal by final solid packing.
    !> @details Erosion is already committed. Deposition removes both mobile
    !> solids and its associated pore carrier; other_loss is independently
    !> limited and therefore included in the fixed depth budget.
    !> @param[in] h Initial mobile thickness [m].
    !> @param[in] alpha Initial total solid volume fraction.
    !> @param[in] solid_erosion Solid-volume influx from erosion [m/s].
    !> @param[in] carrier_erosion Carrier-volume influx from erosion [m/s].
    !> @param[in] solid_deposition Solid-volume sink proposed for deposition [m/s].
    !> @param[in] carrier_deposition Carrier-volume sink proposed by deposition [m/s].
    !> @param[in] shared_loss Additional carrier sink sharing the deposition factor [m/s].
    !> @param[in] other_loss Separately limited carrier sink [m/s].
    !> @param[in] packing Maximum admissible solid volume fraction.
    !> @param[in] dt Physical exchange step duration [s].
    !> @return Maximum permitted common scaling factor in [0,1].
    PURE REAL(wp) FUNCTION packing_withdrawal_factor(h,alpha,solid_erosion, &
              carrier_erosion,solid_deposition,carrier_deposition,shared_loss, &
              other_loss,packing,dt) RESULT(f)
        REAL(wp), INTENT(IN) :: h,alpha,solid_erosion,carrier_erosion
        REAL(wp), INTENT(IN) :: solid_deposition,carrier_deposition,shared_loss
        REAL(wp), INTENT(IN) :: other_loss,packing,dt
        REAL(wp) :: margin,coefficient

        ! Enforce S_final <= packing * H_final without assuming that
        ! deposition leaves the solid volume unchanged.
        margin=packing*(h+dt*(solid_erosion+carrier_erosion-other_loss)) &
                  -(h*alpha+dt*solid_erosion)
        coefficient=packing*(solid_deposition+carrier_deposition+shared_loss) &
                  -solid_deposition
        f=1.0_wp

        IF (coefficient>0.0_wp) THEN

            f=MIN(1.0_wp,MAX(0.0_wp,margin)/(dt*coefficient))
        END IF

    END FUNCTION packing_withdrawal_factor

    !> @brief Add one process contribution to conserved equations and diagnostic rate totals.
    !> @details Sink momentum is advected with the local horizontal velocity.
    !> @param[in] p Rates produced by one mass-exchange process.
    !> @param[in] u Horizontal x velocity of the mobile material [m/s].
    !> @param[in] v Horizontal y velocity of the mobile material [m/s].
    !> @param[in,out] total Accumulated mass-exchange rates for the cell.
    !> @param[in,out] eqns_term Conservative mass, momentum, energy and species source vector.
    SUBROUTINE accumulate_process(p,u,v,total,eqns_term)
        TYPE(mass_exchange_rates_type), INTENT(IN) :: p
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: total
        REAL(wp), INTENT(IN) :: u,v
        REAL(wp), INTENT(INOUT) :: eqns_term(n_eqns)
        eqns_term(1)=eqns_term(1)+p%mixture_mass_rate
        eqns_term(2)=eqns_term(2)+u*MIN(0.0_wp,p%mixture_mass_rate)
        eqns_term(3)=eqns_term(3)+v*MIN(0.0_wp,p%mixture_mass_rate)
        eqns_term(4)=eqns_term(4)+p%energy_rate
        eqns_term(idx_solidEqn_first:idx_solidEqn_last)= &
                  eqns_term(idx_solidEqn_first:idx_solidEqn_last)+p%solid_mass_rate
        IF (gas_flag .AND. liquid_flag) eqns_term(n_eqns)=eqns_term(n_eqns)+p%liquid_mass_rate
        IF (n_add_gas>0) eqns_term(idx_addGasEqn_first:idx_addGasEqn_last)= &
                  eqns_term(idx_addGasEqn_first:idx_addGasEqn_last)+p%add_gas_mass_rate
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
        total%gas_withdrawal_volume_rate= &
                  total%gas_withdrawal_volume_rate+p%gas_withdrawal_volume_rate
        total%liquid_withdrawal_volume_rate= &
                  total%liquid_withdrawal_volume_rate+p%liquid_withdrawal_volume_rate
    END SUBROUTINE accumulate_process

    !> @brief Apply a common resource-limiting factor to every contribution from one process.
    !> @param[in,out] p Rates produced by one mass-exchange process.
    !> @param[in] factor Common dimensionless reduction factor between zero and one.
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
        p%gas_withdrawal_volume_rate=p%gas_withdrawal_volume_rate*factor
        p%liquid_withdrawal_volume_rate=p%liquid_withdrawal_volume_rate*factor
    END SUBROUTINE scale_process

END MODULE mass_exchange_2d
