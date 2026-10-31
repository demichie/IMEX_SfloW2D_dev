!> @brief mass exchange 2d: physical mass-exchange components and rate evaluation.
!> @details Explicit physics is evaluated locally; unsupported treatment entry points
!> fail explicitly. The coordinator alone enforces shared carrier inventories.
MODULE mass_exchange_2d
    USE parameters_2d, ONLY : idx_totMassEqn, idx_uEqn, idx_vEqn, idx_engyEqn
    USE parameters_2d, ONLY : wp, n_vars, n_eqns, n_solid, n_add_gas,                              &
              idx_stochEqn, idx_poreEqn,                                                           &
              idx_solidEqn_first, idx_solidEqn_last,                                               &
              idx_addGasEqn_first, idx_addGasEqn_last, gas_flag, liquid_flag,                      &
              pore_pressure_flag, stoch_transport_flag, mass_exchange_rates_type,                  &
              reset_mass_exchange_rates
    USE diagnostics_2d, ONLY : fatal_error
    USE constitutive_parameters_2d, ONLY : maximum_solid_packing
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
    USE mass_exchange_erosion_2d, ONLY : eval_mass_exchange_erosion_explicit
    USE mass_exchange_deposition_2d, ONLY : eval_mass_exchange_deposition_explicit
    USE mass_exchange_carrier_liquid_loss_2d, ONLY : eval_mass_exchange_carrier_liquid_loss_explicit
    USE mass_exchange_carrier_darcy_loss_2d, ONLY : eval_mass_exchange_carrier_darcy_loss_explicit
    USE mass_exchange_entrainment_2d, ONLY : eval_mass_exchange_entrainment_explicit
    USE mass_exchange_vaporization_2d, ONLY : eval_mass_exchange_vaporization_explicit
    USE mass_exchange_lava_deposition_2d, ONLY : lava_deposition_flag,                             &
        eval_mass_exchange_lava_deposition_explicit

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: eval_mass_exchange_explicit, eval_mass_exchange_semi_implicit
    PUBLIC :: eval_mass_exchange_implicit, solve_mass_exchange_analytic

CONTAINS

    SUBROUTINE eval_mass_exchange_explicit(physical_state,B_zone,erodible,dt,total_rates,eqns_term)

        TYPE(physical_state_type), INTENT(IN) :: physical_state
        REAL(wp), INTENT(IN) :: erodible(n_solid),dt
        INTEGER, INTENT(IN) :: B_zone
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: total_rates
        REAL(wp), INTENT(OUT) :: eqns_term(n_eqns)
        TYPE(mass_exchange_rates_type) :: process_rates,deposition_rates,liquid_loss_rates,        &
            darcy_loss_rates
        REAL(wp) :: alpha
        REAL(wp) :: carrier_available,carrier_factor,carrier_volume,carrier_loss
        REAL(wp) :: carrier_withdrawals(2),solid_erosion,carrier_erosion
        REAL(wp) :: solid_deposition,packing_factor,erosion_factor

        CALL reset_mass_exchange_rates(total_rates)
        eqns_term=0.0_wp
        IF (physical_state%h<=EPSILON(1.0_wp)) RETURN

        ! Processes requiring a gas/liquid carrier are evaluated only when that
        ! physical phase exists. Lava deposition is dispatched independently below.
        IF (gas_flag .OR. liquid_flag) THEN
            alpha=SUM(physical_state%alphas)
            CALL eval_mass_exchange_erosion_explicit(physical_state,erodible, &
                      dt,process_rates)
            ! Erosion can increase packing by importing a solid-rich substrate.
            ! Limit it before accumulating any conservative or bed contribution.
            erosion_factor=erosion_packing_factor(physical_state%h,alpha,                                         &
                      SUM(MAX(0.0_wp,process_rates%solid_erosion_rate)),                           &
                      MAX(0.0_wp,process_rates%carrier_erosion_rate),                              &
                      maximum_solid_packing,dt)
            CALL scale_process(process_rates,erosion_factor)
            CALL accumulate_process(process_rates,physical_state,total_rates,eqns_term)
            CALL eval_mass_exchange_deposition_explicit(physical_state,dt,deposition_rates)
            CALL eval_mass_exchange_carrier_liquid_loss_explicit(physical_state,liquid_loss_rates)
            CALL eval_mass_exchange_carrier_darcy_loss_explicit(physical_state,darcy_loss_rates)

            ! A single carrier inventory is shared by pore filling and its
            ! physical loss law. Gas and liquid have distinct laws, but the
            ! admissibility and maximum-packing limits are identical.
            carrier_volume=MAX(0.0_wp,deposition_rates%carrier_deposition_rate)
            solid_deposition=SUM(MAX(0.0_wp,deposition_rates%solid_deposition_rate))
            solid_erosion=SUM(MAX(0.0_wp,process_rates%solid_erosion_rate))
            carrier_erosion=MAX(0.0_wp,process_rates%carrier_erosion_rate)

            ! The process-specific modules determine which carrier is removed.
            ! Mutually exclusive carrier flags ensure that only one loss is active.
            carrier_loss=0.0_wp

            IF (gas_flag) THEN

                carrier_loss=MAX(0.0_wp,darcy_loss_rates%carrier_withdrawal_volume_rate)

            ELSE

                carrier_loss=MAX(0.0_wp,liquid_loss_rates%carrier_withdrawal_volume_rate)
            END IF

            ! The carrier occupies the volume not filled by dispersed solids.
            ! Limit gross withdrawals independently of positive erosion inputs.
            carrier_available=physical_state%h*MAX(0.0_wp,1.0_wp-alpha)/dt
            carrier_withdrawals=[-carrier_volume,-carrier_loss]
            carrier_factor=shared_withdrawal_factor(carrier_withdrawals,carrier_available)

            ! The packing constraint acts on the coupled changes in solid and
            ! total volumes, including the solid removed by deposition.
            packing_factor=packing_withdrawal_factor(physical_state%h,alpha,solid_erosion,                        &
                      carrier_erosion,solid_deposition,carrier_volume,                             &
                      carrier_loss,0.0_wp,maximum_solid_packing,dt)
            carrier_factor=MIN(carrier_factor,packing_factor)

            CALL scale_process(deposition_rates,carrier_factor)

            IF (gas_flag) THEN

                CALL scale_process(darcy_loss_rates,carrier_factor)

            ELSE

                CALL scale_process(liquid_loss_rates,carrier_factor)
            END IF

            ! Assemble only the rates after all common-inventory limiters have acted.
            CALL accumulate_process(deposition_rates,physical_state,total_rates,eqns_term)
            CALL accumulate_process(liquid_loss_rates,physical_state,total_rates,eqns_term)
            CALL accumulate_process(darcy_loss_rates,physical_state,total_rates,eqns_term)
            CALL eval_mass_exchange_entrainment_explicit(physical_state,process_rates)
            CALL accumulate_process(process_rates,physical_state,total_rates,eqns_term)
            CALL eval_mass_exchange_vaporization_explicit(physical_state,B_zone, &
                deposition_rates%solid_deposition_rate,process_rates)
            CALL accumulate_process(process_rates,physical_state,total_rates,eqns_term)

            ! Scalar transport and pore-pressure sources inherit the net mass exchange.
            IF (stoch_transport_flag) eqns_term(idx_stochEqn)=eqns_term(idx_totMassEqn)*physical_state%Zs
            IF (pore_pressure_flag) eqns_term(idx_poreEqn)=                                        &
                      eqns_term(idx_totMassEqn)*physical_state%exc_pore_pres
        END IF  ! gas_flag .OR. liquid_flag

        ! Lava deposition follows the same process_rates-rate accumulation contract.
        IF (lava_deposition_flag) THEN
            CALL eval_mass_exchange_lava_deposition_explicit(physical_state,dt,process_rates)
            CALL accumulate_process(process_rates,physical_state,total_rates,eqns_term)
        END IF

    END SUBROUTINE eval_mass_exchange_explicit


    !> This family currently supports only explicit coordination.
    SUBROUTINE eval_mass_exchange_semi_implicit(state,total_rates,eqns_term)
        TYPE(physical_state_type), INTENT(IN) :: state
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: total_rates
        REAL(wp), INTENT(OUT) :: eqns_term(n_eqns)
        CALL reset_mass_exchange_rates(total_rates)
        eqns_term=0.0_wp
        CALL fatal_error('mass_exchange: semi-implicit coordination is unsupported')
    END SUBROUTINE eval_mass_exchange_semi_implicit

    SUBROUTINE eval_mass_exchange_implicit(state,total_rates,eqns_term)
        TYPE(complex_physical_state_type), INTENT(IN) :: state
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: total_rates
        COMPLEX(wp), INTENT(OUT) :: eqns_term(n_eqns)
        CALL reset_mass_exchange_rates(total_rates)
        eqns_term=CMPLX(0.0_wp,0.0_wp,wp)
        CALL fatal_error('mass_exchange: implicit coordination is unsupported')
    END SUBROUTINE eval_mass_exchange_implicit

    SUBROUTINE solve_mass_exchange_analytic(state,total_rates,eqns_term)
        TYPE(physical_state_type), INTENT(IN) :: state
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: total_rates
        REAL(wp), INTENT(OUT) :: eqns_term(n_eqns)
        CALL reset_mass_exchange_rates(total_rates)
        eqns_term=0.0_wp
        CALL fatal_error('mass_exchange: analytic coordination is unsupported')
    END SUBROUTINE solve_mass_exchange_analytic


    !> @brief Coordinate cell-local mass exchange, shared phase budgets and conservative source
    !> assembly.
    !> @details A single gas or liquid carrier inventory is constrained before conservative
    !> assembly.
    !> Positive sources never cancel gross withdrawals in a shared inventory.
    !> @param[in] qp Local primitive flow state; component order follows model configuration.
    !> @param[in] B_zone Substrate zone identifier controlling vaporization eligibility.
    !> @param[in] erodible Available erodible solid inventory by species [m].
    !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
    !> @param[in,out] total_rates Accumulated process rates and bed changes for this cell.
    !> @param[out] eqns_term Conservative mass, momentum, energy and species source vector.


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
        IF (total_withdrawal>0.0_wp)                                                               &
                  f=MIN(1.0_wp,MAX(0.0_wp,available)/total_withdrawal)
    END FUNCTION shared_withdrawal_factor

    !> @brief Bound solid-rich substrate erosion by the final packing limit.
    !> @details The eroded carrier and solid volumes are scaled together.
    !> @param[in] h Initial mobile thickness [m].
    !> @param[in] alpha Initial solid volume fraction.
    !> @param[in] solid_erosion Proposed solid erosion flux [m/s].
    !> @param[in] carrier_erosion Proposed carrier erosion flux [m/s].
    !> @param[in] packing Maximum admissible solid volume fraction.
    !> @param[in] dt Physical timestep [s].
    !> @return Common erosion factor within [0,1].
    PURE REAL(wp) FUNCTION erosion_packing_factor(h,alpha,solid_erosion,                           &
              carrier_erosion,packing,dt) RESULT(f)
        REAL(wp), INTENT(IN) :: h,alpha,solid_erosion,carrier_erosion,packing,dt
        REAL(wp) :: deficit,margin

        f=1.0_wp
        IF (dt<=0.0_wp) RETURN
        deficit=(1.0_wp-packing)*solid_erosion-packing*carrier_erosion
        margin=MAX(0.0_wp,h*(packing-alpha))
        IF (deficit>0.0_wp) f=MIN(1.0_wp,margin/(dt*deficit))
    END FUNCTION erosion_packing_factor

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
    PURE REAL(wp) FUNCTION packing_withdrawal_factor(h,alpha,solid_erosion,                        &
              carrier_erosion,solid_deposition,carrier_deposition,shared_loss,                     &
              other_loss,packing,dt) RESULT(f)
        REAL(wp), INTENT(IN) :: h,alpha,solid_erosion,carrier_erosion
        REAL(wp), INTENT(IN) :: solid_deposition,carrier_deposition,shared_loss
        REAL(wp), INTENT(IN) :: other_loss,packing,dt
        REAL(wp) :: margin,coefficient

        ! Enforce S_final <= packing * H_final without assuming that
        ! deposition leaves the solid volume unchanged.
        margin=packing*(h+dt*(solid_erosion+carrier_erosion-other_loss))                           &
                  -(h*alpha+dt*solid_erosion)
        coefficient=packing*(solid_deposition+carrier_deposition+shared_loss)                      &
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
    SUBROUTINE accumulate_process(p,state,total,eqns_term)

        TYPE(mass_exchange_rates_type), INTENT(IN) :: p
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: total
        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(INOUT) :: eqns_term(n_eqns)

        eqns_term(idx_totMassEqn)=eqns_term(idx_totMassEqn)+p%mixture_mass_rate
        eqns_term(idx_uEqn)=eqns_term(idx_uEqn)+state%u*MIN(0.0_wp,p%mixture_mass_rate)
        eqns_term(idx_vEqn)=eqns_term(idx_vEqn)+state%v*MIN(0.0_wp,p%mixture_mass_rate)
        eqns_term(idx_engyEqn)=eqns_term(idx_engyEqn)+p%energy_rate
        eqns_term(idx_solidEqn_first:idx_solidEqn_last)=                                           &
                  eqns_term(idx_solidEqn_first:idx_solidEqn_last)+p%solid_mass_rate
        IF (n_add_gas>0) eqns_term(idx_addGasEqn_first:idx_addGasEqn_last)=                        &
                  eqns_term(idx_addGasEqn_first:idx_addGasEqn_last)+p%add_gas_mass_rate
        total%mixture_mass_rate=total%mixture_mass_rate+p%mixture_mass_rate
        total%energy_rate=total%energy_rate+p%energy_rate
        total%solid_mass_rate=total%solid_mass_rate+p%solid_mass_rate
        total%add_gas_mass_rate=total%add_gas_mass_rate+p%add_gas_mass_rate
        total%solid_erosion_rate=total%solid_erosion_rate+p%solid_erosion_rate
        total%solid_deposition_rate=total%solid_deposition_rate+p%solid_deposition_rate
        total%carrier_erosion_rate=                                                                &
                  total%carrier_erosion_rate+p%carrier_erosion_rate
        total%carrier_deposition_rate=                                                             &
                  total%carrier_deposition_rate+p%carrier_deposition_rate
        total%dB_dt=total%dB_dt+p%dB_dt
        total%lava_deposit_rate=total%lava_deposit_rate+p%lava_deposit_rate
        total%carrier_withdrawal_volume_rate=                                                      &
                  total%carrier_withdrawal_volume_rate+p%carrier_withdrawal_volume_rate
    END SUBROUTINE accumulate_process

    !> @brief Apply a common resource-limiting factor to every contribution from one process.
    !> @param[in,out] p Rates produced by one mass-exchange process.
    !> @param[in] factor Common dimensionless reduction factor between zero and one.
    SUBROUTINE scale_process(p,factor)

        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: p
        REAL(wp), INTENT(IN) :: factor

        p%mixture_mass_rate=p%mixture_mass_rate*factor
        p%energy_rate=p%energy_rate*factor
        p%solid_mass_rate=p%solid_mass_rate*factor
        p%add_gas_mass_rate=p%add_gas_mass_rate*factor
        p%solid_erosion_rate=p%solid_erosion_rate*factor
        p%solid_deposition_rate=p%solid_deposition_rate*factor
        p%carrier_erosion_rate=p%carrier_erosion_rate*factor
        p%carrier_deposition_rate=p%carrier_deposition_rate*factor
        p%dB_dt=p%dB_dt*factor
        p%lava_deposit_rate=p%lava_deposit_rate*factor
        p%carrier_withdrawal_volume_rate=p%carrier_withdrawal_volume_rate*factor
    END SUBROUTINE scale_process

END MODULE mass_exchange_2d
