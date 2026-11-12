!> \brief Common physical source assembly entry points.
!> \details Assembles physical source rates and their coupled limiters.
!> Spatial finite-volume operators belong exclusively to spatial modules.
!> The Lie executor applies the resulting limited rates and publishes
!> fluid/substrate candidates at group boundaries.
MODULE source_terms_2d
    USE parameters_2d, ONLY : idx_h
    USE mass_injection_2d, ONLY : eval_mass_injection_explicit
    USE mass_injection_coefficients_2d, ONLY : injection_material_set_type
    USE mass_injection_coefficients_2d
    USE parameters_2d
    USE constitutive_parameters_2d
    USE state_types, ONLY : physical_state_real, physical_state_complex,                                           &
        initialize_physical_state_real, initialize_physical_state_complex
    USE state_physical, ONLY : build_physical_state_from_q, build_physical_state_from_qp, evaluate_temperature_from_q
    USE state_mixture, ONLY : evaluate_mixture_heat_capacity
    USE thermal_laws_2d, ONLY : eval_thermal_explicit, eval_thermal_implicit,                                          &
            eval_thermal_implicit_real, eval_thermal_analytic, eval_thermal_implicit_temperature
    USE pore_pressure_2d, ONLY : eval_pore_pressure_analytic
    USE mass_exchange_2d, ONLY : eval_mass_exchange_explicit,mass_exchange_scratch_type
    USE rheology_parameters_2d, ONLY : rheology_l1
    USE rheology_2d, ONLY : curvature_acceleration, eval_rheology_analytic,                                            &
            eval_rheology_implicit_id, eval_rheology_semi_implicit_id,                                                 &
            eval_rheology_implicit_real_id, eval_rheology_explicit
    USE diagnostics_2d, ONLY : fatal_error
    USE state_composition, ONLY : conservative_mass_admissible
    USE source_catalog_2d, ONLY : SOURCE_METHOD_NEWTON
    USE source_catalog_2d, ONLY : SOURCE_METHOD_EXPLICIT, SOURCE_METHOD_ANALYTIC
    USE source_registry_2d, ONLY : source_group_type
    USE source_catalog_2d, ONLY : SOURCE_ID_GEOMETRY_CURVATURE, SOURCE_ID_BOTTOM_CIRCULAR_INJECTION,                   &
        SOURCE_ID_BOTTOM_FISSURAL_INJECTION, SOURCE_ID_RADIAL_INJECTION, SOURCE_ID_THERMAL_COOLING
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_source_explicit
    PUBLIC :: limit_explicit_source_rates,assemble_lie_explicit_transaction
    ! Explicit physical-source evaluation provides distinct cell/geometry
    ! and mass-exchange signatures, independent of time-integration stages.
    INTERFACE eval_source_explicit
        MODULE PROCEDURE eval_geometry_cell_source_explicit, eval_mass_exchange_source_explicit
    END INTERFACE
    PUBLIC :: eval_source_semi_implicit, eval_source_implicit, eval_source_implicit_real, eval_source_analytic
    PUBLIC :: source_cell_context_type
    PUBLIC :: eval_geometry_curvature_explicit
    PUBLIC :: bottom_coeff_step, bottom_vel_step, lateral_coeff_step, fissure_coeff_step
    ! Cell-local geometry and stochastic context for analytic source solves.
    TYPE :: source_cell_context_type
        REAL(wp) :: Bprime_x=0.0_wp, Bprime_y=0.0_wp
        REAL(wp) :: Bsecond_xx=0.0_wp, Bsecond_xy=0.0_wp, Bsecond_yy=0.0_wp
        REAL(wp) :: grav_coeff=0.0_wp, Zij=0.0_wp
    END TYPE source_cell_context_type
CONTAINS


    !> Apply a uniform physical factor to all coupled exchange rates.
    SUBROUTINE scale_mass_exchange_transfer(transfer,scale)
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: transfer
        REAL(wp), INTENT(IN) :: scale
        transfer%solid_deposition_rate=scale*transfer%solid_deposition_rate
        transfer%solid_erosion_rate=scale*transfer%solid_erosion_rate
        transfer%solid_mass_rate=scale*transfer%solid_mass_rate
        transfer%add_gas_mass_rate=scale*transfer%add_gas_mass_rate
        transfer%mixture_mass_rate=scale*transfer%mixture_mass_rate
        transfer%energy_rate=scale*transfer%energy_rate
        transfer%carrier_erosion_rate=scale*transfer%carrier_erosion_rate
        transfer%carrier_deposition_rate=scale*transfer%carrier_deposition_rate
        transfer%carrier_withdrawal_volume_rate=scale*transfer%carrier_withdrawal_volume_rate
        transfer%lava_deposit_rate=scale*transfer%lava_deposit_rate
        transfer%dB_dt=scale*transfer%dB_dt
    END SUBROUTINE scale_mass_exchange_transfer

    !> Assemble and limit a complete explicit Lie group, without committing state.
    SUBROUTINE assemble_lie_explicit_transaction(qp_eval,q_recipient,dt,t,group,                                       &
            inventory,rate,transfer,has_exchange,valid,stochastic_value,                                               &
            B_zone_cell,context,cell_fract,cell_fissure_fract,arc_perim,arc_n_x,arc_n_y,cell_area,                     &
            non_exchange_group,exchange_selected,physical_scratch,exchange_scratch,q_preexisting,materials)
        INTEGER, INTENT(IN) :: B_zone_cell
        TYPE(source_cell_context_type), INTENT(IN) :: context
        REAL(wp), INTENT(IN) :: cell_fract,cell_fissure_fract(:)
        REAL(wp), INTENT(IN) :: arc_perim,arc_n_x,arc_n_y,cell_area
        REAL(wp), INTENT(IN) :: qp_eval(n_vars),q_recipient(n_eqns),dt,t
        TYPE(source_group_type), INTENT(IN) :: group
        ! Prepartitioned immutable metadata supplied by the Lie executor.
        TYPE(source_group_type), INTENT(IN) :: non_exchange_group
        LOGICAL, INTENT(IN) :: exchange_selected
        TYPE(physical_state_real), INTENT(INOUT), OPTIONAL :: physical_scratch
        TYPE(mass_exchange_scratch_type), INTENT(INOUT), OPTIONAL :: exchange_scratch
        TYPE(injection_material_set_type), INTENT(IN) :: materials
        ! Inventory at the beginning of the timestep: exchange cannot consume
        ! source mass introduced by preceding groups during this same step.
        REAL(wp), INTENT(IN), OPTIONAL :: q_preexisting(n_eqns)
        REAL(wp), INTENT(IN) :: inventory(n_solid)
        REAL(wp), INTENT(IN) :: stochastic_value
        REAL(wp), INTENT(OUT) :: rate(n_eqns)
        ! Retain allocations in the thread-owned transfer record across cells.
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: transfer
        LOGICAL, INTENT(OUT) :: has_exchange,valid
        REAL(wp) :: physical(n_eqns),exchange(n_eqns),coverage,scale
        LOGICAL :: physical_ready
        has_exchange=exchange_selected
        valid=.FALSE.
        physical=0.0_wp
        exchange=0.0_wp
        physical_ready=.FALSE.
        ! Partitioned immutable metadata comes from the Lie registry/executor.
        ! Never copy allocatable group members or PACK source names inside cells.
        IF (non_exchange_group%source_mask/=0) THEN
            CALL eval_source_explicit(context%Bprime_x,context%Bprime_y,                                               &
                context%Bsecond_xx,context%Bsecond_xy,context%Bsecond_yy,                                              &
                context%grav_coeff,qp_eval,physical,t,cell_fract,cell_fissure_fract,                                   &
                arc_perim,arc_n_x,arc_n_y,cell_area,non_exchange_group,stochastic_value,physical_scratch,              &
                physical_ready,materials)
        END IF

        IF (has_exchange) THEN
            CALL eval_source_explicit(qp_eval,B_zone_cell,                                                             &
                context%Bprime_x,context%Bprime_y,inventory,dt,                                                        &
                group%source_mask,transfer,exchange,physical_scratch,exchange_scratch,physical_ready)

            IF (bottom_radial_source_flag .OR. bottom_fissural_source_flag) THEN
                coverage=1.0_wp-cell_fract
                IF (coverage<0.0_wp .OR. coverage>1.0_wp) RETURN
                exchange=coverage*exchange
                CALL scale_mass_exchange_transfer(transfer,coverage)

            END IF

        END IF

        rate=physical+exchange

        IF (has_exchange) THEN
            CALL limit_explicit_source_rates(q_recipient,inventory,dt,rate,                                            &
                transfer%solid_erosion_rate,scale,valid,q_preexisting)
        ELSE
            CALL limit_explicit_source_rates(q_recipient,[REAL(wp)::],dt,rate,                                         &
                [REAL(wp)::],scale,valid)

        END IF

        IF (.NOT.valid) RETURN
        rate=scale*rate
        IF (has_exchange) CALL scale_mass_exchange_transfer(transfer,scale)
    END SUBROUTINE assemble_lie_explicit_transaction

    !> Scale a complete cell-local physical source transaction only after
    !> assembling all explicit source rates. Reference inventories belong to
    !> the post-spatial Lie state; group zero rate evaluation stays at q^n.
    SUBROUTINE limit_explicit_source_rates(q_reference,inventory,dt,rate,                                              &
            solid_erosion_rate,limiter,valid,q_preexisting)
        REAL(wp), INTENT(IN) :: q_reference(n_eqns),inventory(:),dt
        REAL(wp), INTENT(IN) :: rate(n_eqns),solid_erosion_rate(:)
        REAL(wp), INTENT(IN), OPTIONAL :: q_preexisting(n_eqns)
        REAL(wp), INTENT(OUT) :: limiter
        LOGICAL, INTENT(OUT) :: valid
        REAL(wp) :: available_carrier,carrier_consumption
        REAL(wp) :: subnormal_reserve,withdrawal
        REAL(wp) :: candidate(n_eqns)
        INTEGER :: solid_eq,i,gas_eq
        limiter=1.0_wp
        valid=.FALSE.
        IF (dt<=0.0_wp) RETURN
        IF (SIZE(inventory)/=SIZE(solid_erosion_rate)) RETURN
        IF (rate(idx_totMassEqn)<0.0_wp) limiter=MIN(limiter,                                                          &
            MAX(0.0_wp,q_reference(idx_totMassEqn)/(-dt*rate(idx_totMassEqn))))
        ! Additional gases belong to the gas carrier, whose withdrawal is
        ! limited as a whole. Their deposition rates follow gas composition.
        DO solid_eq=idx_solidEqn_first,idx_solidEqn_last
            IF (rate(solid_eq)<0.0_wp) limiter=MIN(limiter,                                                            &
                MAX(0.0_wp,q_reference(solid_eq)/(-dt*rate(solid_eq))))

        END DO

        ! Gas species also have separate nonnegative conservative budgets.
        DO gas_eq=idx_addGasEqn_first,idx_addGasEqn_last
            IF (rate(gas_eq)<0.0_wp) limiter=MIN(limiter, &
                MAX(0.0_wp,q_reference(gas_eq)/(-dt*rate(gas_eq))))
        END DO

        DO i=1,SIZE(inventory)
            IF (solid_erosion_rate(i)>0.0_wp) limiter=MIN(limiter,                                                     &
                MAX(0.0_wp,inventory(i)/(dt*solid_erosion_rate(i))))

        END DO

        available_carrier=q_reference(idx_totMassEqn)-                                                                 &
            SUM(q_reference(idx_solidEqn_first:idx_solidEqn_last))
        carrier_consumption=SUM(rate(idx_solidEqn_first:idx_solidEqn_last))-rate(idx_totMassEqn)
        IF (carrier_consumption>0.0_wp) limiter=MIN(limiter,                                                           &
            MAX(0.0_wp,available_carrier/(dt*carrier_consumption)))
        ! Preserve the historical source-splitting withdrawal inventory,
        ! but verify the resulting rate against the actual updated recipient.
        ! This prevents settling from consuming newly injected solids within
        ! the same step while avoiding the old stale-recipient D19 failure.
        IF (PRESENT(q_preexisting)) THEN
            IF (rate(idx_totMassEqn)<0.0_wp) limiter=MIN(limiter, &
                MAX(0.0_wp,q_preexisting(idx_totMassEqn)/(-dt*rate(idx_totMassEqn))))
            DO solid_eq=idx_solidEqn_first,idx_solidEqn_last
                IF (rate(solid_eq)<0.0_wp) limiter=MIN(limiter, &
                    MAX(0.0_wp,q_preexisting(solid_eq)/(-dt*rate(solid_eq))))
            END DO
            DO gas_eq=idx_addGasEqn_first,idx_addGasEqn_last
                IF (rate(gas_eq)<0.0_wp) limiter=MIN(limiter, &
                    MAX(0.0_wp,q_preexisting(gas_eq)/(-dt*rate(gas_eq))))
            END DO
            available_carrier=q_preexisting(idx_totMassEqn)- &
                SUM(q_preexisting(idx_solidEqn_first:idx_solidEqn_last))
            IF (carrier_consumption>0.0_wp) limiter=MIN(limiter, &
                MAX(0.0_wp,available_carrier/(dt*carrier_consumption)))
        END IF
        ! When one or more species balances are saturated, reserve a small
        ! relative buffer before the floating-point candidate check. A bound
        ! with zero remaining species is susceptible to reassociation/FMA
        ! differences between this routine and the caller (observed in TAAL
        ! as -1.4293E-20 kg/m2 after the rate was reported admissible).
        ! This is a numerical safety margin on source rates, not clipping q.
        IF (limiter<1.0_wp .AND. limiter>0.0_wp) &
            limiter=limiter*(1.0_wp-4096.0_wp*EPSILON(1.0_wp))
        ! For subnormal inventories, the relative margin above can round
        ! entirely away: 4096*epsilon times ~1E-313 is below the smallest
        ! representable positive real. Reserve a few representable quanta
        ! of a depleted species before computing its withdrawal bound.
        ! This adjusts only the already-existing explicit source limiter,
        ! never the candidate state or the strict admissibility predicate.
        DO i=1,n_eqns
            IF (i/=idx_totMassEqn .AND.                                                              &
                (i<idx_solidEqn_first .OR. i>idx_solidEqn_last) .AND.                                &
                (i<idx_addGasEqn_first .OR. i>idx_addGasEqn_last)) CYCLE
            IF (q_reference(i)<=0.0_wp .OR. q_reference(i)>=TINY(1.0_wp)) CYCLE
            IF (rate(i)>=0.0_wp) CYCLE
            withdrawal=-dt*rate(i)
            IF (withdrawal<=0.0_wp) CYCLE
            subnormal_reserve=MIN(q_reference(i),32.0_wp*SPACING(q_reference(i)))
            limiter=MIN(limiter,MAX(0.0_wp,(q_reference(i)-subnormal_reserve)/withdrawal))
        END DO
        ! Verify the exact floating-point update applied by the Lie executor:
        ! q_candidate=q_reference+dt*(limiter*rate). The analytical ratio can
        ! round up to 1, or leave a negative component at machine precision.
        ! Move the physical rate limiter toward zero by representable steps,
        ! never modifying the conservative candidate or its admissibility test.
        IF (limiter>=0.0_wp .AND. limiter<=1.0_wp) THEN
            DO i=1,256
                candidate=q_reference+dt*(limiter*rate)
                ! A wet conservative state requires strictly positive energy.
                ! Complete carrier depletion can round the energy to exactly
                ! zero while leaving a representable wet mass (Etna lava D19).
                ! Back off the *coupled rate*, never clip q or invent energy.
                IF (conservative_mass_admissible(candidate)) THEN
                    IF (q_reference(idx_engyEqn)<=0.0_wp .OR. rate(idx_engyEqn)>=0.0_wp) EXIT
                    IF (candidate(idx_totMassEqn)<=EPSILON(1.0_wp)) EXIT
                    IF (candidate(idx_engyEqn)>0.0_wp) EXIT
                END IF
                IF (limiter<=0.0_wp) EXIT
                limiter=NEAREST(limiter,-1.0_wp)
            END DO
            candidate=q_reference+dt*(limiter*rate)
        END IF
        valid=limiter>=0.0_wp .AND. limiter<=1.0_wp
        IF (valid) valid=conservative_mass_admissible(candidate)
        IF (valid .AND. q_reference(idx_engyEqn)>0.0_wp .AND. rate(idx_engyEqn)<0.0_wp) THEN
            IF (candidate(idx_totMassEqn)>EPSILON(1.0_wp)) &
                valid=candidate(idx_engyEqn)>0.0_wp
        END IF
    END SUBROUTINE limit_explicit_source_rates

    !> Return the conservative equivalent rate of the exact source substep.
    !> The interval is supplied by the selected Lie group operator.
    SUBROUTINE eval_source_analytic(group,q_before,context,interval,rate,status,physical_scratch)
        TYPE(source_group_type), INTENT(IN) :: group
        REAL(wp), INTENT(IN) :: q_before(n_vars),interval
        TYPE(source_cell_context_type), INTENT(IN) :: context
        REAL(wp), INTENT(OUT) :: rate(n_eqns)
        INTEGER, INTENT(OUT) :: status
        TYPE(physical_state_real), OPTIONAL, TARGET, INTENT(INOUT) :: physical_scratch
        REAL(wp) :: q_after(n_vars)
        rate=0.0_wp

        IF (interval<0.0_wp) THEN
        status=-8
        RETURN

        END IF

        IF (interval==0.0_wp) THEN
        status=0
        RETURN

        END IF

        CALL eval_source_analytic_state(group,q_before,context,interval,q_after,status,physical_scratch)
        IF (status/=0) RETURN
        rate(1:n_vars)=(q_after-q_before)/interval
    END SUBROUTINE eval_source_analytic

    SUBROUTINE eval_source_analytic_state(group,q_before,context,interval,q_after,status,physical_scratch)
        TYPE(source_group_type), INTENT(IN) :: group
        REAL(wp), INTENT(IN) :: q_before(n_vars),interval
        TYPE(source_cell_context_type), INTENT(IN) :: context
        REAL(wp), INTENT(OUT) :: q_after(n_vars)
        INTEGER, INTENT(OUT) :: status
        TYPE(physical_state_real), OPTIONAL, TARGET, INTENT(INOUT) :: physical_scratch
        TYPE(physical_state_real), TARGET :: local_state
        TYPE(physical_state_real), POINTER :: physical_state
        REAL(wp) :: horizontal_speed,speed_ratio,curvature_accel,velocity_after(2)
        REAL(wp) :: heat_capacity,temperature_after,temperature_before,cp_m,cp_c
        REAL(wp) :: xs(n_solid),xg(n_add_gas)
        LOGICAL :: needs_physical_state
        REAL(wp) :: pore_rate(n_eqns)
        LOGICAL :: pore_ok
        LOGICAL :: supported
        CHARACTER(LEN=64) :: selected_model_name
        q_after=q_before
        status=-5
        IF (group%method/=SOURCE_METHOD_ANALYTIC) RETURN
        IF (.NOT.ALLOCATED(group%source_names)) RETURN
        IF (SIZE(group%source_names)/=1) RETURN
        selected_model_name=TRIM(group%source_names(1))
        ! Only rheology and positive-mass Darcy require the complete physical
        ! state. For wet analytic cooling the sole thermodynamic dependency is
        ! cp_m, which determines T from conservative energy. Keep the old dry
        ! branch to preserve its ambient-temperature and rejection semantics.
        needs_physical_state=INDEX(selected_model_name,'rheology_')==1 .OR.                                            &
            (selected_model_name=='pore_pressure_dissipation' .AND. q_before(1)>0.0_wp) .OR.                            &
            (selected_model_name=='thermal_cooling' .AND. q_before(1)<=0.0_wp)
        IF (needs_physical_state) THEN
            physical_state=>local_state
            IF (PRESENT(physical_scratch)) THEN
                physical_state=>physical_scratch
            ELSE
                CALL initialize_physical_state_real(local_state)
            END IF
            CALL build_physical_state_from_q(q_before,physical_state,                                                  &
                context%Bprime_x,context%Bprime_y,context%Bsecond_xx,                                                  &
                context%Bsecond_xy,context%Bsecond_yy)
        END IF

        IF (INDEX(selected_model_name,'rheology_')==1) THEN
            horizontal_speed=physical_state%mod_vel_hor
            speed_ratio=1.0_wp
            IF (horizontal_speed>0.0_wp)                                                                               &
                    speed_ratio=physical_state%mod_vel/horizontal_speed
            curvature_accel=physical_state%curvature_acceleration
            velocity_after=[physical_state%u,physical_state%v]

            IF (LEN_TRIM(selected_model_name)==0) THEN
                status=-7
                RETURN

            END IF

            CALL eval_rheology_analytic(rheology_l1,[selected_model_name], interval, physical_state,                               &
                    speed_ratio, context%grav_coeff, curvature_accel, pore_pressure_flag,                              &
                    stochastic_flag, context%Zij, velocity_after, supported,                                           &
                    context%Bprime_x, context%Bprime_y, context%Bsecond_xx,                                            &
                    context%Bsecond_xy, context%Bsecond_yy)

            IF (.NOT.supported) THEN
                status=-7
                RETURN

            END IF

            q_after(2:3)=q_before(1)*velocity_after
            status=0

        ELSE IF (selected_model_name=='thermal_cooling') THEN
            IF (q_before(1)>0.0_wp) THEN
                xs=q_before(idx_solid_first:idx_solid_last)/q_before(1)
                xg=q_before(idx_add_gas_first:idx_add_gas_last)/q_before(1)
                CALL evaluate_mixture_heat_capacity(xs,xg,cp_m,cp_c)
                temperature_before=q_before(4)/(q_before(1)*cp_m)
            ELSE
                temperature_before=physical_state%T
            END IF
            IF (temperature_before<=0.0_wp) THEN
                status=-3
                RETURN

            END IF

            ! Preserve the original conservative-energy division and the
            ! frozen heat capacity used by the exact cooling integrator.
            heat_capacity=q_before(4)/temperature_before
            CALL eval_thermal_analytic(['thermal_cooling'],temperature_before,interval,                               &
                    heat_capacity,radiative_term_coeff,convective_term_coeff,                                          &
                    T_ambient,temperature_after)
            q_after(4)=heat_capacity*temperature_after
            status=0

        ELSE IF (selected_model_name=='pore_pressure_dissipation') THEN
            IF (q_before(1)>0.0_wp) THEN
                CALL eval_pore_pressure_analytic(['pore_pressure_darcy'],                                              &
                        physical_state,q_before(1:n_eqns),interval,pore_rate,pore_ok)

                IF (.NOT.pore_ok) THEN
                    status=-9
                    RETURN

                END IF

                q_after(1:n_eqns)=q_before(1:n_eqns)+interval*pore_rate

            END IF

            status=0
        ELSE
            RETURN

        END IF
    END SUBROUTINE eval_source_analytic_state


    !> Evaluate the frozen semi-implicit contribution of one Newton Lie group.
    !> All inputs refer to one cell; no domain indices or shared-state writes.
    SUBROUTINE eval_source_semi_implicit(group,q_cell,Bprime_x,Bprime_y,                                               &
            grav_correction,stochastic_value,rate,physical_scratch)
        TYPE(source_group_type), INTENT(IN) :: group
        REAL(wp), INTENT(IN) :: q_cell(n_vars)
        REAL(wp), INTENT(IN) :: Bprime_x,Bprime_y,grav_correction,stochastic_value
        REAL(wp), INTENT(OUT) :: rate(n_eqns)
        TYPE(physical_state_real), TARGET, OPTIONAL, INTENT(INOUT) :: physical_scratch
        TYPE(physical_state_real), TARGET :: local_state
        TYPE(physical_state_real), POINTER :: physical_state

        rate=0.0_wp
        IF (group%method/=SOURCE_METHOD_NEWTON)                                                                        &
            CALL fatal_error('eval_source_semi_implicit requires a Newton Lie group')
        IF (.NOT.group%rheology_implicit) RETURN
        IF (q_cell(idx_totMassEqn)<=0.0_wp) RETURN

        IF (group%rheology_id==0) RETURN
        physical_state=>local_state
        IF (PRESENT(physical_scratch)) THEN
            physical_state=>physical_scratch
        ELSE
            CALL initialize_physical_state_real(local_state)
        END IF
        CALL build_physical_state_from_q(q_cell,physical_state,Bprime_x,Bprime_y,                                      &
            0.0_wp,0.0_wp,0.0_wp)
        CALL eval_rheology_semi_implicit_id(rheology_l1,group%rheology_id,physical_state,                                          &
            grav_correction,pore_pressure_flag,stochastic_flag,stochastic_value,rate(2:3))
    END SUBROUTINE eval_source_semi_implicit

    !> Evaluate the REAL implicit source directly for Newton residuals.
    !> The source formulas are identical to those used by COMPLEX Jacobians.
    SUBROUTINE eval_source_implicit_real(Bprime_x,Bprime_y,Zij,q_cell,rate,group,physical_scratch,accepted)
        REAL(wp), INTENT(IN) :: Bprime_x,Bprime_y,Zij,q_cell(n_vars)
        REAL(wp), INTENT(OUT) :: rate(n_eqns)
        TYPE(source_group_type), INTENT(IN) :: group
        TYPE(physical_state_real), OPTIONAL, TARGET, INTENT(INOUT) :: physical_scratch
        TYPE(physical_state_real), TARGET :: local_state
        TYPE(physical_state_real), POINTER :: physical_state
        REAL(wp) :: temperature
        LOGICAL, INTENT(OUT), OPTIONAL :: accepted
        LOGICAL :: state_ok

        IF (group%method/=SOURCE_METHOD_NEWTON) &
            CALL fatal_error('eval_source_implicit_real requires a Newton Lie group')
        rate=0.0_wp
        IF (PRESENT(accepted)) accepted=.TRUE.
        IF (.NOT.group%rheology_implicit .AND. .NOT.group%thermal_active) RETURN
        IF (group%thermal_active .AND. .NOT.group%rheology_implicit) THEN
            CALL evaluate_temperature_from_q(q_cell,temperature,state_ok)
            IF (.NOT.state_ok) THEN
                IF (PRESENT(accepted)) accepted=.FALSE.
                RETURN
            END IF
            CALL eval_thermal_implicit_temperature(SOURCE_ID_THERMAL_COOLING,temperature,                               &
                radiative_term_coeff,convective_term_coeff,T_ambient,rate(4))
            RETURN
        END IF
        physical_state=>local_state
        IF (PRESENT(physical_scratch)) THEN
            physical_state=>physical_scratch
        ELSE
            CALL initialize_physical_state_real(local_state)
        END IF
        CALL build_physical_state_from_q(q_cell,physical_state,Bprime_x,Bprime_y,0.0_wp,0.0_wp,0.0_wp, &
            accepted=state_ok)
        IF (.NOT.state_ok) THEN
            IF (PRESENT(accepted)) accepted=.FALSE.
            RETURN
        END IF
        IF (group%rheology_implicit) THEN
            CALL eval_rheology_implicit_real_id(rheology_l1,group%rheology_id,physical_state, &
                stochastic_flag,Zij,rate(2:3))
        END IF
        IF (group%thermal_active) THEN
            CALL eval_thermal_implicit_real(SOURCE_ID_THERMAL_COOLING,physical_state, &
                radiative_term_coeff,convective_term_coeff,T_ambient,rate(4))
        END IF
    END SUBROUTINE eval_source_implicit_real

    !> Evaluate the COMPLEX implicit rate for Newton complex-step Jacobians.
    !> The real residual uses the distinct eval_source_implicit_real contract.
    SUBROUTINE eval_source_implicit(Bprimej_x,Bprimej_y,Zij,c_qj,c_nh_term_impl,group,physical_scratch,accepted)
        REAL(wp), INTENT(IN) :: Bprimej_x,Bprimej_y,Zij
        COMPLEX(wp), INTENT(IN) :: c_qj(n_vars)
        COMPLEX(wp), INTENT(OUT) :: c_nh_term_impl(n_eqns)
        TYPE(source_group_type), INTENT(IN) :: group
        TYPE(physical_state_complex), OPTIONAL, TARGET, INTENT(INOUT) :: physical_scratch
        TYPE(physical_state_complex), TARGET :: local_state
        TYPE(physical_state_complex), POINTER :: physical_state
        COMPLEX(wp) :: temperature
        LOGICAL, INTENT(OUT), OPTIONAL :: accepted
        LOGICAL :: state_ok

        IF (group%method/=SOURCE_METHOD_NEWTON) &
            CALL fatal_error('eval_source_implicit requires a Newton Lie group')
        c_nh_term_impl=CMPLX(0.0_wp,0.0_wp,wp)
        IF (PRESENT(accepted)) accepted=.TRUE.
        IF (.NOT.group%rheology_implicit .AND. .NOT.group%thermal_active) RETURN
        IF (group%thermal_active .AND. .NOT.group%rheology_implicit) THEN
            CALL evaluate_temperature_from_q(c_qj,temperature,state_ok)
            IF (.NOT.state_ok) THEN
                IF (PRESENT(accepted)) accepted=.FALSE.
                RETURN
            END IF
            CALL eval_thermal_implicit_temperature(SOURCE_ID_THERMAL_COOLING,temperature,                               &
                radiative_term_coeff,convective_term_coeff,T_ambient,c_nh_term_impl(4))
            RETURN
        END IF

        physical_state=>local_state
        IF (PRESENT(physical_scratch)) THEN
            physical_state=>physical_scratch
        ELSE
            CALL initialize_physical_state_complex(local_state)
        END IF
        CALL build_physical_state_from_q(c_qj,physical_state,Bprimej_x,Bprimej_y,0.0_wp,0.0_wp,0.0_wp, &
            accepted=state_ok)
        IF (.NOT.state_ok) THEN
            IF (PRESENT(accepted)) accepted=.FALSE.
            RETURN
        END IF
        IF (group%rheology_implicit) THEN
            CALL eval_rheology_implicit_id(rheology_l1,group%rheology_id,physical_state, &
                stochastic_flag,Zij,c_nh_term_impl(2:3))
        END IF
        IF (group%thermal_active) THEN
            CALL eval_thermal_implicit(SOURCE_ID_THERMAL_COOLING,physical_state, &
                radiative_term_coeff,convective_term_coeff,T_ambient,c_nh_term_impl(4))
        END IF
    END SUBROUTINE eval_source_implicit




    SUBROUTINE eval_geometry_cell_source_explicit( Bprimej_x, Bprimej_y, Bsecondj_xx , Bsecondj_xy ,                   &
            Bsecondj_yy, grav_coeff, qpj, expl_term, time, cell_fract_jk,                                              &
            cell_fissure_fract,                                                                                        &
            lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, cell_area_jk, group, stochastic_value,physical_scratch,          &
            physical_ready,materials )

    USE parameters_2d, ONLY : pore_pressure_flag



    IMPLICIT NONE

        REAL(wp), INTENT(IN) :: Bprimej_x
        REAL(wp), INTENT(IN) :: Bprimej_y
        REAL(wp), INTENT(IN) :: Bsecondj_xx
        REAL(wp), INTENT(IN) :: Bsecondj_xy
        REAL(wp), INTENT(IN) :: Bsecondj_yy
        REAL(wp), INTENT(IN) :: grav_coeff

        REAL(wp), INTENT(IN) :: qpj(n_vars)
        REAL(wp), INTENT(OUT) :: expl_term(n_eqns)

        REAL(wp), INTENT(IN) :: time
        REAL(wp), INTENT(IN) :: cell_fract_jk
        REAL(wp), INTENT(IN) :: cell_fissure_fract(:)

        ! Lateral radial-source per-cell geometry. Pass 0 to disable the lateral
        ! injection.
        REAL(wp), INTENT(IN) :: lat_arc_perim_jk
        REAL(wp), INTENT(IN) :: lat_n_x_jk
        REAL(wp), INTENT(IN) :: lat_n_y_jk
        REAL(wp), INTENT(IN) :: cell_area_jk
        TYPE(source_group_type), INTENT(IN) :: group
        REAL(wp), INTENT(IN) :: stochastic_value
        TYPE(physical_state_real), TARGET, INTENT(INOUT), OPTIONAL :: physical_scratch
        LOGICAL, INTENT(OUT), OPTIONAL :: physical_ready
        TYPE(injection_material_set_type), INTENT(IN) :: materials

        REAL(wp) :: r_h          !< real-value flow thickness
        REAL(wp) :: r_u          !< real-value x-velocity
        REAL(wp) :: r_v          !< real-value y-velocity
        REAL(wp) :: r_rho_m      !< real-value mixture density [kg/m3]



        TYPE(physical_state_real), TARGET :: local_state
        TYPE(physical_state_real), POINTER :: physical_state


        REAL(wp) :: q1

        REAL(wp) :: rheology_rate(2), thermal_rate
        REAL(wp) :: injection_rate(n_eqns)


        expl_term(1:n_eqns) = 0.0_wp
        IF (PRESENT(physical_ready)) physical_ready=.FALSE.

        q1 = 0.0_wp


        IF ( ( qpj(idx_h) .LE. EPSILON(1.0_wp) ) .AND. ( cell_fract_jk .EQ. 0.0_wp )                                   &
                    .AND. ( lat_arc_perim_jk .EQ. 0.0_wp ) .AND.                                                       &
                    ALL(cell_fissure_fract .EQ. 0.0_wp) ) RETURN

        ! Injection into a dry cell still needs density and pore-pressure state.
        physical_state=>local_state
        IF (PRESENT(physical_scratch)) THEN
            physical_state=>physical_scratch
        ELSE
            CALL initialize_physical_state_real(local_state)
        END IF
        CALL build_physical_state_from_qp(qpj, physical_state, Bprimej_x, Bprimej_y,                                   &
            Bsecondj_xx, Bsecondj_xy, Bsecondj_yy)
        IF (PRESENT(physical_ready)) physical_ready=PRESENT(physical_scratch)

        ! Gravity terms - only meaningful when the cell has mass.
        IF ( qpj(idx_h) .GT. EPSILON(1.0_wp) ) THEN

            r_h = qpj(idx_h)
            r_u = qpj(idx_u)
            r_v = qpj(idx_v)

            r_rho_m = physical_state%rho_m

            q1 = r_h * r_rho_m

            IF ( BTEST(group%source_mask,SOURCE_ID_GEOMETRY_CURVATURE-1) .AND.                                         &
                r_h .GT. dry_thickness_tolerance ) THEN

                ! Hydrostatic pressure and the ordinary bed-slope source are already
                ! represented by the HP path. Curvature remains a separate,
                ! velocity-dependent cell source using the fitted Hessian and G_c.
                ! Below the shared dry threshold the auxiliary velocity is not a
                ! resolved physical velocity: do not square it in this source.
                ! This guard neither discards mass nor disables external injection.
                CALL eval_geometry_curvature_explicit( Bprimej_x, Bprimej_y,                                           &
                    Bsecondj_xx, Bsecondj_xy, Bsecondj_yy, grav_coeff, q1,                                             &
                    r_u, r_v, expl_term(2), expl_term(3) )

            END IF

            ! The hydrostatic path/source has no thermal-energy component.
            expl_term(4) = 0.0_wp

            ! Source rates are evaluated at the same unmodified stage state.
            ! The construction of EXPLICIT stays inside each physical module.
            IF (group%rheology_id>0 .AND. r_h>dry_thickness_tolerance) THEN
                CALL eval_rheology_explicit(rheology_l1,group%rheology_id,physical_state,                                          &
                    grav_coeff, pore_pressure_flag,                                                                    &
                    stochastic_flag, stochastic_value, rheology_rate)
                expl_term(2:3)=expl_term(2:3)+rheology_rate

            END IF

            IF (group%thermal_active) THEN
                CALL eval_thermal_explicit(SOURCE_ID_THERMAL_COOLING,physical_state,                                   &
                    radiative_term_coeff,convective_term_coeff,T_ambient,thermal_rate)
                expl_term(4)=expl_term(4)+thermal_rate

            END IF

        END IF

        CALL eval_mass_injection_explicit(physical_state,time,cell_fract_jk,                                           &
            cell_fissure_fract,lat_arc_perim_jk,lat_n_x_jk,lat_n_y_jk,cell_area_jk,                                    &
            BTEST(group%source_mask,SOURCE_ID_BOTTOM_CIRCULAR_INJECTION-1),                                            &
            BTEST(group%source_mask,SOURCE_ID_BOTTOM_FISSURAL_INJECTION-1),                                            &
            BTEST(group%source_mask,SOURCE_ID_RADIAL_INJECTION-1),injection_rate,materials)
        expl_term=expl_term+injection_rate

        RETURN

    END SUBROUTINE eval_geometry_cell_source_explicit











    PURE SUBROUTINE eval_geometry_curvature_explicit( Bprime_x, Bprime_y,                                              &
            Bsecond_xx, Bsecond_xy, Bsecond_yy, grav_coefficient,                                                      &
            mixture_mass, velocity_x, velocity_y, source_x, source_y )

        REAL(wp), INTENT(IN) :: Bprime_x, Bprime_y
        REAL(wp), INTENT(IN) :: Bsecond_xx, Bsecond_xy, Bsecond_yy
        REAL(wp), INTENT(IN) :: grav_coefficient
        REAL(wp), INTENT(IN) :: mixture_mass, velocity_x, velocity_y
        REAL(wp), INTENT(OUT) :: source_x, source_y
        REAL(wp) :: acceleration

        acceleration = curvature_acceleration( Bsecond_xx, Bsecond_xy,                                                 &
            Bsecond_yy, velocity_x, velocity_y )
        source_x = -grav_coefficient * mixture_mass * acceleration * Bprime_x
        source_y = -grav_coefficient * mixture_mass * acceleration * Bprime_y

    END SUBROUTINE eval_geometry_curvature_explicit




    !> \brief Evaluate the atomic explicit mass-exchange rate.
    !> \details Returns rates only; it never mutates conservative state, bed or inventory.
    !> \param[in] qpj Cell-local primitive state.
    !> \param[in] B_zone Connected-water-region indicator.
    !> \param[in] B_prime_x Filtered bed x derivative.
    !> \param[in] B_prime_y Filtered bed y derivative.
    !> \param[in] erodible Per-class available substrate inventory.
    !> \param[in] dt Time-step length.
    !> \param[out] erosion_term Per-class substrate erosion rate.
    !> \param[out] deposition_term Per-class deposition rate.
    !> \param[out] carrier_erosion_term Carrier entrainment rate.
    !> \param[out] carrier_loss_term Limited carrier loss rate.
    !> \param[out] eqns_term Conservative equation source-rate vector.
    !> \param[out] topo_term Proposed bed-elevation rate.
    SUBROUTINE eval_mass_exchange_source_explicit(qpj,B_zone,B_prime_x,B_prime_y,                                      &
            erodible,dt,source_mask,rates,eqns_term,physical_scratch,exchange_scratch,physical_ready)
        REAL(wp), INTENT(IN) :: qpj(n_vars),B_prime_x,B_prime_y
        REAL(wp), INTENT(IN) :: erodible(n_solid),dt
        INTEGER, INTENT(IN) :: B_zone
        INTEGER, INTENT(IN) :: source_mask
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: rates
        REAL(wp), INTENT(OUT) :: eqns_term(n_eqns)
        TYPE(physical_state_real), TARGET, INTENT(INOUT), OPTIONAL :: physical_scratch
        TYPE(mass_exchange_scratch_type), INTENT(INOUT), OPTIONAL :: exchange_scratch
        LOGICAL, INTENT(IN), OPTIONAL :: physical_ready
        TYPE(physical_state_real), TARGET :: local_state
        TYPE(physical_state_real), POINTER :: exchange_state
        LOGICAL :: reuse_prepared_state

        reuse_prepared_state=.FALSE.
        IF (PRESENT(physical_ready)) THEN
            IF (PRESENT(physical_scratch)) reuse_prepared_state=physical_ready
        END IF
        exchange_state=>local_state
        IF (PRESENT(physical_scratch)) THEN
            exchange_state=>physical_scratch
        ELSE
            CALL initialize_physical_state_real(local_state)
        END IF
        ! The combined explicit group already prepared this same qp in its
        ! worker-local scratch. Standalone exchange calls still build normally.
        IF (.NOT.reuse_prepared_state) &
            CALL build_physical_state_from_qp(qpj,exchange_state,B_prime_x,B_prime_y)
        CALL eval_mass_exchange_explicit(rheology_l1,exchange_state,B_zone,                                                        &
            erodible,dt,source_mask,rates,eqns_term,exchange_scratch)
    END SUBROUTINE eval_mass_exchange_source_explicit
END MODULE source_terms_2d
