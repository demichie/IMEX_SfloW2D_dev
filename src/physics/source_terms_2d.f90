!> \brief Common physical source assembly entry points.
!> \details Assembles physical source rates and their coupled limiters.
!> Spatial finite-volume operators belong exclusively to spatial modules.
!> The Lie executor applies the resulting limited rates and publishes
!> fluid/substrate candidates at group boundaries.
MODULE source_terms_2d
    USE parameters_2d, ONLY : idx_h
    USE mass_injection_2d, ONLY : eval_mass_injection_explicit
    USE mass_injection_coefficients_2d
    USE parameters_2d
    USE constitutive_parameters_2d
    USE state_conversion_2d
    USE thermal_laws_2d, ONLY : eval_thermal_explicit, eval_thermal_implicit,                                          &
            eval_thermal_analytic
    USE pore_pressure_2d, ONLY : eval_pore_pressure_analytic
    USE mass_exchange_2d, ONLY : eval_mass_exchange_explicit,mass_exchange_scratch_type
    USE rheology_2d, ONLY : curvature_acceleration, eval_rheology_analytic,                                            &
            eval_rheology_implicit_id, eval_rheology_semi_implicit_id,                                                 &
            eval_rheology_explicit
    USE diagnostics_2d, ONLY : fatal_error
    USE source_catalog_2d, ONLY : SOURCE_METHOD_NEWTON
    USE source_catalog_2d, ONLY : SOURCE_METHOD_EXPLICIT, SOURCE_METHOD_ANALYTIC
    USE source_registry_2d, ONLY : source_group_type, group_contains_name
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_source_explicit
    PUBLIC :: limit_explicit_source_rates,assemble_lie_explicit_transaction
    ! Explicit physical-source evaluation provides distinct cell/geometry
    ! and mass-exchange signatures, independent of time-integration stages.
    INTERFACE eval_source_explicit
        MODULE PROCEDURE eval_geometry_cell_source_explicit, eval_mass_exchange_source_explicit
    END INTERFACE
    PUBLIC :: eval_source_semi_implicit, eval_source_implicit, eval_source_analytic
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
            B_zone_cell,context,cell_fract,cell_fissure_fract,arc_perim,arc_n_x,arc_n_y,cell_area,                    &
            non_exchange_group,exchange_selected,physical_scratch,exchange_scratch)
        INTEGER, INTENT(IN) :: B_zone_cell
        TYPE(source_cell_context_type), INTENT(IN) :: context
        REAL(wp), INTENT(IN) :: cell_fract,cell_fissure_fract(:)
        REAL(wp), INTENT(IN) :: arc_perim,arc_n_x,arc_n_y,cell_area
        REAL(wp), INTENT(IN) :: qp_eval(n_vars+2),q_recipient(n_eqns),dt,t
        TYPE(source_group_type), INTENT(IN) :: group
        ! Prepartitioned immutable metadata supplied by the Lie executor.
        TYPE(source_group_type), INTENT(IN), OPTIONAL :: non_exchange_group
        LOGICAL, INTENT(IN), OPTIONAL :: exchange_selected
        TYPE(physical_state_type), INTENT(INOUT), OPTIONAL :: physical_scratch
        TYPE(mass_exchange_scratch_type), INTENT(INOUT), OPTIONAL :: exchange_scratch
        REAL(wp), INTENT(IN) :: inventory(n_solid)
        REAL(wp), INTENT(IN) :: stochastic_value
        REAL(wp), INTENT(OUT) :: rate(n_eqns)
        ! Retain allocations in the thread-owned transfer record across cells.
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: transfer
        LOGICAL, INTENT(OUT) :: has_exchange,valid
        TYPE(source_group_type) :: non_exchange
        REAL(wp) :: physical(n_eqns),exchange(n_eqns),coverage,scale
        IF (PRESENT(exchange_selected)) THEN
            has_exchange=exchange_selected
        ELSE
            has_exchange=ANY(INDEX(group%source_names,'mass_exchange_')==1)
        END IF
        valid=.FALSE.
        physical=0.0_wp
        exchange=0.0_wp
        ! Production groups use a one-time partition. Preserve the original
        ! interface for standalone callers, including mixed explicit groups.
        IF (PRESENT(non_exchange_group)) THEN
            IF (SIZE(non_exchange_group%source_names)>0)                                                               &
                CALL eval_source_explicit(context%Bprime_x,context%Bprime_y,                                           &
                    context%Bsecond_xx,context%Bsecond_xy,context%Bsecond_yy,                                          &
                    context%grav_coeff,qp_eval,physical,t,cell_fract,cell_fissure_fract,                               &
                    arc_perim,arc_n_x,arc_n_y,cell_area,non_exchange_group,stochastic_value,physical_scratch)
        ELSEIF (.NOT.has_exchange) THEN
            CALL eval_source_explicit(context%Bprime_x,context%Bprime_y,                                               &
                context%Bsecond_xx,context%Bsecond_xy,context%Bsecond_yy,                                              &
                context%grav_coeff,qp_eval,physical,t,cell_fract,cell_fissure_fract,                                   &
                arc_perim,arc_n_x,arc_n_y,cell_area,group,stochastic_value,physical_scratch)
        ELSEIF (ANY(INDEX(group%source_names,'mass_exchange_')/=1)) THEN
            non_exchange=group
            non_exchange%source_names=PACK(group%source_names,                                                         &
                INDEX(group%source_names,'mass_exchange_')/=1)
            CALL eval_source_explicit(context%Bprime_x,context%Bprime_y,                                               &
                context%Bsecond_xx,context%Bsecond_xy,context%Bsecond_yy,                                              &
                context%grav_coeff,qp_eval,physical,t,cell_fract,cell_fissure_fract,                                   &
                arc_perim,arc_n_x,arc_n_y,cell_area,non_exchange,stochastic_value,physical_scratch)
        END IF

        IF (has_exchange) THEN
            CALL eval_source_explicit(qp_eval,B_zone_cell,                                                             &
                context%Bprime_x,context%Bprime_y,inventory,dt,                                                        &
                group%source_names,transfer,exchange,physical_scratch,exchange_scratch)

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
                transfer%solid_erosion_rate,scale,valid)
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
            solid_erosion_rate,limiter,valid)
        REAL(wp), INTENT(IN) :: q_reference(n_eqns),inventory(:),dt
        REAL(wp), INTENT(IN) :: rate(n_eqns),solid_erosion_rate(:)
        REAL(wp), INTENT(OUT) :: limiter
        LOGICAL, INTENT(OUT) :: valid
        REAL(wp) :: available_carrier,carrier_consumption
        INTEGER :: solid_eq,i
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

        DO i=1,SIZE(inventory)
            IF (solid_erosion_rate(i)>0.0_wp) limiter=MIN(limiter,                                                     &
                MAX(0.0_wp,inventory(i)/(dt*solid_erosion_rate(i))))

        END DO

        available_carrier=q_reference(idx_totMassEqn)-                                                                 &
            SUM(q_reference(idx_solidEqn_first:idx_solidEqn_last))
        carrier_consumption=SUM(rate(idx_solidEqn_first:idx_solidEqn_last))-rate(idx_totMassEqn)
        IF (carrier_consumption>0.0_wp) limiter=MIN(limiter,                                                           &
            MAX(0.0_wp,available_carrier/(dt*carrier_consumption)))
        ! The limiting ratio places a depleted component exactly at zero in
        ! real arithmetic. Two rounded products (dt*rate and limiter*rate)
        ! can otherwise leave a tiny negative mass (~1e-24) at a dry front.
        ! Stay a few representable steps inside the admissible interval,
        ! only when an actual physical transfer is limited.
        IF (limiter>0.0_wp .AND. limiter<1.0_wp) THEN

            DO i=1,8
                limiter=NEAREST(limiter,-1.0_wp)

            END DO

        END IF

        valid=limiter>=0.0_wp .AND. limiter<=1.0_wp
    END SUBROUTINE limit_explicit_source_rates

    !> Return the conservative equivalent rate of the exact source substep.
    !> The interval is supplied by the selected Lie group operator.
    SUBROUTINE eval_source_analytic(group,q_before,context,interval,rate,status)
        TYPE(source_group_type), INTENT(IN) :: group
        REAL(wp), INTENT(IN) :: q_before(n_vars),interval
        TYPE(source_cell_context_type), INTENT(IN) :: context
        REAL(wp), INTENT(OUT) :: rate(n_eqns)
        INTEGER, INTENT(OUT) :: status
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

        CALL eval_source_analytic_state(group,q_before,context,interval,q_after,status)
        IF (status/=0) RETURN
        rate(1:n_vars)=(q_after-q_before)/interval
    END SUBROUTINE eval_source_analytic

    RECURSIVE SUBROUTINE eval_source_analytic_state(group,q_before,context,interval,q_after,status)
        TYPE(source_group_type), INTENT(IN) :: group
        REAL(wp), INTENT(IN) :: q_before(n_vars),interval
        TYPE(source_cell_context_type), INTENT(IN) :: context
        REAL(wp), INTENT(OUT) :: q_after(n_vars)
        INTEGER, INTENT(OUT) :: status
        TYPE(physical_state_type) :: physical_state
        REAL(wp) :: horizontal_speed,speed_ratio,curvature_accel,velocity_after(2)
        REAL(wp) :: heat_capacity,temperature_after
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

        IF (INDEX(selected_model_name,'rheology_')==1) THEN

            CALL real_physical_state_from_q(q_before, physical_state,                                                  &
                    context%Bprime_x,context%Bprime_y,context%Bsecond_xx,                                              &
                    context%Bsecond_xy,context%Bsecond_yy)
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

            CALL eval_rheology_analytic([selected_model_name], interval, physical_state,                               &
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
            CALL real_physical_state_from_q(q_before, physical_state,                                                  &
                    context%Bprime_x,context%Bprime_y,context%Bsecond_xx,                                              &
                    context%Bsecond_xy,context%Bsecond_yy)

            IF (physical_state%T<=0.0_wp) THEN
                status=-3
                RETURN

            END IF

            heat_capacity=q_before(4)/physical_state%T
            CALL eval_thermal_analytic(['thermal_cooling'],physical_state,interval,                                    &
                    heat_capacity,radiative_term_coeff,convective_term_coeff,                                          &
                    T_env,temperature_after)
            q_after(4)=heat_capacity*temperature_after
            status=0

        ELSE IF (selected_model_name=='pore_pressure_dissipation') THEN
            IF (q_before(1)>0.0_wp) THEN

                CALL real_physical_state_from_q(q_before, physical_state,                                              &
                    context%Bprime_x,context%Bprime_y,context%Bsecond_xx,                                              &
                    context%Bsecond_xy,context%Bsecond_yy)
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
            grav_correction,stochastic_value,rate)
        TYPE(source_group_type), INTENT(IN) :: group
        REAL(wp), INTENT(IN) :: q_cell(n_vars)
        REAL(wp), INTENT(IN) :: Bprime_x,Bprime_y,grav_correction,stochastic_value
        REAL(wp), INTENT(OUT) :: rate(n_eqns)
        TYPE(physical_state_type) :: physical_state

        rate=0.0_wp
        IF (group%method/=SOURCE_METHOD_NEWTON)                                                                        &
            CALL fatal_error('eval_source_semi_implicit requires a Newton Lie group')
        IF (.NOT.group%rheology_implicit) RETURN
        IF (q_cell(idx_totMassEqn)<=0.0_wp) RETURN

        IF (group%rheology_id==0) RETURN
        CALL real_physical_state_from_q(q_cell,physical_state,Bprime_x,Bprime_y,                                       &
            0.0_wp,0.0_wp,0.0_wp)
        CALL eval_rheology_semi_implicit_id(group%rheology_id,physical_state,                                          &
            grav_correction,pore_pressure_flag,stochastic_flag,stochastic_value,rate(2:3))
    END SUBROUTINE eval_source_semi_implicit

    SUBROUTINE eval_source_implicit( Bprimej_x, Bprimej_y, Zij, c_qj,                                                  &
            c_nh_term_impl, r_qj , r_nh_term_impl, group, physical_scratch )

    USE COMPLEXIFY

    IMPLICIT NONE

        REAL(wp), INTENT(IN) :: Bprimej_x
        REAL(wp), INTENT(IN) :: Bprimej_y
        REAL(wp), INTENT(IN):: Zij
        COMPLEX(wp), INTENT(IN), OPTIONAL :: c_qj(n_vars)
        COMPLEX(wp), INTENT(OUT), OPTIONAL :: c_nh_term_impl(n_eqns)
        REAL(wp), INTENT(IN), OPTIONAL :: r_qj(n_vars)
        REAL(wp), INTENT(OUT), OPTIONAL :: r_nh_term_impl(n_eqns)
        TYPE(source_group_type), INTENT(IN) :: group
        TYPE(complex_physical_state_type), OPTIONAL, TARGET, INTENT(INOUT) :: physical_scratch


        COMPLEX(wp) :: qj(n_vars)
        COMPLEX(wp) :: source_term(n_eqns)

        INTEGER :: i

        TYPE(complex_physical_state_type), TARGET :: physical_state_local
        TYPE(complex_physical_state_type), POINTER :: physical_state

        LOGICAL :: selected_rheology, selected_thermal

        IF (PRESENT(c_qj) .AND. PRESENT(c_nh_term_impl) .AND.                                                          &
            .NOT.PRESENT(r_qj) .AND. .NOT.PRESENT(r_nh_term_impl)) THEN

            qj = c_qj

        ELSEIF (PRESENT(r_qj) .AND. PRESENT(r_nh_term_impl) .AND.                                                      &
            .NOT.PRESENT(c_qj) .AND. .NOT.PRESENT(c_nh_term_impl)) THEN

            DO i = 1,n_vars

                qj(i) = CMPLX( r_qj(i),0.0_wp,wp )

            END DO

        ELSE

            CALL fatal_error('eval_source_implicit: provide matching real or complex input/output arguments')

        END IF

        IF (group%method/=SOURCE_METHOD_NEWTON)                                                                        &
            CALL fatal_error('eval_source_implicit requires a Newton Lie group')
        selected_rheology=group%rheology_implicit
        selected_thermal=group%thermal_active

        ! No cross-cell or module-shared state: use the caller's optional scratch.
        IF (PRESENT(physical_scratch)) THEN
            physical_state=>physical_scratch
        ELSE
            physical_state=>physical_state_local
        END IF

        ! initialize the source terms
        source_term(1:n_eqns) = CMPLX(0.0_wp,0.0_wp,wp)

        IF (selected_rheology .OR. selected_thermal) THEN
            CALL physical_state_from_q(qj, physical_state, Bprimej_x, Bprimej_y,                                       &
                0.0_wp, 0.0_wp, 0.0_wp)

        END IF

        IF (selected_rheology) THEN

            CALL eval_rheology_implicit_id(group%rheology_id,physical_state,                                           &
                stochastic_flag,Zij,source_term(2:3))

        END IF

        IF (selected_thermal) THEN
            CALL eval_thermal_implicit(['thermal_cooling'],physical_state,                                             &
                radiative_term_coeff,convective_term_coeff,T_env,source_term(4))

        END IF

        IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

            c_nh_term_impl = source_term

        ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

            r_nh_term_impl = REAL(source_term)

        END IF

        RETURN

    END SUBROUTINE eval_source_implicit




    SUBROUTINE eval_geometry_cell_source_explicit( Bprimej_x, Bprimej_y, Bsecondj_xx , Bsecondj_xy ,                   &
            Bsecondj_yy, grav_coeff, qpj, expl_term, time, cell_fract_jk,                                              &
            cell_fissure_fract,                                                                                        &
            lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, cell_area_jk, group, stochastic_value,physical_scratch )

    USE parameters_2d, ONLY : vel_source , T_source , xs_source , xg_source,                                           &
            time_param , bottom_radial_source_flag,                                                                    &
            pore_pressure_flag , pore_pres_fract ,                                                                     &
            n_intervals , t_intervals , vel_intervals ,                                                                &
            radial_source_flag , h_source, bottom_fissural_source_flag,                                                &
            n_fissures, linear_vel_fissures, T_fissures, time_param_fissures,                                          &
            xs_fissures, xg_fissures



    IMPLICIT NONE

        REAL(wp), INTENT(IN) :: Bprimej_x
        REAL(wp), INTENT(IN) :: Bprimej_y
        REAL(wp), INTENT(IN) :: Bsecondj_xx
        REAL(wp), INTENT(IN) :: Bsecondj_xy
        REAL(wp), INTENT(IN) :: Bsecondj_yy
        REAL(wp), INTENT(IN) :: grav_coeff

        REAL(wp), INTENT(IN) :: qpj(n_vars+2)
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
        TYPE(physical_state_type), TARGET, INTENT(INOUT), OPTIONAL :: physical_scratch

        REAL(wp) :: r_h          !< real-value flow thickness
        REAL(wp) :: r_u          !< real-value x-velocity
        REAL(wp) :: r_v          !< real-value y-velocity
        REAL(wp) :: r_rho_m      !< real-value mixture density [kg/m3]



        TYPE(physical_state_type), TARGET :: local_state
        TYPE(physical_state_type), POINTER :: physical_state


        REAL(wp) :: q1

        REAL(wp) :: rheology_rate(2), thermal_rate, pressure_rate
        REAL(wp) :: injection_rate(n_eqns)
        CHARACTER(LEN=64) :: active_rheology


        expl_term(1:n_eqns) = 0.0_wp

        q1 = 0.0_wp


        IF ( ( qpj(idx_h) .LE. EPSILON(1.0_wp) ) .AND. ( cell_fract_jk .EQ. 0.0_wp )                                   &
                    .AND. ( lat_arc_perim_jk .EQ. 0.0_wp ) .AND.                                                       &
                    ALL(cell_fissure_fract .EQ. 0.0_wp) ) RETURN

        ! Injection into a dry cell still needs density and pore-pressure state.
        physical_state=>local_state
        IF (PRESENT(physical_scratch)) physical_state=>physical_scratch
        CALL physical_state_from_qp(qpj, physical_state, Bprimej_x, Bprimej_y,                                         &
            Bsecondj_xx, Bsecondj_xy, Bsecondj_yy)

        ! Gravity terms - only meaningful when the cell has mass.
        IF ( qpj(idx_h) .GT. EPSILON(1.0_wp) ) THEN

            r_h = qpj(idx_h)
            r_u = qpj(idx_u)
            r_v = qpj(idx_v)

            r_rho_m = physical_state%rho_m

            q1 = r_h * r_rho_m

            IF ( group_contains_name(group,'geometry_curvature') .AND.                                                 &
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
            active_rheology=group%rheology_name
            IF (LEN_TRIM(active_rheology)>0 .AND.                                                                      &
                r_h>dry_thickness_tolerance) THEN
                CALL eval_rheology_explicit([active_rheology], physical_state,                                         &
                    grav_coeff, pore_pressure_flag,                                                                    &
                    stochastic_flag, stochastic_value, rheology_rate)
                expl_term(2:3)=expl_term(2:3)+rheology_rate

            END IF

            IF (group%thermal_active) THEN
                CALL eval_thermal_explicit(['thermal_cooling'],physical_state,                                         &
                    radiative_term_coeff,convective_term_coeff,T_env,thermal_rate)
                expl_term(4)=expl_term(4)+thermal_rate

            END IF

        END IF

        CALL eval_mass_injection_explicit(physical_state,time,cell_fract_jk,                                           &
            cell_fissure_fract,lat_arc_perim_jk,lat_n_x_jk,lat_n_y_jk,cell_area_jk,                                    &
            group_contains_name(group,'bottom_circular_injection'),                                                    &
            group_contains_name(group,'bottom_fissural_injection'),                                                    &
            group_contains_name(group,'radial_injection'),injection_rate)
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
            erodible,dt,selected_names,rates,eqns_term,physical_scratch,exchange_scratch)
        REAL(wp), INTENT(IN) :: qpj(n_vars+2),B_prime_x,B_prime_y
        REAL(wp), INTENT(IN) :: erodible(n_solid),dt
        INTEGER, INTENT(IN) :: B_zone
        CHARACTER(*), INTENT(IN) :: selected_names(:)
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: rates
        REAL(wp), INTENT(OUT) :: eqns_term(n_eqns)
        TYPE(physical_state_type), TARGET, INTENT(INOUT), OPTIONAL :: physical_scratch
        TYPE(mass_exchange_scratch_type), INTENT(INOUT), OPTIONAL :: exchange_scratch
        TYPE(physical_state_type), TARGET :: local_state
        TYPE(physical_state_type), POINTER :: exchange_state

        exchange_state=>local_state
        IF (PRESENT(physical_scratch)) exchange_state=>physical_scratch
        CALL physical_state_from_qp(qpj,exchange_state,B_prime_x,B_prime_y)
        CALL eval_mass_exchange_explicit(exchange_state,B_zone,                                                        &
            erodible,dt,selected_names,rates,eqns_term,exchange_scratch)
    END SUBROUTINE eval_mass_exchange_source_explicit
END MODULE source_terms_2d
