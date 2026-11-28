!> @brief Cell-local simultaneous inventory limiter for paired L1/L2 mass exchange.
!> @details Settling and elutriation supply independent solid-volume fluxes.
!>          The coordinator limits both donors, carrier co-lift and final packing
!>          before returning any accepted flux. It does not update conservative
!>          states, own grid indices or perform a Lie timestep.
MODULE interlayer_mass_exchange_2d
    USE parameters_2d, ONLY : wp,n_solid,n_add_gas,n_vars,n_eqns
    USE parameters_2d, ONLY : idx_totMassEqn,idx_uEqn,idx_vEqn,idx_engyEqn
    USE parameters_2d, ONLY : idx_solidEqn_first,idx_addGasEqn_first
    USE parameters_2d, ONLY : idx_T,idx_solid_first,idx_solid_last,idx_add_gas_first,idx_add_gas_last
    USE diagnostics_2d, ONLY : fatal_error
    USE state_types, ONLY : physical_state_real,physical_state_complex
    USE state_physical, ONLY : build_physical_state_from_qp
    USE state_conversions, ONLY : conservative_to_primitive
    USE state_mixture, ONLY : evaluate_mixture_phase_properties
    USE state_composition, ONLY : conservative_mass_admissible
    USE state_fields, ONLY : initialize_dry_physical_cell
    USE rheology_parameters_2d, ONLY : rheology_config_type
    USE mass_exchange_parameters_2d, ONLY : mass_exchange_config_type
    USE solid_transport_parameters_2d, ONLY : rho_s,sp_heat_s,maximum_solid_packing
    USE interlayer_settling_2d, ONLY : eval_interlayer_settling_explicit
    USE interlayer_elutriation_2d, ONLY : eval_interlayer_elutriation_explicit
    USE interlayer_mass_exchange_parameters_2d, ONLY : interlayer_mass_exchange_config_type
    USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: limit_interlayer_mass_transfer, eval_interlayer_mass_exchange_explicit
    PUBLIC :: apply_interlayer_lie_group

    ! Four physical treatment entry points; unsupported methods fail explicitly.
    PUBLIC :: eval_interlayer_mass_exchange_semi_implicit
    PUBLIC :: eval_interlayer_mass_exchange_implicit
    PUBLIC :: solve_interlayer_mass_exchange_analytic

CONTAINS

    !> @brief Limit gross two-way solid and carrier fluxes against both layer inventories.
    !> @details All volumes are phase volumes per horizontal area [m]; proposed
    !>          and accepted fluxes are [m/s]. The solid rates for all classes in
    !>          each direction share a donor scale, but the two process scales
    !>          are independent when constrained by packing or carrier supply.
    !>          Gross withdrawals never spend same-substep receipts; final
    !>          packing accounts for both gross processes and carrier streams.
    !>          LIMIT_SOLIDS preserves the parcel solid fraction; ALLOW_CONCENTRATED
    !>          may reduce co-lift independently. With the default -1 fraction,
    !>          the lifted parcel takes its solid fraction from the current L1
    !>          donor inventory. Downward carrier is the minimum
    !>          required to keep L1 admissible after both solid fluxes and co-lift.
    !> @param[in] solid_l1 Existing L1 per-class solid volume inventories [m].
    !> @param[in] solid_l2 Existing L2 per-class solid volume inventories [m].
    !> @param[in] carrier_l1 Existing L1 carrier volume inventory [m].
    !> @param[in] carrier_l2 Existing L2 carrier volume inventory [m].
    !> @param[in] settling Proposed per-class downward solid rates [m/s].
    !> @param[in] elutriation Proposed per-class upward solid rates [m/s].
    !> @param[in] dt Common timestep [s].
    !> @param[in] packing Maximum admissible solid volume fraction of either layer.
    !> @param[in] config Canonical pair-wide co-lift configuration.
    !> @param[out] accepted_settling Limited downward solid fluxes [m/s].
    !> @param[out] accepted_elutriation Limited upward solid fluxes [m/s].
    !> @param[out] carrier_down Limited downward carrier volume flux [m/s].
    !> @param[out] carrier_up Limited upward carrier volume flux [m/s].
    !> @param[out] valid False on invalid input; outputs are then zero.
    !> @param[in] rho_carrier_l1 Optional L1 carrier donor density [kg/m3].
    !> @param[in] rho_carrier_l2 Optional L2 carrier donor density [kg/m3].
    !> @param[in] candidate_rho_l1 Optional L1 recipient carrier density [kg/m3].
    !> @param[in] candidate_rho_l2 Optional L2 recipient carrier density [kg/m3].
    !>          Omitted densities default to unity or their corresponding donors.
    PURE SUBROUTINE limit_interlayer_mass_transfer(solid_l1,solid_l2,carrier_l1,carrier_l2, &
            settling,elutriation,dt,packing,config,accepted_settling,accepted_elutriation, &
            carrier_down,carrier_up,valid,rho_carrier_l1,rho_carrier_l2,candidate_rho_l1,candidate_rho_l2)
        REAL(wp), INTENT(IN) :: solid_l1(:),solid_l2(:),carrier_l1,carrier_l2
        REAL(wp), INTENT(IN) :: settling(:),elutriation(:),dt,packing
        TYPE(interlayer_mass_exchange_config_type), INTENT(IN) :: config
        REAL(wp), INTENT(OUT) :: accepted_settling(:),accepted_elutriation(:),carrier_down,carrier_up
        LOGICAL, INTENT(OUT) :: valid
        REAL(wp), INTENT(IN), OPTIONAL :: rho_carrier_l1,rho_carrier_l2
        REAL(wp), INTENT(IN), OPTIONAL :: candidate_rho_l1,candidate_rho_l2
        REAL(wp) :: donor_rho1,donor_rho2,receiver_rho1,receiver_rho2
        REAL(wp) :: up_to_l1,down_to_l1,up_to_l2,down_to_l2,available_up,available_down
        REAL(wp) :: f_settle,f_elut,co_lift,requested,alpha_lift,beta,lo,hi,mid,scale
        REAL(wp) :: initial_scale,down_candidate,f_carrier
        REAL(wp) :: settle_total,elut_total,s1,s2,c1,c2
        INTEGER :: i,k,n

        accepted_settling=0.0_wp
        accepted_elutriation=0.0_wp
        carrier_down=0.0_wp
        carrier_up=0.0_wp
        valid=.FALSE.
        n=SIZE(solid_l1)
        IF (SIZE(solid_l2)/=n .OR. SIZE(settling)/=n .OR. SIZE(elutriation)/=n) RETURN
        IF (SIZE(accepted_settling)/=n .OR. SIZE(accepted_elutriation)/=n) RETURN
        IF (.NOT.ieee_is_finite(dt) .OR. .NOT.ieee_is_finite(packing)) RETURN
        IF (dt<=0.0_wp .OR. packing<=0.0_wp .OR. packing>=1.0_wp) RETURN
        IF (.NOT.ALL(ieee_is_finite(solid_l1)) .OR. .NOT.ALL(ieee_is_finite(solid_l2))) RETURN
        IF (.NOT.ALL(ieee_is_finite(settling)) .OR. .NOT.ALL(ieee_is_finite(elutriation))) RETURN
        IF (.NOT.ieee_is_finite(carrier_l1) .OR. .NOT.ieee_is_finite(carrier_l2)) RETURN
        IF (ANY(solid_l1<0.0_wp) .OR. ANY(solid_l2<0.0_wp) .OR. &
            ANY(settling<0.0_wp) .OR. ANY(elutriation<0.0_wp)) RETURN
        IF (carrier_l1<0.0_wp .OR. carrier_l2<0.0_wp) RETURN
        donor_rho1=1.0_wp
        donor_rho2=1.0_wp
        IF (PRESENT(rho_carrier_l1)) donor_rho1=rho_carrier_l1
        IF (PRESENT(rho_carrier_l2)) donor_rho2=rho_carrier_l2
        receiver_rho1=donor_rho1
        receiver_rho2=donor_rho2
        IF (PRESENT(candidate_rho_l1)) receiver_rho1=candidate_rho_l1
        IF (PRESENT(candidate_rho_l2)) receiver_rho2=candidate_rho_l2
        IF (.NOT.ieee_is_finite(donor_rho1) .OR. .NOT.ieee_is_finite(donor_rho2)) RETURN
        IF (.NOT.ieee_is_finite(receiver_rho1) .OR. .NOT.ieee_is_finite(receiver_rho2)) RETURN
        IF (MIN(donor_rho1,donor_rho2,receiver_rho1,receiver_rho2)<=0.0_wp) RETURN
        up_to_l1=donor_rho1/receiver_rho1
        down_to_l1=donor_rho2/receiver_rho1
        up_to_l2=donor_rho1/receiver_rho2
        down_to_l2=donor_rho2/receiver_rho2
        available_up=carrier_l1/(dt*up_to_l1)
        available_down=carrier_l2/(dt*down_to_l2)
        IF (.NOT.config%settling_flag .AND. ANY(settling>0.0_wp)) RETURN
        IF (.NOT.config%elutriation_flag .AND. ANY(elutriation>0.0_wp)) RETURN
        s1=SUM(solid_l1)
        s2=SUM(solid_l2)
        c1=carrier_l1
        c2=carrier_l2
        IF (config%elutriation_flag) THEN
            alpha_lift=config%elutriation_solid_volume_fraction
            IF (alpha_lift==-1.0_wp) THEN
                ! An omitted fraction uses the current L1 donor inventory,
                ! not the pre-spatial state used to evaluate the solid rate.
                ! With no donor solids, no solid co-lift is possible.
                alpha_lift=1.0_wp
                IF (s1>0.0_wp) alpha_lift=s1/(s1+c1)
            END IF
            IF (.NOT.ieee_is_finite(alpha_lift)) RETURN
            IF (alpha_lift<=0.0_wp .OR. alpha_lift>1.0_wp) RETURN
            SELECT CASE(TRIM(config%elutriation_carrier_mode))
            CASE('LIMIT_SOLIDS','ALLOW_CONCENTRATED')
            CASE DEFAULT
                RETURN
            END SELECT
        ELSE
            alpha_lift=1.0_wp
        END IF

        beta=(1.0_wp-packing)/packing
        IF (s1>packing*(s1+c1)+128.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,s1+c1)) RETURN
        IF (s2>packing*(s2+c2)+128.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,s2+c2)) RETURN

        ! Gross per-species donor bounds: no reliance on simultaneous receipts.
        f_settle=1.0_wp
        f_elut=1.0_wp
        DO i=1,n
            IF (settling(i)>0.0_wp) f_settle=MIN(f_settle,solid_l2(i)/(dt*settling(i)))
            IF (elutriation(i)>0.0_wp) f_elut=MIN(f_elut,solid_l1(i)/(dt*elutriation(i)))
        END DO
        accepted_settling=MIN(1.0_wp,f_settle)*settling
        accepted_elutriation=MIN(1.0_wp,f_elut)*elutriation
        settle_total=SUM(accepted_settling)
        elut_total=SUM(accepted_elutriation)

        co_lift=0.0_wp
        IF (elut_total>0.0_wp .AND. alpha_lift<1.0_wp) THEN
            requested=(1.0_wp-alpha_lift)/alpha_lift*elut_total
            IF (.NOT.ieee_is_finite(requested)) THEN
                valid=.FALSE.
                accepted_settling=0.0_wp
                accepted_elutriation=0.0_wp
                RETURN
            END IF
            IF (config%elutriation_carrier_mode=='LIMIT_SOLIDS' .AND. requested>available_up) THEN
                f_elut=available_up/requested
                accepted_elutriation=f_elut*accepted_elutriation
                elut_total=SUM(accepted_elutriation)
                co_lift=(1.0_wp-alpha_lift)/alpha_lift*elut_total
            ELSE
                co_lift=MIN(requested,available_up)
            END IF
        END IF

        ! Start with a jointly feasible scale. The old common-scale
        ! solution remains a conservative admissible point; it is not the
        ! final physical limit when the two gross transfers can differ.
        ! Packing and carrier caps define a convex feasible region in
        ! the solid-process and carrier scales. The searches start from
        ! an admissible point and cannot spend donor receipts.
        ! In ALLOW_CONCENTRATED, reduce upward carrier independently
        ! before sacrificing either solid flux. This never applies to
        ! LIMIT_SOLIDS, whose lifted parcel composition is prescribed.
        IF (config%elutriation_carrier_mode=='ALLOW_CONCENTRATED' .AND. &
            .NOT.feasible(1.0_wp,1.0_wp,1.0_wp) .AND. feasible(1.0_wp,1.0_wp,0.0_wp)) THEN
            f_settle=1.0_wp
            f_elut=1.0_wp
            lo=0.0_wp
            hi=1.0_wp
            DO k=1,56
                mid=0.5_wp*(lo+hi)
                IF (feasible(1.0_wp,1.0_wp,mid)) THEN
                    lo=mid
                ELSE
                    hi=mid
                END IF
            END DO
            f_carrier=lo*(1.0_wp-128.0_wp*EPSILON(1.0_wp))
        ELSE
            IF (.NOT.feasible(1.0_wp,1.0_wp,1.0_wp)) THEN
                lo=0.0_wp
                hi=1.0_wp
                DO k=1,56
                    mid=0.5_wp*(lo+hi)
                    IF (feasible(mid,mid,mid)) THEN
                        lo=mid
                    ELSE
                        hi=mid
                    END IF
                END DO
                scale=lo*(1.0_wp-128.0_wp*EPSILON(1.0_wp))
            ELSE
                scale=1.0_wp
            END IF

            initial_scale=scale
            f_settle=scale
            f_elut=scale
            IF (scale<1.0_wp) THEN
                ! Relax the descending process first, then ascending co-lift.
                ! Each search starts from the accepted candidate and changes
                ! only its own physical rate; no new limiter representation.
                IF (feasible(1.0_wp,f_elut,f_elut)) THEN
                    f_settle=1.0_wp
                ELSE
                    lo=f_settle
                    hi=1.0_wp
                    DO k=1,56
                        mid=0.5_wp*(lo+hi)
                        IF (feasible(mid,f_elut,f_elut)) THEN
                            lo=mid
                        ELSE
                            hi=mid
                        END IF
                    END DO
                    f_settle=MAX(initial_scale,lo*(1.0_wp-128.0_wp*EPSILON(1.0_wp)))
                END IF
                IF (feasible(f_settle,1.0_wp,1.0_wp)) THEN
                    f_elut=1.0_wp
                ELSE
                    lo=f_elut
                    hi=1.0_wp
                    DO k=1,56
                        mid=0.5_wp*(lo+hi)
                        IF (feasible(f_settle,mid,mid)) THEN
                            lo=mid
                        ELSE
                            hi=mid
                        END IF
                    END DO
                    f_elut=MAX(initial_scale,lo*(1.0_wp-128.0_wp*EPSILON(1.0_wp)))
                END IF
            END IF

            f_carrier=f_elut
        END IF

        accepted_settling=f_settle*accepted_settling
        accepted_elutriation=f_elut*accepted_elutriation
        CALL candidate_flux(f_settle,f_elut,f_carrier,down_candidate,carrier_up)
        carrier_down=down_candidate
        ! Recheck the final independently scaled candidate before exposing
        ! any accepted rate. Failure must never leak partially accepted fluxes.
        IF (ieee_is_finite(carrier_down)) THEN
            IF (carrier_down<=available_down+128.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,available_down)) &
                valid=feasible(f_settle,f_elut,f_carrier)
        END IF
        IF (.NOT.valid) THEN
            accepted_settling=0.0_wp
            accepted_elutriation=0.0_wp
            carrier_down=0.0_wp
            carrier_up=0.0_wp
        END IF

    CONTAINS

        !> @brief Evaluate simultaneous upward/downward carrier rates at independent process scales.
        !> @param[in] fs Settling reduction factor in [0,1].
        !> @param[in] fe Elutriation factor in [0,1].
        !> @param[in] fc Upward carrier factor in [0,fe]; equal to fe in LIMIT_SOLIDS.
        !> @param[out] down Required downward carrier rate [m/s].
        !> @param[out] up Prescribed upward carrier rate [m/s].
        PURE SUBROUTINE candidate_flux(fs,fe,fc,down,up)
            REAL(wp), INTENT(IN) :: fs,fe,fc
            REAL(wp), INTENT(OUT) :: down,up
            up=fc*co_lift
            down=MAX(0.0_wp,(beta*(s1+dt*(fs*settle_total-fe*elut_total))-c1+dt*up_to_l1*up)/ &
                (dt*down_to_l1))
        END SUBROUTINE candidate_flux

        !> @brief Check L2 carrier inventory and packing at the proposed process scales.
        !> @param[in] fs Settling reduction factor in [0,1].
        !> @param[in] fe Elutriation factor in [0,1].
        !> @param[in] fc Upward carrier factor in [0,fe]; equal to fe in LIMIT_SOLIDS.
        !> @return Whether the candidate carrier exchange is admissible in L2.
        PURE LOGICAL FUNCTION feasible(fs,fe,fc) RESULT(ok)
            REAL(wp), INTENT(IN) :: fs,fe,fc
            REAL(wp) :: down,up,solid2_after,carrier2_after,tol
            CALL candidate_flux(fs,fe,fc,down,up)
            solid2_after=s2+dt*(fe*elut_total-fs*settle_total)
            carrier2_after=c2+dt*(up_to_l2*up-down_to_l2*down)
            tol=128.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,s2+c2)
            ok=ieee_is_finite(down) .AND. ieee_is_finite(carrier2_after)
            IF (.NOT.ok) RETURN
            ok=down<=available_down+tol/dt .AND. &
                carrier2_after>=-tol .AND. beta*solid2_after<=carrier2_after+tol
        END FUNCTION feasible
    END SUBROUTINE limit_interlayer_mass_transfer

    !> @brief Evaluate physical rates and assemble a conservative paired exchange.
    !> @details Law evaluations use one common source snapshot. A separate pair of
    !>          post-spatial candidate states supplies inventories and packing bounds,
    !>          preserving the Lie group-0 evaluation/recipient distinction. The
    !>          two source vectors are antisymmetric for mass, species, momentum
    !>          and donor-state thermal content; no terrain or atmosphere is touched.
    !> @param[in] config Pair-wide physical source configuration.
    !> @param[in] rheology_l2 Donor L2 rheology closure used by settling.
    !> @param[in] mass_l2 L2 particle settling drag and carrier properties.
    !> @param[in] eval_l1 L1 physical state for source evaluation.
    !> @param[in] eval_l2 L2 physical state for source evaluation.
    !> @param[in] limit_l1 L1 recipient state for inventory checks.
    !> @param[in] limit_l2 L2 recipient state for inventory checks.
    !> @param[in] dt Physical Lie-substep duration [s].
    !> @param[out] source_l1 Signed L1 conservative rate vector.
    !> @param[out] source_l2 Signed L2 conservative rate vector.
    !> @param[out] valid Whether an admissible paired rate has been assembled.
    !> @param[in] rate_scale Optional common reduction of the already limited transfer rates.
    !> @param[in] carrier_down_extra Optional additional L2-to-L1 carrier volume rate [m/s].
    !> @param[out] carrier_down_used Optional original limited downward carrier rate [m/s].
    SUBROUTINE eval_interlayer_mass_exchange_explicit(config,rheology_l2,mass_l2,eval_l1,eval_l2, &
            limit_l1,limit_l2,dt,source_l1,source_l2,valid,rate_scale,carrier_down_extra,carrier_down_used)
        TYPE(interlayer_mass_exchange_config_type), INTENT(IN) :: config
        TYPE(rheology_config_type), INTENT(IN) :: rheology_l2
        TYPE(mass_exchange_config_type), INTENT(IN) :: mass_l2
        TYPE(physical_state_real), INTENT(IN) :: eval_l1,eval_l2,limit_l1,limit_l2
        REAL(wp), INTENT(IN) :: dt
        REAL(wp), INTENT(OUT) :: source_l1(n_eqns),source_l2(n_eqns)
        LOGICAL, INTENT(OUT) :: valid
        REAL(wp), INTENT(IN), OPTIONAL :: rate_scale,carrier_down_extra
        REAL(wp), INTENT(OUT), OPTIONAL :: carrier_down_used
        REAL(wp) :: proposed_s(n_solid),proposed_e(n_solid)
        REAL(wp) :: limited_s(n_solid),limited_e(n_solid)
        REAL(wp) :: solid1(n_solid),solid2(n_solid),carrier1,carrier2
        REAL(wp) :: down,up,chi,interfacial_u,interfacial_v,solid_mass,carrier_mass
        REAL(wp) :: donor_xg(n_add_gas),gas_mass_fraction
        REAL(wp) :: rho_eval_l1,rho_eval_l2,rho_limit_l1,rho_limit_l2
        INTEGER :: i

        source_l1=0.0_wp
        source_l2=0.0_wp
        valid=.FALSE.
        IF (PRESENT(carrier_down_used)) carrier_down_used=0.0_wp
        IF (.NOT.ieee_is_finite(dt) .OR. dt<=0.0_wp) RETURN

        ! The dry-state carrier density is an ambient placeholder, not the
        ! density of material first received from the other layer. Use the
        ! donor density to convert transferred carrier mass into recipient
        ! volume when that recipient has no material yet. The evaluation
        ! snapshot and the post-spatial inventory may have different wetness.
        rho_eval_l1=eval_l1%rho_c
        rho_eval_l2=eval_l2%rho_c
        rho_limit_l1=limit_l1%rho_c
        rho_limit_l2=limit_l2%rho_c
        IF (eval_l1%h<=0.0_wp) rho_eval_l1=eval_l2%rho_c
        IF (eval_l2%h<=0.0_wp) rho_eval_l2=eval_l1%rho_c
        IF (limit_l1%h<=0.0_wp) rho_limit_l1=eval_l2%rho_c
        IF (limit_l2%h<=0.0_wp) rho_limit_l2=eval_l1%rho_c
        IF (rho_limit_l1<=0.0_wp .OR. rho_limit_l2<=0.0_wp) RETURN

        IF (config%settling_flag .AND. eval_l2%h>0.0_wp) THEN
            CALL eval_interlayer_settling_explicit(config,rheology_l2,mass_l2,eval_l2,dt,proposed_s)
        ELSE
            proposed_s=0.0_wp
        END IF
        IF (config%elutriation_flag .AND. eval_l1%h>0.0_wp) THEN
            CALL eval_interlayer_elutriation_explicit(config,mass_l2,eval_l1,eval_l2,dt,proposed_e)
        ELSE
            proposed_e=0.0_wp
        END IF
        solid1=limit_l1%h*limit_l1%alphas
        solid2=limit_l2%h*limit_l2%alphas
        carrier1=limit_l1%h*MAX(0.0_wp,1.0_wp-SUM(limit_l1%alphas))
        carrier2=limit_l2%h*MAX(0.0_wp,1.0_wp-SUM(limit_l2%alphas))
        CALL limit_interlayer_mass_transfer(solid1,solid2,carrier1,carrier2,proposed_s,proposed_e, &
            dt,maximum_solid_packing,config,limited_s,limited_e,down,up,valid, &
            rho_carrier_l1=rho_eval_l1,rho_carrier_l2=rho_eval_l2, &
            candidate_rho_l1=rho_limit_l1,candidate_rho_l2=rho_limit_l2)
        IF (.NOT.valid) RETURN
        IF (PRESENT(carrier_down_used)) carrier_down_used=down

        ! The optional adjustment is used only after testing the reconstructed
        ! conservative candidate. It preserves gross donor budgets: even a
        ! compensating carrier transfer cannot spend simultaneous receipts.
        IF (PRESENT(rate_scale)) THEN
            IF (.NOT.ieee_is_finite(rate_scale)) THEN
                valid=.FALSE.
                RETURN
            END IF
            IF (rate_scale<0.0_wp .OR. rate_scale>1.0_wp) THEN
                valid=.FALSE.
                RETURN
            END IF
            limited_s=rate_scale*limited_s
            limited_e=rate_scale*limited_e
            down=rate_scale*down
            up=rate_scale*up
        END IF
        IF (PRESENT(carrier_down_extra)) THEN
            IF (.NOT.ieee_is_finite(carrier_down_extra)) THEN
                valid=.FALSE.
                RETURN
            END IF
            IF (carrier_down_extra<0.0_wp) THEN
                valid=.FALSE.
                RETURN
            END IF
            down=down+carrier_down_extra
        END IF
        IF (dt*down*eval_l2%rho_c>carrier2*limit_l2%rho_c+ &
            128.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,carrier2*limit_l2%rho_c)) THEN
            valid=.FALSE.
            RETURN
        END IF

        ! The same interface velocity is debited from L2 and credited to L1.
        chi=config%settling_velocity_weight+(1.0_wp-config%settling_velocity_weight)* &
            config%settling_plug_thickness/(MAX(0.0_wp,eval_l2%h)+config%settling_plug_thickness)
        interfacial_u=eval_l1%u+chi*(eval_l2%u-eval_l1%u)
        interfacial_v=eval_l1%v+chi*(eval_l2%v-eval_l1%v)
        DO i=1,n_solid
            solid_mass=rho_s(i)*limited_s(i)
            CALL add_pair(solid_mass,interfacial_u,interfacial_v, &
                solid_mass*sp_heat_s(i)*eval_l2%T,1.0_wp,idx_solidEqn_first+i-1)
            solid_mass=rho_s(i)*limited_e(i)
            CALL add_pair(solid_mass,eval_l1%u,eval_l1%v, &
                solid_mass*sp_heat_s(i)*eval_l1%T,-1.0_wp,idx_solidEqn_first+i-1)
        END DO

        carrier_mass=eval_l2%rho_c*down
        IF (carrier_mass>0.0_wp) THEN
            donor_xg=eval_l2%xg
            gas_mass_fraction=1.0_wp-SUM(eval_l2%xs)
            IF (gas_mass_fraction<=0.0_wp) THEN
                valid=.FALSE.
                source_l1=0.0_wp
                source_l2=0.0_wp
                RETURN
            END IF
            CALL add_pair(carrier_mass,eval_l2%u,eval_l2%v, &
                carrier_mass*eval_l2%cp_c*eval_l2%T,1.0_wp,0)
            DO i=1,n_add_gas
                CALL add_gas(carrier_mass*donor_xg(i)/gas_mass_fraction,1.0_wp,i)
            END DO
        END IF
        carrier_mass=eval_l1%rho_c*up
        IF (carrier_mass>0.0_wp) THEN
            donor_xg=eval_l1%xg
            gas_mass_fraction=1.0_wp-SUM(eval_l1%xs)
            IF (gas_mass_fraction<=0.0_wp) THEN
                valid=.FALSE.
                source_l1=0.0_wp
                source_l2=0.0_wp
                RETURN
            END IF
            CALL add_pair(carrier_mass,eval_l1%u,eval_l1%v, &
                carrier_mass*eval_l1%cp_c*eval_l1%T,-1.0_wp,0)
            DO i=1,n_add_gas
                CALL add_gas(carrier_mass*donor_xg(i)/gas_mass_fraction,-1.0_wp,i)
            END DO
        END IF

    CONTAINS
        !> @brief Add a signed, paired mass, momentum and thermal energy transfer.
        !> @param[in] mass Positive transferred phase mass rate [kg/m2/s].
        !> @param[in] ux Transfer horizontal x velocity [m/s].
        !> @param[in] uy Transfer horizontal y velocity [m/s].
        !> @param[in] thermal Donor-state sensible thermal rate [J/m2/s].
        !> @param[in] direction Positive L2-to-L1, negative L1-to-L2.
        !> @param[in] species_eq Species equation index, zero for carrier.
        SUBROUTINE add_pair(mass,ux,uy,thermal,direction,species_eq)
            REAL(wp), INTENT(IN) :: mass,ux,uy,thermal,direction
            INTEGER, INTENT(IN) :: species_eq
            source_l1(idx_totMassEqn)=source_l1(idx_totMassEqn)+direction*mass
            source_l1(idx_uEqn)=source_l1(idx_uEqn)+direction*mass*ux
            source_l1(idx_vEqn)=source_l1(idx_vEqn)+direction*mass*uy
            source_l1(idx_engyEqn)=source_l1(idx_engyEqn)+direction*thermal
            IF (species_eq>0) source_l1(species_eq)=source_l1(species_eq)+direction*mass
            source_l2(idx_totMassEqn)=source_l2(idx_totMassEqn)-direction*mass
            source_l2(idx_uEqn)=source_l2(idx_uEqn)-direction*mass*ux
            source_l2(idx_vEqn)=source_l2(idx_vEqn)-direction*mass*uy
            source_l2(idx_engyEqn)=source_l2(idx_engyEqn)-direction*thermal
            IF (species_eq>0) source_l2(species_eq)=source_l2(species_eq)-direction*mass
        END SUBROUTINE add_pair

        !> @brief Add equal and opposite transported additional-carrier species mass.
        !> @param[in] mass Additional carrier constituent mass rate [kg/m2/s].
        !> @param[in] direction Positive L2-to-L1, negative L1-to-L2.
        !> @param[in] index Additional gas species index (1-based).
        SUBROUTINE add_gas(mass,direction,index)
            REAL(wp), INTENT(IN) :: mass,direction
            INTEGER, INTENT(IN) :: index
            source_l1(idx_addGasEqn_first+index-1)=source_l1(idx_addGasEqn_first+index-1)+direction*mass
            source_l2(idx_addGasEqn_first+index-1)=source_l2(idx_addGasEqn_first+index-1)-direction*mass
        END SUBROUTINE add_gas
    END SUBROUTINE eval_interlayer_mass_exchange_explicit

    !> @brief Perform one paired cell-local Lie transfer as an atomic transaction.
    !> @details The physical law uses the group-zero pre-spatial snapshot when
    !>          requested, but always withdraws from the post-spatial/current
    !>          recipients. All four physical workspaces are owned by the caller
    !>          (one quartet per OpenMP worker). No grid indices, allocations or
    !>          additional time integrator are introduced.
    !>          A solid withdrawal that exceeds the conservative donor by only
    !>          floating-point roundoff receives a common paired rate margin
    !>          before committing; any other invalid candidate remains fatal.
    !> @param[in] config Canonical interlayer configuration.
    !> @param[in] rheology_l2 Settling rheology of the donor upper layer.
    !> @param[in] mass_l2 Upper-layer deposition/settling closure.
    !> @param[in] pre_spatial_qp Pre-spatial primitives for both layers.
    !> @param[in] group_zero True to evaluate on the pre-spatial snapshot.
    !> @param[in] dt Positive shared Lie timestep [s].
    !> @param[in,out] q Both conservative candidates, unchanged on failure.
    !> @param[in,out] qp Both primitive candidates, unchanged on failure.
    !> @param[in,out] physical Four caller-owned physical-state workspaces.
    !> @param[out] status Zero on success; 18-21 identify rejected transactions.
    SUBROUTINE apply_interlayer_lie_group(config,rheology_l2,mass_l2,pre_spatial_qp, &
            group_zero,dt,q,qp,physical,status)
        TYPE(interlayer_mass_exchange_config_type), INTENT(IN) :: config
        TYPE(rheology_config_type), INTENT(IN) :: rheology_l2
        TYPE(mass_exchange_config_type), INTENT(IN) :: mass_l2
        REAL(wp), INTENT(IN) :: pre_spatial_qp(n_vars,2)
        LOGICAL, INTENT(IN) :: group_zero
        REAL(wp), INTENT(IN) :: dt
        REAL(wp), INTENT(INOUT) :: q(n_vars,2),qp(n_vars,2)
        TYPE(physical_state_real), INTENT(INOUT) :: physical(4)
        INTEGER, INTENT(OUT) :: status
        REAL(wp) :: candidate_q(n_vars,2),candidate_qp(n_vars,2)
        REAL(wp) :: rate_l1(n_eqns),rate_l2(n_eqns),rate_margin,down_baseline
        REAL(wp) :: best_q(n_vars,2),best_qp(n_vars,2)
        LOGICAL :: ok,roundoff_depletion,has_negative_solid,packed(2)
        LOGICAL :: trial_ok
        INTEGER :: layer,i,eq

        status=18
        IF (group_zero) THEN
            CALL build_physical_state_from_qp(pre_spatial_qp(:,1),physical(1),accepted=ok)
        ELSE
            CALL build_physical_state_from_qp(qp(:,1),physical(1),accepted=ok)
        END IF
        IF (.NOT.ok) RETURN
        IF (group_zero) THEN
            CALL build_physical_state_from_qp(pre_spatial_qp(:,2),physical(2),accepted=ok)
        ELSE
            CALL build_physical_state_from_qp(qp(:,2),physical(2),accepted=ok)
        END IF
        IF (.NOT.ok) RETURN
        CALL build_physical_state_from_qp(qp(:,1),physical(3),accepted=ok)
        IF (.NOT.ok) RETURN
        CALL build_physical_state_from_qp(qp(:,2),physical(4),accepted=ok)
        IF (.NOT.ok) RETURN

        status=19
        CALL eval_interlayer_mass_exchange_explicit(config,rheology_l2,mass_l2, &
            physical(1),physical(2),physical(3),physical(4),dt,rate_l1,rate_l2,ok, &
            carrier_down_used=down_baseline)
        IF (.NOT.ok) RETURN
        candidate_q=q
        candidate_q(1:n_eqns,1)=candidate_q(1:n_eqns,1)+dt*rate_l1
        candidate_q(1:n_eqns,2)=candidate_q(1:n_eqns,2)+dt*rate_l2
        status=20
        IF (.NOT.conservative_mass_admissible(candidate_q(1:n_eqns,1)) .OR. &
            .NOT.conservative_mass_admissible(candidate_q(1:n_eqns,2))) THEN
            ! The physical limiter uses primitive phase volumes, whereas the
            ! recipient is a rounded conservative mass. An exhausted donor
            ! solid may differ by one ULP after the rate is assembled. Reserve
            ! a relative roundoff margin on BOTH paired source vectors, never
            ! clip the conservative state or relax strict admissibility.
            roundoff_depletion=.TRUE.
            has_negative_solid=.FALSE.
            DO layer=1,2
                DO i=1,n_solid
                    eq=idx_solidEqn_first+i-1
                    IF (candidate_q(eq,layer)>=0.0_wp) CYCLE
                    has_negative_solid=.TRUE.
                    IF (q(eq,layer)<=0.0_wp) roundoff_depletion=.FALSE.
                    IF (ABS(candidate_q(eq,layer))>16.0_wp*SPACING(q(eq,layer))) &
                        roundoff_depletion=.FALSE.
                END DO
            END DO
            IF (has_negative_solid .AND. roundoff_depletion) THEN
                rate_margin=1.0_wp-64.0_wp*EPSILON(1.0_wp)
                candidate_q(1:n_eqns,1)=q(1:n_eqns,1)+dt*(rate_margin*rate_l1)
                candidate_q(1:n_eqns,2)=q(1:n_eqns,2)+dt*(rate_margin*rate_l2)
            END IF
            IF (.NOT.conservative_mass_admissible(candidate_q(1:n_eqns,1)) .OR. &
                .NOT.conservative_mass_admissible(candidate_q(1:n_eqns,2))) RETURN
        END IF

        candidate_qp=qp
        status=21
        DO layer=1,2
            IF (.NOT.ANY(candidate_q(:,layer)/=q(:,layer))) CYCLE
            IF (candidate_q(idx_totMassEqn,layer)>0.0_wp) THEN
                CALL conservative_to_primitive(candidate_q(:,layer),candidate_qp(:,layer),accepted=ok)
                IF (.NOT.ok) RETURN
            ELSE
                CALL initialize_dry_physical_cell(candidate_qp(:,layer))
            END IF
        END DO
        CALL evaluate_candidate_packing(candidate_q,candidate_qp,packed,ok)
        IF (.NOT.ok) RETURN
        IF (.NOT.ALL(packed)) THEN
            ! Keep the transaction atomic. First seek extra descending carrier,
            ! then reduce the original coupled source rates when insufficient.
            CALL recover_thermodynamic_packing(best_q,best_qp,trial_ok)
            IF (.NOT.trial_ok) RETURN
            candidate_q=best_q
            candidate_qp=best_qp
        END IF
        q=candidate_q
        qp=candidate_qp
        status=0

    CONTAINS

        !> @brief Check the true post-transfer solid fractions from the canonical EOS.
        !> @param[in] candidate Conservative paired state after a proposed transfer.
        !> @param[in] primitive Reconstructed primitive paired state.
        !> @param[out] fits Packing admissibility for each layer.
        !> @param[out] accepted Whether thermodynamic reconstruction was valid.
        SUBROUTINE evaluate_candidate_packing(candidate,primitive,fits,accepted)
            REAL(wp), INTENT(IN) :: candidate(n_vars,2),primitive(n_vars,2)
            LOGICAL, INTENT(OUT) :: fits(2),accepted
            REAL(wp) :: rho_m,inv_rho_m,rho_c,inv_rho_c,alpha(n_solid),alphag(n_add_gas)
            INTEGER :: k

            fits=.FALSE.
            accepted=.TRUE.
            DO k=1,2
                IF (candidate(idx_totMassEqn,k)<=0.0_wp) THEN
                    fits(k)=.TRUE.
                    CYCLE
                END IF
                IF (.NOT.ieee_is_finite(primitive(idx_T,k)) .OR. primitive(idx_T,k)<=0.0_wp) THEN
                    accepted=.FALSE.
                    RETURN
                END IF
                CALL evaluate_mixture_phase_properties(primitive(idx_T,k), &
                    primitive(idx_solid_first:idx_solid_last,k), &
                    primitive(idx_add_gas_first:idx_add_gas_last,k), &
                    rho_m,inv_rho_m,rho_c,inv_rho_c,alpha,alphag)
                accepted=ieee_is_finite(rho_m) .AND. ieee_is_finite(rho_c) .AND. &
                    ALL(ieee_is_finite(alpha))
                IF (.NOT.accepted) RETURN
                fits(k)=SUM(alpha)<=maximum_solid_packing+ &
                    128.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,maximum_solid_packing)
            END DO
        END SUBROUTINE evaluate_candidate_packing

        !> @brief Assemble and validate one candidate without committing either layer.
        !> @param[in] scale Common factor on the previously accepted physical exchange.
        !> @param[in] extra Additional L2-to-L1 carrier volume rate [m/s].
        !> @param[out] trial_q Conservative candidate, unchanged outside this call.
        !> @param[out] trial_qp Primitive candidate, unchanged outside this call.
        !> @param[out] fits Layer-wise thermodynamic packing validity.
        !> @param[out] accepted True when conservative and primitive states are valid.
        SUBROUTINE trial_transfer(scale,extra,trial_q,trial_qp,fits,accepted)
            REAL(wp), INTENT(IN) :: scale,extra
            REAL(wp), INTENT(OUT) :: trial_q(n_vars,2),trial_qp(n_vars,2)
            LOGICAL, INTENT(OUT) :: fits(2),accepted
            REAL(wp) :: trial_l1(n_eqns),trial_l2(n_eqns)
            INTEGER :: k

            fits=.FALSE.
            CALL eval_interlayer_mass_exchange_explicit(config,rheology_l2,mass_l2, &
                physical(1),physical(2),physical(3),physical(4),dt,trial_l1,trial_l2,accepted, &
                rate_scale=scale,carrier_down_extra=extra)
            IF (.NOT.accepted) RETURN
            trial_q=q
            trial_q(1:n_eqns,1)=q(1:n_eqns,1)+dt*trial_l1
            trial_q(1:n_eqns,2)=q(1:n_eqns,2)+dt*trial_l2
            DO k=1,2
                IF (.NOT.conservative_mass_admissible(trial_q(1:n_eqns,k))) THEN
                    accepted=.FALSE.
                    RETURN
                END IF
            END DO
            trial_qp=qp
            DO k=1,2
                IF (.NOT.ANY(trial_q(:,k)/=q(:,k))) CYCLE
                IF (trial_q(idx_totMassEqn,k)>0.0_wp) THEN
                    CALL conservative_to_primitive(trial_q(:,k),trial_qp(:,k),accepted=accepted)
                    IF (.NOT.accepted) RETURN
                ELSE
                    CALL initialize_dry_physical_cell(trial_qp(:,k))
                END IF
            END DO
            CALL evaluate_candidate_packing(trial_q,trial_qp,fits,accepted)
        END SUBROUTINE trial_transfer

        !> @brief Find a physically admissible correction without sequencing source laws.
        !> @details Sample carrier increments and refine a feasible transition.
        !>          No global monotonicity is assumed across different source scales.
        !> @param[out] final_q Conservative pair for the accepted adjustment.
        !> @param[out] final_qp Primitive pair for the accepted adjustment.
        !> @param[out] found Whether a thermodynamically admissible pair was found.
        SUBROUTINE recover_thermodynamic_packing(final_q,final_qp,found)
            REAL(wp), INTENT(OUT) :: final_q(n_vars,2),final_qp(n_vars,2)
            LOGICAL, INTENT(OUT) :: found
            REAL(wp) :: scale,maximum_extra,extra,previous_extra,lo,hi,mid
            REAL(wp) :: trial_q(n_vars,2),trial_qp(n_vars,2)
            LOGICAL :: fits(2),accepted
            INTEGER :: attempt,j,k

            found=.FALSE.
            DO attempt=1,13
                scale=2.0_wp**(1-attempt)
                ! An extra transfer can use only the post-spatial L2 carrier
                ! inventory remaining after the original gross downward flux.
                maximum_extra=physical(4)%h*MAX(0.0_wp,1.0_wp-SUM(physical(4)%alphas))* &
                    physical(4)%rho_c/(dt*physical(2)%rho_c)-scale*down_baseline
                maximum_extra=MAX(0.0_wp,maximum_extra*(1.0_wp-64.0_wp*EPSILON(1.0_wp)))
                previous_extra=0.0_wp
                DO j=0,16
                    extra=maximum_extra*REAL(j,wp)/16.0_wp
                    CALL trial_transfer(scale,extra,trial_q,trial_qp,fits,accepted)
                    IF (accepted .AND. ALL(fits)) THEN
                        final_q=trial_q
                        final_qp=trial_qp
                        found=.TRUE.
                        IF (j==0) RETURN
                        ! Refine the first sampled admissible transition. The
                        ! complete pair is checked at each trial; failure does
                        ! not mutate either physical recipient.
                        lo=previous_extra
                        hi=extra
                        DO k=1,56
                            mid=0.5_wp*(lo+hi)
                            CALL trial_transfer(scale,mid,trial_q,trial_qp,fits,accepted)
                            IF (accepted .AND. ALL(fits)) THEN
                                hi=mid
                                final_q=trial_q
                                final_qp=trial_qp
                            ELSE
                                lo=mid
                            END IF
                        END DO
                        RETURN
                    END IF
                    previous_extra=extra
                END DO
            END DO

            ! A zero source is admissible whenever the pre-transfer pair is.
            CALL trial_transfer(0.0_wp,0.0_wp,trial_q,trial_qp,fits,accepted)
            IF (accepted .AND. ALL(fits)) THEN
                final_q=trial_q
                final_qp=trial_qp
                found=.TRUE.
            END IF
        END SUBROUTINE recover_thermodynamic_packing
    END SUBROUTINE apply_interlayer_lie_group
    !> @brief Reject unsupported semi-implicit treatment for interlayer mass exchange.
    !> @details No candidate state is committed; report the unsupported method.
    !> @param[in] config Pair-wide solid transfer and carrier co-lift closure.
    !> @param[in] rheology_l2 L2 rheology configuration used by interlayer settling.
    !> @param[in] mass_l2 L2 particle drag and carrier configuration.
    !> @param[in] eval_l1 Lower-layer physical snapshot for law evaluation.
    !> @param[in] eval_l2 Upper-layer physical snapshot for law evaluation.
    !> @param[in] limit_l1 Lower-layer recipient state for donor and packing limits.
    !> @param[in] limit_l2 Upper-layer recipient state for donor and packing limits.
    !> @param[in] dt Proposed paired-source duration [s].
    !> @param[out] source_l1 Zeroed lower-layer conservative source vector before failure.
    !> @param[out] source_l2 Zeroed upper-layer conservative source vector before failure.
    !> @param[out] valid False because this source method is unsupported.
    SUBROUTINE eval_interlayer_mass_exchange_semi_implicit(config,rheology_l2,mass_l2,eval_l1,eval_l2, &
        limit_l1,limit_l2,dt,source_l1,source_l2,valid)
        TYPE(interlayer_mass_exchange_config_type), INTENT(IN) :: config
        TYPE(rheology_config_type), INTENT(IN) :: rheology_l2
        TYPE(mass_exchange_config_type), INTENT(IN) :: mass_l2
        TYPE(physical_state_real), INTENT(IN) :: eval_l1,eval_l2,limit_l1,limit_l2
        REAL(wp), INTENT(IN) :: dt
        REAL(wp), INTENT(OUT) :: source_l1(n_eqns),source_l2(n_eqns)
        LOGICAL, INTENT(OUT) :: valid

        source_l1=0.0_wp
        source_l2=0.0_wp
        valid=.FALSE.
        CALL fatal_error('interlayer_mass_exchange: unsupported semi-implicit physical treatment')
    END SUBROUTINE eval_interlayer_mass_exchange_semi_implicit

    !> @brief Reject unsupported implicit treatment for interlayer mass exchange.
    !> @details No candidate state is committed; report the unsupported method.
    !> @param[in] config Pair-wide solid transfer and carrier co-lift closure.
    !> @param[in] rheology_l2 L2 rheology configuration used by interlayer settling.
    !> @param[in] mass_l2 L2 particle drag and carrier configuration.
    !> @param[in] eval_l1 Lower-layer physical snapshot for law evaluation.
    !> @param[in] eval_l2 Upper-layer physical snapshot for law evaluation.
    !> @param[in] limit_l1 Lower-layer recipient state for donor and packing limits.
    !> @param[in] limit_l2 Upper-layer recipient state for donor and packing limits.
    !> @param[in] dt Proposed paired-source duration [s].
    !> @param[out] source_l1 Zeroed lower-layer conservative source vector before failure.
    !> @param[out] source_l2 Zeroed upper-layer conservative source vector before failure.
    !> @param[out] valid False because this source method is unsupported.
    SUBROUTINE eval_interlayer_mass_exchange_implicit(config,rheology_l2,mass_l2,eval_l1,eval_l2,limit_l1,limit_l2,dt, &
        source_l1,source_l2,valid)
        TYPE(interlayer_mass_exchange_config_type), INTENT(IN) :: config
        TYPE(rheology_config_type), INTENT(IN) :: rheology_l2
        TYPE(mass_exchange_config_type), INTENT(IN) :: mass_l2
        TYPE(physical_state_complex), INTENT(IN) :: eval_l1,eval_l2,limit_l1,limit_l2
        REAL(wp), INTENT(IN) :: dt
        REAL(wp), INTENT(OUT) :: source_l1(n_eqns),source_l2(n_eqns)
        LOGICAL, INTENT(OUT) :: valid

        source_l1=0.0_wp
        source_l2=0.0_wp
        valid=.FALSE.
        CALL fatal_error('interlayer_mass_exchange: unsupported implicit physical treatment')
    END SUBROUTINE eval_interlayer_mass_exchange_implicit

    !> @brief Reject unsupported analytic treatment for interlayer mass exchange.
    !> @details No candidate state is committed; report the unsupported method.
    !> @param[in] config Pair-wide solid transfer and carrier co-lift closure.
    !> @param[in] rheology_l2 L2 rheology configuration used by interlayer settling.
    !> @param[in] mass_l2 L2 particle drag and carrier configuration.
    !> @param[in] eval_l1 Lower-layer physical snapshot for law evaluation.
    !> @param[in] eval_l2 Upper-layer physical snapshot for law evaluation.
    !> @param[in] limit_l1 Lower-layer recipient state for donor and packing limits.
    !> @param[in] limit_l2 Upper-layer recipient state for donor and packing limits.
    !> @param[in] dt Proposed paired-source duration [s].
    !> @param[out] source_l1 Zeroed lower-layer conservative source vector before failure.
    !> @param[out] source_l2 Zeroed upper-layer conservative source vector before failure.
    !> @param[out] valid False because this source method is unsupported.
    SUBROUTINE solve_interlayer_mass_exchange_analytic(config,rheology_l2,mass_l2,eval_l1,eval_l2,limit_l1, &
            limit_l2,dt, &
        source_l1,source_l2,valid)
        TYPE(interlayer_mass_exchange_config_type), INTENT(IN) :: config
        TYPE(rheology_config_type), INTENT(IN) :: rheology_l2
        TYPE(mass_exchange_config_type), INTENT(IN) :: mass_l2
        TYPE(physical_state_real), INTENT(IN) :: eval_l1,eval_l2,limit_l1,limit_l2
        REAL(wp), INTENT(IN) :: dt
        REAL(wp), INTENT(OUT) :: source_l1(n_eqns),source_l2(n_eqns)
        LOGICAL, INTENT(OUT) :: valid

        source_l1=0.0_wp
        source_l2=0.0_wp
        valid=.FALSE.
        CALL fatal_error('interlayer_mass_exchange: unsupported analytic physical treatment')
    END SUBROUTINE solve_interlayer_mass_exchange_analytic

END MODULE interlayer_mass_exchange_2d
