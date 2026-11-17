!> @brief Cell-local simultaneous inventory limiter for paired L1/L2 mass exchange.
!> @details Settling and elutriation supply independent solid-volume fluxes.
!>          The coordinator limits both donors, carrier co-lift and final packing
!>          before returning any accepted flux. It does not update conservative
!>          states, own grid indices or perform a Lie timestep.
MODULE interlayer_exchange_2d
    USE parameters_2d, ONLY : wp
    USE interlayer_exchange_parameters_2d, ONLY : interlayer_exchange_config_type
    USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: limit_interlayer_mass_transfer

CONTAINS

    !> @brief Limit gross two-way solid and carrier fluxes against both layer inventories.
    !> @details All volumes are phase volumes per horizontal area [m]; proposed
    !>          and accepted fluxes are [m/s]. The solid rates for all classes in
    !>          each direction share a donor scale. Gross withdrawals never spend
    !>          same-substep receipts, whereas final packing uses both processes.
    !>          LIMIT_SOLIDS preserves the parcel solid fraction; ALLOW_CONCENTRATED
    !>          may reduce co-lift independently. Downward carrier is the minimum
    !>          required to keep L1 admissible after both solid fluxes and co-lift.
    !> @param[in] solid_l1,solid_l2 Existing per-class solid volume inventories [m].
    !> @param[in] carrier_l1,carrier_l2 Existing carrier volume inventories [m].
    !> @param[in] settling,elutriation Proposed per-class positive solid rates [m/s].
    !> @param[in] dt Common timestep [s].
    !> @param[in] packing Maximum admissible solid volume fraction of either layer.
    !> @param[in] config Canonical pair-wide co-lift configuration.
    !> @param[out] accepted_settling,accepted_elutriation Limited solid fluxes [m/s].
    !> @param[out] carrier_down,carrier_up Limited carrier volume fluxes [m/s].
    !> @param[out] valid False on invalid input; outputs are then zero.
    PURE SUBROUTINE limit_interlayer_mass_transfer(solid_l1,solid_l2,carrier_l1,carrier_l2, &
            settling,elutriation,dt,packing,config,accepted_settling,accepted_elutriation, &
            carrier_down,carrier_up,valid)
        REAL(wp), INTENT(IN) :: solid_l1(:),solid_l2(:),carrier_l1,carrier_l2
        REAL(wp), INTENT(IN) :: settling(:),elutriation(:),dt,packing
        TYPE(interlayer_exchange_config_type), INTENT(IN) :: config
        REAL(wp), INTENT(OUT) :: accepted_settling(:),accepted_elutriation(:),carrier_down,carrier_up
        LOGICAL, INTENT(OUT) :: valid
        REAL(wp) :: f_settle,f_elut,co_lift,requested,alpha_lift,beta,lo,hi,mid,scale
        REAL(wp) :: settle_total,elut_total,s1,s2,c1,c2,unused_down,unused_up
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
        IF (.NOT.config%settling_flag .AND. ANY(settling>0.0_wp)) RETURN
        IF (.NOT.config%elutriation_flag .AND. ANY(elutriation>0.0_wp)) RETURN
        IF (config%elutriation_flag) THEN
            alpha_lift=config%elutriation_solid_volume_fraction
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
        s1=SUM(solid_l1)
        s2=SUM(solid_l2)
        c1=carrier_l1
        c2=carrier_l2
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
            IF (config%elutriation_carrier_mode=='LIMIT_SOLIDS' .AND. requested>c1/dt) THEN
                f_elut=(c1/dt)/requested
                accepted_elutriation=f_elut*accepted_elutriation
                elut_total=SUM(accepted_elutriation)
                co_lift=(1.0_wp-alpha_lift)/alpha_lift*elut_total
            ELSE
                co_lift=MIN(requested,c1/dt)
            END IF
        END IF

        ! Final packing uses the combined changes, but no gross donor withdrawal
        ! can use a simultaneous credit. The feasible shared scale interval
        ! contains zero; bisection is entered only for a binding packing bound.
        CALL candidate_flux(1.0_wp,unused_down,unused_up)
        IF (.NOT.feasible(1.0_wp)) THEN
            lo=0.0_wp
            hi=1.0_wp
            DO k=1,56
                mid=0.5_wp*(lo+hi)
                IF (feasible(mid)) THEN
                    lo=mid
                ELSE
                    hi=mid
                END IF
            END DO
            scale=lo*(1.0_wp-128.0_wp*EPSILON(1.0_wp))
        ELSE
            scale=1.0_wp
        END IF

        accepted_settling=scale*accepted_settling
        accepted_elutriation=scale*accepted_elutriation
        carrier_up=scale*co_lift
        carrier_down=MAX(0.0_wp,(beta*(s1+dt* &
            (SUM(accepted_settling)-SUM(accepted_elutriation)))-c1+dt*carrier_up)/dt)
        IF (.NOT.ieee_is_finite(carrier_down)) RETURN
        IF (carrier_down>c2/dt+128.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,c2/dt)) RETURN
        valid=.TRUE.

    CONTAINS

        PURE SUBROUTINE candidate_flux(f,down,up)
            REAL(wp), INTENT(IN) :: f
            REAL(wp), INTENT(OUT) :: down,up
            up=f*co_lift
            down=MAX(0.0_wp,(beta*(s1+dt*f*(settle_total-elut_total))-c1+dt*up)/dt)
        END SUBROUTINE candidate_flux

        PURE LOGICAL FUNCTION feasible(f) RESULT(ok)
            REAL(wp), INTENT(IN) :: f
            REAL(wp) :: down,up,solid2_after,carrier2_after,tol
            CALL candidate_flux(f,down,up)
            solid2_after=s2+dt*f*(elut_total-settle_total)
            carrier2_after=c2+dt*(up-down)
            tol=128.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,s2+c2)
            ok=ieee_is_finite(down) .AND. ieee_is_finite(carrier2_after)
            IF (.NOT.ok) RETURN
            ok=down<=c2/dt+tol/dt .AND. &
                carrier2_after>=-tol .AND. beta*solid2_after<=carrier2_after+tol
        END FUNCTION feasible
    END SUBROUTINE limit_interlayer_mass_transfer
END MODULE interlayer_exchange_2d
