!> @brief Shear-driven and Rouse-filtered L1 to L2 solid elutriation.
!> @details Implements the LINEAR_USTAR and SHEAR_POWER families of IMEX_v3,
!>          Eqs. (54)-(60). Carrier co-lift and inventory/packing bounds belong
!>          to the shared interlayer coordinator, not to this physical law.
MODULE interlayer_elutriation_2d
    USE parameters_2d, ONLY : wp,n_solid
    USE constitutive_parameters_2d, ONLY : grav,kin_visc_c,rho_a_amb,vonK
    USE solid_transport_parameters_2d, ONLY : rho_s,diam_s,maximum_solid_packing
    USE state_types, ONLY : physical_state_real
    USE interlayer_mass_exchange_parameters_2d, ONLY : interlayer_mass_exchange_config_type
    USE mass_exchange_parameters_2d, ONLY : mass_exchange_config_type
    USE mass_exchange_deposition_2d, ONLY : settling_velocity
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_interlayer_elutriation_explicit
CONTAINS
    !> @brief Return per-class solid volume fluxes upward, without carrier or state updates.
    !> @param[in] config Pair-wide elutriation closure and its input coefficients.
    !> @param[in] mass_l2 Donor-independent particle-drag parameters used in L2.
    !> @param[in] state_l1 Local lower-layer physical state (donor).
    !> @param[in] state_l2 Local upper-layer physical state (may be dry).
    !> @param[in] dt Common source timestep [s].
    !> @param[out] elutriation Per-class upward solid volume rates [m/s].
    SUBROUTINE eval_interlayer_elutriation_explicit(config,mass_l2,state_l1,state_l2,dt,elutriation)
        TYPE(interlayer_mass_exchange_config_type), INTENT(IN) :: config
        TYPE(mass_exchange_config_type), INTENT(IN) :: mass_l2
        TYPE(physical_state_real), INTENT(IN) :: state_l1,state_l2
        REAL(wp), INTENT(IN) :: dt
        REAL(wp), INTENT(OUT) :: elutriation(n_solid)
        REAL(wp), PARAMETER :: hind_exp=4.65_wp, ustar_min=1.0E-6_wp
        REAL(wp) :: gate,slip,speed1,tau12,tausurf,taueff,rho_top,rho_carrier
        REAL(wp) :: ustar,alpha1,alpha2,availability,ws,rouse,phi,shields,phi_shields
        REAL(wp) :: work,shields_denominator,grain_buoyancy,rate
        INTEGER :: i

        elutriation=0.0_wp
        IF (.NOT.config%elutriation_flag) RETURN
        IF (dt<=0.0_wp .OR. state_l1%h<=0.0_wp) RETURN
        alpha1=SUM(state_l1%alphas)
        IF (alpha1<=0.0_wp) RETURN
        alpha2=0.0_wp
        IF (state_l2%h>0.0_wp) alpha2=SUM(state_l2%alphas)
        gate=MAX(0.0_wp,state_l2%h)/(MAX(0.0_wp,state_l2%h)+config%upper_depth_gate)
        speed1=SQRT(state_l1%u**2+state_l1%v**2)
        slip=SQRT((state_l1%u-gate*state_l2%u)**2+ &
            (state_l1%v-gate*state_l2%v)**2)
        rho_top=rho_a_amb
        IF (state_l2%h>0.0_wp) rho_top=state_l2%rho_m
        IF (rho_top<=0.0_wp) RETURN
        rho_carrier=state_l1%rho_c
        IF (state_l2%h>0.0_wp) rho_carrier=state_l2%rho_c
        IF (rho_carrier<=0.0_wp .OR. kin_visc_c<=0.0_wp) RETURN
        tau12=config%shear_interlayer*rho_top*slip**2
        tausurf=config%shear_surface*state_l1%rho_m*speed1**2
        taueff=tau12+tausurf
        ustar=SQRT(MAX(0.0_wp,taueff/rho_top))
        IF (ustar<=ustar_min) RETURN
        availability=MAX(0.0_wp,1.0_wp-MIN(1.0_wp,alpha1/maximum_solid_packing) &
            **config%bed_availability_exponent)

        DO i=1,n_solid
            IF (state_l1%alphas(i)<=0.0_wp) CYCLE
            grain_buoyancy=rho_s(i)-rho_carrier
            IF (grain_buoyancy<=0.0_wp .OR. diam_s(i)<=0.0_wp) CYCLE
            ws=settling_velocity(mass_l2,diam_s(i),rho_s(i),rho_carrier,1.0_wp/kin_visc_c)
            rouse=ws*MAX(0.0_wp,1.0_wp-alpha2)**hind_exp/ &
                (vonK*config%sediment_diffusivity_ratio*ustar)
            phi=MAX(0.0_wp,1.0_wp- &
                (rouse/config%rouse_critical)**config%rouse_exponent)
            SELECT CASE(TRIM(config%elutriation_law))
            CASE('LINEAR_USTAR')
                rate=config%elutriation_efficiency*state_l1%alphas(i)*ustar*phi*availability
            CASE('SHEAR_POWER')
                shields_denominator=grain_buoyancy*grav*diam_s(i)
                shields=taueff/shields_denominator
                phi_shields=0.0_wp
                IF (shields>0.0_wp) phi_shields= &
                    MAX(0.0_wp,1.0_wp-config%shields_critical/shields)**config%shields_exponent
                work=grain_buoyancy*grav*config%lift_length_diameters*diam_s(i)
                rate=config%elutriation_efficiency*state_l1%alphas(i)* &
                    (tau12*slip+tausurf*speed1)/work*phi*phi_shields*availability
            CASE DEFAULT
                CYCLE
            END SELECT
            elutriation(i)=MIN(MAX(0.0_wp,rate), &
                config%elutriation_donor_cap*state_l1%h*state_l1%alphas(i)/dt)
        END DO
    END SUBROUTINE eval_interlayer_elutriation_explicit
END MODULE interlayer_elutriation_2d
