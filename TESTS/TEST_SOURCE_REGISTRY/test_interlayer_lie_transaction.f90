!> @brief Exercise the exact paired Lie transaction used by the production executor.
PROGRAM test_interlayer_lie_transaction
    USE parameters_2d, ONLY : wp,n_layers,n_solid,n_add_gas,n_vars,n_eqns,n_stoch_vars,n_pore_vars
    USE parameters_2d, ONLY : gas_flag,sutherland_flag,configure_canonical_state_indices,idx_h,idx_T,idx_u,idx_v
    USE parameters_2d, ONLY : idx_solid_first,idx_solid_last,idx_stoch,idx_pore
    USE parameters_2d, ONLY : idx_totMassEqn,idx_uEqn,idx_vEqn,idx_engyEqn
    USE parameters_2d, ONLY : idx_solidEqn_first,idx_solidEqn_last
    USE constitutive_parameters_2d, ONLY : sp_heat_a,sp_gas_const_a,sp_heat_g,sp_gas_const_g
    USE constitutive_parameters_2d, ONLY : pres,inv_pres,T_ambient,grav,kin_visc_c,rho_a_amb,vonK
    USE solid_transport_parameters_2d, ONLY : rho_s,inv_rho_s,sp_heat_s,diam_s,maximum_solid_packing
    USE state_types, ONLY : physical_state_real,initialize_physical_state_real
    USE state_physical, ONLY : build_physical_state_from_qp
    USE state_mixture, ONLY : evaluate_mixture_from_volume
    USE state_conversions, ONLY : primitive_to_conservative,conservative_to_primitive
    USE interlayer_mass_exchange_2d, ONLY : apply_interlayer_lie_group
    USE mass_exchange_2d, ONLY : eval_mass_exchange_explicit
    USE source_catalog_2d, ONLY : SOURCE_ID_EXCHANGE_DEPOSITION, SOURCE_ID_EXCHANGE_VAPORIZATION
    USE source_catalog_2d, ONLY : SOURCE_ID_EXCHANGE_ENTRAINMENT
    USE parameters_2d, ONLY : mass_exchange_rates_type
    USE interlayer_momentum_exchange_2d, ONLY : apply_interlayer_momentum_exchange_analytic, &
        interlayer_momentum_exchange_config_type
    USE interlayer_mass_exchange_parameters_2d, ONLY : interlayer_mass_exchange_config_type
    USE mass_exchange_parameters_2d, ONLY : mass_exchange_config_type
    USE rheology_parameters_2d, ONLY : rheology_config_type
    IMPLICIT NONE

    TYPE(physical_state_real) :: scratch(4),check1,check2
    TYPE(mass_exchange_rates_type) :: deposition_rates
    TYPE(interlayer_mass_exchange_config_type) :: config
    TYPE(mass_exchange_config_type) :: mass_l2
    TYPE(rheology_config_type) :: rheology_l2
    TYPE(interlayer_momentum_exchange_config_type) :: drag
    REAL(wp) :: q(6,2),qp(6,2),before(6,2),primitives(6,2),snapshot(6,2)
    REAL(wp) :: after_group0(6,2),after_later(6,2)
    REAL(wp) :: after_transfer(6,2),slip_before,slip_after
    REAL(wp) :: dt,exchange_source(6)
    REAL(wp) :: alpha_test(2),alphag_test(0),xs_test(2),xg_test(0)
    REAL(wp) :: rho_test,inv_rho_test,rho_c_test,xc_test,cp_c_test,cp_m_test
    LOGICAL :: ok
    INTEGER :: i,status

    n_layers=2; n_solid=2; n_add_gas=0; n_vars=6; n_eqns=6
    n_stoch_vars=0; n_pore_vars=0
    CALL configure_canonical_state_indices()
    gas_flag=.TRUE.; sutherland_flag=.FALSE.
    pres=1.0E5_wp; inv_pres=1.0_wp/pres
    sp_gas_const_a=287.0_wp; sp_heat_a=1000.0_wp
    T_ambient=300.0_wp; rho_a_amb=pres/(sp_gas_const_a*T_ambient)
    grav=9.81_wp; kin_visc_c=1.0E-5_wp; vonK=0.41_wp
    maximum_solid_packing=0.6_wp
    ALLOCATE(rho_s(2),inv_rho_s(2),sp_heat_s(2),diam_s(2))
    ALLOCATE(sp_heat_g(0),sp_gas_const_g(0))
    rho_s=[2200.0_wp,2500.0_wp]
    inv_rho_s=1.0_wp/rho_s
    sp_heat_s=[800.0_wp,900.0_wp]
    diam_s=[2.0E-5_wp,5.0E-5_wp]
    DO i=1,4
        CALL initialize_physical_state_real(scratch(i))
    END DO
    CALL initialize_physical_state_real(check1)
    CALL initialize_physical_state_real(check2)
    dt=0.2_wp
    config%elutriation_flag=.TRUE.
    config%settling_flag=.TRUE.
    config%elutriation_solid_volume_fraction=0.1_wp
    config%rouse_critical=100.0_wp

    ! The four physical column configurations exercise the same paired call.
    CALL make_column(0.0_wp,0.0_wp,20.0_wp,2.0_wp,q,qp)
    before=q; snapshot=qp
    CALL apply_interlayer_lie_group(config,rheology_l2,mass_l2,snapshot,.FALSE.,dt,q,qp,scratch,status)
    CALL require(status==0,'fully dry paired Lie group is valid')
    CALL require(ALL(q==before),'fully dry paired group remains a no-op')

    CALL make_column(1.0_wp,0.0_wp,20.0_wp,0.0_wp,q,qp)
    before=q; snapshot=qp
    CALL apply_interlayer_lie_group(config,rheology_l2,mass_l2,snapshot,.FALSE.,dt,q,qp,scratch,status)
    CALL require(status==0,'elutriation creates a formerly dry L2')
    CALL require(q(1,2)>0.0_wp .AND. qp(idx_h,2)>0.0_wp,'new L2 has positive mass and thickness')
    CALL check_conservation(before,q,'cold-start L2')
    CALL check_admissibility(qp)

    CALL make_column(0.0_wp,1.0_wp,0.0_wp,2.0_wp,q,qp)
    before=q; snapshot=qp
    CALL apply_interlayer_lie_group(config,rheology_l2,mass_l2,snapshot,.FALSE.,dt,q,qp,scratch,status)
    CALL require(status==0,'settling creates a formerly dry L1')
    CALL require(q(1,1)>0.0_wp .AND. qp(idx_h,1)>0.0_wp,'new L1 has positive mass and thickness')
    CALL check_conservation(before,q,'cold-start L1')
    CALL check_admissibility(qp)

    ! The canonical dry-state density is ambient, but the first material
    ! transferred into the empty layer has the hot donor's carrier density.
    CALL make_column(0.0_wp,1.0_wp,0.0_wp,2.0_wp,q,qp)
    qp(idx_T,2)=650.0_wp
    CALL primitive_to_conservative(qp(:,2),q(:,2),accepted=ok)
    CALL require(ok,'hot L2 donor conservative conversion')
    CALL build_physical_state_from_qp(qp(:,1),check1,accepted=ok)
    CALL require(ok,'dry L1 physical reconstruction')
    CALL build_physical_state_from_qp(qp(:,2),check2,accepted=ok)
    CALL require(ok,'hot L2 physical reconstruction')
    CALL require(ABS(check1%rho_c-check2%rho_c)>0.1_wp*check1%rho_c, &
        'cold-start test distinguishes ambient and donor carrier densities')
    before=q; snapshot=qp
    CALL apply_interlayer_lie_group(config,rheology_l2,mass_l2,snapshot,.FALSE.,dt,q,qp,scratch,status)
    CALL require(status==0,'hot L2 settles into initially dry L1')
    CALL require(q(idx_totMassEqn,1)>0.0_wp,'hot donor creates lower-layer mass')
    CALL build_physical_state_from_qp(qp(:,1),check1,accepted=ok)
    CALL require(ok,'hot donor produces a valid lower-layer physical state')
    IF (ABS(SUM(check1%alphas)-maximum_solid_packing)>=1.0E-8_wp) THEN
        WRITE(*,'(A,ES24.16)') 'DIAG dry L1 solid fraction: ',SUM(check1%alphas)
        WRITE(*,'(A,ES24.16)') 'DIAG maximum packing: ',maximum_solid_packing
        WRITE(*,'(A,ES24.16)') 'DIAG final L1 temperature: ',check1%T
        WRITE(*,'(A,ES24.16)') 'DIAG final L1 carrier density: ',check1%rho_c
        WRITE(*,'(A,ES24.16)') 'DIAG donor L2 carrier density: ',check2%rho_c
        WRITE(*,'(A,ES24.16)') 'DIAG final L1 thickness: ',qp(idx_h,1)
    END IF
    CALL require(ABS(SUM(check1%alphas)-maximum_solid_packing)<1.0E-8_wp, &
        'first settled L1 deposit uses donor carrier density at packing')
    CALL check_conservation(before,q,'hot donor cold-start L1')
    CALL check_admissibility(qp)

    ! The reverse transfer must likewise initialize an empty L2 using the
    ! L1 donor properties rather than its ambient dry-state placeholder.
    CALL make_column(1.0_wp,0.0_wp,20.0_wp,0.0_wp,q,qp)
    qp(idx_T,1)=650.0_wp
    CALL primitive_to_conservative(qp(:,1),q(:,1),accepted=ok)
    CALL require(ok,'hot L1 donor conservative conversion')
    before=q; snapshot=qp
    CALL apply_interlayer_lie_group(config,rheology_l2,mass_l2,snapshot,.FALSE.,dt,q,qp,scratch,status)
    CALL require(status==0,'hot L1 elutriates into initially dry L2')
    CALL require(q(idx_totMassEqn,2)>0.0_wp,'hot donor creates upper-layer mass')
    CALL build_physical_state_from_qp(qp(:,2),check2,accepted=ok)
    CALL require(ok,'hot donor produces a valid upper-layer physical state')
    CALL require(ABS(SUM(check2%alphas)-config%elutriation_solid_volume_fraction)<1.0E-8_wp, &
        'first elutriated L2 parcel retains prescribed solid fraction')
    CALL check_conservation(before,q,'hot donor cold-start L2')
    CALL check_admissibility(qp)

    ! A concentrated hot L1 cools under settling from cold L2. The ideal-gas
    ! carrier contracts, so packing must be checked after energy reconstruction
    ! rather than with the carrier densities at the start of the group.
    config%elutriation_flag=.FALSE.
    CALL make_column(0.5_wp,10.0_wp,0.0_wp,2.0_wp,q,qp)
    alpha_test=[0.35_wp,0.2499_wp]
    CALL evaluate_mixture_from_volume(900.0_wp,alphag_test,alpha_test, &
        rho_test,inv_rho_test,rho_c_test,xs_test,xg_test,xc_test,cp_c_test,cp_m_test)
    qp(idx_solid_first:idx_solid_last,1)=xs_test
    qp(idx_T,1)=900.0_wp
    qp(idx_T,2)=280.0_wp
    DO i=1,2
        CALL primitive_to_conservative(qp(:,i),q(:,i),accepted=ok)
        CALL require(ok,'concentrated thermal donor conservative conversion')
    END DO
    CALL build_physical_state_from_qp(qp(:,1),check1,accepted=ok)
    CALL require(ok,'concentrated thermal recipient initialized')
    CALL require(ABS(SUM(check1%alphas)-SUM(alpha_test))<1.0E-10_wp, &
        'concentrated thermal recipient initial packing')
    before=q; snapshot=qp
    CALL apply_interlayer_lie_group(config,rheology_l2,mass_l2,snapshot,.FALSE.,dt,q,qp,scratch,status)
    CALL require(status==0,'thermodynamic packing transaction accepted')
    CALL check_conservation(before,q,'thermodynamic packing correction')
    CALL check_admissibility(qp)
    config%elutriation_flag=.TRUE.

    ! A nearly packed hot L2 receives cold solids from L1. Even with a
    ! downward-free transaction, cooling can require additional upward carrier.
    config%settling_flag=.FALSE.
    config%elutriation_solid_volume_fraction=1.0_wp
    CALL make_column(1.0_wp,0.2_wp,20.0_wp,2.0_wp,q,qp)
    alpha_test=[0.35_wp,0.2499_wp]
    CALL evaluate_mixture_from_volume(900.0_wp,alphag_test,alpha_test, &
        rho_test,inv_rho_test,rho_c_test,xs_test,xg_test,xc_test,cp_c_test,cp_m_test)
    qp(idx_solid_first:idx_solid_last,2)=xs_test
    qp(idx_T,1)=280.0_wp
    qp(idx_T,2)=900.0_wp
    DO i=1,2
        CALL primitive_to_conservative(qp(:,i),q(:,i),accepted=ok)
        CALL require(ok,'upper-layer thermal packing fixture conversion')
    END DO
    before=q; snapshot=qp
    CALL apply_interlayer_lie_group(config,rheology_l2,mass_l2,snapshot,.FALSE.,dt,q,qp,scratch,status)
    CALL require(status==0,'hot packed L2 accepts cold elutriated solids and carrier correction')
    CALL require(q(idx_totMassEqn,2)>before(idx_totMassEqn,2),'cold solids reach packed L2')
    CALL check_conservation(before,q,'upper-layer thermodynamic packing correction')
    CALL check_admissibility(qp)
    config%settling_flag=.TRUE.
    config%elutriation_solid_volume_fraction=0.1_wp

    CALL make_column(1.0_wp,1.0_wp,20.0_wp,2.0_wp,q,qp)
    before=q; snapshot=qp
    CALL apply_interlayer_lie_group(config,rheology_l2,mass_l2,snapshot,.FALSE.,dt,q,qp,scratch,status)
    CALL require(status==0,'both wet layers exchange conservatively')
    CALL require(ANY(q(:,1)/=before(:,1)),'both wet layers change state')
    CALL check_conservation(before,q,'wet L1 and L2')
    CALL check_admissibility(qp)

    ! Validate the omitted-fraction mode with simultaneous bidirectional
    ! transfer on complete conservative/primitive recipient states.
    config%elutriation_solid_volume_fraction=-1.0_wp
    CALL make_column(1.0_wp,1.0_wp,20.0_wp,2.0_wp,q,qp)
    before=q; snapshot=qp
    CALL apply_interlayer_lie_group(config,rheology_l2,mass_l2,snapshot,.FALSE.,dt,q,qp,scratch,status)
    CALL require(status==0,'dynamic co-lift with settling accepts both wet layers')
    CALL check_conservation(before,q,'dynamic two-way Lie mass exchange')
    CALL check_admissibility(qp)
    config%elutriation_solid_volume_fraction=0.1_wp

    ! Geometric contact must not let a sole wet L2 deposit directly to bed.
    ! The same physical deposition law is allowed for L1, including after
    ! interlayer settling has created that formerly dry lower layer.
    mass_l2%deposition%settling_flag=.TRUE.
    CALL make_column(0.0_wp,1.0_wp,0.0_wp,2.0_wp,q,qp)
    CALL build_physical_state_from_qp(qp(:,2),check2,accepted=ok)
    CALL require(ok,'wet L2 physical state for substrate-deposition test')
    CALL eval_mass_exchange_explicit(rheology_l2,mass_l2,check2,0,[0.0_wp,0.0_wp],dt, &
        IBSET(0,SOURCE_ID_EXCHANGE_DEPOSITION-1),deposition_rates,exchange_source, &
        transport_stochastic=.FALSE.,pore_enabled=.FALSE.,terrain_exposed=.TRUE., &
        atmosphere_exposed=.TRUE.,static_deposition_allowed=.FALSE.)
    CALL require(ALL(deposition_rates%solid_deposition_rate==0.0_wp) .AND. &
        ALL(exchange_source==0.0_wp),'isolated L2 cannot deposit directly into substrate')
    ! L2 vaporization must also remain disabled: it derives from the
    ! substrate deposition rather than being an independent mass source.
    CALL eval_mass_exchange_explicit(rheology_l2,mass_l2,check2,0,[0.0_wp,0.0_wp],dt, &
        IBSET(IBSET(0,SOURCE_ID_EXCHANGE_DEPOSITION-1),SOURCE_ID_EXCHANGE_VAPORIZATION-1), &
        deposition_rates,exchange_source,transport_stochastic=.FALSE.,pore_enabled=.FALSE., &
        terrain_exposed=.TRUE.,atmosphere_exposed=.TRUE.,static_deposition_allowed=.FALSE.)
    CALL require(ALL(deposition_rates%solid_deposition_rate==0.0_wp) .AND. &
        ALL(exchange_source==0.0_wp),'L2 cannot drive deposition-coupled vaporization')
    CALL build_physical_state_from_qp(qp(:,2),check1,accepted=ok)
    CALL require(ok,'equivalent L1 physical state for deposition test')
    CALL eval_mass_exchange_explicit(rheology_l2,mass_l2,check1,0,[0.0_wp,0.0_wp],dt, &
        IBSET(0,SOURCE_ID_EXCHANGE_DEPOSITION-1),deposition_rates,exchange_source, &
        transport_stochastic=.FALSE.,pore_enabled=.FALSE.,terrain_exposed=.TRUE., &
        atmosphere_exposed=.TRUE.,static_deposition_allowed=.TRUE.)
    CALL require(ANY(deposition_rates%solid_deposition_rate>0.0_wp), &
        'L1 is allowed to deposit after settling creates it')
    CALL eval_mass_exchange_explicit(rheology_l2,mass_l2,check1,0,[0.0_wp,0.0_wp],dt, &
        IBSET(0,SOURCE_ID_EXCHANGE_DEPOSITION-1),deposition_rates,exchange_source, &
        transport_stochastic=.FALSE.,pore_enabled=.FALSE.,terrain_exposed=.FALSE., &
        atmosphere_exposed=.TRUE.,static_deposition_allowed=.TRUE.)
    CALL require(ALL(deposition_rates%solid_deposition_rate==0.0_wp), &
        'L1 static deposition also requires physical terrain contact')

    ! Group zero uses the immutable pre-spatial snapshot, but applies and
    ! limits on the spatial candidate. A later group reads that candidate.
    config%settling_flag=.FALSE.
    CALL make_column(1.0_wp,0.0_wp,20.0_wp,0.0_wp,before,snapshot)
    CALL make_column(0.7_wp,0.0_wp,0.0_wp,0.0_wp,q,qp)
    primitives=qp
    before=q
    CALL apply_interlayer_lie_group(config,rheology_l2,mass_l2,snapshot,.TRUE.,dt,q,qp,scratch,status)
    CALL require(status==0,'group zero evaluates pre-spatial source snapshot')
    CALL require(q(1,2)>0.0_wp,'group zero lifts solids despite stationary spatial recipient')
    CALL check_conservation(before,q,'group-zero snapshot/recipient separation')
    after_group0=q
    CALL make_column(0.7_wp,0.0_wp,0.0_wp,0.0_wp,q,qp)
    CALL apply_interlayer_lie_group(config,rheology_l2,mass_l2,snapshot,.FALSE.,dt,q,qp,scratch,status)
    CALL require(status==0,'later Lie group accepts stationary donor')
    after_later=q
    CALL require(ALL(after_later==before),'later group reads stationary candidate, not group zero snapshot')
    CALL require(ANY(after_group0/=after_later),'group-zero and later-group evaluation differ')
    CALL check_admissibility(primitives)

    ! An external/spatial step can exhaust the lower layer while group zero
    ! still sees its wet pre-spatial donor; the inventory limiter must veto
    ! the requested flux instead of allowing a transient negative donor.
    CALL make_column(0.0_wp,0.0_wp,0.0_wp,0.0_wp,q,qp)
    before=q; primitives=qp
    CALL apply_interlayer_lie_group(config,rheology_l2,mass_l2,snapshot,.TRUE.,dt,q,qp,scratch,status)
    CALL require(status==0,'group zero with exhausted post-spatial donor remains valid')
    CALL require(ALL(q==before) .AND. ALL(qp==primitives), &
        'pre-spatial mass cannot be spent after post-spatial depletion')

    ! A recipient violating the gross donor budget must reject BOTH layers,
    ! preserving the caller's state exactly (no partial commit).
    CALL make_column(1.0_wp,0.0_wp,20.0_wp,0.0_wp,q,qp)
    snapshot=qp
    q(idx_solid_first,1)=0.0_wp
    before=q; primitives=qp
    CALL apply_interlayer_lie_group(config,rheology_l2,mass_l2,snapshot,.TRUE.,dt,q,qp,scratch,status)
    CALL require(status==20,'negative candidate species mass rejects the entire pair')
    CALL require(ALL(q==before) .AND. ALL(qp==primitives),'paired failure does not commit either layer')

    ! Sequence the actual paired physical source and analytic Chezy drag on
    ! the same accepted candidates. Mass exchange may change both inertias;
    ! the following drag must read those updated masses and velocities.
    config%settling_flag=.TRUE.
    config%elutriation_flag=.TRUE.
    config%rouse_critical=100.0_wp
    CALL make_column(1.0_wp,1.0_wp,20.0_wp,2.0_wp,q,qp)
    qp(idx_T,1)=420.0_wp
    qp(idx_T,2)=275.0_wp
    DO i=1,2
        CALL primitive_to_conservative(qp(:,i),q(:,i),accepted=ok)
        CALL require(ok,'thermally distinct donor candidate conversion')
    END DO
    before=q
    snapshot=qp
    CALL apply_interlayer_lie_group(config,rheology_l2,mass_l2,snapshot,.FALSE.,dt,q,qp,scratch,status)
    CALL require(status==0,'two-way physical mass group accepted before drag')
    CALL check_conservation(before,q,'two-way thermal mass exchange')
    CALL check_admissibility(qp)
    after_transfer=q
    slip_before=ABS(qp(idx_u,2)-qp(idx_u,1))
    drag%model='CHEZY'
    drag%chezy_drag_coeff=0.03_wp
    CALL apply_interlayer_momentum_exchange_analytic(drag,dt,q,qp,ok)
    CALL require(ok,'drag after mass group accepts updated layer inertia')
    CALL check_conservation(before,q,'two-way mass and subsequent analytic drag')
    CALL require(ALL(q(1,:)==after_transfer(1,:)), &
        'drag leaves post-transfer layer masses unchanged')
    CALL require(ALL(q(4:6,:)==after_transfer(4:6,:)), &
        'drag leaves thermal energy and solid-species inventories unchanged')
    slip_after=ABS(qp(idx_u,2)-qp(idx_u,1))
    CALL require(slip_after<=slip_before,'drag relaxes relative velocity after mass transfer')
    CALL check_admissibility(qp)

    CALL test_external_and_interlayer_budget()

    ! Optional transported scalars belong to their configured layer even
    ! when the interlayer transaction changes mass and recomputes primitives.
    ! No stochastic/Darcy state is silently transferred with solid or carrier.
    BLOCK
        TYPE(physical_state_real) :: optional_physical(4)
        REAL(wp) :: optional_q(8,2),optional_qp(8,2),optional_before(8,2),optional_snapshot(8,2)
        INTEGER :: layer

        n_stoch_vars=1
        n_pore_vars=1
        n_vars=8
        n_eqns=8
        CALL configure_canonical_state_indices()
        optional_q=0.0_wp
        optional_qp=0.0_wp
        DO layer=1,2
            optional_qp(idx_h,layer)=1.0_wp
            optional_qp(idx_T,layer)=300.0_wp
            optional_qp(idx_u,layer)=REAL(10*layer,wp)
            optional_qp(idx_solid_first:idx_solid_last,layer)=[0.6_wp,0.2_wp]
        END DO
        optional_qp(idx_pore,1)=0.15_wp
        optional_qp(idx_stoch,2)=0.35_wp
        DO layer=1,2
            CALL primitive_to_conservative(optional_qp(:,layer),optional_q(:,layer),accepted=ok)
            CALL require(ok,'optional layer candidate conversion')
        END DO
        optional_before=optional_q
        optional_snapshot=optional_qp
        DO layer=1,4
            CALL initialize_physical_state_real(optional_physical(layer))
        END DO
        CALL apply_interlayer_lie_group(config,rheology_l2,mass_l2,optional_snapshot,.FALSE.,dt, &
            optional_q,optional_qp,optional_physical,status)
        CALL require(status==0,'paired transfer with independent optional layers accepted')
        CALL require(ALL(optional_q(idx_stoch,:)==optional_before(idx_stoch,:)), &
            'interlayer transfer does not invent or move stochastic conservative scalar')
        CALL require(ALL(optional_q(idx_pore,:)==optional_before(idx_pore,:)), &
            'interlayer transfer does not invent or move Darcy conservative scalar')
        CALL require(optional_q(idx_pore,2)==0.0_wp .AND. optional_q(idx_stoch,1)==0.0_wp, &
            'inactive optional layer components remain exactly zero')
        DO layer=1,8
            CALL require(ABS(SUM(optional_q(layer,:))-SUM(optional_before(layer,:)))<= &
                1.0E-10_wp*MAX(1.0_wp,ABS(SUM(optional_before(layer,:)))), &
                'optional pair preserves every conservative component')
        END DO
    END BLOCK

    PRINT *, 'PASS paired Lie transaction four wet/dry cases, group-zero snapshot, conservation, atomic rejection'

CONTAINS

    !> @brief Close mass, species, momentum and thermal budgets across Lie group-4 sources.
    !> @details Reproduces the existing L1 deposition, L2 air entrainment,
    !>          then dynamic two-way interlayer sequence using independent
    !>          conservative candidates. Internal transfers must cancel.
    SUBROUTINE test_external_and_interlayer_budget()
        TYPE(mass_exchange_config_type) :: basal,atmospheric
        TYPE(interlayer_mass_exchange_config_type) :: pair
        TYPE(mass_exchange_rates_type) :: bed_rates,air_rates
        REAL(wp) :: mobile(6,2),primitive(6,2),initial(6,2),after_external(6,2)
        REAL(wp) :: source_bed(6),source_air(6),snapshot(6,2)
        REAL(wp) :: bed_mass,bed_solids,air_mass,before_gas,after_gas,external_gas
        REAL(wp) :: residual,scale
        LOGICAL :: accepted
        INTEGER :: component,status

        CALL make_column(1.0_wp,1.0_wp,20.0_wp,2.0_wp,mobile,primitive)
        initial=mobile
        snapshot=primitive
        basal%deposition%settling_flag=.TRUE.
        atmospheric%entrainment%flag=.TRUE.
        pair%settling_flag=.TRUE.
        pair%elutriation_flag=.TRUE.
        pair%elutriation_solid_volume_fraction=-1.0_wp
        pair%rouse_critical=100.0_wp

        CALL build_physical_state_from_qp(primitive(:,1),check1,accepted=accepted)
        CALL require(accepted,'L1 deposition donor state')
        CALL eval_mass_exchange_explicit(rheology_l2,basal,check1,0,[0.0_wp,0.0_wp],dt, &
            IBSET(0,SOURCE_ID_EXCHANGE_DEPOSITION-1),bed_rates,source_bed, &
            transport_stochastic=.FALSE.,pore_enabled=.FALSE.,terrain_exposed=.TRUE., &
            atmosphere_exposed=.FALSE.,static_deposition_allowed=.TRUE.)
        CALL require(ANY(bed_rates%solid_deposition_rate>0.0_wp), &
            'active basal solid deposition before interlayer exchange')
        mobile(:,1)=mobile(:,1)+dt*source_bed
        CALL conservative_to_primitive(mobile(:,1),primitive(:,1),accepted=accepted)
        CALL require(accepted,'L1 after substrate deposition')

        CALL build_physical_state_from_qp(primitive(:,2),check2,accepted=accepted)
        CALL require(accepted,'L2 entrainment donor state')
        ! Select a positive-Ri, finite-speed physical fixture to exercise
        ! atmospheric entrainment without changing the production law.
        check2%Ri=1.0_wp
        check2%mod_vel2=4.0_wp
        CALL eval_mass_exchange_explicit(rheology_l2,atmospheric,check2,0,[0.0_wp,0.0_wp],dt, &
            IBSET(0,SOURCE_ID_EXCHANGE_ENTRAINMENT-1),air_rates,source_air, &
            transport_stochastic=.FALSE.,pore_enabled=.FALSE.,terrain_exposed=.FALSE., &
            atmosphere_exposed=.TRUE.,static_deposition_allowed=.FALSE.)
        CALL require(source_air(idx_totMassEqn)>0.0_wp,'ambient air adds mass')
        CALL require(ALL(source_air(idx_solidEqn_first:idx_solidEqn_last)==0.0_wp), &
            'air entrainment carries no solid species')
        CALL require(source_air(idx_uEqn)==0.0_wp .AND. source_air(idx_vEqn)==0.0_wp, &
            'air entrainment carries no horizontal momentum')
        mobile(:,2)=mobile(:,2)+dt*source_air
        CALL conservative_to_primitive(mobile(:,2),primitive(:,2),accepted=accepted)
        CALL require(accepted,'L2 after atmospheric entrainment')

        after_external=mobile
        CALL apply_interlayer_lie_group(pair,rheology_l2,atmospheric,snapshot,.FALSE., &
            dt,mobile,primitive,scratch,status)
        CALL require(status==0,'interlayer source accepts external recipients')
        CALL require(ANY(mobile(:,1)/=after_external(:,1)), &
            'both-way interlayer source updates L1')
        CALL check_conservation(after_external,mobile,'mass exchange after external sources')
        CALL check_admissibility(primitive)

        ! Each component in the combined mobile layers must change only by
        ! the independently calculated terrain and atmospheric sources.
        DO component=1,n_eqns
            residual=SUM(mobile(component,:))-SUM(initial(component,:))- &
                dt*(source_bed(component)+source_air(component))
            scale=MAX(1.0_wp,ABS(SUM(initial(component,:))), &
                ABS(dt*(source_bed(component)+source_air(component))))
            CALL require(ABS(residual)<=1.0E-10_wp*scale, &
                'combined external-only conservative source budget')
        END DO

        bed_mass=-dt*source_bed(idx_totMassEqn)
        bed_solids=dt*SUM(rho_s*bed_rates%solid_deposition_rate)
        air_mass=dt*source_air(idx_totMassEqn)
        CALL require(bed_mass>0.0_wp .AND. bed_solids>0.0_wp .AND. air_mass>0.0_wp, &
            'terrain and atmospheric mass sources both active')
        residual=SUM(initial(idx_solidEqn_first:idx_solidEqn_last,:))- &
            SUM(mobile(idx_solidEqn_first:idx_solidEqn_last,:))-bed_solids
        CALL require(ABS(residual)<=1.0E-10_wp*MAX(1.0_wp,bed_solids), &
            'mobile and deposited substrate solid mass closure')
        residual=SUM(mobile(idx_totMassEqn,:))-SUM(initial(idx_totMassEqn,:))+ &
            bed_mass-air_mass
        CALL require(ABS(residual)<=1.0E-10_wp*MAX(1.0_wp,bed_mass,air_mass), &
            'total mass balance with terrain and atmosphere')

        before_gas=SUM(initial(idx_totMassEqn,:))- &
            SUM(initial(idx_solidEqn_first:idx_solidEqn_last,:))
        after_gas=SUM(mobile(idx_totMassEqn,:))- &
            SUM(mobile(idx_solidEqn_first:idx_solidEqn_last,:))
        external_gas=dt*(source_bed(idx_totMassEqn)- &
            SUM(source_bed(idx_solidEqn_first:idx_solidEqn_last))+source_air(idx_totMassEqn))
        CALL require(ABS(after_gas-before_gas-external_gas)<= &
            1.0E-10_wp*MAX(1.0_wp,ABS(external_gas)), &
            'carrier mass changes only through external sources')
        CALL require(ABS(SUM(mobile(idx_engyEqn,:))-SUM(initial(idx_engyEqn,:))- &
            dt*(source_bed(idx_engyEqn)+source_air(idx_engyEqn)))<= &
            1.0E-10_wp*MAX(1.0_wp,ABS(SUM(initial(idx_engyEqn,:)))), &
            'sensible thermal budget closes after source sequence')
    END SUBROUTINE test_external_and_interlayer_budget

    !> @brief Build consistent conservative and primitive layer states for a column.
    SUBROUTINE make_column(h1,h2,u1,u2,conservative,primitive)
        REAL(wp), INTENT(IN) :: h1,h2,u1,u2
        REAL(wp), INTENT(OUT) :: conservative(6,2),primitive(6,2)
        REAL(wp) :: h(2),velocity(2)
        INTEGER :: k
        h=[h1,h2]; velocity=[u1,u2]
        primitive=0.0_wp
        DO k=1,2
            primitive(idx_h,k)=h(k)
            primitive(idx_T,k)=300.0_wp
            primitive(idx_u,k)=velocity(k)
            primitive(idx_v,k)=0.0_wp
            IF (h(k)>0.0_wp) primitive(idx_solid_first:idx_solid_last,k)=[0.6_wp,0.2_wp]
            CALL primitive_to_conservative(primitive(:,k),conservative(:,k),accepted=ok)
            CALL require(ok,'fixture conversion to conservative state')
        END DO
    END SUBROUTINE make_column

    !> @brief Check paired conservation for every transported conservative component.
    SUBROUTINE check_conservation(original,candidate,description)
        REAL(wp), INTENT(IN) :: original(6,2),candidate(6,2)
        CHARACTER(*), INTENT(IN) :: description
        INTEGER :: k
        DO k=1,6
            CALL require(ABS(SUM(candidate(k,:))-SUM(original(k,:)))<= &
                1.0E-10_wp*MAX(1.0_wp,ABS(SUM(original(k,:)))), &
                'component conservation: '//description)
        END DO
    END SUBROUTINE check_conservation

    !> @brief Confirm post-transaction mixture packing and finite physical states.
    SUBROUTINE check_admissibility(primitive)
        REAL(wp), INTENT(IN) :: primitive(6,2)
        CALL build_physical_state_from_qp(primitive(:,1),check1,accepted=ok)
        CALL require(ok,'reconstruct final physical L1')
        CALL build_physical_state_from_qp(primitive(:,2),check2,accepted=ok)
        CALL require(ok,'reconstruct final physical L2')
        CALL require(SUM(check1%alphas)<=maximum_solid_packing+1.0E-10_wp,'L1 final packing')
        CALL require(SUM(check2%alphas)<=maximum_solid_packing+1.0E-10_wp,'L2 final packing')
    END SUBROUTINE check_admissibility

    !> @brief Stop the numerical regression on a failed physical invariant.
    SUBROUTINE require(condition,description)
        LOGICAL, INTENT(IN) :: condition
        CHARACTER(*), INTENT(IN) :: description
        IF (.NOT.condition) THEN
            PRINT *, 'FAIL interlayer Lie: ',TRIM(description)
            ERROR STOP 1
        END IF
    END SUBROUTINE require
END PROGRAM test_interlayer_lie_transaction
