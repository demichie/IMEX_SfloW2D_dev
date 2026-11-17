!> @brief Exercise the exact paired Lie transaction used by the production executor.
PROGRAM test_interlayer_lie_transaction
    USE parameters_2d, ONLY : wp,n_layers,n_solid,n_add_gas,n_vars,n_eqns,n_stoch_vars,n_pore_vars
    USE parameters_2d, ONLY : gas_flag,sutherland_flag,configure_canonical_state_indices,idx_h,idx_T,idx_u,idx_v
    USE parameters_2d, ONLY : idx_solid_first,idx_solid_last
    USE constitutive_parameters_2d, ONLY : sp_heat_a,sp_gas_const_a,sp_heat_g,sp_gas_const_g
    USE constitutive_parameters_2d, ONLY : pres,inv_pres,T_ambient,grav,kin_visc_c,rho_a_amb,vonK
    USE solid_transport_parameters_2d, ONLY : rho_s,inv_rho_s,sp_heat_s,diam_s,maximum_solid_packing
    USE state_types, ONLY : physical_state_real,initialize_physical_state_real
    USE state_physical, ONLY : build_physical_state_from_qp
    USE state_conversions, ONLY : primitive_to_conservative
    USE interlayer_exchange_2d, ONLY : apply_interlayer_lie_group
    USE interlayer_exchange_parameters_2d, ONLY : interlayer_exchange_config_type
    USE mass_exchange_parameters_2d, ONLY : mass_exchange_config_type
    USE rheology_parameters_2d, ONLY : rheology_config_type
    IMPLICIT NONE

    TYPE(physical_state_real) :: scratch(4),check1,check2
    TYPE(interlayer_exchange_config_type) :: config
    TYPE(mass_exchange_config_type) :: mass_l2
    TYPE(rheology_config_type) :: rheology_l2
    REAL(wp) :: q(6,2),qp(6,2),before(6,2),primitives(6,2),snapshot(6,2)
    REAL(wp) :: after_group0(6,2),after_later(6,2)
    REAL(wp) :: dt
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

    CALL make_column(1.0_wp,1.0_wp,20.0_wp,2.0_wp,q,qp)
    before=q; snapshot=qp
    CALL apply_interlayer_lie_group(config,rheology_l2,mass_l2,snapshot,.FALSE.,dt,q,qp,scratch,status)
    CALL require(status==0,'both wet layers exchange conservatively')
    CALL require(ANY(q(:,1)/=before(:,1)),'both wet layers change state')
    CALL check_conservation(before,q,'wet L1 and L2')
    CALL check_admissibility(qp)

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

    PRINT *, 'PASS paired Lie transaction four wet/dry cases, group-zero snapshot, conservation, atomic rejection'

CONTAINS

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
