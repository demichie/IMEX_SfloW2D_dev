PROGRAM test_state_equivalence
    USE parameters_2d
    USE constitutive_parameters_2d
    USE state_conversion_2d, ONLY : old_qc_to_qp=>qc_to_qp, old_qp_to_qc=>qp_to_qc, &
        old_state_qp=>physical_state_from_qp, old_state_cq=>physical_state_from_q, &
        physical_state_type, complex_physical_state_type, old_mixture=>mixt_var, &
        old_from_volume=>eval_mixture_properties_from_volume_fractions, &
        old_primitive_volume=>primitive_to_volume_fractions
    USE equation_terms_2d, ONLY : eval_hydrostatic_coefficient, eval_local_speeds_direction
    USE state, ONLY : conservative_to_primitive, primitive_to_conservative, &
        new_state_qp=>build_physical_state_from_qp, new_state_q=>build_physical_state_from_q, &
        evaluate_mixture_density, evaluate_mixture_heat_capacity, &
        evaluate_mixture_from_volume, primitive_to_volume_fractions, &
        physical_state_real, physical_state_complex, &
        initialize_physical_state_real, initialize_physical_state_complex
    IMPLICIT NONE
    REAL(wp), ALLOCATABLE :: q(:),qp(:),q_old(:),qp_old(:),qp_new(:)
    COMPLEX(wp), ALLOCATABLE :: cq(:)
    TYPE(physical_state_type) :: old_real
    TYPE(complex_physical_state_type) :: old_complex
    TYPE(physical_state_real) :: new_real
    TYPE(physical_state_complex) :: new_complex
    LOGICAL :: ok

    CALL initialize_case
    ALLOCATE(q(n_vars),qp(n_vars),q_old(n_vars),qp_old(n_vars),qp_new(n_vars),cq(n_vars))
    CALL initialize_physical_state_real(new_real)
    CALL initialize_physical_state_complex(new_complex)
    CALL compare_case(2.0_wp,300.0_wp,1.5_wp,-0.25_wp,.15_wp,.10_wp,.05_wp,'wet gas-solid')
    liquid_flag=.TRUE.; gas_flag=.FALSE.
    CALL compare_case(2.0_wp,300.0_wp,1.5_wp,-0.25_wp,.15_wp,.10_wp,0.0_wp,'wet liquid-solid')
    liquid_flag=.FALSE.; gas_flag=.TRUE.
    CALL compare_case(0.0_wp,300.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,'dry')
    CALL compare_case(1.0E-12_wp,300.0_wp,100.0_wp,-50.0_wp,.15_wp,.10_wp,.05_wp,'near-dry')
    ! Unified strict-policy trial: a substantially overfilled trace is
    ! rejected rather than silently normalized.
    CALL make_qp(1.0_wp,300.0_wp,1.5_wp,-0.25_wp,0.8_wp,0.8_wp,0.2_wp,qp)
    CALL primitive_to_conservative(qp,q,ok)
    IF (ok) ERROR STOP 'overfilled numerical face accepted under strict trial'
    CALL make_qp(1.0_wp,300.0_wp,1.5_wp,-0.25_wp, &
        -32.0_wp*EPSILON(1.0_wp),0.15_wp,0.05_wp,qp)
    CALL check_spatial_face_properties(qp,'roundoff face')
    CALL check_dry_decoder_limits
    CALL check_settling_trace_equivalence
    CALL check_volume_conversion_equivalence
    WRITE(*,'(A)') 'PASS: state equivalence cases verified'
CONTAINS
    SUBROUTINE compare_case(h,T,u,v,xs1,xs2,xg1,label)
        REAL(wp), INTENT(IN) :: h,T,u,v,xs1,xs2,xg1
        CHARACTER(*), INTENT(IN) :: label
        REAL(wp) :: pd_old,pd_new
        REAL(wp) :: material_qp(n_vars),old_ri,old_rho,old_carrier,old_grav
        REAL(wp) :: old_cp_c,old_cp_m,new_rho,new_inv_rho,new_cp_c,new_cp_m,new_grav
        CALL make_qp(h,T,u,v,xs1,xs2,xg1,qp)
        CALL old_qp_to_qc(qp,q_old); CALL primitive_to_conservative(qp,q)
        CALL old_qc_to_qp(q_old,qp_old,pd_old); CALL conservative_to_primitive(q,qp_new,ok,pd_new)
        IF (.NOT.ok) ERROR STOP 'new conversion rejected valid state'
        CALL assert_scalar(label//' dynamic pressure',pd_new,pd_old,1.0E-9_wp)
        CALL assert_vector(label//' q',q,q_old,1.0E-10_wp)
        CALL assert_vector(label//' qp',qp_new,qp_old,1.0E-9_wp)
        CALL old_state_qp(qp,old_real); CALL new_state_qp(qp,new_real,accepted=ok)
        IF (.NOT.ok) ERROR STOP 'new REAL state rejected valid state'
        CALL assert_scalar(label//' density',new_real%rho_m,old_real%rho_m,1.0E-9_wp)
        CALL assert_scalar(label//' gravity',new_real%red_grav,old_real%red_grav,1.0E-9_wp)
        CALL assert_scalar(label//' Richardson',new_real%Ri,old_real%Ri,1.0E-8_wp)
        cq=CMPLX(q_old,0.0_wp,wp); CALL old_state_cq(cq,old_complex); CALL new_state_q(cq,new_complex,accepted=ok)
        IF (.NOT.ok) ERROR STOP 'new COMPLEX state rejected valid state'
        CALL assert_scalar(label//' complex density',REAL(new_complex%rho_m,wp),REAL(old_complex%rho_m,wp),1.0E-9_wp)
        ! The source material has its own composition, independently of the recipient cell.
        material_qp=qp
        material_qp(idx_h)=1.0_wp
        material_qp(idx_u)=0.0_wp; material_qp(idx_v)=0.0_wp
        CALL check_spatial_face_properties(qp,label)
        CALL old_mixture(material_qp,old_ri,old_rho,old_carrier,old_grav,old_cp_c,old_cp_m)
        CALL evaluate_mixture_density(material_qp(idx_T), &
            material_qp(idx_solid_first:idx_solid_last), &
            material_qp(idx_add_gas_first:idx_add_gas_last),new_rho,new_inv_rho)
        CALL evaluate_mixture_heat_capacity(material_qp(idx_solid_first:idx_solid_last), &
            material_qp(idx_add_gas_first:idx_add_gas_last),new_cp_m,new_cp_c)
        new_grav=(new_rho-rho_a_amb)*new_inv_rho*grav
        CALL assert_scalar(label//' injection rho',new_rho,old_rho,1.0E-9_wp)
        CALL assert_scalar(label//' injection cp',new_cp_m,old_cp_m,1.0E-9_wp)
        CALL assert_scalar(label//' injection red grav',new_grav,old_grav,1.0E-9_wp)

    END SUBROUTINE compare_case
    SUBROUTINE check_spatial_face_properties(face,label)
        REAL(wp), INTENT(IN) :: face(n_vars)
        CHARACTER(*), INTENT(IN) :: label
        REAL(wp) :: old_ri,old_rho,old_rho_c,old_grav,old_cp_c,old_cp_m
        REAL(wp) :: new_grav,gamma,minimum_speed,maximum_speed
        REAL(wp) :: expected_min,expected_max,celerity,normal_speed
        REAL(wp) :: reference_q(n_vars),new_q(n_vars)
        LOGICAL :: face_ok
        INTEGER :: direction,velocity_index

        CALL old_qp_to_qc(face,reference_q)
        CALL primitive_to_conservative(face,new_q,face_ok)
        IF (.NOT.face_ok) ERROR STOP 'valid reconstructed face conversion rejected'
        CALL assert_vector(label//' face conservative',new_q,reference_q,1.0E-10_wp)
        CALL old_mixture(face,old_ri,old_rho,old_rho_c,old_grav,old_cp_c,old_cp_m)
        CALL eval_hydrostatic_coefficient(face,new_grav,gamma)
        CALL assert_scalar(label//' face gravity',new_grav,old_grav,1.0E-10_wp)
        CALL assert_scalar(label//' face gamma',gamma,old_rho*old_grav,1.0E-10_wp)
        DO direction=1,2
            IF (direction==1) THEN
                velocity_index=idx_u
            ELSE
                velocity_index=idx_v
            END IF
            normal_speed=face(velocity_index)
            IF (old_grav*face(idx_h)<0.0_wp) THEN
                expected_min=normal_speed
                expected_max=normal_speed
            ELSE
                celerity=SQRT(old_grav*face(idx_h))
                expected_min=normal_speed-celerity
                expected_max=normal_speed+celerity
            END IF
            CALL eval_local_speeds_direction(face,1.0_wp,velocity_index,minimum_speed,maximum_speed)
            CALL assert_scalar(label//' face speed min',minimum_speed,expected_min,1.0E-10_wp)
            CALL assert_scalar(label//' face speed max',maximum_speed,expected_max,1.0E-10_wp)
        END DO
    END SUBROUTINE check_spatial_face_properties

    !> Regression: F5 source-state closure must retain tiny positive solids.
    !> The historical constructor never zeroed a positive mass fraction,
    !> even when its magnitude was within the negative-roundoff tolerance.
    SUBROUTINE check_settling_trace_equivalence
        REAL(wp), PARAMETER :: tiny_fraction=0.25_wp*4096.0_wp*EPSILON(1.0_wp)
        REAL(wp), PARAMETER :: velocity=0.01_wp,dt=1.0_wp
        REAL(wp) :: old_rate,new_rate
        COMPLEX(wp) :: ctrace(n_vars)
        INTEGER :: mode
        LOGICAL :: accepted

        DO mode=1,3
            SELECT CASE(mode)
            CASE(1)
                CALL make_qp(2.0_wp,300.0_wp,0.0_wp,0.0_wp, &
                    tiny_fraction,0.0_wp,0.0_wp,qp)
            CASE(2)
                CALL make_qp(2.0_wp,300.0_wp,0.0_wp,0.0_wp, &
                    0.0_wp,tiny_fraction,0.0_wp,qp)
            CASE(3)
                CALL make_qp(2.0_wp,300.0_wp,0.0_wp,0.0_wp, &
                    tiny_fraction,0.0_wp,tiny_fraction,qp)
            END SELECT
            CALL old_state_qp(qp,old_real)
            CALL new_state_qp(qp,new_real,accepted=accepted)
            IF (.NOT.accepted) ERROR STOP 'tiny-positive source state rejected'
            IF (ANY(new_real%xs/=old_real%xs)) &
                ERROR STOP 'F5 tiny-positive solid fraction differs from legacy'
            IF (ANY(new_real%xg/=old_real%xg)) &
                ERROR STOP 'F5 tiny-positive gas fraction differs from legacy'
            IF (ANY(new_real%alphas/=old_real%alphas)) &
                ERROR STOP 'F5 settling volume fraction differs from legacy'
            IF (ANY(qp(idx_solid_first:idx_solid_last)>0.0_wp) .AND. &
                .NOT.ANY(new_real%alphas>0.0_wp)) &
                ERROR STOP 'F5 source closure removed positive solid volume'
            old_rate=MIN(old_real%alphas(1)*velocity, &
                old_real%h*old_real%alphas(1)/dt)
            new_rate=MIN(new_real%alphas(1)*velocity, &
                new_real%h*new_real%alphas(1)/dt)
            IF (new_rate/=old_rate) &
                ERROR STOP 'F5 near-depleted settling input rate differs'
            ctrace=CMPLX(qp,0.0_wp,wp)
            ctrace(idx_solid_first)=ctrace(idx_solid_first)+ &
                CMPLX(0.0_wp,1.0E-20_wp,wp)
            CALL new_state_qp(ctrace,new_complex,accepted=accepted)
            IF (.NOT.accepted) ERROR STOP 'tiny-positive complex source rejected'
            IF (AIMAG(new_complex%xs(1))/=1.0E-20_wp) &
                ERROR STOP 'tiny-positive complex-step perturbation lost'
        END DO
    END SUBROUTINE check_settling_trace_equivalence

    !> Test the single machine-epsilon dry boundary; there is no separate
    !> numerical reconstruction or near-dry fallback path.
    SUBROUTINE check_dry_decoder_limits
        REAL(wp) :: decoded(n_vars),mass,p_dyn
        LOGICAL :: accepted
        INTEGER :: case_id
        DO case_id=1,3
            mass=0.5_wp*EPSILON(1.0_wp)
            IF (case_id==2) mass=EPSILON(1.0_wp)
            IF (case_id==3) mass=2.0_wp*EPSILON(1.0_wp)
            q=0.0_wp
            q(idx_totMassEqn)=mass
            q(idx_engyEqn)=mass*sp_heat_a*900.0_wp
            q(idx_uEqn)=mass*20.0_wp
            q(idx_vEqn)=-mass*5.0_wp
            CALL conservative_to_primitive(q,decoded,accepted,p_dyn)
            IF (.NOT.accepted) ERROR STOP 'epsilon-boundary decoder rejected valid mass'
            IF (mass<=EPSILON(1.0_wp)) THEN
                CALL assert_scalar('epsilon-dry temperature',decoded(idx_T),T_ambient,0.0_wp)
                CALL assert_scalar('epsilon-dry depth',decoded(idx_h),0.0_wp,0.0_wp)
                CALL assert_scalar('epsilon-dry pressure',p_dyn,0.0_wp,0.0_wp)
            ELSE
                CALL assert_scalar('positive-mass temperature',decoded(idx_T),900.0_wp,1.0E-12_wp)
            END IF
        END DO
        q(idx_totMassEqn)=1.0_wp
        q(idx_solid_first)=-1.0E-3_wp
        CALL conservative_to_primitive(q,decoded,accepted,p_dyn)
        IF (accepted) ERROR STOP 'decoder accepted materially negative solids'
    END SUBROUTINE check_dry_decoder_limits

    !> I/O and initialization use the inverse volume->mass composition
    !> independently of the conservative/primitive numerical conversions.
    SUBROUTINE check_volume_conversion_equivalence
        REAL(wp) :: volume_old(n_solid),volume_new(n_solid)
        REAL(wp) :: gas_volume(n_add_gas),xs_old(n_solid),xs_new(n_solid)
        REAL(wp) :: xg_old(n_add_gas),xg_new(n_add_gas)
        REAL(wp) :: rho_old,inv_old,carrier_old,xc_old,heat_c_old,heat_m_old
        REAL(wp) :: rho_new,inv_new,carrier_new,xc_new,heat_c_new,heat_m_new
        REAL(wp) :: out_solids_old(n_solid),out_solids_new(n_solid)
        REAL(wp) :: out_gas_old(n_add_gas),out_gas_new(n_add_gas)
        INTEGER :: configuration

        DO configuration=1,2
            gas_flag=configuration==1
            liquid_flag=.NOT.gas_flag
            volume_old=[0.02_wp,0.03_wp]
            volume_new=volume_old
            gas_volume=0.0_wp
            IF (gas_flag) gas_volume=0.01_wp
            CALL old_from_volume(350.0_wp,gas_volume,volume_old,rho_old,inv_old, &
                carrier_old,xs_old,xg_old,xc_old,heat_c_old,heat_m_old)
            CALL evaluate_mixture_from_volume(350.0_wp,gas_volume,volume_new, &
                rho_new,inv_new,carrier_new,xs_new,xg_new,xc_new,heat_c_new,heat_m_new)
            CALL assert_vector('inverse volume solids',volume_new,volume_old,1.0E-12_wp)
            CALL assert_vector('inverse volume mass solids',xs_new,xs_old,1.0E-12_wp)
            CALL assert_vector('inverse volume mass gas',xg_new,xg_old,1.0E-12_wp)
            CALL assert_scalar('inverse volume density',rho_new,rho_old,1.0E-12_wp)
            CALL assert_scalar('inverse volume inverse density',inv_new,inv_old,1.0E-12_wp)
            CALL assert_scalar('inverse volume carrier',carrier_new,carrier_old,1.0E-12_wp)
            CALL assert_scalar('inverse volume heat carrier',heat_c_new,heat_c_old,1.0E-12_wp)
            CALL assert_scalar('inverse volume heat mixture',heat_m_new,heat_m_old,1.0E-12_wp)
            CALL make_qp(0.25_wp,350.0_wp,1.0_wp,0.0_wp, &
                xs_new(1),xs_new(2),SUM(xg_new),qp)
            CALL old_primitive_volume(qp,out_solids_old,out_gas_old)
            CALL primitive_to_volume_fractions(qp,out_solids_new,out_gas_new)
            CALL assert_vector('primitive solids volume',out_solids_new,out_solids_old,1.0E-12_wp)
            CALL assert_vector('primitive gas volume',out_gas_new,out_gas_old,1.0E-12_wp)
        END DO
        gas_flag=.TRUE.
        liquid_flag=.FALSE.
    END SUBROUTINE check_volume_conversion_equivalence

    SUBROUTINE make_qp(h,T,u,v,xs1,xs2,xg1,out)
        REAL(wp), INTENT(IN) :: h,T,u,v,xs1,xs2,xg1
        REAL(wp), INTENT(OUT) :: out(n_vars)
        out=0.0_wp; out(idx_h)=h; out(idx_T)=T; out(idx_u)=u; out(idx_v)=v
        out(idx_solid_first:idx_solid_last)=[xs1,xs2]; out(idx_add_gas_first:idx_add_gas_last)=xg1
        out(idx_stoch)=0.2_wp; out(idx_pore)=1000.0_wp
    END SUBROUTINE make_qp
    SUBROUTINE initialize_case
        n_solid=2; n_add_gas=1; n_stoch_vars=1; n_pore_vars=1; n_vars=9; n_eqns=n_vars
        CALL configure_canonical_state_indices
        ALLOCATE(rho_s(2),inv_rho_s(2),sp_heat_s(2),sp_heat_g(1),sp_gas_const_g(1))
        rho_s=[2500.0_wp,1800.0_wp]; inv_rho_s=1.0_wp/rho_s; sp_heat_s=[900.0_wp,1100.0_wp]
        sp_heat_a=1005.0_wp; sp_gas_const_a=287.05_wp; sp_heat_g=[1850.0_wp]; sp_gas_const_g=[461.5_wp]
        sp_heat_l=4180.0_wp; rho_l=1000.0_wp; inv_rho_l=1.0_wp/rho_l; pres=101325.0_wp; inv_pres=1.0_wp/pres
        T_ambient=300.0_wp; rho_a_amb=pres/(sp_gas_const_a*T_ambient); grav=9.81_wp; eps_sing=1.0E-8_wp; eps_sing4=eps_sing**4
        gas_flag=.TRUE.; liquid_flag=.FALSE.; stoch_transport_flag=.TRUE.
        pore_pressure_flag=.TRUE.; slope_correction_flag=.FALSE.; curvature_term_flag=.FALSE.
    END SUBROUTINE initialize_case
    SUBROUTINE assert_scalar(label,a,b,tol)
        CHARACTER(*), INTENT(IN) :: label; REAL(wp), INTENT(IN) :: a,b,tol
        IF (ABS(a-b)>tol*MAX(1.0_wp,ABS(a),ABS(b))) ERROR STOP label
    END SUBROUTINE assert_scalar
    SUBROUTINE assert_vector(label,a,b,tol)
        CHARACTER(*), INTENT(IN) :: label; REAL(wp), INTENT(IN) :: a(:),b(:),tol
        IF (ANY(ABS(a-b)>tol*MAX(1.0_wp,MAXVAL(ABS(a)),MAXVAL(ABS(b))))) ERROR STOP label
    END SUBROUTINE assert_vector
END PROGRAM test_state_equivalence
