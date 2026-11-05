PROGRAM test_state_equivalence
    USE parameters_2d
    USE constitutive_parameters_2d
    USE state_conversion_2d, ONLY : old_qc_to_qp=>qc_to_qp, old_qp_to_qc=>qp_to_qc, &
        old_state_qp=>physical_state_from_qp, old_state_cq=>physical_state_from_q, &
        physical_state_type, complex_physical_state_type
    USE state_new, ONLY : conservative_to_primitive, primitive_to_conservative, &
        new_state_qp=>build_physical_state_from_qp, new_state_q=>build_physical_state_from_q, &
        state_workspace_real, state_workspace_complex, &
        initialize_state_workspace_real, initialize_state_workspace_complex
    IMPLICIT NONE
    REAL(wp), ALLOCATABLE :: q(:),qp(:),q_old(:),qp_old(:),qp_new(:)
    COMPLEX(wp), ALLOCATABLE :: cq(:)
    TYPE(physical_state_type) :: old_real
    TYPE(complex_physical_state_type) :: old_complex
    TYPE(state_workspace_real) :: new_real
    TYPE(state_workspace_complex) :: new_complex
    LOGICAL :: ok

    CALL initialize_case
    ALLOCATE(q(n_vars),qp(n_vars),q_old(n_vars),qp_old(n_vars),qp_new(n_vars),cq(n_vars))
    CALL initialize_state_workspace_real(new_real)
    CALL initialize_state_workspace_complex(new_complex)
    CALL compare_case(2.0_wp,300.0_wp,1.5_wp,-0.25_wp,.15_wp,.10_wp,.05_wp,'wet gas-solid')
    liquid_flag=.TRUE.; gas_flag=.FALSE.
    CALL compare_case(2.0_wp,300.0_wp,1.5_wp,-0.25_wp,.15_wp,.10_wp,0.0_wp,'wet liquid-solid')
    liquid_flag=.FALSE.; gas_flag=.TRUE.
    CALL compare_case(0.0_wp,300.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,'dry')
    CALL compare_case(1.0E-12_wp,300.0_wp,100.0_wp,-50.0_wp,.15_wp,.10_wp,.05_wp,'near-dry')
    WRITE(*,'(A)') 'PASS: state_new equivalence cases verified'
CONTAINS
    SUBROUTINE compare_case(h,T,u,v,xs1,xs2,xg1,label)
        REAL(wp), INTENT(IN) :: h,T,u,v,xs1,xs2,xg1
        CHARACTER(*), INTENT(IN) :: label
        REAL(wp) :: pd_old,pd_new
        CALL make_qp(h,T,u,v,xs1,xs2,xg1,qp)
        CALL old_qp_to_qc(qp,q_old); CALL primitive_to_conservative(qp,q)
        CALL old_qc_to_qp(q_old,qp_old,pd_old); CALL conservative_to_primitive(q,qp_new,ok,pd_new)
        IF (.NOT.ok) ERROR STOP 'new conversion rejected valid state'
        CALL assert_vector(label//' q',q,q_old,1.0E-10_wp)
        CALL assert_vector(label//' qp',qp_new,qp_old,1.0E-9_wp)
        CALL old_state_qp(qp,old_real); CALL new_state_qp(qp,new_real%state,accepted=ok)
        IF (.NOT.ok) ERROR STOP 'new REAL state rejected valid state'
        CALL assert_scalar(label//' density',new_real%state%rho_m,old_real%rho_m,1.0E-9_wp)
        CALL assert_scalar(label//' gravity',new_real%state%red_grav,old_real%red_grav,1.0E-9_wp)
        CALL assert_scalar(label//' Richardson',new_real%state%Ri,old_real%Ri,1.0E-8_wp)
        cq=CMPLX(q_old,0.0_wp,wp); CALL old_state_cq(cq,old_complex); CALL new_state_q(cq,new_complex%state,accepted=ok)
        IF (.NOT.ok) ERROR STOP 'new COMPLEX state rejected valid state'
        CALL assert_scalar(label//' complex density',REAL(new_complex%state%rho_m,wp),REAL(old_complex%rho_m,wp),1.0E-9_wp)
    END SUBROUTINE compare_case
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
