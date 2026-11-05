PROGRAM test_state_new
    USE parameters_2d
    USE constitutive_parameters_2d
    USE state_new_types
    USE state_new_conversions
    USE state_new_physical
    IMPLICIT NONE

    REAL(wp), ALLOCATABLE :: qp(:), qp_roundtrip(:), q(:)
    COMPLEX(wp), ALLOCATABLE :: cqp(:), cq(:)
    TYPE(state_workspace_real) :: real_workspace
    TYPE(state_workspace_complex) :: complex_workspace
    LOGICAL :: accepted

    CALL initialize_case
    ALLOCATE(qp(n_vars),qp_roundtrip(n_vars),q(n_vars),cqp(n_vars),cq(n_vars))
    qp=0.0_wp
    qp(idx_h)=2.0_wp; qp(idx_T)=300.0_wp; qp(idx_u)=1.5_wp; qp(idx_v)=-0.25_wp
    qp(idx_solid_first:idx_solid_last)=[0.15_wp,0.10_wp]
    qp(idx_add_gas_first:idx_add_gas_last)=0.05_wp
    qp(idx_stoch)=0.2_wp; qp(idx_pore)=1000.0_wp

    CALL initialize_state_workspace_real(real_workspace)
    CALL initialize_state_workspace_complex(complex_workspace)
    CALL primitive_to_conservative(qp,q)
    CALL conservative_to_primitive(q,qp_roundtrip)
    CALL assert_close('round-trip h',qp_roundtrip(idx_h),qp(idx_h),1.0E-11_wp)
    CALL assert_close('round-trip T',qp_roundtrip(idx_T),qp(idx_T),1.0E-11_wp)
    CALL assert_close('round-trip u',qp_roundtrip(idx_u),qp(idx_u),1.0E-11_wp)
    CALL check_admissibility(qp,q)

    CALL build_physical_state_from_qp(qp,real_workspace%state,0.1_wp,-0.2_wp,0.0_wp,0.0_wp,0.0_wp)
    cqp=CMPLX(qp,0.0_wp,wp)
    CALL build_physical_state_from_qp(cqp,complex_workspace%state,0.1_wp,-0.2_wp,0.0_wp,0.0_wp,0.0_wp)
    CALL assert_close('REAL/COMPLEX density',real_workspace%state%rho_m,REAL(complex_workspace%state%rho_m,wp),1.0E-11_wp)
    CALL assert_close('REAL/COMPLEX Richardson',real_workspace%state%Ri,REAL(complex_workspace%state%Ri,wp),1.0E-11_wp)
    CALL assert_close('REAL/COMPLEX kinematics',real_workspace%state%mod_vel,REAL(complex_workspace%state%mod_vel,wp),1.0E-11_wp)

    ! A velocity change must refresh both kinematics and Richardson without
    ! recomputing mixture properties, for both REAL and COMPLEX workspaces.
    real_workspace%state%u=3.0_wp
    real_workspace%state%v=4.0_wp
    complex_workspace%state%u=CMPLX(3.0_wp,0.0_wp,wp)
    complex_workspace%state%v=CMPLX(4.0_wp,0.0_wp,wp)
    CALL update_physical_state_kinematics(real_workspace%state,0.1_wp,-0.2_wp,0.0_wp,0.0_wp,0.0_wp)
    CALL update_physical_state_kinematics(complex_workspace%state,0.1_wp,-0.2_wp,0.0_wp,0.0_wp,0.0_wp)
    CALL assert_close('velocity update Richardson',real_workspace%state%Ri, &
        real_workspace%state%red_grav*real_workspace%state%h/25.0_wp,1.0E-11_wp)
    CALL assert_close('velocity update REAL/COMPLEX Richardson',real_workspace%state%Ri, &
        REAL(complex_workspace%state%Ri,wp),1.0E-11_wp)
    CALL assert_close('velocity update REAL/COMPLEX speed',real_workspace%state%mod_vel, &
        REAL(complex_workspace%state%mod_vel,wp),1.0E-11_wp)
    CALL update_physical_state_kinematics(real_workspace%state,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp)
    IF (real_workspace%state%layout_valid .EQV. .FALSE.) ERROR STOP 'workspace layout lost on velocity update'

    liquid_flag=.TRUE.; gas_flag=.FALSE.
    CALL build_physical_state_from_qp(qp,real_workspace%state,accepted=accepted)
    IF (.NOT.accepted) ERROR STOP 'valid liquid/solid state rejected'
    CALL assert_close('liquid carrier heat capacity',real_workspace%state%cp_c,sp_heat_l,1.0E-12_wp)
    liquid_flag=.FALSE.; gas_flag=.TRUE.

    cqp=CMPLX(qp,0.0_wp,wp)
    cqp(idx_T)=CMPLX(qp(idx_T),1.0E-12_wp,wp)
    CALL build_physical_state_from_qp(cqp,complex_workspace%state,accepted=accepted)
    IF (.NOT.accepted) ERROR STOP 'valid complex state rejected'
    IF (AIMAG(complex_workspace%state%T)==0.0_wp) ERROR STOP 'complex temperature perturbation lost'

    qp(idx_h)=1.0E-12_wp
    qp(idx_u)=100.0_wp; qp(idx_v)=-50.0_wp
    CALL build_physical_state_from_qp(qp,real_workspace%state,accepted=accepted)
    IF (.NOT.accepted) ERROR STOP 'valid near-dry state rejected'
    IF (.NOT.(real_workspace%state%mod_vel>=0.0_wp)) ERROR STOP 'near-dry velocity invalid'

    qp=0.0_wp
    qp(idx_h)=2.0_wp; qp(idx_T)=300.0_wp
    qp(idx_solid_first)=1.2_wp
    CALL build_physical_state_from_qp(qp,real_workspace%state,accepted=accepted)
    IF (accepted) ERROR STOP 'invalid physical state accepted'

    qp=0.0_wp
    CALL build_physical_state_from_qp(qp,real_workspace%state)
    IF (real_workspace%state%Ri /= 0.0_wp) ERROR STOP 'dry Richardson state is not zero'
    WRITE(*,'(A)') 'PASS: state_new kernels verified'

CONTAINS

    SUBROUTINE check_admissibility(reference_qp,reference_q)
        REAL(wp), INTENT(IN) :: reference_qp(n_vars), reference_q(n_vars)
        REAL(wp) :: bad_qp(n_vars), bad_q(n_vars), p_dyn
        LOGICAL :: accepted
        bad_qp=reference_qp
        bad_qp(idx_solid_first)=-1.0E-3_wp
        CALL primitive_to_conservative(bad_qp,bad_q,accepted)
        IF (accepted) ERROR STOP 'negative composition accepted'
        bad_qp=reference_qp
        bad_qp(idx_solid_first:idx_solid_last)=0.8_wp
        bad_qp(idx_add_gas_first:idx_add_gas_last)=0.4_wp
        CALL primitive_to_conservative(bad_qp,bad_q,accepted)
        IF (accepted) ERROR STOP 'overfilled composition accepted'
        bad_q=reference_q
        bad_q(idx_solid_first)= -1.0E-3_wp
        CALL conservative_to_primitive(bad_q,bad_qp,accepted,p_dyn)
        IF (accepted) ERROR STOP 'negative conservative component accepted'
        CALL conservative_to_primitive(reference_q,bad_qp,accepted,p_dyn)
        IF (.NOT.accepted) ERROR STOP 'valid conservative state rejected'
        IF (p_dyn<=0.0_wp) ERROR STOP 'dynamic pressure was not computed'
    END SUBROUTINE check_admissibility

    SUBROUTINE initialize_case
        n_solid=2; n_add_gas=1; n_stoch_vars=1; n_pore_vars=1
        n_vars=4+n_solid+n_add_gas+n_stoch_vars+n_pore_vars
        n_eqns=n_vars
        CALL configure_canonical_state_indices
        ALLOCATE(rho_s(n_solid),inv_rho_s(n_solid),sp_heat_s(n_solid))
        ALLOCATE(sp_heat_g(n_add_gas),sp_gas_const_g(n_add_gas))
        rho_s=[2500.0_wp,1800.0_wp]; inv_rho_s=1.0_wp/rho_s; sp_heat_s=[900.0_wp,1100.0_wp]
        sp_heat_a=1005.0_wp; sp_gas_const_a=287.05_wp; sp_heat_g=[1850.0_wp]; sp_gas_const_g=[461.5_wp]
        sp_heat_l=4180.0_wp; rho_l=1000.0_wp; inv_rho_l=1.0_wp/rho_l
        pres=101325.0_wp; inv_pres=1.0_wp/pres; T_ambient=300.0_wp
        rho_a_amb=pres/(sp_gas_const_a*T_ambient); grav=9.81_wp
        gas_flag=.TRUE.; liquid_flag=.FALSE.; stoch_transport_flag=.TRUE.; pore_pressure_flag=.TRUE.
        slope_correction_flag=.TRUE.; curvature_term_flag=.TRUE.
    END SUBROUTINE initialize_case

    SUBROUTINE assert_close(label,a,b,tol)
        CHARACTER(*), INTENT(IN) :: label
        REAL(wp), INTENT(IN) :: a,b,tol
        IF (ABS(a-b)>tol*MAX(1.0_wp,ABS(a),ABS(b))) ERROR STOP label
    END SUBROUTINE assert_close
END PROGRAM test_state_new
