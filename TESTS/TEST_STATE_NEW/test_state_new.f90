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

    CALL build_physical_state_from_qp(qp,real_workspace%state,0.1_wp,-0.2_wp,0.0_wp,0.0_wp,0.0_wp)
    cqp=CMPLX(qp,0.0_wp,wp)
    CALL build_physical_state_from_qp(cqp,complex_workspace%state,0.1_wp,-0.2_wp,0.0_wp,0.0_wp,0.0_wp)
    CALL assert_close('REAL/COMPLEX density',real_workspace%state%rho_m,REAL(complex_workspace%state%rho_m,wp),1.0E-11_wp)
    CALL assert_close('REAL/COMPLEX Richardson',real_workspace%state%Ri,REAL(complex_workspace%state%Ri,wp),1.0E-11_wp)
    CALL assert_close('REAL/COMPLEX kinematics',real_workspace%state%mod_vel,REAL(complex_workspace%state%mod_vel,wp),1.0E-11_wp)

    qp=0.0_wp
    CALL build_physical_state_from_qp(qp,real_workspace%state)
    IF (real_workspace%state%Ri /= 0.0_wp) ERROR STOP 'dry Richardson state is not zero'
    WRITE(*,'(A)') 'PASS: state_new kernels verified'

CONTAINS

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
