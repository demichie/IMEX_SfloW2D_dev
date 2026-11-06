PROGRAM test_state
    USE parameters_2d
    USE constitutive_parameters_2d
    USE state_types
    USE state_mixture, ONLY : sauter_diameter, average_density_solids, &
        evaluate_mixture_density, evaluate_mixture_phase_properties, evaluate_reduced_gravity, &
        evaluate_gas_specific_volume, evaluate_mixture_from_volume, evaluate_mixture
    USE state_composition, ONLY : close_composition
    USE state_conversions
    USE state_physical
    IMPLICIT NONE

    REAL(wp), ALLOCATABLE :: qp(:), qp_roundtrip(:), q(:)
    COMPLEX(wp), ALLOCATABLE :: cqp(:), cq(:)
    TYPE(physical_state_real) :: real_workspace
    TYPE(physical_state_complex) :: complex_workspace
    LOGICAL :: accepted

    CALL initialize_case
    CALL check_particle_descriptors
    CALL check_shared_density_kernel
    CALL check_shared_conservative_math
    CALL check_gas_eos
    CALL check_inverse_volume_conversion
    CALL check_complex_closure_limits
    CALL check_real_complex_invalid_candidates
    CALL check_dry_refresh_consistency
    CALL check_canonical_dry_kinematics
    CALL check_epsilon_dry_real_complex
    CALL check_nonpositive_heat_capacity
    ALLOCATE(qp(n_vars),qp_roundtrip(n_vars),q(n_vars),cqp(n_vars),cq(n_vars))
    qp=0.0_wp
    qp(idx_h)=2.0_wp; qp(idx_T)=300.0_wp; qp(idx_u)=1.5_wp; qp(idx_v)=-0.25_wp
    qp(idx_solid_first:idx_solid_last)=[0.15_wp,0.10_wp]
    qp(idx_add_gas_first:idx_add_gas_last)=0.05_wp
    qp(idx_stoch)=0.2_wp; qp(idx_pore)=1000.0_wp

    CALL initialize_physical_state_real(real_workspace)
    CALL initialize_physical_state_complex(complex_workspace)
    CALL primitive_to_conservative(qp,q)
    CALL conservative_to_primitive(q,qp_roundtrip)
    CALL assert_close('round-trip h',qp_roundtrip(idx_h),qp(idx_h),1.0E-11_wp)
    CALL assert_close('round-trip T',qp_roundtrip(idx_T),qp(idx_T),1.0E-11_wp)
    CALL assert_close('round-trip u',qp_roundtrip(idx_u),qp(idx_u),1.0E-11_wp)
    CALL check_admissibility(qp,q)

    CALL build_physical_state_from_qp(qp,real_workspace,0.1_wp,-0.2_wp,0.0_wp,0.0_wp,0.0_wp)
    cqp=CMPLX(qp,0.0_wp,wp)
    CALL build_physical_state_from_qp(cqp,complex_workspace,0.1_wp,-0.2_wp,0.0_wp,0.0_wp,0.0_wp)
    CALL assert_close('REAL/COMPLEX density',real_workspace%rho_m,REAL(complex_workspace%rho_m,wp),1.0E-11_wp)
    CALL assert_close('REAL/COMPLEX Richardson',real_workspace%Ri,REAL(complex_workspace%Ri,wp),1.0E-11_wp)
    CALL assert_close('REAL/COMPLEX kinematics',real_workspace%mod_vel,REAL(complex_workspace%mod_vel,wp),1.0E-11_wp)

    ! A velocity change must refresh both kinematics and Richardson without
    ! recomputing mixture properties, for both REAL and COMPLEX workspaces.
    real_workspace%u=3.0_wp
    real_workspace%v=4.0_wp
    complex_workspace%u=CMPLX(3.0_wp,0.0_wp,wp)
    complex_workspace%v=CMPLX(4.0_wp,0.0_wp,wp)
    CALL update_physical_state_kinematics(real_workspace,0.1_wp,-0.2_wp,0.0_wp,0.0_wp,0.0_wp)
    CALL update_physical_state_kinematics(complex_workspace,0.1_wp,-0.2_wp,0.0_wp,0.0_wp,0.0_wp)
    CALL assert_close('velocity update Richardson',real_workspace%Ri, &
        real_workspace%red_grav*real_workspace%h/25.0_wp,1.0E-11_wp)
    CALL assert_close('velocity update REAL/COMPLEX Richardson',real_workspace%Ri, &
        REAL(complex_workspace%Ri,wp),1.0E-11_wp)
    CALL assert_close('velocity update REAL/COMPLEX speed',real_workspace%mod_vel, &
        REAL(complex_workspace%mod_vel,wp),1.0E-11_wp)
    CALL update_physical_state_kinematics(real_workspace,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp)
    IF (real_workspace%layout_valid .EQV. .FALSE.) ERROR STOP 'workspace layout lost on velocity update'

    liquid_flag=.TRUE.; gas_flag=.FALSE.
    CALL build_physical_state_from_qp(qp,real_workspace,accepted=accepted)
    IF (.NOT.accepted) ERROR STOP 'valid liquid/solid state rejected'
    CALL assert_close('liquid carrier heat capacity',real_workspace%cp_c,sp_heat_l,1.0E-12_wp)
    liquid_flag=.FALSE.; gas_flag=.TRUE.

    cqp=CMPLX(qp,0.0_wp,wp)
    cqp(idx_T)=CMPLX(qp(idx_T),1.0E-12_wp,wp)
    CALL build_physical_state_from_qp(cqp,complex_workspace,accepted=accepted)
    IF (.NOT.accepted) ERROR STOP 'valid complex state rejected'
    IF (AIMAG(complex_workspace%T)==0.0_wp) ERROR STOP 'complex temperature perturbation lost'

    qp(idx_h)=1.0E-12_wp
    qp(idx_u)=100.0_wp; qp(idx_v)=-50.0_wp
    CALL build_physical_state_from_qp(qp,real_workspace,accepted=accepted)
    IF (.NOT.accepted) ERROR STOP 'valid near-dry state rejected'
    IF (.NOT.(real_workspace%mod_vel>=0.0_wp)) ERROR STOP 'near-dry velocity invalid'

    qp=0.0_wp
    qp(idx_h)=2.0_wp; qp(idx_T)=300.0_wp
    qp(idx_solid_first)=1.2_wp
    CALL build_physical_state_from_qp(qp,real_workspace,accepted=accepted)
    IF (accepted) ERROR STOP 'invalid physical state accepted'

    qp=0.0_wp
    CALL build_physical_state_from_qp(qp,real_workspace)
    IF (real_workspace%Ri /= 0.0_wp) ERROR STOP 'dry Richardson state is not zero'
    WRITE(*,'(A)') 'PASS: state kernels verified'

CONTAINS

    !> Both APIs must use the same density closure for REAL and COMPLEX,
    !> including the gas/liquid switch and the zero-carrier limit.
    SUBROUTINE check_shared_density_kernel
        REAL(wp) :: xs(2),xg(1),rho,inv_rho
        REAL(wp) :: phase_rho,phase_inv_rho,rho_c,inv_rho_c
        REAL(wp) :: alphas(2),alphag(1)
        COMPLEX(wp) :: cxs(2),cxg(1),crho,cinv_rho
        COMPLEX(wp) :: cphase_rho,cphase_inv_rho,crho_c,cinv_rho_c
        COMPLEX(wp) :: calphas(2),calphag(1),temperature
        INTEGER :: carrier_case,composition_case

        DO carrier_case=1,2
            gas_flag=carrier_case==1
            liquid_flag=.NOT.gas_flag
            DO composition_case=1,3
                xs=[0.15_wp,0.10_wp]
                xg=[0.05_wp]
                IF (composition_case==2) xs=[0.0_wp,0.0_wp]
                IF (composition_case==3) THEN
                    xs=[0.6_wp,0.4_wp]
                    xg=[0.0_wp]
                END IF
                CALL evaluate_mixture_density(350.0_wp,xs,xg,rho,inv_rho)
                CALL evaluate_mixture_phase_properties(350.0_wp,xs,xg, &
                    phase_rho,phase_inv_rho,rho_c,inv_rho_c,alphas,alphag)
                CALL assert_close('shared REAL mixture rho',rho,phase_rho,0.0_wp)
                CALL assert_close('shared REAL mixture inverse rho',inv_rho,phase_inv_rho,0.0_wp)

                cxs=CMPLX(xs,0.0_wp,wp)
                cxg=CMPLX(xg,0.0_wp,wp)
                temperature=CMPLX(350.0_wp,1.0E-20_wp,wp)
                CALL evaluate_mixture_density(temperature,cxs,cxg,crho,cinv_rho)
                CALL evaluate_mixture_phase_properties(temperature,cxs,cxg, &
                    cphase_rho,cphase_inv_rho,crho_c,cinv_rho_c,calphas,calphag)
                IF (crho/=cphase_rho .OR. cinv_rho/=cphase_inv_rho) &
                    ERROR STOP 'shared COMPLEX density kernel mismatch'
            END DO
        END DO
        gas_flag=.TRUE.
        liquid_flag=.FALSE.
    END SUBROUTINE check_shared_density_kernel

    !> Two frontends consume the same thermodynamic and momentum formulas.
    SUBROUTINE check_shared_conservative_math
        REAL(wp) :: qp_case(n_vars),q_case(n_vars),qp_check(n_vars), gravity_value
        COMPLEX(wp) :: cq_case(n_vars), cqp_check(n_vars)
        TYPE(physical_state_real) :: physical
        TYPE(physical_state_complex) :: cphysical
        REAL(wp) :: tolerance
        INTEGER :: sample
        LOGICAL :: ok

        CALL initialize_physical_state_real(physical)
        CALL initialize_physical_state_complex(cphysical)
        tolerance=5.0E-12_wp
        DO sample=1,3
            qp_case=0.0_wp
            qp_case(idx_h)=2.0_wp
            qp_case(idx_T)=350.0_wp
            qp_case(idx_u)=2.1_wp
            qp_case(idx_v)=-0.7_wp
            qp_case(idx_solid_first:idx_solid_last)=[0.2_wp,0.1_wp]
            qp_case(idx_add_gas_first:idx_add_gas_last)=0.05_wp
            IF (sample==2) THEN
                qp_case(idx_h)=1.0E-12_wp
                qp_case(idx_u)=0.02_wp
            ELSE IF (sample==3) THEN
                qp_case(idx_h)=1.0E-8_wp
            END IF
            CALL primitive_to_conservative(qp_case,q_case,ok)
            IF (.NOT.ok) ERROR STOP 'shared math: failed to create q'
            CALL conservative_to_primitive(q_case,qp_check,ok)
            IF (.NOT.ok) ERROR STOP 'shared math: q->qp rejected valid q'
            CALL build_physical_state_from_q(q_case,physical,accepted=ok)
            IF (.NOT.ok) ERROR STOP 'shared math: q->state rejected valid q'
            CALL assert_close('shared q->qp T',qp_check(idx_T),physical%T,tolerance)
            CALL assert_close('shared q->qp u',qp_check(idx_u),physical%u,tolerance)
            CALL assert_close('shared q->qp v',qp_check(idx_v),physical%v,tolerance)
            CALL assert_close('shared q->qp h',qp_check(idx_h),physical%h,tolerance)
            gravity_value=evaluate_reduced_gravity(physical%rho_m,physical%inv_rhom)
            CALL assert_close('shared reduced gravity',gravity_value,physical%red_grav,0.0_wp)

            cq_case=CMPLX(q_case,0.0_wp,wp)
            CALL conservative_to_primitive(cq_case,cqp_check,ok)
            IF (.NOT.ok) ERROR STOP 'shared complex q->qp rejected'
            CALL build_physical_state_from_q(cq_case,cphysical,accepted=ok)
            IF (.NOT.ok) ERROR STOP 'shared complex q->state rejected'
            CALL assert_close('shared complex T',REAL(cqp_check(idx_T),wp),REAL(cphysical%T,wp),tolerance)
            CALL assert_close('shared complex u',REAL(cqp_check(idx_u),wp),REAL(cphysical%u,wp),tolerance)
            CALL assert_close('shared complex v',REAL(cqp_check(idx_v),wp),REAL(cphysical%v,wp),tolerance)
        END DO
    END SUBROUTINE check_shared_conservative_math

    !> Gas law must be consistent for constituent and carrier evaluations.
    SUBROUTINE check_gas_eos
        REAL(wp) :: T,expected,actual
        COMPLEX(wp) :: complex_temperature,complex_specific_volume
        T=325.0_wp
        actual=evaluate_gas_specific_volume(T,sp_gas_const_a)
        expected=sp_gas_const_a*T*inv_pres
        CALL assert_close('ideal-gas specific volume',actual,expected,0.0_wp)
        CALL assert_close('ideal-gas density reciprocal',1.0_wp/actual, &
            pres/(sp_gas_const_a*T),5.0E-15_wp)
        complex_temperature=CMPLX(T,1.0E-24_wp,wp)
        complex_specific_volume=evaluate_gas_specific_volume(complex_temperature,sp_gas_const_a)
        CALL assert_close('ideal-gas complex derivative', &
            AIMAG(complex_specific_volume)/1.0E-24_wp,sp_gas_const_a*inv_pres,5.0E-15_wp)
    END SUBROUTINE check_gas_eos

    !> Initial and source materials are supplied as volume fractions. The
    !> inverse conversion must agree with the canonical mass-based EOS.
    SUBROUTINE check_inverse_volume_conversion
        REAL(wp) :: alphas(2),alphag(1),given_alphas(2),given_alphag(1)
        REAL(wp) :: xs(2),xg(1),xc,rho_m,inv_rhom,rho_c,cp_m,cp_c
        REAL(wp) :: rho_ref,inv_rhom_ref,rho_c_ref,inv_rho_c_ref
        REAL(wp) :: cp_m_ref,cp_c_ref,alphas_ref(2),alphag_ref(1)
        REAL(wp) :: expected_density
        INTEGER :: carrier_case,composition_case
        DO carrier_case=1,2
            gas_flag=carrier_case==1
            liquid_flag=.NOT.gas_flag
            DO composition_case=1,4
                SELECT CASE(composition_case)
                CASE(1)
                    given_alphas=[0.08_wp,0.12_wp]
                    given_alphag=[0.03_wp]
                CASE(2)
                    given_alphas=0.0_wp
                    given_alphag=[0.02_wp]
                CASE(3)
                    given_alphas=[0.45_wp,0.55_wp]
                    given_alphag=0.0_wp
                CASE(4)
                    ! Preserve the historical normalization of excess solids.
                    given_alphas=[0.75_wp,0.75_wp]
                    given_alphag=0.0_wp
                END SELECT
                alphas=given_alphas
                alphag=given_alphag
                CALL evaluate_mixture_from_volume(350.0_wp,alphag,alphas, &
                    rho_m,inv_rhom,rho_c,xs,xg,xc,cp_c,cp_m)
                CALL evaluate_mixture(350.0_wp,xs,xg,rho_ref,inv_rhom_ref, &
                    rho_c_ref,inv_rho_c_ref,cp_m_ref,cp_c_ref,alphas_ref,alphag_ref)
                CALL assert_close('volume inverse mixture density',rho_m,rho_ref,0.0_wp)
                CALL assert_close('volume inverse inverse density',inv_rhom,inv_rhom_ref,0.0_wp)
                CALL assert_close('volume inverse carrier density',rho_c,rho_c_ref,0.0_wp)
                CALL assert_close('volume inverse cp_m',cp_m,cp_m_ref,0.0_wp)
                CALL assert_close('volume inverse cp_c',cp_c,cp_c_ref,0.0_wp)
                CALL assert_close('volume inverse xc',xc,1.0_wp-SUM(xs),0.0_wp)
                IF (composition_case==4) THEN
                    CALL assert_close('volume inverse normalized solid 1',alphas(1),0.5_wp,0.0_wp)
                    CALL assert_close('volume inverse normalized solid 2',alphas(2),0.5_wp,0.0_wp)
                END IF
                IF (gas_flag) THEN
                    expected_density=SUM(alphas*rho_s)+ &
                        (1.0_wp-SUM(alphas)-SUM(alphag))* &
                        (1.0_wp/evaluate_gas_specific_volume(350.0_wp,sp_gas_const_a))+ &
                        SUM(alphag*(1.0_wp/evaluate_gas_specific_volume(350.0_wp,sp_gas_const_g)))
                ELSE
                    expected_density=SUM(alphas*rho_s)+(1.0_wp-SUM(alphas))*rho_l
                END IF
                CALL assert_close('volume inverse mass budget',rho_m,expected_density,1.0E-12_wp)
                CALL assert_close('volume inverse solid mass 1',xs(1)*rho_m,alphas(1)*rho_s(1),1.0E-12_wp)
                CALL assert_close('volume inverse solid mass 2',xs(2)*rho_m,alphas(2)*rho_s(2),1.0E-12_wp)
                IF (gas_flag) THEN
                    CALL assert_close('volume inverse gas mass',xg(1)*rho_m, &
                        alphag(1)/evaluate_gas_specific_volume(350.0_wp,sp_gas_const_g(1)),1.0E-12_wp)
                END IF
            END DO
        END DO
        gas_flag=.TRUE.
        liquid_flag=.FALSE.
    END SUBROUTINE check_inverse_volume_conversion

    !> Complex closure must reproduce REAL clipping and its local derivative.
    SUBROUTINE check_complex_closure_limits
        REAL(wp) :: xs(2),xg(1),base,total,expected_derivative
        COMPLEX(wp) :: cxs(2),cxg(1)
        REAL(wp), PARAMETER :: delta=1.0E-25_wp
        LOGICAL :: real_ok,complex_ok

        xs=[-1.0E-14_wp,0.2_wp]; xg=0.1_wp
        cxs=CMPLX(xs,0.0_wp,wp); cxg=CMPLX(xg,0.0_wp,wp)
        cxs(1)=CMPLX(xs(1),delta,wp)
        real_ok=close_composition(xs,xg)
        complex_ok=close_composition(cxs,cxg)
        IF (.NOT.real_ok .OR. .NOT.complex_ok) ERROR STOP 'roundoff closure rejected'
        CALL assert_close('clipped REAL species',xs(1),0.0_wp,0.0_wp)
        CALL assert_close('clipped COMPLEX species',REAL(cxs(1),wp),0.0_wp,0.0_wp)
        CALL assert_close('clipped derivative',AIMAG(cxs(1)),0.0_wp,0.0_wp)

        xs=[0.7_wp,0.3_wp+1.0E-13_wp]; xg=0.0_wp
        base=xs(1); total=SUM(xs)+SUM(xg)
        cxs=CMPLX(xs,0.0_wp,wp); cxg=CMPLX(xg,0.0_wp,wp)
        cxs(1)=CMPLX(xs(1),delta,wp)
        real_ok=close_composition(xs,xg)
        complex_ok=close_composition(cxs,cxg)
        IF (.NOT.real_ok .OR. .NOT.complex_ok) ERROR STOP 'tolerated overfill rejected'
        CALL assert_close('overfill REAL/COMPLEX baseline',REAL(cxs(1),wp),xs(1),1.0E-14_wp)
        expected_derivative=(total-base)/(total*total)
        CALL assert_close('complex overfill normalization derivative', &
            AIMAG(cxs(1))/delta,expected_derivative,1.0E-12_wp)

        xs=[0.2_wp,0.3_wp]; xg=0.0_wp
        cxs=CMPLX(xs,0.0_wp,wp); cxg=CMPLX(xg,0.0_wp,wp)
        cxs(1)=CMPLX(xs(1),delta,wp)
        IF (.NOT.close_composition(cxs,cxg)) ERROR STOP 'positive complex fraction rejected'
        CALL assert_close('positive fraction derivative retained',AIMAG(cxs(1))/delta,1.0_wp,1.0E-14_wp)

        xs=[0.8_wp,0.5_wp]; xg=0.0_wp
        cxs=CMPLX(xs,0.0_wp,wp); cxg=CMPLX(xg,0.0_wp,wp)
        IF (close_composition(xs,xg)) ERROR STOP 'overfilled REAL composition accepted'
        IF (close_composition(cxs,cxg)) ERROR STOP 'overfilled COMPLEX composition accepted'
    END SUBROUTINE check_complex_closure_limits

    SUBROUTINE check_real_complex_invalid_candidates
        REAL(wp) :: qp_case(n_vars),q_case(n_vars),qp_result(n_vars),temperature
        COMPLEX(wp) :: cq_case(n_vars),cqp_result(n_vars),ctemperature
        TYPE(physical_state_real) :: rs
        TYPE(physical_state_complex) :: cs
        LOGICAL :: ok
        CALL initialize_physical_state_real(rs)
        CALL initialize_physical_state_complex(cs)
        qp_case=0.0_wp
        qp_case(idx_h)=1.0_wp; qp_case(idx_T)=350.0_wp
        qp_case(idx_solid_first:idx_solid_last)=[0.2_wp,0.1_wp]
        qp_case(idx_add_gas_first:idx_add_gas_last)=0.05_wp
        CALL primitive_to_conservative(qp_case,q_case,ok)
        IF (.NOT.ok) ERROR STOP 'invalid candidate fixture creation failed'
        q_case(idx_solid_first)=1.3_wp*q_case(idx_totMassEqn)
        cq_case=CMPLX(q_case,0.0_wp,wp)

        CALL conservative_to_primitive(q_case,qp_result,ok)
        IF (ok) ERROR STOP 'invalid REAL q->qp composition accepted'
        CALL conservative_to_primitive(cq_case,cqp_result,ok)
        IF (ok) ERROR STOP 'invalid COMPLEX q->qp composition accepted'
        CALL build_physical_state_from_q(q_case,rs,accepted=ok)
        IF (ok) ERROR STOP 'invalid REAL q->state composition accepted'
        CALL build_physical_state_from_q(cq_case,cs,accepted=ok)
        IF (ok) ERROR STOP 'invalid COMPLEX q->state composition accepted'
        CALL evaluate_temperature_from_q(q_case,temperature,ok)
        IF (ok) ERROR STOP 'invalid REAL q temperature unexpectedly accepted'
        CALL evaluate_temperature_from_q(cq_case,ctemperature,ok)
        IF (ok) ERROR STOP 'invalid COMPLEX q temperature unexpectedly accepted'
        CALL assert_close('invalid real temperature sentinel',temperature,0.0_wp,0.0_wp)
        CALL assert_close('invalid complex temperature sentinel',REAL(ctemperature,wp),0.0_wp,0.0_wp)
    END SUBROUTINE check_real_complex_invalid_candidates

    SUBROUTINE check_dry_refresh_consistency
        REAL(wp) :: qp_case(n_vars),before_rho,expected_cp
        COMPLEX(wp) :: cqp_case(n_vars)
        TYPE(physical_state_real) :: rs
        TYPE(physical_state_complex) :: cs
        INTEGER :: carrier_case
        LOGICAL :: ok
        CALL initialize_physical_state_real(rs)
        CALL initialize_physical_state_complex(cs)
        DO carrier_case=1,2
            gas_flag=carrier_case==1; liquid_flag=.NOT.gas_flag
            expected_cp=sp_heat_a
            IF (liquid_flag) expected_cp=sp_heat_l
            qp_case=0.0_wp; qp_case(idx_h)=0.0_wp
            qp_case(idx_T)=350.0_wp
            qp_case(idx_solid_first:idx_solid_last)=[0.2_wp,0.1_wp]
            qp_case(idx_add_gas_first:idx_add_gas_last)=0.05_wp
            cqp_case=CMPLX(qp_case,0.0_wp,wp)
            CALL build_physical_state_from_qp(qp_case,rs,accepted=ok)
            IF (.NOT.ok) ERROR STOP 'dry REAL state rejected'
            CALL build_physical_state_from_qp(cqp_case,cs,accepted=ok)
            IF (.NOT.ok) ERROR STOP 'dry COMPLEX state rejected'
            before_rho=rs%rho_m
            rs%T=500.0_wp; cs%T=CMPLX(500.0_wp,1.0E-20_wp,wp)
            CALL refresh_physical_state(rs,accepted=ok)
            IF (.NOT.ok) ERROR STOP 'dry REAL refresh rejected'
            CALL refresh_physical_state(cs,accepted=ok)
            IF (.NOT.ok) ERROR STOP 'dry COMPLEX refresh rejected'
            CALL assert_close('dry REAL density unchanged',rs%rho_m,before_rho,0.0_wp)
            CALL assert_close('dry density ambient',rs%rho_m,rho_a_amb,0.0_wp)
            CALL assert_close('dry COMPLEX density',REAL(cs%rho_m,wp),rho_a_amb,0.0_wp)
            CALL assert_close('dry REAL temperature ambient',rs%T,T_ambient,0.0_wp)
            CALL assert_close('dry COMPLEX temperature ambient',REAL(cs%T,wp),T_ambient,0.0_wp)
            CALL assert_close('dry REAL heat capacity',rs%cp_m,expected_cp,0.0_wp)
            CALL assert_close('dry COMPLEX heat capacity',REAL(cs%cp_m,wp),expected_cp,0.0_wp)
            CALL assert_close('dry REAL reduced gravity',rs%red_grav,0.0_wp,0.0_wp)
            CALL assert_close('dry COMPLEX reduced gravity',REAL(cs%red_grav,wp),0.0_wp,0.0_wp)
            CALL assert_close('dry REAL zero phase volume',SUM(rs%alphas),0.0_wp,0.0_wp)
            CALL assert_close('dry COMPLEX zero phase volume',REAL(SUM(cs%alphas),wp),0.0_wp,0.0_wp)
        END DO
        gas_flag=.TRUE.; liquid_flag=.FALSE.
    END SUBROUTINE check_dry_refresh_consistency

    SUBROUTINE check_canonical_dry_kinematics
        REAL(wp) :: qp_case(n_vars),q_case(n_vars),rho_mix,red_grav
        COMPLEX(wp) :: cqp_case(n_vars),cq_case(n_vars)
        TYPE(physical_state_real) :: rs
        TYPE(physical_state_complex) :: cs
        LOGICAL :: ok
        CALL initialize_physical_state_real(rs)
        CALL initialize_physical_state_complex(cs)
        qp_case=0.0_wp
        qp_case(idx_h)=0.5_wp*EPSILON(1.0_wp)
        qp_case(idx_T)=650.0_wp; qp_case(idx_u)=17.0_wp; qp_case(idx_v)=-9.0_wp
        qp_case(idx_solid_first)=0.2_wp
        qp_case(idx_stoch)=2.0_wp; qp_case(idx_pore)=5.0_wp
        cqp_case=CMPLX(qp_case,0.0_wp,wp)
        CALL build_physical_state_from_qp(qp_case,rs,accepted=ok)
        IF (.NOT.ok) ERROR STOP 'canonical REAL dry qp rejected'
        CALL build_physical_state_from_qp(cqp_case,cs,accepted=ok)
        IF (.NOT.ok) ERROR STOP 'canonical COMPLEX dry qp rejected'
        CALL assert_close('dry REAL velocity',rs%u,0.0_wp,0.0_wp)
        CALL assert_close('dry COMPLEX velocity',REAL(cs%u,wp),0.0_wp,0.0_wp)
        CALL assert_close('dry REAL thickness',rs%h,0.0_wp,0.0_wp)
        CALL assert_close('dry COMPLEX thickness',REAL(cs%h,wp),0.0_wp,0.0_wp)
        CALL assert_close('dry REAL species',SUM(rs%xs),0.0_wp,0.0_wp)
        CALL assert_close('dry COMPLEX species',REAL(SUM(cs%xs),wp),0.0_wp,0.0_wp)
        CALL assert_close('dry REAL scalar',rs%Zs,0.0_wp,0.0_wp)
        rs%u=11.0_wp; rs%h=0.5_wp*EPSILON(1.0_wp)
        cs%u=CMPLX(11.0_wp,1.0E-20_wp,wp)
        cs%h=CMPLX(0.5_wp*EPSILON(1.0_wp),0.0_wp,wp)
        CALL refresh_physical_state(rs,accepted=ok)
        IF (.NOT.ok) ERROR STOP 'canonical REAL dry refresh rejected'
        CALL refresh_physical_state(cs,accepted=ok)
        IF (.NOT.ok) ERROR STOP 'canonical COMPLEX dry refresh rejected'
        CALL assert_close('dry REAL refreshed velocity',rs%u,0.0_wp,0.0_wp)
        CALL assert_close('dry COMPLEX refreshed velocity',REAL(cs%u,wp),0.0_wp,0.0_wp)
        CALL assert_close('dry REAL refreshed depth',rs%h,0.0_wp,0.0_wp)
        qp_case(idx_h)=-1.0_wp
        CALL evaluate_density_from_qp(qp_case,rho_mix,red_grav,ok)
        IF (ok) ERROR STOP 'negative-depth density evaluation accepted'
        qp_case(idx_h)=1.0_wp; qp_case(idx_T)=-10.0_wp
        CALL primitive_to_conservative(qp_case,q_case,ok)
        IF (ok) ERROR STOP 'negative-temperature qp accepted'
        cq_case=CMPLX(qp_case,0.0_wp,wp)
        CALL primitive_to_conservative(cq_case,cqp_case,ok)
        IF (ok) ERROR STOP 'negative-temperature COMPLEX qp accepted'
        CALL build_physical_state_from_qp(qp_case,rs,accepted=ok)
        IF (ok) ERROR STOP 'negative-temperature physical state accepted'
        cqp_case=CMPLX(qp_case,0.0_wp,wp)
        CALL build_physical_state_from_qp(cqp_case,cs,accepted=ok)
        IF (ok) ERROR STOP 'negative-temperature COMPLEX state accepted'
    END SUBROUTINE check_canonical_dry_kinematics

    SUBROUTINE check_nonpositive_heat_capacity
        REAL(wp) :: qp_case(n_vars),q_case(n_vars),qp_out(n_vars),temperature,old_cp_s(n_solid)
        COMPLEX(wp) :: cq_case(n_vars),cqp_out(n_vars),ctemperature
        TYPE(physical_state_real) :: rs
        TYPE(physical_state_complex) :: cs
        LOGICAL :: ok
        CALL initialize_physical_state_real(rs)
        CALL initialize_physical_state_complex(cs)
        qp_case=0.0_wp; qp_case(idx_h)=1.0_wp; qp_case(idx_T)=350.0_wp
        qp_case(idx_solid_first:idx_solid_last)=[0.5_wp,0.4_wp]
        qp_case(idx_add_gas_first:idx_add_gas_last)=0.05_wp
        CALL primitive_to_conservative(qp_case,q_case,ok)
        IF (.NOT.ok) ERROR STOP 'cp test initial state creation failed'
        cq_case=CMPLX(q_case,0.0_wp,wp)
        old_cp_s=sp_heat_s
        sp_heat_s=-1.0E6_wp
        CALL conservative_to_primitive(q_case,qp_out,ok)
        IF (ok) ERROR STOP 'negative cp REAL conversion accepted'
        CALL conservative_to_primitive(cq_case,cqp_out,ok)
        IF (ok) ERROR STOP 'negative cp COMPLEX conversion accepted'
        CALL build_physical_state_from_q(q_case,rs,accepted=ok)
        IF (ok) ERROR STOP 'negative cp REAL state accepted'
        CALL build_physical_state_from_q(cq_case,cs,accepted=ok)
        IF (ok) ERROR STOP 'negative cp COMPLEX state accepted'
        CALL evaluate_temperature_from_q(q_case,temperature,ok)
        IF (ok) ERROR STOP 'negative cp REAL temperature accepted'
        CALL evaluate_temperature_from_q(cq_case,ctemperature,ok)
        IF (ok) ERROR STOP 'negative cp COMPLEX temperature accepted'
        CALL assert_close('negative cp REAL temperature',temperature,0.0_wp,0.0_wp)
        CALL assert_close('negative cp COMPLEX temperature',REAL(ctemperature,wp),0.0_wp,0.0_wp)
        sp_heat_s=old_cp_s
    END SUBROUTINE check_nonpositive_heat_capacity

    SUBROUTINE check_epsilon_dry_real_complex
        REAL(wp) :: q_case(n_vars),qp_case(n_vars),decoded(n_vars),temperature
        COMPLEX(wp) :: cq_case(n_vars),cqp_case(n_vars),cdecoded(n_vars),ctemperature
        TYPE(physical_state_real) :: real_state
        TYPE(physical_state_complex) :: complex_state
        LOGICAL :: ok
        INTEGER :: case_id
        CALL initialize_physical_state_real(real_state)
        CALL initialize_physical_state_complex(complex_state)
        DO case_id=1,3
            q_case=0.0_wp
            SELECT CASE(case_id)
            CASE(1)
                q_case(idx_totMassEqn)=0.5_wp*EPSILON(1.0_wp)
            CASE(2)
                q_case(idx_totMassEqn)=EPSILON(1.0_wp)
            CASE(3)
                q_case(idx_totMassEqn)=2.0_wp*EPSILON(1.0_wp)
            END SELECT
            q_case(idx_engyEqn)=q_case(idx_totMassEqn)*sp_heat_a*400.0_wp
            q_case(idx_stoch)=q_case(idx_totMassEqn)*0.2_wp
            q_case(idx_pore)=q_case(idx_totMassEqn)*1000.0_wp
            cq_case=CMPLX(q_case,0.0_wp,wp)
            CALL conservative_to_primitive(q_case,decoded,ok)
            IF (.NOT.ok) ERROR STOP 'epsilon REAL q conversion rejected'
            CALL conservative_to_primitive(cq_case,cdecoded,ok)
            IF (.NOT.ok) ERROR STOP 'epsilon COMPLEX q conversion rejected'
            CALL build_physical_state_from_q(q_case,real_state,accepted=ok)
            IF (.NOT.ok) ERROR STOP 'epsilon REAL state rejected'
            CALL build_physical_state_from_q(cq_case,complex_state,accepted=ok)
            IF (.NOT.ok) ERROR STOP 'epsilon COMPLEX state rejected'
            CALL evaluate_temperature_from_q(q_case,temperature)
            CALL evaluate_temperature_from_q(cq_case,ctemperature)
            CALL assert_close('epsilon REAL/COMPLEX q temperature',decoded(idx_T), &
                REAL(cdecoded(idx_T),wp),0.0_wp)
            CALL assert_close('epsilon REAL/COMPLEX state temperature',real_state%T, &
                REAL(complex_state%T,wp),0.0_wp)
            IF (case_id<=2) THEN
                CALL assert_close('epsilon dry q temperature',decoded(idx_T),T_ambient,0.0_wp)
                CALL assert_close('epsilon dry state h',real_state%h,0.0_wp,0.0_wp)
                CALL assert_close('epsilon dry scalar',REAL(complex_state%Zs,wp),0.0_wp,0.0_wp)
                CALL assert_close('epsilon dry temperature only',temperature,T_ambient,0.0_wp)
            ELSE
                CALL assert_close('epsilon wet q temperature',decoded(idx_T),400.0_wp,1.0E-12_wp)
                CALL assert_close('epsilon wet temperature only',temperature,400.0_wp,1.0E-12_wp)
            END IF
        END DO
        DO case_id=1,3
            qp_case=0.0_wp
            SELECT CASE(case_id)
            CASE(1)
                qp_case(idx_h)=0.5_wp*EPSILON(1.0_wp)
            CASE(2)
                qp_case(idx_h)=EPSILON(1.0_wp)
            CASE(3)
                qp_case(idx_h)=2.0_wp*EPSILON(1.0_wp)
            END SELECT
            qp_case(idx_T)=400.0_wp
            cqp_case=CMPLX(qp_case,0.0_wp,wp)
            CALL primitive_to_conservative(qp_case,q_case,ok)
            IF (.NOT.ok) ERROR STOP 'epsilon REAL qp conversion rejected'
            CALL primitive_to_conservative(cqp_case,cq_case,ok)
            IF (.NOT.ok) ERROR STOP 'epsilon COMPLEX qp conversion rejected'
            CALL build_physical_state_from_qp(qp_case,real_state,accepted=ok)
            IF (.NOT.ok) ERROR STOP 'epsilon REAL qp state rejected'
            CALL build_physical_state_from_qp(cqp_case,complex_state,accepted=ok)
            IF (.NOT.ok) ERROR STOP 'epsilon COMPLEX qp state rejected'
            IF (case_id<=2) THEN
                CALL assert_close('epsilon dry q mass',q_case(idx_totMassEqn),0.0_wp,0.0_wp)
                CALL assert_close('epsilon dry state temperature',real_state%T,T_ambient,0.0_wp)
                CALL assert_close('epsilon dry complex state temperature',REAL(complex_state%T,wp), &
                    T_ambient,0.0_wp)
            ELSE
                IF (q_case(idx_totMassEqn)<=0.0_wp) ERROR STOP 'epsilon wet mass vanished'
                CALL assert_close('epsilon wet state temperature',real_state%T,400.0_wp,0.0_wp)
            END IF
        END DO
    END SUBROUTINE check_epsilon_dry_real_complex

    !> Preserve the granular descriptors from the original mixture formulas.
    SUBROUTINE check_particle_descriptors
        REAL(wp) :: fractions(2),expected_d,expected_rho
        fractions=[0.2_wp,0.3_wp]
        expected_d=SUM(fractions)/SUM(fractions/diam_s)
        expected_rho=SUM(fractions*rho_s)/SUM(fractions)
        CALL assert_close('weighted Sauter diameter',sauter_diameter(fractions),expected_d,1.0E-14_wp)
        CALL assert_close('weighted solid density',average_density_solids(fractions),expected_rho,1.0E-14_wp)
        fractions=0.0_wp
        CALL assert_close('zero-solid diameter',sauter_diameter(fractions),SUM(diam_s)/n_solid,1.0E-14_wp)
        CALL assert_close('zero-solid density',average_density_solids(fractions),SUM(rho_s)/n_solid,1.0E-14_wp)
        fractions=[0.25_wp*EPSILON(1.0_wp),0.25_wp*EPSILON(1.0_wp)]
        CALL assert_close('near-zero solid diameter',sauter_diameter(fractions),SUM(diam_s)/n_solid,1.0E-14_wp)
        CALL assert_close('near-zero solid density',average_density_solids(fractions),SUM(rho_s)/n_solid,1.0E-14_wp)
    END SUBROUTINE check_particle_descriptors

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
        ALLOCATE(rho_s(n_solid),inv_rho_s(n_solid),sp_heat_s(n_solid),diam_s(n_solid))
        ALLOCATE(sp_heat_g(n_add_gas),sp_gas_const_g(n_add_gas))
        rho_s=[2500.0_wp,1800.0_wp]; inv_rho_s=1.0_wp/rho_s; sp_heat_s=[900.0_wp,1100.0_wp]
        diam_s=[0.0001_wp,0.001_wp]
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
END PROGRAM test_state
