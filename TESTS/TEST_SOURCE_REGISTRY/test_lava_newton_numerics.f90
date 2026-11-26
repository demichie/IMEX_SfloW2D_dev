!> @brief Check Newton, complex-step and exact linear-drag updates for the lava viscous law.
!> @details With tau0=0, fixed temperature and h=1, the existing lava viscous law
!> is linear in velocity. No separate test-only constitutive law is registered.
PROGRAM test_lava_newton_numerics
    USE parameters_2d, ONLY : wp
    USE rheology_parameters_2d, ONLY : rheology_config_type
    USE source_catalog_2d, ONLY : source_term_descriptor_by_name, canonical_rheology_source_name,                 &
        source_term_type, rheology_name_has_role, SOURCE_ROLE_IMPLICIT, SOURCE_ROLE_SEMI_IMPLICIT
    USE state_types, ONLY : physical_state_real, physical_state_complex
    USE rheology_lava_2d, ONLY : eval_rheology_lava_implicit, eval_rheology_lava_implicit_real,                      &
        eval_rheology_lava_semi_implicit
    IMPLICIT NONE

    TYPE(source_term_type) :: term
    TYPE(rheology_config_type) :: config
    TYPE(physical_state_real) :: rs
    TYPE(physical_state_complex) :: cs
    REAL(wp), PARAMETER :: k=0.1_wp, eps=1.0E-24_wp
    REAL(wp) :: rho, velocity(2), q_before(2), q_after(2), expected(2), actual_real(2)
    REAL(wp) :: dt, analytic_jac(2,2), cs_jac(2,2), residual(2)
    REAL(wp) :: si_x,si_y,det,inverse_jac(2,2)
    COMPLEX(wp) :: complex_velocity(2), rate(2)
    INTEGER :: i,j,case_number

    term=source_term_descriptor_by_name('rheology_lava')
    IF (TRIM(term%name)/='rheology_lava') ERROR STOP 'Lava rheology missing from source catalog'
    IF (TRIM(canonical_rheology_source_name('lava'))/=TRIM(term%name)) &
        ERROR STOP 'Noncanonical lava rheology lookup'
    IF (.NOT.term%has_newton .OR. .NOT.term%has_explicit .OR. .NOT.term%has_analytic) &
        ERROR STOP 'Incorrect lava rheology capabilities'
    IF (.NOT.rheology_name_has_role(term%name,SOURCE_ROLE_IMPLICIT)) &
        ERROR STOP 'Lava implicit contribution not selected'
    IF (.NOT.rheology_name_has_role(term%name,SOURCE_ROLE_SEMI_IMPLICIT)) &
        ERROR STOP 'Lava semi-implicit contribution not selected'

    ! With h=1, T=T_ref and tau0=0, the existing viscous rate is -k*rho*u.
    config%nu_ref=k/3.0_wp
    config%visc_par=0.0_wp
    config%T_ref=300.0_wp
    config%tau0=0.0_wp
    rs%h=1.0_wp
    rs%T=config%T_ref
    cs%h=CMPLX(rs%h,0.0_wp,wp)
    cs%T=CMPLX(rs%T,0.0_wp,wp)

    DO case_number=1,4
        SELECT CASE(case_number)
        CASE(1)
            rho=1000.0_wp
            velocity=[3.0_wp,-4.0_wp]
            dt=0.75_wp
        CASE(2)
            rho=1800.0_wp
            velocity=[0.0_wp,0.0_wp]
            dt=0.25_wp
        CASE(3)
            rho=850.0_wp
            velocity=[-1.0E-5_wp,2.0E-5_wp]
            dt=10.0_wp
        CASE(4)
            rho=2750.0_wp
            velocity=[-8.0_wp,-6.0_wp]
            dt=0.001_wp
        END SELECT

        rs%rho_m=rho
        rs%u=velocity(1)
        rs%v=velocity(2)
        rs%mod_vel=SQRT(SUM(velocity**2))
        cs%rho_m=CMPLX(rho,0.0_wp,wp)
        complex_velocity=CMPLX(velocity,0.0_wp,wp)
        CALL set_complex_velocity(cs,complex_velocity)

        CALL eval_rheology_lava_semi_implicit(config,rs,si_x,si_y)
        IF (si_x/=0.0_wp .OR. si_y/=0.0_wp) ERROR STOP 'Zero-yield lava semi-implicit rate is nonzero'

        CALL eval_rheology_lava_implicit(config,cs,rate(1),rate(2))
        CALL eval_rheology_lava_implicit_real(config,rs,actual_real(1),actual_real(2))
        expected=-k*rho*velocity
        CALL close_vector(REAL(rate,wp),expected,'implicit rate')
        CALL close_vector(actual_real,expected,'real implicit rate')

        analytic_jac=0.0_wp
        analytic_jac(1,1)=-k*rho
        analytic_jac(2,2)=-k*rho

        DO j=1,2
            complex_velocity=CMPLX(velocity,0.0_wp,wp)
            complex_velocity(j)=complex_velocity(j)+CMPLX(0.0_wp,eps,wp)
            CALL set_complex_velocity(cs,complex_velocity)
            CALL eval_rheology_lava_implicit(config,cs,rate(1),rate(2))
            cs_jac(:,j)=AIMAG(rate)/eps
            CALL close_vector(cs_jac(:,j),analytic_jac(:,j),'complex-step column')
        END DO

        ! For fixed density, (I-dt*dS/dq)*q_after=q_before with dS/dq=-k*I.
        q_before=rho*velocity
        q_after=q_before
        DO i=1,5
            complex_velocity=CMPLX(q_after/rho,0.0_wp,wp)
            CALL set_complex_velocity(cs,complex_velocity)
            CALL eval_rheology_lava_implicit(config,cs,rate(1),rate(2))
            residual=q_after-q_before-dt*REAL(rate,wp)
            inverse_jac=0.0_wp
            det=(1.0_wp+k*dt)**2
            inverse_jac(1,1)=(1.0_wp+k*dt)/det
            inverse_jac(2,2)=(1.0_wp+k*dt)/det
            q_after=q_after-MATMUL(inverse_jac,residual)
        END DO

        expected=q_before/(1.0_wp+k*dt)
        CALL close_vector(q_after,expected,'Newton 2x2 versus exact step')
        complex_velocity=CMPLX(q_after/rho,0.0_wp,wp)
        CALL set_complex_velocity(cs,complex_velocity)
        CALL eval_rheology_lava_implicit(config,cs,rate(1),rate(2))
        residual=q_after-q_before-dt*REAL(rate,wp)
        CALL close_vector(residual,[0.0_wp,0.0_wp],'Newton final residual')
    END DO

    WRITE(*,'(A)') 'PASS: lava viscous Newton, REAL/COMPLEX rate, complex-step Jacobian, exact map (4 cases)'

CONTAINS

    !> @brief Update complex velocity and its derived norms for a local physical state.
    !> @param[in,out] state Local complex physical state to update.
    !> @param[in] velocity Two horizontal velocity components.
    SUBROUTINE set_complex_velocity(state,velocity)
        TYPE(physical_state_complex), INTENT(INOUT) :: state
        COMPLEX(wp), INTENT(IN) :: velocity(2)

        state%u=velocity(1)
        state%v=velocity(2)
        state%mod_vel2=SUM(velocity**2)
        state%mod_vel_hor=SQRT(state%mod_vel2)
        state%mod_vel=state%mod_vel_hor
    END SUBROUTINE set_complex_velocity

    !> @brief Assert that a two-component rate matches the independent expected value.
    !> @param[in] actual Computed rate or Newton residual.
    !> @param[in] reference Independently derived rate or residual.
    !> @param[in] label Name of the verified numerical invariant.
    SUBROUTINE close_vector(actual,reference,label)
        REAL(wp), INTENT(IN) :: actual(2),reference(2)
        CHARACTER(*), INTENT(IN) :: label
        REAL(wp) :: tolerance

        tolerance=3.0E-11_wp*MAX(1.0_wp,MAXVAL(ABS(reference)))
        IF (ANY(ABS(actual-reference)>tolerance)) THEN
            WRITE(*,*) 'FAIL ',label,actual,reference
            ERROR STOP 'Lava Newton regression failed'
        END IF
    END SUBROUTINE close_vector
END PROGRAM test_lava_newton_numerics
