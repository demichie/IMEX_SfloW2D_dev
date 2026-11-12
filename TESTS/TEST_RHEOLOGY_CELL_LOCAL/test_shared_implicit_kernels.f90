!> Verify the single REAL/COMPLEX implicit formula used by seven rheologies.
PROGRAM test_shared_implicit_kernels
    USE parameters_2d, ONLY : wp, n_solid, n_add_gas
    USE constitutive_parameters_2d, ONLY : diam_s, rho_s, grav, inv_grav
    USE rheology_parameters_2d, ONLY : rheology_config_type, xi, I_0, mu_2, mu_s, muI_inf, T_ref
    USE rheology_2d, ONLY : eval_rheology_implicit_real_id
    USE rheology_newmodel_2d, ONLY : eval_rheology_newmodel_implicit, eval_rheology_newmodel_implicit_real
    USE state_types, ONLY : physical_state_real, physical_state_complex, &
        initialize_physical_state_real, promote_physical_state
    USE rheology_voellmy_2d, ONLY : voellmy_mu=>mu, &
        eval_rheology_voellmy_explicit, eval_rheology_voellmy_semi_implicit, eval_rheology_voellmy_implicit
    USE rheology_mui_voellmy_2d, ONLY : eval_rheology_mui_voellmy_explicit, &
        eval_rheology_mui_voellmy_semi_implicit, eval_rheology_mui_voellmy_implicit
    USE rheology_lava_2d, ONLY : tau0, nu_ref, visc_par, &
        eval_rheology_lava_explicit, eval_rheology_lava_semi_implicit, eval_rheology_lava_implicit
    USE rheology_plastic_2d, ONLY : &
        eval_rheology_plastic_explicit, eval_rheology_plastic_semi_implicit, eval_rheology_plastic_implicit
    USE rheology_quadratic_2d, ONLY : &
        eval_rheology_quadratic_explicit, eval_rheology_quadratic_semi_implicit, eval_rheology_quadratic_implicit
    USE rheology_empirical_2d, ONLY : eval_rheology_empirical_explicit, &
        eval_rheology_empirical_semi_implicit, eval_rheology_empirical_implicit
    USE rheology_lahar_2d, ONLY : alpha2, beta2, alpha1_coeff, beta1, Kappa, n_td2, &
        eval_rheology_lahar_explicit, eval_rheology_lahar_semi_implicit, eval_rheology_lahar_implicit
    IMPLICIT NONE

    TYPE(rheology_config_type) :: quadratic_config
    REAL(wp) :: friction_factor

    REAL(wp), PARAMETER :: eps_cs=1.0E-20_wp, fd_step=1.0E-6_wp
    TYPE(physical_state_real) :: state,shifted
    TYPE(physical_state_complex) :: complex_state
    REAL(wp) :: explicit_rate(2),semi_rate(2),plus_exp(2),plus_semi(2)
    REAL(wp) :: minus_exp(2),minus_semi(2),jacobian_fd(2),jacobian_cs(2),implicit_real(2)
    COMPLEX(wp) :: implicit_rate(2)
    INTEGER, PARAMETER :: model_ids(7)=[1,10,3,2,6,5,4]
    INTEGER :: model,component

    n_solid=1
    n_add_gas=0
    ALLOCATE(diam_s(1),rho_s(1))
    diam_s=0.001_wp
    rho_s=2500.0_wp
    grav=9.81_wp
    inv_grav=1.0_wp/grav
    xi=100.0_wp
    I_0=0.3_wp
    mu_s=0.35_wp
    mu_2=0.6_wp
    muI_inf=0.0_wp
    T_ref=300.0_wp
    voellmy_mu=0.2_wp
    tau0=150.0_wp
    nu_ref=1.0E-4_wp
    visc_par=0.002_wp
    quadratic_config%tau=0.15_wp
    friction_factor=0.01_wp
    quadratic_config%friction_factor=friction_factor
    alpha2=150.0_wp
    beta2=2.0_wp
    alpha1_coeff=1.0_wp
    beta1=0.2_wp
    Kappa=2.0_wp
    n_td2=0.005_wp

    CALL initialize_physical_state_real(state)
    state%h=2.0_wp
    state%T=320.0_wp
    state%rho_m=1800.0_wp
    state%inv_rhom=1.0_wp/state%rho_m
    state%red_grav=5.0_wp
    state%alphas=0.15_wp
    state%curvature_acceleration=0.0_wp
    state%exc_pore_pres=0.0_wp
    CALL set_real_velocity(state,[3.0_wp,4.0_wp])

    DO model=1,7
        CALL evaluate_real(model,state,explicit_rate,semi_rate)
        CALL promote_physical_state(state,complex_state)
        CALL evaluate_complex(model,complex_state,implicit_rate)
        CALL eval_rheology_implicit_real_id(quadratic_config,model_ids(model),state,.FALSE.,0.0_wp,implicit_real)
        CALL assert_rate(implicit_real,REAL(implicit_rate,wp),1.0E-12_wp,model,-3)
        CALL assert_rate(explicit_rate,semi_rate+REAL(implicit_rate,wp),1.0E-12_wp,model,0)

        DO component=1,2
            shifted=state
            IF (component==1) THEN
                CALL set_real_velocity(shifted,[state%u+fd_step,state%v])
            ELSE
                CALL set_real_velocity(shifted,[state%u,state%v+fd_step])
            END IF
            CALL evaluate_real(model,shifted,plus_exp,plus_semi)

            IF (component==1) THEN
                CALL set_real_velocity(shifted,[state%u-fd_step,state%v])
            ELSE
                CALL set_real_velocity(shifted,[state%u,state%v-fd_step])
            END IF
            CALL evaluate_real(model,shifted,minus_exp,minus_semi)
            jacobian_fd=((plus_exp-plus_semi)-(minus_exp-minus_semi))/(2.0_wp*fd_step)

            CALL promote_physical_state(state,complex_state)
            IF (component==1) THEN
                complex_state%u=complex_state%u+CMPLX(0.0_wp,eps_cs,wp)
            ELSE
                complex_state%v=complex_state%v+CMPLX(0.0_wp,eps_cs,wp)
            END IF
            complex_state%mod_vel2=complex_state%u**2+complex_state%v**2
            complex_state%mod_vel_hor=SQRT(complex_state%mod_vel2)
            complex_state%mod_vel=complex_state%mod_vel_hor
            CALL evaluate_complex(model,complex_state,implicit_rate)
            jacobian_cs=AIMAG(implicit_rate)/eps_cs
            CALL assert_rate(jacobian_cs,jacobian_fd,2.0E-5_wp,model,component)
        END DO

        CALL set_real_velocity(state,[0.0_wp,0.0_wp])
        CALL evaluate_real(model,state,explicit_rate,semi_rate)
        CALL promote_physical_state(state,complex_state)
        CALL evaluate_complex(model,complex_state,implicit_rate)
        CALL assert_rate(explicit_rate,[0.0_wp,0.0_wp],0.0_wp,model,-1)
        CALL assert_rate(REAL(implicit_rate,wp),[0.0_wp,0.0_wp],0.0_wp,model,-2)
        CALL set_real_velocity(state,[3.0_wp,4.0_wp])
    END DO

    ! The Newton-only linear fixture also has identical REAL/COMPLEX kernels.
    CALL promote_physical_state(state,complex_state)
    CALL eval_rheology_newmodel_implicit_real(state,implicit_real(1),implicit_real(2))
    CALL eval_rheology_newmodel_implicit(complex_state,implicit_rate(1),implicit_rate(2))
    CALL assert_rate(implicit_real,REAL(implicit_rate,wp),1.0E-12_wp,8,0)
    PRINT *, 'PASS: shared implicit rheology REAL/COMPLEX, complex-step and zero speed'
CONTAINS
    SUBROUTINE set_real_velocity(s,v)
        TYPE(physical_state_real), INTENT(INOUT) :: s
        REAL(wp), INTENT(IN) :: v(2)
        s%u=v(1)
        s%v=v(2)
        s%mod_vel2=SUM(v**2)
        s%mod_vel_hor=SQRT(s%mod_vel2)
        s%mod_vel=s%mod_vel_hor
    END SUBROUTINE set_real_velocity

    SUBROUTINE evaluate_real(model,s,explicit_part,semi_part)
        INTEGER, INTENT(IN) :: model
        TYPE(physical_state_real), INTENT(IN) :: s
        REAL(wp), INTENT(OUT) :: explicit_part(2),semi_part(2)
        REAL(wp) :: sx,sy
        SELECT CASE(model)
        CASE(1)
            CALL eval_rheology_voellmy_explicit(s,1.0_wp,.FALSE.,.FALSE.,0.0_wp,explicit_part)
            CALL eval_rheology_voellmy_semi_implicit(s,1.0_wp,.FALSE.,sx,sy)
        CASE(2)
            CALL eval_rheology_mui_voellmy_explicit(s,1.0_wp,.FALSE.,.FALSE.,0.0_wp,explicit_part)
            CALL eval_rheology_mui_voellmy_semi_implicit(s,1.0_wp,.FALSE.,sx,sy)
        CASE(3)
            CALL eval_rheology_lava_explicit(s,explicit_part)
            CALL eval_rheology_lava_semi_implicit(s,sx,sy)
        CASE(4)
            CALL eval_rheology_plastic_explicit(quadratic_config,s,explicit_part)
            CALL eval_rheology_plastic_semi_implicit(quadratic_config,s,sx,sy)
        CASE(5)
            CALL eval_rheology_quadratic_explicit(quadratic_config,s,explicit_part)
            CALL eval_rheology_quadratic_semi_implicit(quadratic_config,s,sx,sy)
        CASE(6)
            CALL eval_rheology_empirical_explicit(s,explicit_part)
            CALL eval_rheology_empirical_semi_implicit(s,sx,sy)
        CASE(7)
            CALL eval_rheology_lahar_explicit(s,explicit_part)
            CALL eval_rheology_lahar_semi_implicit(s,sx,sy)
        CASE DEFAULT
            ERROR STOP 'Unknown test rheology'
        END SELECT
        semi_part=[sx,sy]
    END SUBROUTINE evaluate_real

    SUBROUTINE evaluate_complex(model,s,implicit_part)
        INTEGER, INTENT(IN) :: model
        TYPE(physical_state_complex), INTENT(IN) :: s
        COMPLEX(wp), INTENT(OUT) :: implicit_part(2)
        SELECT CASE(model)
        CASE(1)
            CALL eval_rheology_voellmy_implicit(s,.FALSE.,0.0_wp,implicit_part(1),implicit_part(2))
        CASE(2)
            CALL eval_rheology_mui_voellmy_implicit(s,.FALSE.,0.0_wp,implicit_part(1),implicit_part(2))
        CASE(3)
            CALL eval_rheology_lava_implicit(s,implicit_part(1),implicit_part(2))
        CASE(4)
            CALL eval_rheology_plastic_implicit(quadratic_config,s,implicit_part(1),implicit_part(2))
        CASE(5)
            CALL eval_rheology_quadratic_implicit(quadratic_config,s,implicit_part(1),implicit_part(2))
        CASE(6)
            CALL eval_rheology_empirical_implicit(s,implicit_part(1),implicit_part(2))
        CASE(7)
            CALL eval_rheology_lahar_implicit(s,implicit_part(1),implicit_part(2))
        CASE DEFAULT
            ERROR STOP 'Unknown test rheology'
        END SELECT
    END SUBROUTINE evaluate_complex

    SUBROUTINE assert_rate(actual,expected,tolerance,model,component)
        REAL(wp), INTENT(IN) :: actual(2),expected(2),tolerance
        INTEGER, INTENT(IN) :: model,component
        REAL(wp) :: limit
        limit=tolerance*MAX(1.0_wp,MAXVAL(ABS(expected)))
        IF (ANY(ABS(actual-expected)>limit)) THEN
            PRINT *, 'FAIL shared implicit rheology model/component:',model,component
            PRINT *, 'Actual and expected:',actual,expected
            ERROR STOP 'Shared implicit rheology regression failed'
        END IF
    END SUBROUTINE assert_rate
END PROGRAM test_shared_implicit_kernels
