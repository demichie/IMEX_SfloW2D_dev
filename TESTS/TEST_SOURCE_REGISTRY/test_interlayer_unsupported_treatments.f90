!> @brief Verify that all unsupported interlayer physical treatments fail explicitly.
PROGRAM test_interlayer_unsupported_treatments
    USE parameters_2d, ONLY : wp,n_solid,n_vars,n_eqns
    USE state_types, ONLY : physical_state_real,physical_state_complex
    USE rheology_parameters_2d, ONLY : rheology_config_type
    USE mass_exchange_parameters_2d, ONLY : mass_exchange_config_type
    USE interlayer_mass_exchange_parameters_2d, ONLY : interlayer_mass_exchange_config_type
    USE interlayer_momentum_exchange_2d, ONLY : interlayer_momentum_exchange_config_type
    USE interlayer_settling_2d, ONLY : &
        eval_interlayer_settling_semi_implicit,eval_interlayer_settling_implicit, &
        solve_interlayer_settling_analytic
    USE interlayer_elutriation_2d, ONLY : &
        eval_interlayer_elutriation_semi_implicit,eval_interlayer_elutriation_implicit, &
        solve_interlayer_elutriation_analytic
    USE interlayer_mass_exchange_2d, ONLY : &
        eval_interlayer_mass_exchange_semi_implicit,eval_interlayer_mass_exchange_implicit, &
        solve_interlayer_mass_exchange_analytic
    USE interlayer_drag_chezy_2d, ONLY : &
        evaluate_interlayer_drag_chezy_explicit,evaluate_interlayer_drag_chezy_semi_implicit, &
        evaluate_interlayer_drag_chezy_implicit
    USE interlayer_momentum_exchange_2d, ONLY : &
        apply_interlayer_momentum_exchange_explicit,apply_interlayer_momentum_exchange_semi_implicit, &
        apply_interlayer_momentum_exchange_implicit
    IMPLICIT NONE

    TYPE(interlayer_mass_exchange_config_type) :: mass_pair
    TYPE(interlayer_momentum_exchange_config_type) :: momentum_pair
    TYPE(rheology_config_type) :: rheology
    TYPE(mass_exchange_config_type) :: mass
    TYPE(physical_state_real) :: physical
    TYPE(physical_state_complex) :: physical_complex
    REAL(wp) :: settling(1),elutriation(1),source_l1(6),source_l2(6)
    REAL(wp) :: impulse(2),q(6,2),qp(6,2),u1(2),u2(2)
    COMPLEX(wp) :: impulse_complex(2),q_complex(6,2),qp_complex(6,2)
    COMPLEX(wp) :: u1_complex(2),u2_complex(2)
    INTEGER :: index,status
    LOGICAL :: accepted
    CHARACTER(LEN=16) :: argument

    CALL GET_COMMAND_ARGUMENT(1,argument,STATUS=status)
    IF (status/=0) ERROR STOP 'Missing unsupported interlayer treatment index'
    READ(argument,*,IOSTAT=status) index
    IF (status/=0) ERROR STOP 'Invalid unsupported interlayer treatment index'
    n_solid=1
    n_vars=6
    n_eqns=6
    q=0.0_wp
    qp=0.0_wp
    q_complex=CMPLX(0.0_wp,0.0_wp,KIND=wp)
    qp_complex=CMPLX(0.0_wp,0.0_wp,KIND=wp)
    u1=0.0_wp
    u2=0.0_wp
    u1_complex=CMPLX(0.0_wp,0.0_wp,KIND=wp)
    u2_complex=CMPLX(0.0_wp,0.0_wp,KIND=wp)

    SELECT CASE(index)
    CASE(1)
        CALL eval_interlayer_settling_semi_implicit(mass_pair,rheology,mass,physical,1.0_wp,settling)
    CASE(2)
        CALL eval_interlayer_settling_implicit(mass_pair,rheology,mass,physical_complex,1.0_wp,settling)
    CASE(3)
        CALL solve_interlayer_settling_analytic(mass_pair,rheology,mass,physical,1.0_wp,settling)
    CASE(4)
        CALL eval_interlayer_elutriation_semi_implicit(mass_pair,mass,physical,physical,1.0_wp,elutriation)
    CASE(5)
        CALL eval_interlayer_elutriation_implicit(mass_pair,mass,physical_complex,physical_complex,1.0_wp,elutriation)
    CASE(6)
        CALL solve_interlayer_elutriation_analytic(mass_pair,mass,physical,physical,1.0_wp,elutriation)
    CASE(7)
        CALL eval_interlayer_mass_exchange_semi_implicit(mass_pair,rheology,mass,physical,physical, &
            physical,physical,1.0_wp,source_l1,source_l2,accepted)
    CASE(8)
        CALL eval_interlayer_mass_exchange_implicit(mass_pair,rheology,mass,physical_complex,physical_complex, &
            physical_complex,physical_complex,1.0_wp,source_l1,source_l2,accepted)
    CASE(9)
        CALL solve_interlayer_mass_exchange_analytic(mass_pair,rheology,mass,physical,physical, &
            physical,physical,1.0_wp,source_l1,source_l2,accepted)
    CASE(10)
        CALL evaluate_interlayer_drag_chezy_explicit(1.0_wp,1.0_wp,1.0_wp,u1,u2,1.0_wp,1.0_wp,impulse,accepted)
    CASE(11)
        CALL evaluate_interlayer_drag_chezy_semi_implicit(1.0_wp,1.0_wp,1.0_wp,u1,u2,1.0_wp,1.0_wp,impulse,accepted)
    CASE(12)
        CALL evaluate_interlayer_drag_chezy_implicit(1.0_wp,1.0_wp,1.0_wp,u1_complex,u2_complex, &
            1.0_wp,1.0_wp,impulse_complex,accepted)
    CASE(13)
        CALL apply_interlayer_momentum_exchange_explicit(momentum_pair,1.0_wp,q,qp,accepted)
    CASE(14)
        CALL apply_interlayer_momentum_exchange_semi_implicit(momentum_pair,1.0_wp,q,qp,accepted)
    CASE(15)
        CALL apply_interlayer_momentum_exchange_implicit(momentum_pair,1.0_wp,q_complex,qp_complex,accepted)
    CASE DEFAULT
        ERROR STOP 'Invalid unsupported interlayer treatment index'
    END SELECT

    ERROR STOP 'Unsupported interlayer treatment returned without a fatal error'
END PROGRAM test_interlayer_unsupported_treatments
