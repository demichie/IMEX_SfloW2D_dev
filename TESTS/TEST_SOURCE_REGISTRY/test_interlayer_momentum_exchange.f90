!> @brief Conservative analytic Chezy drag regressions for the two-layer momentum pair.
PROGRAM test_interlayer_momentum_exchange
    USE parameters_2d, ONLY : wp,n_layers,n_vars,n_eqns,n_solid,n_add_gas,n_stoch_vars,n_pore_vars
    USE parameters_2d, ONLY : configure_canonical_state_indices,idx_h,idx_T,idx_u,idx_v
    USE parameters_2d, ONLY : idx_totMassEqn,idx_uEqn,idx_vEqn,idx_engyEqn
    USE interlayer_momentum_exchange_2d, ONLY : interlayer_momentum_exchange_config_type, &
        apply_interlayer_momentum_exchange_analytic
    IMPLICIT NONE
    TYPE(interlayer_momentum_exchange_config_type) :: config
    REAL(wp) :: q(4,2),qp(4,2),old_q(4,2),old_qp(4,2)
    REAL(wp) :: momentum_before(2),momentum_after(2),kinetic_before,kinetic_after
    REAL(wp) :: old_slip,expected_slip,mu,lambda,dt
    LOGICAL :: accepted
    INTEGER :: axis

    n_layers=2; n_solid=0; n_add_gas=0; n_stoch_vars=0; n_pore_vars=0
    n_vars=4; n_eqns=4
    CALL configure_canonical_state_indices()
    config%chezy_drag_coeff=0.03_wp
    dt=2.0_wp

    ! Unequal masses, oblique velocities, and exact two-dimensional relaxation.
    CALL init_pair(2.0_wp,3.0_wp,[5.0_wp,-2.0_wp],[-1.0_wp,4.0_wp])
    old_q=q; old_qp=qp
    momentum_before=total_momentum(q)
    kinetic_before=kinetic_energy(q)
    old_slip=SQRT(SUM((old_qp(idx_u:idx_v,2)-old_qp(idx_u:idx_v,1))**2))
    lambda=config%chezy_drag_coeff*(q(1,2)/qp(idx_h,2))*(1.0_wp/q(1,1)+1.0_wp/q(1,2))
    expected_slip=old_slip/(1.0_wp+lambda*old_slip*dt)
    CALL apply_interlayer_momentum_exchange_analytic(config,dt,q,qp,accepted)
    CALL require(accepted,'oblique Chezy update accepted')
    momentum_after=total_momentum(q)
    DO axis=1,2
        CALL close(momentum_after(axis),momentum_before(axis),'total momentum')
    END DO
    CALL close(SQRT(SUM((qp(idx_u:idx_v,2)-qp(idx_u:idx_v,1))**2)), &
        expected_slip,'exact quadratic relative-velocity decay')
    kinetic_after=kinetic_energy(q)
    CALL require(kinetic_after<=kinetic_before,'Chezy drag dissipates mechanical kinetic energy')
    mu=old_q(1,1)*old_q(1,2)/(old_q(1,1)+old_q(1,2))
    CALL close(kinetic_before-kinetic_after,0.5_wp*mu*(old_slip**2-expected_slip**2), &
        'exact loss of relative kinetic energy')
    CALL require(ALL(q(idx_totMassEqn,:)==old_q(idx_totMassEqn,:)), &
        'each layer preserves its own mass')
    CALL require(ALL(q(idx_engyEqn,:)==old_q(idx_engyEqn,:)), &
        'each layer preserves its own thermal energy')
    CALL require(ALL(qp(idx_h,:)==old_qp(idx_h,:)) .AND. &
        ALL(qp(idx_T,:)==old_qp(idx_T,:)), &
        'thickness and temperature remain unchanged')

    ! A very stiff substep must approach barycentric velocity without overshoot.
    CALL init_pair(0.5_wp,200.0_wp,[100.0_wp,30.0_wp],[-50.0_wp,-10.0_wp])
    old_q=q
    CALL apply_interlayer_momentum_exchange_analytic(config,1.0E6_wp,q,qp,accepted)
    CALL require(accepted,'stiff analytic drag is stable')
    CALL require(kinetic_energy(q)<=kinetic_energy(old_q),'stiff drag does not add kinetic energy')
    CALL close(SUM(q(idx_uEqn,:)),SUM(old_q(idx_uEqn,:)),'stiff drag x momentum')
    CALL close(SUM(q(idx_vEqn,:)),SUM(old_q(idx_vEqn,:)),'stiff drag y momentum')

    ! All four wet/dry configurations and the 1 mm physical interface gate.
    CALL init_pair(0.0_wp,0.0_wp,[0.0_wp,0.0_wp],[0.0_wp,0.0_wp])
    CALL assert_noop('both layers dry')
    CALL init_pair(1.0_wp,0.0_wp,[5.0_wp,0.0_wp],[0.0_wp,0.0_wp])
    CALL assert_noop('L1 alone')
    CALL init_pair(0.0_wp,1.0_wp,[0.0_wp,0.0_wp],[5.0_wp,0.0_wp])
    CALL assert_noop('L2 alone')
    CALL init_pair(1.0E-3_wp,1.0_wp,[5.0_wp,0.0_wp],[0.0_wp,0.0_wp])
    CALL assert_noop('L1 at the interface threshold')
    CALL init_pair(1.0_wp,1.0E-3_wp,[5.0_wp,0.0_wp],[0.0_wp,0.0_wp])
    CALL assert_noop('L2 at the interface threshold')
    CALL init_pair(1.0_wp,1.0_wp,[5.0_wp,0.0_wp],[5.0_wp,0.0_wp])
    CALL assert_noop('both wet with identical velocities')
    CALL init_pair(1.0_wp,1.0_wp,[5.0_wp,0.0_wp],[-5.0_wp,0.0_wp])
    config%chezy_drag_coeff=0.0_wp
    CALL assert_noop('zero drag coefficient')

    ! Invalid inputs may not commit one recipient before the other.
    config%chezy_drag_coeff=-1.0_wp
    old_q=q; old_qp=qp
    CALL apply_interlayer_momentum_exchange_analytic(config,dt,q,qp,accepted)
    CALL require(.NOT.accepted,'negative Chezy coefficient rejected')
    CALL require(ALL(q==old_q) .AND. ALL(qp==old_qp),'failed pair leaves both candidates untouched')
    config%chezy_drag_coeff=0.03_wp
    CALL apply_interlayer_momentum_exchange_analytic(config,-dt,q,qp,accepted)
    CALL require(.NOT.accepted,'negative timestep rejected')
    CALL require(ALL(q==old_q) .AND. ALL(qp==old_qp),'failed step leaves both candidates untouched')

    PRINT *, 'PASS interlayer Chezy analytic momentum conservation, exact decay, energy and dry gates'
CONTAINS

    !> @brief Initialize a consistent test pair with unit mixture density.
    !> @param[in] h1 L1 thickness [m].
    !> @param[in] h2 L2 thickness [m].
    !> @param[in] v1 L1 velocity components [m/s].
    !> @param[in] v2 L2 velocity components [m/s].
    SUBROUTINE init_pair(h1,h2,v1,v2)
        REAL(wp), INTENT(IN) :: h1,h2,v1(2),v2(2)
        REAL(wp) :: h(2),vel(2,2)
        INTEGER :: layer
        h=[h1,h2]
        vel(:,1)=v1
        vel(:,2)=v2
        q=0.0_wp
        qp=0.0_wp
        DO layer=1,2
            q(idx_totMassEqn,layer)=h(layer)
            q(idx_uEqn,layer)=h(layer)*vel(1,layer)
            q(idx_vEqn,layer)=h(layer)*vel(2,layer)
            q(idx_engyEqn,layer)=h(layer)*300.0_wp
            qp(idx_h,layer)=h(layer)
            qp(idx_u,layer)=vel(1,layer)
            qp(idx_v,layer)=vel(2,layer)
            qp(idx_T,layer)=300.0_wp
        END DO
    END SUBROUTINE init_pair

    !> @brief Assert unchanged conservative and primitive states.
    !> @param[in] message Name of this zero-drag case.
    SUBROUTINE assert_noop(message)
        CHARACTER(*), INTENT(IN) :: message
        old_q=q; old_qp=qp
        CALL apply_interlayer_momentum_exchange_analytic(config,dt,q,qp,accepted)
        CALL require(accepted,'zero-drag accepted: '//message)
        CALL require(ALL(q==old_q) .AND. ALL(qp==old_qp),'zero-drag no-op: '//message)
    END SUBROUTINE assert_noop

    !> @brief Return total horizontal momentum of the two layers.
    !> @param[in] state Conservative state pair.
    !> @return Sum of x- and y-momenta.
    PURE FUNCTION total_momentum(state) RESULT(total)
        REAL(wp), INTENT(IN) :: state(4,2)
        REAL(wp) :: total(2)
        total=[SUM(state(idx_uEqn,:)),SUM(state(idx_vEqn,:))]
    END FUNCTION total_momentum

    !> @brief Calculate total kinetic energy per horizontal area.
    !> @param[in] state Conservative state pair.
    !> @return Kinetic energy [J/m2].
    PURE REAL(wp) FUNCTION kinetic_energy(state) RESULT(energy)
        REAL(wp), INTENT(IN) :: state(4,2)
        INTEGER :: layer
        energy=0.0_wp
        DO layer=1,2
            IF (state(idx_totMassEqn,layer)<=0.0_wp) CYCLE
            energy=energy+0.5_wp* &
                (state(idx_uEqn,layer)**2+state(idx_vEqn,layer)**2)/state(idx_totMassEqn,layer)
        END DO
    END FUNCTION kinetic_energy

    !> @brief Require a floating-point invariant with scaled tolerance.
    !> @param[in] actual Observed real quantity.
    !> @param[in] expected Analytic expected value.
    !> @param[in] message Name of the invariant.
    SUBROUTINE close(actual,expected,message)
        REAL(wp), INTENT(IN) :: actual,expected
        CHARACTER(*), INTENT(IN) :: message
        CALL require(ABS(actual-expected)<=1.0E-10_wp*MAX(1.0_wp,ABS(expected)),message)
    END SUBROUTINE close

    !> @brief Fail the regression when a physical invariant is violated.
    !> @param[in] condition Required condition.
    !> @param[in] message Invariant description.
    SUBROUTINE require(condition,message)
        LOGICAL, INTENT(IN) :: condition
        CHARACTER(*), INTENT(IN) :: message
        IF (.NOT.condition) THEN
            PRINT *, 'FAIL interlayer Chezy drag: ',TRIM(message)
            ERROR STOP 1
        END IF
    END SUBROUTINE require
END PROGRAM test_interlayer_momentum_exchange
