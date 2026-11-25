!> @brief Coordinate pair-wide interlayer momentum exchange without mass transfer.
!> @details Owns the single drag configuration and dispatches physical drag laws.
!>          No cell indices, grid storage, or alternative time integrator are introduced.
MODULE interlayer_momentum_exchange_2d
    USE parameters_2d, ONLY : wp,n_vars,idx_h,idx_u,idx_v,idx_totMassEqn,idx_uEqn,idx_vEqn
    USE diagnostics_2d, ONLY : fatal_error
    USE interlayer_drag_chezy_2d, ONLY : evaluate_interlayer_drag_chezy_analytic
    USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: interlayer_momentum_exchange_config_type,interlayer_momentum_exchange
    PUBLIC :: apply_interlayer_momentum_exchange_analytic

    !> @brief Canonical configuration for the pair-wide interlayer drag source.
    TYPE :: interlayer_momentum_exchange_config_type
        CHARACTER(LEN=24) :: model='CHEZY'         !< Physical interfacial drag law.
        REAL(wp) :: chezy_drag_coeff=0.0_wp        !< Dimensionless Chezy coefficient.
        CHARACTER(LEN=16) :: method='ANALYTIC'     !< Requested coupled Lie method.
        INTEGER :: group=5                         !< Positive Lie group; group zero is not supported.
    END TYPE interlayer_momentum_exchange_config_type

    TYPE(interlayer_momentum_exchange_config_type) :: interlayer_momentum_exchange

    ! Four physical treatment entry points; unsupported methods fail explicitly.
    PUBLIC :: apply_interlayer_momentum_exchange_explicit
    PUBLIC :: apply_interlayer_momentum_exchange_semi_implicit
    PUBLIC :: apply_interlayer_momentum_exchange_implicit

CONTAINS

    !> @brief Apply an exact drag-only substep to the two current Lie candidates.
    !> @details Active only for h1,h2>1 mm, as in IMEX_v3 Eq. (63).
    !>          Both momenta are staged before committing; every other component,
    !>          including each layer's mass, composition and thermal energy,
    !>          remains bitwise unchanged. Use only after the spatial/Lie recipient
    !>          state is accepted; no group-zero snapshot approximation is made.
    !> @param[in] config Canonical pair-wide drag law and coefficient.
    !> @param[in] dt Positive duration of the drag-only Lie substep [s].
    !> @param[in,out] q Conservative candidates, component by layer.
    !> @param[in,out] qp Primitive candidates, component by layer.
    !> @param[out] accepted True if the entire paired update is admissible.
    SUBROUTINE apply_interlayer_momentum_exchange_analytic(config,dt,q,qp,accepted)
        TYPE(interlayer_momentum_exchange_config_type), INTENT(IN) :: config
        REAL(wp), INTENT(IN) :: dt
        REAL(wp), INTENT(INOUT) :: q(n_vars,2),qp(n_vars,2)
        LOGICAL, INTENT(OUT) :: accepted
        REAL(wp), PARAMETER :: min_contact_depth=1.0E-3_wp
        REAL(wp) :: mass1,mass2,rho2,velocity1(2),velocity2(2),impulse(2)
        REAL(wp) :: momenta1(2),momenta2(2),updated_velocity1(2),updated_velocity2(2)
        LOGICAL :: valid

        accepted=.FALSE.
        IF (config%model/='CHEZY' .OR. config%method/='ANALYTIC' .OR. config%group<=0) RETURN
        IF (.NOT.ieee_is_finite(config%chezy_drag_coeff) .OR. config%chezy_drag_coeff<0.0_wp) RETURN
        IF (.NOT.ieee_is_finite(dt) .OR. dt<=0.0_wp) RETURN
        IF (.NOT.ALL(ieee_is_finite(q)) .OR. .NOT.ALL(ieee_is_finite(qp))) RETURN
        IF (ANY(qp(idx_h,:)<=min_contact_depth) .OR. config%chezy_drag_coeff==0.0_wp) THEN
            accepted=.TRUE.
            RETURN
        END IF

        mass1=q(idx_totMassEqn,1)
        mass2=q(idx_totMassEqn,2)
        IF (mass1<=0.0_wp .OR. mass2<=0.0_wp) RETURN
        rho2=mass2/qp(idx_h,2)
        velocity1=[q(idx_uEqn,1),q(idx_vEqn,1)]/mass1
        velocity2=[q(idx_uEqn,2),q(idx_vEqn,2)]/mass2
        CALL evaluate_interlayer_drag_chezy_analytic(mass1,mass2,rho2,velocity1,velocity2, &
            config%chezy_drag_coeff,dt,impulse,valid)
        IF (.NOT.valid) RETURN

        momenta1=[q(idx_uEqn,1),q(idx_vEqn,1)]+impulse
        momenta2=[q(idx_uEqn,2),q(idx_vEqn,2)]-impulse
        updated_velocity1=momenta1/mass1
        updated_velocity2=momenta2/mass2
        IF (.NOT.ALL(ieee_is_finite(momenta1)) .OR. .NOT.ALL(ieee_is_finite(momenta2)) .OR. &
            .NOT.ALL(ieee_is_finite(updated_velocity1)) .OR. &
            .NOT.ALL(ieee_is_finite(updated_velocity2))) RETURN

        q(idx_uEqn,1)=momenta1(1)
        q(idx_vEqn,1)=momenta1(2)
        q(idx_uEqn,2)=momenta2(1)
        q(idx_vEqn,2)=momenta2(2)
        qp(idx_u,1)=updated_velocity1(1)
        qp(idx_v,1)=updated_velocity1(2)
        qp(idx_u,2)=updated_velocity2(1)
        qp(idx_v,2)=updated_velocity2(2)
        accepted=.TRUE.
    END SUBROUTINE apply_interlayer_momentum_exchange_analytic
    !> @brief Reject unsupported explicit treatment for interlayer momentum exchange.
    !> @details No candidate state is committed; report the unsupported method.
    !> @param[in] config Pair-wide drag model and coefficient.
    !> @param[in] dt Proposed drag-only Lie duration [s].
    !> @param[in,out] q Conservative L1/L2 candidates, unchanged on fatal failure.
    !> @param[in,out] qp Primitive L1/L2 candidates, unchanged on fatal failure.
    !> @param[out] accepted False because the requested source method is unavailable.
    SUBROUTINE apply_interlayer_momentum_exchange_explicit(config,dt,q,qp,accepted)
        TYPE(interlayer_momentum_exchange_config_type), INTENT(IN) :: config
        REAL(wp), INTENT(IN) :: dt
        REAL(wp), INTENT(INOUT) :: q(n_vars,2),qp(n_vars,2)
        LOGICAL, INTENT(OUT) :: accepted

        accepted=.FALSE.
        CALL fatal_error('interlayer_momentum_exchange: unsupported explicit physical treatment')
    END SUBROUTINE apply_interlayer_momentum_exchange_explicit

    !> @brief Reject unsupported semi-implicit treatment for interlayer momentum exchange.
    !> @details No candidate state is committed; report the unsupported method.
    !> @param[in] config Pair-wide drag model and coefficient.
    !> @param[in] dt Proposed drag-only Lie duration [s].
    !> @param[in,out] q Conservative L1/L2 candidates, unchanged on fatal failure.
    !> @param[in,out] qp Primitive L1/L2 candidates, unchanged on fatal failure.
    !> @param[out] accepted False because the requested source method is unavailable.
    SUBROUTINE apply_interlayer_momentum_exchange_semi_implicit(config,dt,q,qp,accepted)
        TYPE(interlayer_momentum_exchange_config_type), INTENT(IN) :: config
        REAL(wp), INTENT(IN) :: dt
        REAL(wp), INTENT(INOUT) :: q(n_vars,2),qp(n_vars,2)
        LOGICAL, INTENT(OUT) :: accepted

        accepted=.FALSE.
        CALL fatal_error('interlayer_momentum_exchange: unsupported semi-implicit physical treatment')
    END SUBROUTINE apply_interlayer_momentum_exchange_semi_implicit

    !> @brief Reject unsupported implicit treatment for interlayer momentum exchange.
    !> @details No candidate state is committed; report the unsupported method.
    !> @param[in] config Pair-wide drag model and coefficient.
    !> @param[in] dt Proposed drag-only Lie duration [s].
    !> @param[in,out] q Conservative L1/L2 candidates, unchanged on fatal failure.
    !> @param[in,out] qp Primitive L1/L2 candidates, unchanged on fatal failure.
    !> @param[out] accepted False because the requested source method is unavailable.
    SUBROUTINE apply_interlayer_momentum_exchange_implicit(config,dt,q,qp,accepted)
        TYPE(interlayer_momentum_exchange_config_type), INTENT(IN) :: config
        REAL(wp), INTENT(IN) :: dt
        COMPLEX(wp), INTENT(INOUT) :: q(n_vars,2),qp(n_vars,2)
        LOGICAL, INTENT(OUT) :: accepted

        accepted=.FALSE.
        CALL fatal_error('interlayer_momentum_exchange: unsupported implicit physical treatment')
    END SUBROUTINE apply_interlayer_momentum_exchange_implicit

END MODULE interlayer_momentum_exchange_2d
