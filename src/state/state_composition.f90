MODULE state_composition
    USE parameters_2d, ONLY : wp,n_solid,n_add_gas,n_eqns
    USE parameters_2d, ONLY : idx_totMassEqn,idx_solidEqn_first,idx_solidEqn_last, &
        idx_add_gas_first,idx_add_gas_last
    IMPLICIT NONE
    PRIVATE
    REAL(wp), PARAMETER :: tol=4096.0_wp*EPSILON(1.0_wp)
    PUBLIC :: close_composition, normalize_face_composition, conservative_mass_admissible
    INTERFACE close_composition
        MODULE PROCEDURE close_composition_real, close_composition_complex
    END INTERFACE
CONTAINS
    !> Strict conservative mass budget shared by the solver and source limiter.
    !> Never silently clip a negative component or accept a tolerance violation.
    PURE LOGICAL FUNCTION conservative_mass_admissible(q_cell) RESULT(valid)
        REAL(wp), INTENT(IN) :: q_cell(n_eqns)
        valid=.FALSE.
        IF (q_cell(idx_totMassEqn)<0.0_wp) RETURN
        IF (ANY(q_cell(idx_solidEqn_first:idx_solidEqn_last)<0.0_wp)) RETURN
        IF (ANY(q_cell(idx_add_gas_first:idx_add_gas_last)<0.0_wp)) RETURN
        IF (SUM(q_cell(idx_solidEqn_first:idx_solidEqn_last))+ &
            SUM(q_cell(idx_add_gas_first:idx_add_gas_last))>q_cell(idx_totMassEqn)) RETURN
        valid=.TRUE.
    END FUNCTION conservative_mass_admissible

    !> Closure of reconstructed face traces (distinct from strict state conversion).
    !> Retains legacy HP-PCCU semantics: clip roundoff-negative species and
    !> normalize any positive overfill, while rejecting materially negative input.
    LOGICAL FUNCTION normalize_face_composition(xs,xg)
        REAL(wp), INTENT(INOUT) :: xs(n_solid), xg(n_add_gas)
        REAL(wp) :: total
        normalize_face_composition=.FALSE.
        IF (ANY(xs<-tol) .OR. ANY(xg<-tol)) RETURN
        xs=MAX(xs,0.0_wp)
        xg=MAX(xg,0.0_wp)
        total=SUM(xs)+SUM(xg)
        IF (total>1.0_wp) THEN
            xs=xs/total
            xg=xg/total
        END IF
        normalize_face_composition=.TRUE.
    END FUNCTION normalize_face_composition

    LOGICAL FUNCTION close_composition_real(xs,xg)
        REAL(wp), INTENT(INOUT) :: xs(n_solid), xg(n_add_gas)
        REAL(wp) :: total
        close_composition_real=.FALSE.
        IF (ANY(xs < -tol) .OR. ANY(xg < -tol)) RETURN
        ! Retain every positive species fraction, including sub-roundoff values.
        ! The historical source-state closure only clipped negative roundoff;
        ! zeroing positive solids would suppress near-depleted settling rates.
        xs=MAX(xs,0.0_wp)
        xg=MAX(xg,0.0_wp)
        total=SUM(xs)+SUM(xg)
        IF (total>1.0_wp+tol) RETURN
        IF (total>1.0_wp) THEN
            xs=xs/total; xg=xg/total
        END IF
        close_composition_real=.TRUE.
    END FUNCTION close_composition_real
    LOGICAL FUNCTION close_composition_complex(xs,xg)
        COMPLEX(wp), INTENT(INOUT) :: xs(n_solid), xg(n_add_gas)
        REAL(wp) :: total
        COMPLEX(wp) :: total_complex
        close_composition_complex=.FALSE.
        IF (ANY(REAL(xs,wp)<-tol) .OR. ANY(REAL(xg,wp)<-tol)) RETURN
        ! As in REAL, clip only roundoff-negative species. A clipped
        ! fraction has a zero local derivative; others retain complex steps.
        WHERE (REAL(xs,wp)<0.0_wp) xs=(0.0_wp,0.0_wp)
        WHERE (REAL(xg,wp)<0.0_wp) xg=(0.0_wp,0.0_wp)
        total=SUM(REAL(xs,wp))+SUM(REAL(xg,wp))
        IF (total>1.0_wp+tol) RETURN
        IF (total>1.0_wp) THEN
            ! Select the branch by real composition but differentiate
            ! the complete normalization including the denominator.
            total_complex=SUM(xs)+SUM(xg)
            xs=xs/total_complex; xg=xg/total_complex
        END IF
        close_composition_complex=.TRUE.
    END FUNCTION close_composition_complex
END MODULE state_composition
