MODULE state_new_composition
    USE parameters_2d, ONLY : wp, n_solid, n_add_gas
    IMPLICIT NONE
    PRIVATE
    REAL(wp), PARAMETER :: tol=4096.0_wp*EPSILON(1.0_wp)
    PUBLIC :: close_composition
    INTERFACE close_composition
        MODULE PROCEDURE close_composition_real, close_composition_complex
    END INTERFACE
CONTAINS
    LOGICAL FUNCTION close_composition_real(xs,xg)
        REAL(wp), INTENT(INOUT) :: xs(n_solid), xg(n_add_gas)
        REAL(wp) :: total
        close_composition_real=.FALSE.
        IF (ANY(xs < -tol) .OR. ANY(xg < -tol)) RETURN
        WHERE (ABS(xs)<=tol) xs=0.0_wp
        WHERE (ABS(xg)<=tol) xg=0.0_wp
        WHERE (xs<0.0_wp) xs=0.0_wp
        WHERE (xg<0.0_wp) xg=0.0_wp
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
        close_composition_complex=.FALSE.
        IF (ANY(REAL(xs,wp)<-tol) .OR. ANY(REAL(xg,wp)<-tol)) RETURN
        ! Match REAL composition closure without destroying imaginary
        ! complex-step perturbations of a nearly zero species.
        WHERE (ABS(REAL(xs,wp))<=tol) xs=CMPLX(0.0_wp,AIMAG(xs),wp)
        WHERE (ABS(REAL(xg,wp))<=tol) xg=CMPLX(0.0_wp,AIMAG(xg),wp)
        total=SUM(REAL(xs,wp))+SUM(REAL(xg,wp))
        IF (total>1.0_wp+tol) RETURN
        IF (total>1.0_wp) THEN
            xs=xs/total; xg=xg/total
        END IF
        close_composition_complex=.TRUE.
    END FUNCTION close_composition_complex
END MODULE state_new_composition
