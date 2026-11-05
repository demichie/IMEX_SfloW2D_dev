!*******************************************************************************
!> Conservative/primitive conversions for the new state architecture.
!*******************************************************************************
MODULE state_new_conversions

    USE parameters_2d, ONLY : wp, n_vars, n_solid, n_add_gas
    USE parameters_2d, ONLY : idx_totMassEqn, idx_uEqn, idx_vEqn, idx_engyEqn
    USE parameters_2d, ONLY : idx_h, idx_T, idx_solid_first, idx_solid_last
    USE parameters_2d, ONLY : idx_add_gas_first, idx_add_gas_last, idx_u, idx_v
    USE parameters_2d, ONLY : idx_stoch, idx_pore
    USE parameters_2d, ONLY : eps_sing, eps_sing4
    USE parameters_2d, ONLY : stoch_transport_flag, pore_pressure_flag
    USE state_new_mixture, ONLY : evaluate_mixture, evaluate_mixture_heat_capacity, evaluate_mixture_density
    USE constitutive_parameters_2d, ONLY : T_ambient

    IMPLICIT NONE
    PRIVATE

    REAL(wp), PARAMETER :: closure_tolerance = 4096.0_wp*EPSILON(1.0_wp)

    PUBLIC :: conservative_to_primitive
    PUBLIC :: primitive_to_conservative

    INTERFACE conservative_to_primitive
        MODULE PROCEDURE conservative_to_primitive_real
        MODULE PROCEDURE conservative_to_primitive_complex
    END INTERFACE conservative_to_primitive

    INTERFACE primitive_to_conservative
        MODULE PROCEDURE primitive_to_conservative_real
        MODULE PROCEDURE primitive_to_conservative_complex
    END INTERFACE primitive_to_conservative

CONTAINS

    LOGICAL FUNCTION admissible_real_fractions(xs,xg)
        REAL(wp), INTENT(INOUT) :: xs(n_solid), xg(n_add_gas)
        REAL(wp) :: total
        admissible_real_fractions=.FALSE.
        IF (ANY(xs < -closure_tolerance) .OR. ANY(xg < -closure_tolerance)) RETURN
        WHERE (ABS(xs)<=closure_tolerance) xs=0.0_wp
        WHERE (ABS(xg)<=closure_tolerance) xg=0.0_wp
        WHERE (xs<0.0_wp) xs=0.0_wp
        WHERE (xg<0.0_wp) xg=0.0_wp
        total=SUM(xs)+SUM(xg)
        IF (total>1.0_wp+closure_tolerance) RETURN
        IF (total>1.0_wp) THEN
            IF (total<=0.0_wp) RETURN
            xs=xs/total; xg=xg/total
        END IF
        admissible_real_fractions=.TRUE.
    END FUNCTION admissible_real_fractions

    LOGICAL FUNCTION admissible_complex_fractions(xs,xg)
        COMPLEX(wp), INTENT(INOUT) :: xs(n_solid), xg(n_add_gas)
        REAL(wp) :: total
        admissible_complex_fractions=.FALSE.
        IF (ANY(REAL(xs,wp)<-closure_tolerance) .OR. ANY(REAL(xg,wp)<-closure_tolerance)) RETURN
        total=SUM(REAL(xs,wp))+SUM(REAL(xg,wp))
        IF (total>1.0_wp+closure_tolerance) RETURN
        IF (total>1.0_wp) THEN
            xs=xs/total; xg=xg/total
        END IF
        admissible_complex_fractions=.TRUE.
    END FUNCTION admissible_complex_fractions

    SUBROUTINE conservative_to_primitive_real(q, qp, accepted, p_dyn)
        REAL(wp), INTENT(IN) :: q(n_vars)
        REAL(wp), INTENT(OUT) :: qp(n_vars)
        LOGICAL, INTENT(OUT), OPTIONAL :: accepted
        REAL(wp), INTENT(OUT), OPTIONAL :: p_dyn
        REAL(wp) :: mass, rho_m, inv_rho_m, rho_c, inv_rho_c, cp_m, cp_c
        REAL(wp) :: T, u, v, alphas(n_solid), alphag(n_add_gas)
        REAL(wp) :: xs(n_solid), xg(n_add_gas)

        qp=0.0_wp
        IF (PRESENT(accepted)) accepted=.FALSE.
        IF (PRESENT(p_dyn)) p_dyn=0.0_wp
        mass=q(idx_totMassEqn)
        IF (mass < 0.0_wp) RETURN
        IF (mass == 0.0_wp) THEN
            qp(idx_T)=T_ambient
            IF (PRESENT(accepted)) accepted=.TRUE.
            RETURN
        END IF

        xs=q(idx_solid_first:idx_solid_last)/mass
        xg=q(idx_add_gas_first:idx_add_gas_last)/mass
        IF (.NOT.admissible_real_fractions(xs,xg)) RETURN
        CALL evaluate_mixture_heat_capacity(xs,xg,cp_m,cp_c)
        IF (cp_m <= 0.0_wp) RETURN
        T=q(idx_engyEqn)/(mass*cp_m)
        CALL evaluate_mixture_density(T,xs,xg,rho_m,inv_rho_m)
        IF (mass > eps_sing) THEN
            u=q(idx_uEqn)/mass
            v=q(idx_vEqn)/mass
        ELSE
            u=SQRT(2.0_wp)*mass*q(idx_uEqn)/SQRT(mass**4+eps_sing4)
            v=SQRT(2.0_wp)*mass*q(idx_vEqn)/SQRT(mass**4+eps_sing4)
        END IF

        qp(idx_h)=mass*inv_rho_m
        qp(idx_T)=T
        qp(idx_solid_first:idx_solid_last)=xs
        qp(idx_add_gas_first:idx_add_gas_last)=xg
        IF (stoch_transport_flag) qp(idx_stoch)=q(idx_stoch)/mass
        IF (pore_pressure_flag) qp(idx_pore)=q(idx_pore)/mass
        qp(idx_u)=u
        qp(idx_v)=v
        IF (PRESENT(p_dyn)) p_dyn=0.5_wp*rho_m*(u*u+v*v)
        IF (PRESENT(accepted)) accepted=.TRUE.
    END SUBROUTINE conservative_to_primitive_real

    SUBROUTINE conservative_to_primitive_complex(q, qp, accepted, p_dyn)
        COMPLEX(wp), INTENT(IN) :: q(n_vars)
        COMPLEX(wp), INTENT(OUT) :: qp(n_vars)
        LOGICAL, INTENT(OUT), OPTIONAL :: accepted
        COMPLEX(wp), INTENT(OUT), OPTIONAL :: p_dyn
        COMPLEX(wp) :: mass, rho_m, inv_rho_m, rho_c, inv_rho_c, cp_m, cp_c
        COMPLEX(wp) :: T, u, v, alphas(n_solid), alphag(n_add_gas)
        COMPLEX(wp) :: xs(n_solid), xg(n_add_gas)

        qp=(0.0_wp,0.0_wp)
        IF (PRESENT(accepted)) accepted=.FALSE.
        IF (PRESENT(p_dyn)) p_dyn=(0.0_wp,0.0_wp)
        mass=q(idx_totMassEqn)
        IF (REAL(mass,wp) < 0.0_wp) RETURN
        IF (REAL(mass,wp) == 0.0_wp) THEN
            qp(idx_T)=T_ambient
            IF (PRESENT(accepted)) accepted=.TRUE.
            RETURN
        END IF

        xs=q(idx_solid_first:idx_solid_last)/mass
        xg=q(idx_add_gas_first:idx_add_gas_last)/mass
        IF (.NOT.admissible_complex_fractions(xs,xg)) RETURN
        CALL evaluate_mixture_heat_capacity(xs,xg,cp_m,cp_c)
        T=q(idx_engyEqn)/(mass*cp_m)
        CALL evaluate_mixture_density(T,xs,xg,rho_m,inv_rho_m)
        IF (REAL(mass,wp) > eps_sing) THEN
            u=q(idx_uEqn)/mass
            v=q(idx_vEqn)/mass
        ELSE
            u=SQRT(2.0_wp)*mass*q(idx_uEqn)/SQRT(mass**4+eps_sing4)
            v=SQRT(2.0_wp)*mass*q(idx_vEqn)/SQRT(mass**4+eps_sing4)
        END IF
        qp(idx_h)=mass*inv_rho_m
        qp(idx_T)=T
        qp(idx_solid_first:idx_solid_last)=xs
        qp(idx_add_gas_first:idx_add_gas_last)=xg
        IF (stoch_transport_flag) qp(idx_stoch)=q(idx_stoch)/mass
        IF (pore_pressure_flag) qp(idx_pore)=q(idx_pore)/mass
        qp(idx_u)=u
        qp(idx_v)=v
        IF (PRESENT(p_dyn)) p_dyn=0.5_wp*rho_m*(u*u+v*v)
        IF (PRESENT(accepted)) accepted=.TRUE.
    END SUBROUTINE conservative_to_primitive_complex

    SUBROUTINE primitive_to_conservative_real(qp, q, accepted)
        REAL(wp), INTENT(IN) :: qp(n_vars)
        REAL(wp), INTENT(OUT) :: q(n_vars)
        LOGICAL, INTENT(OUT), OPTIONAL :: accepted
        REAL(wp) :: rho_m, inv_rho_m, rho_c, inv_rho_c, cp_m, cp_c
        REAL(wp) :: alphas(n_solid), alphag(n_add_gas)
        REAL(wp) :: xs(n_solid), xg(n_add_gas), h, mass
        xs=qp(idx_solid_first:idx_solid_last)
        xg=qp(idx_add_gas_first:idx_add_gas_last)
        h=qp(idx_h)
        q=0.0_wp
        IF (PRESENT(accepted)) accepted=.FALSE.
        IF (h < 0.0_wp) RETURN
        IF (h == 0.0_wp) THEN
            IF (PRESENT(accepted)) accepted=.TRUE.
            RETURN
        END IF
        IF (.NOT.admissible_real_fractions(xs,xg)) RETURN
        CALL evaluate_mixture(qp(idx_T),xs,xg,rho_m,inv_rho_m,rho_c,inv_rho_c,cp_m,cp_c,alphas,alphag)
        mass=rho_m*h
        q(idx_totMassEqn)=mass
        q(idx_uEqn)=mass*qp(idx_u)
        q(idx_vEqn)=mass*qp(idx_v)
        q(idx_engyEqn)=mass*cp_m*qp(idx_T)
        q(idx_solid_first:idx_solid_last)=mass*xs
        q(idx_add_gas_first:idx_add_gas_last)=mass*xg
        IF (stoch_transport_flag) q(idx_stoch)=mass*qp(idx_stoch)
        IF (pore_pressure_flag) q(idx_pore)=mass*qp(idx_pore)
        IF (PRESENT(accepted)) accepted=.TRUE.
    END SUBROUTINE primitive_to_conservative_real

    SUBROUTINE primitive_to_conservative_complex(qp, q, accepted)
        COMPLEX(wp), INTENT(IN) :: qp(n_vars)
        COMPLEX(wp), INTENT(OUT) :: q(n_vars)
        LOGICAL, INTENT(OUT), OPTIONAL :: accepted
        COMPLEX(wp) :: rho_m, inv_rho_m, rho_c, inv_rho_c, cp_m, cp_c
        COMPLEX(wp) :: alphas(n_solid), alphag(n_add_gas)
        COMPLEX(wp) :: xs(n_solid), xg(n_add_gas), mass
        xs=qp(idx_solid_first:idx_solid_last)
        xg=qp(idx_add_gas_first:idx_add_gas_last)
        q=(0.0_wp,0.0_wp)
        IF (PRESENT(accepted)) accepted=.FALSE.
        IF (REAL(qp(idx_h),wp) < 0.0_wp) RETURN
        IF (REAL(qp(idx_h),wp) == 0.0_wp) THEN
            IF (PRESENT(accepted)) accepted=.TRUE.
            RETURN
        END IF
        IF (.NOT.admissible_complex_fractions(xs,xg)) RETURN
        CALL evaluate_mixture(qp(idx_T),xs,xg,rho_m,inv_rho_m,rho_c,inv_rho_c,cp_m,cp_c,alphas,alphag)
        mass=rho_m*qp(idx_h)
        q(idx_totMassEqn)=mass
        q(idx_uEqn)=mass*qp(idx_u)
        q(idx_vEqn)=mass*qp(idx_v)
        q(idx_engyEqn)=mass*cp_m*qp(idx_T)
        q(idx_solid_first:idx_solid_last)=mass*xs
        q(idx_add_gas_first:idx_add_gas_last)=mass*xg
        IF (stoch_transport_flag) q(idx_stoch)=mass*qp(idx_stoch)
        IF (pore_pressure_flag) q(idx_pore)=mass*qp(idx_pore)
        IF (PRESENT(accepted)) accepted=.TRUE.
    END SUBROUTINE primitive_to_conservative_complex

END MODULE state_new_conversions
