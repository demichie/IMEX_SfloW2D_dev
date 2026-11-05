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
    USE state_new_composition, ONLY : close_composition, normalize_face_composition

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: conservative_to_primitive
    PUBLIC :: primitive_to_conservative, numerical_primitive_to_conservative

    INTERFACE conservative_to_primitive
        MODULE PROCEDURE conservative_to_primitive_real
        MODULE PROCEDURE conservative_to_primitive_complex
    END INTERFACE conservative_to_primitive

    INTERFACE primitive_to_conservative
        MODULE PROCEDURE primitive_to_conservative_real
        MODULE PROCEDURE primitive_to_conservative_complex
    END INTERFACE primitive_to_conservative

CONTAINS

    !> REAL state conversion, optionally retaining numerical reconstruction
    !> semantics (dry <= epsilon, ambient near-dry T, trace normalization).
    !> The default strict branch is unchanged for physical-state consumers.
    SUBROUTINE conservative_to_primitive_real(q,qp,accepted,p_dyn,numerical_reconstruction)
        REAL(wp), INTENT(IN) :: q(n_vars)
        REAL(wp), INTENT(OUT) :: qp(n_vars)
        LOGICAL, INTENT(OUT), OPTIONAL :: accepted
        REAL(wp), INTENT(OUT), OPTIONAL :: p_dyn
        LOGICAL, INTENT(IN), OPTIONAL :: numerical_reconstruction
        REAL(wp) :: mass,rho_m,inv_rho_m,cp_m,cp_c,T,u,v
        REAL(wp) :: xs(n_solid),xg(n_add_gas)
        REAL(wp) :: closed_xs(n_solid),closed_xg(n_add_gas)
        REAL(wp) :: inv_mass
        LOGICAL :: numerical

        numerical=.FALSE.
        IF (PRESENT(numerical_reconstruction)) numerical=numerical_reconstruction
        qp=0.0_wp
        IF (PRESENT(accepted)) accepted=.FALSE.
        IF (PRESENT(p_dyn)) p_dyn=0.0_wp
        mass=q(idx_totMassEqn)
        IF (mass<0.0_wp) RETURN
        IF ((numerical .AND. mass<=EPSILON(1.0_wp)) .OR. mass==0.0_wp) THEN
            qp(idx_T)=T_ambient
            IF (PRESENT(accepted)) accepted=.TRUE.
            RETURN
        END IF

        IF (numerical) THEN
            inv_mass=1.0_wp/mass
            xs=q(idx_solid_first:idx_solid_last)*inv_mass
            xg=q(idx_add_gas_first:idx_add_gas_last)*inv_mass
            ! Legacy numerical traces evaluate thermodynamics on the original
            ! fractions, then close the output fractions only.
            closed_xs=xs
            closed_xg=xg
            IF (.NOT.normalize_face_composition(closed_xs,closed_xg)) RETURN
        ELSE
            xs=q(idx_solid_first:idx_solid_last)/mass
            xg=q(idx_add_gas_first:idx_add_gas_last)/mass
            IF (.NOT.close_composition(xs,xg)) RETURN
            closed_xs=xs
            closed_xg=xg
        END IF
        CALL evaluate_mixture_heat_capacity(xs,xg,cp_m,cp_c)
        IF (cp_m<=0.0_wp) RETURN
        IF (numerical .AND. mass<=eps_sing) THEN
            T=T_ambient
        ELSE
            T=q(idx_engyEqn)/(mass*cp_m)
            IF (numerical .AND. T<=0.0_wp) T=T_ambient
        END IF
        CALL evaluate_mixture_density(T,xs,xg,rho_m,inv_rho_m)
        IF (mass>eps_sing) THEN
            u=q(idx_uEqn)/mass
            v=q(idx_vEqn)/mass
        ELSE
            u=SQRT(2.0_wp)*mass*q(idx_uEqn)/SQRT(mass**4+eps_sing4)
            v=SQRT(2.0_wp)*mass*q(idx_vEqn)/SQRT(mass**4+eps_sing4)
        END IF
        qp(idx_h)=mass*inv_rho_m
        qp(idx_T)=T
        qp(idx_solid_first:idx_solid_last)=closed_xs
        qp(idx_add_gas_first:idx_add_gas_last)=closed_xg
        IF (stoch_transport_flag) THEN
            IF (numerical) THEN
                qp(idx_stoch)=q(idx_stoch)*inv_mass
            ELSE
                qp(idx_stoch)=q(idx_stoch)/mass
            END IF
        END IF
        IF (pore_pressure_flag) THEN
            IF (numerical) THEN
                qp(idx_pore)=q(idx_pore)*inv_mass
            ELSE
                qp(idx_pore)=q(idx_pore)/mass
            END IF
        END IF
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
        IF (.NOT.close_composition(xs,xg)) RETURN
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
        IF (.NOT.close_composition(xs,xg)) RETURN
        INCLUDE 'conservative_from_primitive.inc'
        IF (PRESENT(accepted)) accepted=.TRUE.
    END SUBROUTINE primitive_to_conservative_real

    !> Convert a numerical primitive state (face trace or initial condition)
    !> to conservative variables with legacy numerical closure and dry threshold.
    !> Unlike strict bulk conversion, positive overfill is normalized.
    SUBROUTINE numerical_primitive_to_conservative(qp,q,accepted)
        REAL(wp), INTENT(IN) :: qp(n_vars)
        REAL(wp), INTENT(OUT) :: q(n_vars)
        LOGICAL, INTENT(OUT) :: accepted
        REAL(wp) :: rho_m,inv_rho_m,rho_c,inv_rho_c,cp_m,cp_c
        REAL(wp) :: alphas(n_solid),alphag(n_add_gas)
        REAL(wp) :: xs(n_solid),xg(n_add_gas),mass

        q=0.0_wp
        accepted=.FALSE.
        IF (qp(idx_h)<=EPSILON(1.0_wp)) THEN
            accepted=.TRUE.
            RETURN
        END IF
        xs=qp(idx_solid_first:idx_solid_last)
        xg=qp(idx_add_gas_first:idx_add_gas_last)
        IF (.NOT.normalize_face_composition(xs,xg)) RETURN
        INCLUDE 'conservative_from_primitive.inc'
        accepted=.TRUE.
    END SUBROUTINE numerical_primitive_to_conservative

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
        IF (.NOT.close_composition(xs,xg)) RETURN
        INCLUDE 'conservative_from_primitive.inc'
        IF (PRESENT(accepted)) accepted=.TRUE.
    END SUBROUTINE primitive_to_conservative_complex

END MODULE state_new_conversions
