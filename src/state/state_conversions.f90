!*******************************************************************************
!> Conservative/primitive conversions for the new state architecture.
!*******************************************************************************
MODULE state_conversions

    USE parameters_2d, ONLY : wp, n_vars, n_solid, n_add_gas
    USE parameters_2d, ONLY : idx_totMassEqn, idx_uEqn, idx_vEqn, idx_engyEqn
    USE parameters_2d, ONLY : idx_h, idx_T, idx_solid_first, idx_solid_last
    USE parameters_2d, ONLY : idx_add_gas_first, idx_add_gas_last, idx_u, idx_v
    USE parameters_2d, ONLY : idx_stoch, idx_pore
    USE parameters_2d, ONLY : eps_sing, eps_sing4
    USE parameters_2d, ONLY : stoch_transport_flag, pore_pressure_flag
    USE state_mixture, ONLY : evaluate_mixture, evaluate_mixture_heat_capacity, evaluate_mixture_density
    USE constitutive_parameters_2d, ONLY : T_ambient
    USE state_composition, ONLY : close_composition

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: conservative_to_primitive
    PUBLIC :: primitive_to_conservative
    PUBLIC :: dynamic_pressure_from_synchronized_state

    INTERFACE conservative_to_primitive
        MODULE PROCEDURE conservative_to_primitive_real
        MODULE PROCEDURE conservative_to_primitive_complex
    END INTERFACE conservative_to_primitive

    INTERFACE primitive_to_conservative
        MODULE PROCEDURE primitive_to_conservative_real
        MODULE PROCEDURE primitive_to_conservative_complex
    END INTERFACE primitive_to_conservative

CONTAINS

    !> Diagnostic dynamic pressure from already synchronized q and qp.
    !> The numerical q->qp conversion computed h=mass/rho_m and the final
    !> primitive velocities. Reconstruct rho_m from that identity instead of
    !> repeating heat capacity, density, composition and the whole conversion.
    !> The reconstructed density can differ in last-bit rounding from rho_m.
    PURE REAL(wp) FUNCTION dynamic_pressure_from_synchronized_state(q_mass,qp) RESULT(p_dyn)
        REAL(wp), INTENT(IN) :: q_mass,qp(n_vars)

        p_dyn=0.0_wp
        IF (q_mass<=EPSILON(1.0_wp) .OR. qp(idx_h)<=0.0_wp) RETURN
        p_dyn=0.5_wp*(q_mass/qp(idx_h))*(qp(idx_u)*qp(idx_u)+qp(idx_v)*qp(idx_v))
    END FUNCTION dynamic_pressure_from_synchronized_state


    !> REAL conservative-to-primitive conversion shared by numerical
    !> reconstruction and physical-state consumers.
    SUBROUTINE conservative_to_primitive_real(q,qp,accepted,p_dyn)
        REAL(wp), INTENT(IN) :: q(n_vars)
        REAL(wp), INTENT(OUT) :: qp(n_vars)
        LOGICAL, INTENT(OUT), OPTIONAL :: accepted
        REAL(wp), INTENT(OUT), OPTIONAL :: p_dyn
        REAL(wp) :: mass,rho_m,inv_rho_m,cp_m,cp_c,T,u,v
        REAL(wp) :: xs(n_solid),xg(n_add_gas)
        qp=0.0_wp
        IF (PRESENT(accepted)) accepted=.FALSE.
        IF (PRESENT(p_dyn)) p_dyn=0.0_wp
        mass=q(idx_totMassEqn)
        IF (mass<0.0_wp) RETURN
        IF (mass<=EPSILON(1.0_wp)) THEN
            qp(idx_T)=T_ambient
            IF (PRESENT(accepted)) accepted=.TRUE.
            RETURN
        END IF

        IF (REAL(q(idx_engyEqn),wp)<=0.0_wp) RETURN
        xs=q(idx_solid_first:idx_solid_last)/mass
        xg=q(idx_add_gas_first:idx_add_gas_last)/mass
        IF (.NOT.close_composition(xs,xg)) RETURN
        CALL evaluate_mixture_heat_capacity(xs,xg,cp_m,cp_c)
        IF (REAL(cp_m,wp)<=0.0_wp) RETURN
        INCLUDE 'inc/temperature_from_energy.inc'
        IF (REAL(T,wp)<=0.0_wp) RETURN
        CALL evaluate_mixture_density(T,xs,xg,rho_m,inv_rho_m)
        INCLUDE 'inc/velocity_from_momentum.inc'
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
        IF (REAL(mass,wp) <= EPSILON(1.0_wp)) THEN
            qp(idx_T)=T_ambient
            IF (PRESENT(accepted)) accepted=.TRUE.
            RETURN
        END IF

        IF (REAL(q(idx_engyEqn),wp)<=0.0_wp) RETURN
        xs=q(idx_solid_first:idx_solid_last)/mass
        xg=q(idx_add_gas_first:idx_add_gas_last)/mass
        IF (.NOT.close_composition(xs,xg)) RETURN
        CALL evaluate_mixture_heat_capacity(xs,xg,cp_m,cp_c)
        IF (REAL(cp_m,wp)<=0.0_wp) RETURN
        INCLUDE 'inc/temperature_from_energy.inc'
        IF (REAL(T,wp)<=0.0_wp) RETURN
        CALL evaluate_mixture_density(T,xs,xg,rho_m,inv_rho_m)
        INCLUDE 'inc/velocity_from_momentum.inc'
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
        IF (h <= EPSILON(1.0_wp)) THEN
            IF (PRESENT(accepted)) accepted=.TRUE.
            RETURN
        END IF
        IF (.NOT.close_composition(xs,xg)) RETURN
        IF (REAL(qp(idx_T),wp)<=0.0_wp) RETURN
        INCLUDE 'inc/conservative_from_primitive.inc'
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
        IF (REAL(qp(idx_h),wp) <= EPSILON(1.0_wp)) THEN
            IF (PRESENT(accepted)) accepted=.TRUE.
            RETURN
        END IF
        IF (.NOT.close_composition(xs,xg)) RETURN
        IF (REAL(qp(idx_T),wp)<=0.0_wp) RETURN
        INCLUDE 'inc/conservative_from_primitive.inc'
        IF (PRESENT(accepted)) accepted=.TRUE.
    END SUBROUTINE primitive_to_conservative_complex

END MODULE state_conversions
