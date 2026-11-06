!*******************************************************************************
!> Complete physical-state construction for the new state architecture.
!*******************************************************************************
MODULE state_physical

    USE parameters_2d, ONLY : wp, n_vars, n_solid, n_add_gas
    USE parameters_2d, ONLY : idx_h, idx_T, idx_solid_first, idx_solid_last
    USE parameters_2d, ONLY : idx_add_gas_first, idx_add_gas_last, idx_u, idx_v
    USE parameters_2d, ONLY : idx_stoch, idx_pore
    USE parameters_2d, ONLY : stoch_transport_flag, pore_pressure_flag, liquid_flag
    USE parameters_2d, ONLY : slope_correction_flag, curvature_term_flag
    USE parameters_2d, ONLY : eps_sing, eps_sing4, idx_totMassEqn, idx_uEqn, idx_vEqn, idx_engyEqn
    USE constitutive_parameters_2d, ONLY : rho_a_amb, sp_heat_a, sp_heat_l, grav, T_ambient
    USE state_types, ONLY : physical_state_real, physical_state_complex
    USE state_types, ONLY : reset_physical_state_real, reset_physical_state_complex
    USE state_composition, ONLY : close_composition
    USE state_mixture, ONLY : evaluate_mixture, evaluate_mixture_heat_capacity, &
        evaluate_mixture_phase_properties, evaluate_mixture_density, evaluate_reduced_gravity

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: build_physical_state_from_q
    PUBLIC :: build_physical_state_from_qp
    PUBLIC :: update_physical_state_kinematics
    PUBLIC :: refresh_physical_state
    PUBLIC :: evaluate_density_from_qp
    PUBLIC :: evaluate_temperature_from_q

    INTERFACE evaluate_temperature_from_q
        MODULE PROCEDURE evaluate_temperature_from_q_real
        MODULE PROCEDURE evaluate_temperature_from_q_complex
    END INTERFACE evaluate_temperature_from_q

    INTERFACE build_physical_state_from_q
        MODULE PROCEDURE build_physical_state_from_q_real
        MODULE PROCEDURE build_physical_state_from_q_complex
    END INTERFACE build_physical_state_from_q

    INTERFACE build_physical_state_from_qp
        MODULE PROCEDURE build_physical_state_from_qp_real
        MODULE PROCEDURE build_physical_state_from_qp_complex
    END INTERFACE build_physical_state_from_qp

    INTERFACE update_physical_state_kinematics
        MODULE PROCEDURE update_kinematics_real
        MODULE PROCEDURE update_kinematics_complex
    END INTERFACE update_physical_state_kinematics

    INTERFACE refresh_physical_state
        MODULE PROCEDURE refresh_physical_state_real
        MODULE PROCEDURE refresh_physical_state_complex
    END INTERFACE refresh_physical_state

CONTAINS

    !> Evaluate only mixture density and reduced gravity from a numerical
    !> primitive cell, without constructing the full physical state.
    !> The same composition closure and dry branch are used by
    !> build_physical_state_from_qp_real. The optional below_ambient output
    !> classifies the mixture against atmospheric density as a state property.
    SUBROUTINE evaluate_density_from_qp(qp,rho_m,red_grav,accepted,below_ambient)
        REAL(wp), INTENT(IN) :: qp(n_vars)
        REAL(wp), INTENT(OUT) :: rho_m,red_grav
        LOGICAL, INTENT(OUT) :: accepted
        LOGICAL, INTENT(OUT), OPTIONAL :: below_ambient
        REAL(wp) :: inv_rho_m,xs(n_solid),xg(n_add_gas)

        accepted=.FALSE.
        rho_m=rho_a_amb
        red_grav=0.0_wp
        IF (PRESENT(below_ambient)) below_ambient=.FALSE.
        xs=qp(idx_solid_first:idx_solid_last)
        xg=qp(idx_add_gas_first:idx_add_gas_last)
        IF (.NOT.close_composition(xs,xg)) RETURN

        IF (qp(idx_h)>EPSILON(1.0_wp)) THEN
            CALL evaluate_mixture_density(qp(idx_T),xs,xg,rho_m,inv_rho_m)
            red_grav=evaluate_reduced_gravity(rho_m,inv_rho_m)
        END IF
        ! Neutral-buoyancy cells are removed as well. This includes
        ! numerical near-dry states that return atmospheric density.
        IF (PRESENT(below_ambient)) below_ambient=rho_m<=rho_a_amb
        accepted=.TRUE.
    END SUBROUTINE evaluate_density_from_qp


    !> Temperature-only primitive from a conservative REAL candidate.
    !> Mirrors the thermodynamic part of build_physical_state_from_q_real
    !> without computing density, phase volumes or kinematics.
    SUBROUTINE evaluate_temperature_from_q_real(q,T)
        REAL(wp), INTENT(IN) :: q(n_vars)
        REAL(wp), INTENT(OUT) :: T
        REAL(wp) :: mass,cp_m,cp_c
        REAL(wp) :: xs(n_solid),xg(n_add_gas)

        T=0.0_wp
        mass=q(1)
        IF (mass<0.0_wp) RETURN
        IF (mass==0.0_wp) THEN
            T=T_ambient
            RETURN
        END IF
        xs=q(idx_solid_first:idx_solid_last)/mass
        xg=q(idx_add_gas_first:idx_add_gas_last)/mass
        IF (.NOT.close_composition(xs,xg)) RETURN
        CALL evaluate_mixture_heat_capacity(xs,xg,cp_m,cp_c)
        IF (REAL(cp_m,wp)<=0.0_wp) RETURN
        INCLUDE 'inc/temperature_from_energy.inc'
    END SUBROUTINE evaluate_temperature_from_q_real

    !> Complex-step temperature-only primitive. Preserve the same composition
    !> closure as build_physical_state_from_q_complex, including early rejection.
    SUBROUTINE evaluate_temperature_from_q_complex(q,T)
        COMPLEX(wp), INTENT(IN) :: q(n_vars)
        COMPLEX(wp), INTENT(OUT) :: T
        COMPLEX(wp) :: mass,cp_m,cp_c
        COMPLEX(wp) :: xs(n_solid),xg(n_add_gas)

        T=(0.0_wp,0.0_wp)
        mass=q(1)
        IF (REAL(mass,wp)<0.0_wp) RETURN
        IF (REAL(mass,wp)==0.0_wp) THEN
            T=CMPLX(T_ambient,0.0_wp,wp)
            RETURN
        END IF
        xs=q(idx_solid_first:idx_solid_last)/mass
        xg=q(idx_add_gas_first:idx_add_gas_last)/mass
        IF (.NOT.close_composition(xs,xg)) RETURN
        CALL evaluate_mixture_heat_capacity(xs,xg,cp_m,cp_c)
        IF (REAL(cp_m,wp)<=0.0_wp) RETURN
        INCLUDE 'inc/temperature_from_energy.inc'
    END SUBROUTINE evaluate_temperature_from_q_complex

    SUBROUTINE refresh_physical_state_real(state,bx,by,bxx,bxy,byy,accepted)
        TYPE(physical_state_real), INTENT(INOUT) :: state
        REAL(wp), INTENT(IN), OPTIONAL :: bx,by,bxx,bxy,byy
        LOGICAL, INTENT(OUT), OPTIONAL :: accepted
        REAL(wp) :: geometry(5)
        CALL require_real_storage(state)
        IF (PRESENT(accepted)) accepted=.FALSE.
        IF (.NOT.close_composition(state%xs,state%xg)) RETURN
        IF (REAL(state%h,wp)<0.0_wp) RETURN
        IF (REAL(state%h,wp)<=EPSILON(1.0_wp)) THEN
            INCLUDE 'inc/dry_thermodynamics.inc'
        ELSE
            IF (REAL(state%T,wp)<0.0_wp) RETURN
            CALL evaluate_mixture(state%T,state%xs,state%xg,state%rho_m,state%inv_rhom, &
                state%rho_c,state%inv_rho_c,state%cp_m,state%cp_c,state%alphas,state%alphag)
            IF (REAL(state%cp_m,wp)<=0.0_wp) RETURN
            state%red_grav=evaluate_reduced_gravity(state%rho_m,state%inv_rhom)
        END IF
        geometry=0.0_wp
        IF (PRESENT(bx)) geometry(1)=bx; IF (PRESENT(by)) geometry(2)=by
        IF (PRESENT(bxx)) geometry(3)=bxx; IF (PRESENT(bxy)) geometry(4)=bxy; IF (PRESENT(byy)) geometry(5)=byy
        CALL update_kinematics_real(state,geometry(1),geometry(2),geometry(3),geometry(4),geometry(5))
        IF (PRESENT(accepted)) accepted=.TRUE.
    END SUBROUTINE refresh_physical_state_real

    SUBROUTINE refresh_physical_state_complex(state,bx,by,bxx,bxy,byy,accepted)
        TYPE(physical_state_complex), INTENT(INOUT) :: state
        REAL(wp), INTENT(IN), OPTIONAL :: bx,by,bxx,bxy,byy
        LOGICAL, INTENT(OUT), OPTIONAL :: accepted
        REAL(wp) :: geometry(5)
        CALL require_complex_storage(state)
        IF (PRESENT(accepted)) accepted=.FALSE.
        IF (.NOT.close_composition(state%xs,state%xg)) RETURN
        IF (REAL(state%h,wp)<0.0_wp) RETURN
        IF (REAL(state%h,wp)<=EPSILON(1.0_wp)) THEN
            INCLUDE 'inc/dry_thermodynamics.inc'
        ELSE
            IF (REAL(state%T,wp)<0.0_wp) RETURN
            CALL evaluate_mixture(state%T,state%xs,state%xg,state%rho_m,state%inv_rhom, &
                state%rho_c,state%inv_rho_c,state%cp_m,state%cp_c,state%alphas,state%alphag)
            IF (REAL(state%cp_m,wp)<=0.0_wp) RETURN
            state%red_grav=evaluate_reduced_gravity(state%rho_m,state%inv_rhom)
        END IF
        geometry=0.0_wp
        IF (PRESENT(bx)) geometry(1)=bx; IF (PRESENT(by)) geometry(2)=by
        IF (PRESENT(bxx)) geometry(3)=bxx; IF (PRESENT(bxy)) geometry(4)=bxy; IF (PRESENT(byy)) geometry(5)=byy
        CALL update_kinematics_complex(state,geometry(1),geometry(2),geometry(3),geometry(4),geometry(5))
        IF (PRESENT(accepted)) accepted=.TRUE.
    END SUBROUTINE refresh_physical_state_complex

    SUBROUTINE require_real_storage(state)
        TYPE(physical_state_real), INTENT(IN) :: state
        IF (.NOT.state%layout_valid) ERROR STOP 'state REAL workspace was not validated at setup'
    END SUBROUTINE require_real_storage

    SUBROUTINE require_complex_storage(state)
        TYPE(physical_state_complex), INTENT(IN) :: state
        IF (.NOT.state%layout_valid) ERROR STOP 'state COMPLEX workspace was not validated at setup'
    END SUBROUTINE require_complex_storage

    SUBROUTINE build_physical_state_from_qp_real(qp, state, bx, by, bxx, bxy, byy, accepted)
        REAL(wp), INTENT(IN) :: qp(n_vars)
        TYPE(physical_state_real), INTENT(INOUT) :: state
        REAL(wp), INTENT(IN), OPTIONAL :: bx, by, bxx, bxy, byy
        LOGICAL, INTENT(OUT), OPTIONAL :: accepted
        REAL(wp) :: geometry(5)

        CALL require_real_storage(state)
        CALL reset_physical_state_real(state)
        IF (PRESENT(accepted)) accepted=.FALSE.
        state%h=qp(idx_h); state%u=qp(idx_u); state%v=qp(idx_v); state%T=qp(idx_T)
        state%xs=qp(idx_solid_first:idx_solid_last)
        state%xg=qp(idx_add_gas_first:idx_add_gas_last)
        IF (state%h < 0.0_wp) RETURN
        IF (.NOT.close_composition(state%xs,state%xg)) RETURN
        state%Zs=0.0_wp; state%exc_pore_pres=0.0_wp
        IF (stoch_transport_flag) state%Zs=qp(idx_stoch)
        IF (pore_pressure_flag) state%exc_pore_pres=qp(idx_pore)

        IF (state%h <= EPSILON(1.0_wp)) THEN
            INCLUDE 'inc/dry_thermodynamics.inc'
        ELSE
            CALL evaluate_mixture(state%T,state%xs,state%xg,state%rho_m,state%inv_rhom, &
                state%rho_c,state%inv_rho_c,state%cp_m,state%cp_c,state%alphas,state%alphag)
            IF (REAL(state%cp_m,wp)<=0.0_wp) RETURN
            state%red_grav=evaluate_reduced_gravity(state%rho_m,state%inv_rhom)
        END IF

        geometry=0.0_wp
        IF (PRESENT(bx)) geometry(1)=bx
        IF (PRESENT(by)) geometry(2)=by
        IF (PRESENT(bxx)) geometry(3)=bxx
        IF (PRESENT(bxy)) geometry(4)=bxy
        IF (PRESENT(byy)) geometry(5)=byy
        CALL update_kinematics_real(state,geometry(1),geometry(2),geometry(3),geometry(4),geometry(5))
        IF (PRESENT(accepted)) accepted=.TRUE.
    END SUBROUTINE build_physical_state_from_qp_real

    SUBROUTINE build_physical_state_from_q_real(q, state, bx, by, bxx, bxy, byy, accepted)
        REAL(wp), INTENT(IN) :: q(n_vars)
        TYPE(physical_state_real), INTENT(INOUT) :: state
        REAL(wp), INTENT(IN), OPTIONAL :: bx, by, bxx, bxy, byy
        LOGICAL, INTENT(OUT), OPTIONAL :: accepted
        REAL(wp) :: mass, rho_m, inv_rho_m, rho_c, inv_rho_c, cp_m, cp_c
        REAL(wp) :: xs(n_solid), xg(n_add_gas), geometry(5), T, u, v

        CALL require_real_storage(state)
        CALL reset_physical_state_real(state)
        IF (PRESENT(accepted)) accepted=.FALSE.
        mass=q(1)
        IF (mass < 0.0_wp) RETURN
        IF (mass <= 0.0_wp) THEN
            state%h=0.0_wp; state%u=0.0_wp; state%v=0.0_wp
            state%xs=0.0_wp; state%xg=0.0_wp
            state%Zs=0.0_wp; state%exc_pore_pres=0.0_wp
            INCLUDE 'inc/dry_thermodynamics.inc'
            geometry=0.0_wp
            IF (PRESENT(bx)) geometry(1)=bx; IF (PRESENT(by)) geometry(2)=by
            IF (PRESENT(bxx)) geometry(3)=bxx; IF (PRESENT(bxy)) geometry(4)=bxy; IF (PRESENT(byy)) geometry(5)=byy
            CALL update_kinematics_real(state,geometry(1),geometry(2),geometry(3),geometry(4),geometry(5))
            IF (PRESENT(accepted)) accepted=.TRUE.
            RETURN
        END IF
        xs=q(idx_solid_first:idx_solid_last)/mass
        xg=q(idx_add_gas_first:idx_add_gas_last)/mass
        IF (.NOT.close_composition(xs,xg)) RETURN
        CALL evaluate_mixture_heat_capacity(xs,xg,cp_m,cp_c)
        IF (REAL(cp_m,wp)<=0.0_wp) RETURN
        INCLUDE 'inc/temperature_from_energy.inc'
        state%T=T
        CALL evaluate_mixture_phase_properties(state%T,xs,xg,rho_m,inv_rho_m,rho_c,inv_rho_c, &
            state%alphas,state%alphag)
        state%h=mass*inv_rho_m
        INCLUDE 'inc/velocity_from_momentum.inc'
        state%u=u; state%v=v
        state%xs=xs; state%xg=xg; state%rho_m=rho_m; state%inv_rhom=inv_rho_m
        state%rho_c=rho_c; state%inv_rho_c=inv_rho_c; state%cp_m=cp_m; state%cp_c=cp_c
        state%Zs=0.0_wp; state%exc_pore_pres=0.0_wp
        IF (stoch_transport_flag) state%Zs=q(idx_stoch)/mass
        IF (pore_pressure_flag) state%exc_pore_pres=q(idx_pore)/mass
        state%red_grav=evaluate_reduced_gravity(rho_m,inv_rho_m)
        geometry=0.0_wp
        IF (PRESENT(bx)) geometry(1)=bx; IF (PRESENT(by)) geometry(2)=by
        IF (PRESENT(bxx)) geometry(3)=bxx; IF (PRESENT(bxy)) geometry(4)=bxy; IF (PRESENT(byy)) geometry(5)=byy
        CALL update_kinematics_real(state,geometry(1),geometry(2),geometry(3),geometry(4),geometry(5))
        IF (PRESENT(accepted)) accepted=.TRUE.
    END SUBROUTINE build_physical_state_from_q_real

    SUBROUTINE build_physical_state_from_qp_complex(qp, state, bx, by, bxx, bxy, byy, accepted)
        COMPLEX(wp), INTENT(IN) :: qp(n_vars)
        TYPE(physical_state_complex), INTENT(INOUT) :: state
        REAL(wp), INTENT(IN), OPTIONAL :: bx, by, bxx, bxy, byy
        LOGICAL, INTENT(OUT), OPTIONAL :: accepted
        REAL(wp) :: geometry(5)

        CALL require_complex_storage(state)
        CALL reset_physical_state_complex(state)
        IF (PRESENT(accepted)) accepted=.FALSE.
        state%h=qp(idx_h); state%u=qp(idx_u); state%v=qp(idx_v); state%T=qp(idx_T)
        state%xs=qp(idx_solid_first:idx_solid_last)
        state%xg=qp(idx_add_gas_first:idx_add_gas_last)
        IF (REAL(state%h,wp) < 0.0_wp) RETURN
        IF (.NOT.close_composition(state%xs,state%xg)) RETURN
        state%Zs=0.0_wp; state%exc_pore_pres=0.0_wp
        IF (stoch_transport_flag) state%Zs=qp(idx_stoch)
        IF (pore_pressure_flag) state%exc_pore_pres=qp(idx_pore)

        IF (REAL(state%h,wp) <= EPSILON(1.0_wp)) THEN
            INCLUDE 'inc/dry_thermodynamics.inc'
        ELSE
            CALL evaluate_mixture(state%T,state%xs,state%xg,state%rho_m,state%inv_rhom, &
                state%rho_c,state%inv_rho_c,state%cp_m,state%cp_c,state%alphas,state%alphag)
            IF (REAL(state%cp_m,wp)<=0.0_wp) RETURN
            state%red_grav=evaluate_reduced_gravity(state%rho_m,state%inv_rhom)
        END IF
        geometry=0.0_wp
        IF (PRESENT(bx)) geometry(1)=bx; IF (PRESENT(by)) geometry(2)=by
        IF (PRESENT(bxx)) geometry(3)=bxx; IF (PRESENT(bxy)) geometry(4)=bxy; IF (PRESENT(byy)) geometry(5)=byy
        CALL update_kinematics_complex(state,geometry(1),geometry(2),geometry(3),geometry(4),geometry(5))
        IF (PRESENT(accepted)) accepted=.TRUE.
    END SUBROUTINE build_physical_state_from_qp_complex

    SUBROUTINE build_physical_state_from_q_complex(q, state, bx, by, bxx, bxy, byy, accepted)
        COMPLEX(wp), INTENT(IN) :: q(n_vars)
        TYPE(physical_state_complex), INTENT(INOUT) :: state
        REAL(wp), INTENT(IN), OPTIONAL :: bx, by, bxx, bxy, byy
        LOGICAL, INTENT(OUT), OPTIONAL :: accepted
        COMPLEX(wp) :: mass, T, u, v, cp_m, cp_c, rho_m, inv_rho_m, rho_c, inv_rho_c
        COMPLEX(wp) :: xs(n_solid), xg(n_add_gas), alphas(n_solid), alphag(n_add_gas)
        REAL(wp) :: geometry(5)

        CALL require_complex_storage(state)
        CALL reset_physical_state_complex(state)
        IF (PRESENT(accepted)) accepted=.FALSE.
        mass=q(1); state%xs=0.0_wp; state%xg=0.0_wp
        state%Zs=0.0_wp; state%exc_pore_pres=0.0_wp
        IF (REAL(mass,wp) < 0.0_wp) RETURN
        IF (stoch_transport_flag .AND. REAL(mass,wp)>0.0_wp) state%Zs=q(idx_stoch)/mass
        IF (pore_pressure_flag .AND. REAL(mass,wp)>0.0_wp) state%exc_pore_pres=q(idx_pore)/mass
        IF (REAL(mass,wp)>0.0_wp) THEN
            xs=q(idx_solid_first:idx_solid_last)/mass; xg=q(idx_add_gas_first:idx_add_gas_last)/mass
            IF (.NOT.close_composition(xs,xg)) RETURN
            CALL evaluate_mixture_heat_capacity(xs,xg,cp_m,cp_c)
            IF (REAL(cp_m,wp)<=0.0_wp) RETURN
            INCLUDE 'inc/temperature_from_energy.inc'
            CALL evaluate_mixture_phase_properties(T,xs,xg,rho_m,inv_rho_m,rho_c,inv_rho_c, &
                alphas,alphag)
            state%h=mass*inv_rho_m
            INCLUDE 'inc/velocity_from_momentum.inc'
            state%u=u; state%v=v
            state%T=T
            state%xs=xs; state%xg=xg; state%alphas=alphas; state%alphag=alphag
            state%rho_m=rho_m; state%inv_rhom=inv_rho_m; state%rho_c=rho_c; state%inv_rho_c=inv_rho_c
            state%cp_m=cp_m; state%cp_c=cp_c; state%red_grav=evaluate_reduced_gravity(rho_m,inv_rho_m)
        ELSE
            state%h=0.0_wp; state%u=0.0_wp; state%v=0.0_wp
            INCLUDE 'inc/dry_thermodynamics.inc'
        END IF
        geometry=0.0_wp
        IF (PRESENT(bx)) geometry(1)=bx; IF (PRESENT(by)) geometry(2)=by
        IF (PRESENT(bxx)) geometry(3)=bxx; IF (PRESENT(bxy)) geometry(4)=bxy; IF (PRESENT(byy)) geometry(5)=byy
        CALL update_kinematics_complex(state,geometry(1),geometry(2),geometry(3),geometry(4),geometry(5))
        IF (PRESENT(accepted)) accepted=.TRUE.
    END SUBROUTINE build_physical_state_from_q_complex

    SUBROUTINE update_kinematics_real(state,bx,by,bxx,bxy,byy)
        TYPE(physical_state_real), INTENT(INOUT) :: state
        REAL(wp), INTENT(IN) :: bx,by,bxx,bxy,byy
        INCLUDE 'inc/kinematics.inc'
        state%Ri=richardson_real(state%h,state%u,state%v,state%red_grav)
    END SUBROUTINE update_kinematics_real

    SUBROUTINE update_kinematics_complex(state,bx,by,bxx,bxy,byy)
        TYPE(physical_state_complex), INTENT(INOUT) :: state
        REAL(wp), INTENT(IN) :: bx,by,bxx,bxy,byy
        INCLUDE 'inc/kinematics.inc'
        state%Ri=richardson_complex(state%h,state%u,state%v,state%red_grav)
    END SUBROUTINE update_kinematics_complex

    REAL(wp) FUNCTION richardson_real(h,u,v,red_grav)
        REAL(wp), INTENT(IN) :: h,u,v,red_grav
        REAL(wp) :: speed2, ri
        INCLUDE 'inc/richardson.inc'
        richardson_real=ri
    END FUNCTION richardson_real

    COMPLEX(wp) FUNCTION richardson_complex(h,u,v,red_grav)
        COMPLEX(wp), INTENT(IN) :: h,u,v,red_grav
        COMPLEX(wp) :: speed2, ri
        INCLUDE 'inc/richardson.inc'
        richardson_complex=ri
    END FUNCTION richardson_complex

END MODULE state_physical
