!*******************************************************************************
!> Complete physical-state construction for the new state architecture.
!*******************************************************************************
MODULE state_new_physical

    USE parameters_2d, ONLY : wp, n_vars, n_solid, n_add_gas
    USE parameters_2d, ONLY : idx_h, idx_T, idx_solid_first, idx_solid_last
    USE parameters_2d, ONLY : idx_add_gas_first, idx_add_gas_last, idx_u, idx_v
    USE parameters_2d, ONLY : idx_stoch, idx_pore
    USE parameters_2d, ONLY : stoch_transport_flag, pore_pressure_flag
    USE parameters_2d, ONLY : slope_correction_flag, curvature_term_flag
    USE constitutive_parameters_2d, ONLY : rho_a_amb, sp_heat_a, sp_heat_l, grav
    USE state_new_types, ONLY : physical_state_real, physical_state_complex
    USE state_new_mixture, ONLY : evaluate_mixture

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: build_physical_state_from_q
    PUBLIC :: build_physical_state_from_qp
    PUBLIC :: update_physical_state_kinematics

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

CONTAINS

    SUBROUTINE require_real_storage(state)
        TYPE(physical_state_real), INTENT(IN) :: state
        IF (.NOT.ALLOCATED(state%xs) .OR. SIZE(state%xs)/=n_solid .OR. &
            .NOT.ALLOCATED(state%xg) .OR. SIZE(state%xg)/=n_add_gas .OR. &
            .NOT.ALLOCATED(state%alphas) .OR. SIZE(state%alphas)/=n_solid .OR. &
            .NOT.ALLOCATED(state%alphag) .OR. SIZE(state%alphag)/=n_add_gas) &
            ERROR STOP 'state_new REAL state must be initialized before evaluation'
    END SUBROUTINE require_real_storage

    SUBROUTINE require_complex_storage(state)
        TYPE(physical_state_complex), INTENT(IN) :: state
        IF (.NOT.ALLOCATED(state%xs) .OR. SIZE(state%xs)/=n_solid .OR. &
            .NOT.ALLOCATED(state%xg) .OR. SIZE(state%xg)/=n_add_gas .OR. &
            .NOT.ALLOCATED(state%alphas) .OR. SIZE(state%alphas)/=n_solid .OR. &
            .NOT.ALLOCATED(state%alphag) .OR. SIZE(state%alphag)/=n_add_gas) &
            ERROR STOP 'state_new COMPLEX state must be initialized before evaluation'
    END SUBROUTINE require_complex_storage

    SUBROUTINE build_physical_state_from_qp_real(qp, state, bx, by, bxx, bxy, byy)
        REAL(wp), INTENT(IN) :: qp(n_vars)
        TYPE(physical_state_real), INTENT(INOUT) :: state
        REAL(wp), INTENT(IN), OPTIONAL :: bx, by, bxx, bxy, byy
        REAL(wp) :: cp_m, cp_c, xc, geometry(5)

        CALL require_real_storage(state)
        state%h=qp(idx_h); state%u=qp(idx_u); state%v=qp(idx_v); state%T=qp(idx_T)
        state%xs=qp(idx_solid_first:idx_solid_last)
        state%xg=qp(idx_add_gas_first:idx_add_gas_last)
        state%Zs=0.0_wp; state%exc_pore_pres=0.0_wp
        IF (stoch_transport_flag) state%Zs=qp(idx_stoch)
        IF (pore_pressure_flag) state%exc_pore_pres=qp(idx_pore)

        IF (state%h <= EPSILON(1.0_wp)) THEN
            state%rho_m=rho_a_amb; state%inv_rhom=1.0_wp/rho_a_amb
            state%rho_c=rho_a_amb; state%inv_rho_c=1.0_wp/rho_a_amb
            IF (sp_heat_l > 0.0_wp) THEN
                state%cp_c=sp_heat_l; state%cp_m=sp_heat_l
            ELSE
                state%cp_c=sp_heat_a; state%cp_m=sp_heat_a
            END IF
            state%alphas=0.0_wp; state%alphag=0.0_wp; state%red_grav=0.0_wp; state%Ri=0.0_wp
        ELSE
            CALL evaluate_mixture(state%T,state%xs,state%xg,state%rho_m,state%inv_rhom, &
                state%rho_c,state%inv_rho_c,state%cp_m,state%cp_c,state%alphas,state%alphag)
            state%red_grav=(state%rho_m-rho_a_amb)*state%inv_rhom*grav
            state%Ri=richardson_real(state%h,state%u,state%v,state%red_grav)
        END IF

        geometry=0.0_wp
        IF (PRESENT(bx)) geometry(1)=bx
        IF (PRESENT(by)) geometry(2)=by
        IF (PRESENT(bxx)) geometry(3)=bxx
        IF (PRESENT(bxy)) geometry(4)=bxy
        IF (PRESENT(byy)) geometry(5)=byy
        CALL update_kinematics_real(state,geometry(1),geometry(2),geometry(3),geometry(4),geometry(5))
    END SUBROUTINE build_physical_state_from_qp_real

    SUBROUTINE build_physical_state_from_q_real(q, state, bx, by, bxx, bxy, byy)
        REAL(wp), INTENT(IN) :: q(n_vars)
        TYPE(physical_state_real), INTENT(INOUT) :: state
        REAL(wp), INTENT(IN), OPTIONAL :: bx, by, bxx, bxy, byy
        REAL(wp) :: qp(n_vars), mass, rho_m, inv_rho_m, rho_c, inv_rho_c, cp_m, cp_c
        REAL(wp) :: alphas(n_solid), alphag(n_add_gas), xs(n_solid), xg(n_add_gas)
        LOGICAL :: accepted

        CALL require_real_storage(state)
        mass=q(1); qp=0.0_wp
        IF (mass <= 0.0_wp) THEN
            CALL build_physical_state_from_qp_real(qp,state,bx,by,bxx,bxy,byy)
            RETURN
        END IF
        xs=q(idx_solid_first:idx_solid_last)/mass
        xg=q(idx_add_gas_first:idx_add_gas_last)/mass
        CALL evaluate_mixture(1.0_wp,xs,xg,rho_m,inv_rho_m,rho_c,inv_rho_c,cp_m,cp_c,alphas,alphag)
        qp(idx_h)=q(1)*0.0_wp
        qp(idx_T)=q(4)/(mass*cp_m)
        CALL evaluate_mixture(qp(idx_T),xs,xg,rho_m,inv_rho_m,rho_c,inv_rho_c,cp_m,cp_c,alphas,alphag)
        qp(idx_h)=mass*inv_rho_m
        qp(idx_u)=q(2)/mass; qp(idx_v)=q(3)/mass
        qp(idx_solid_first:idx_solid_last)=xs; qp(idx_add_gas_first:idx_add_gas_last)=xg
        IF (stoch_transport_flag) qp(idx_stoch)=q(idx_stoch)/mass
        IF (pore_pressure_flag) qp(idx_pore)=q(idx_pore)/mass
        accepted=.TRUE.
        CALL build_physical_state_from_qp_real(qp,state,bx,by,bxx,bxy,byy)
    END SUBROUTINE build_physical_state_from_q_real

    SUBROUTINE build_physical_state_from_qp_complex(qp, state, bx, by, bxx, bxy, byy)
        COMPLEX(wp), INTENT(IN) :: qp(n_vars)
        TYPE(physical_state_complex), INTENT(INOUT) :: state
        REAL(wp), INTENT(IN), OPTIONAL :: bx, by, bxx, bxy, byy
        REAL(wp) :: geometry(5)

        CALL require_complex_storage(state)
        state%h=qp(idx_h); state%u=qp(idx_u); state%v=qp(idx_v); state%T=qp(idx_T)
        state%xs=qp(idx_solid_first:idx_solid_last)
        state%xg=qp(idx_add_gas_first:idx_add_gas_last)
        state%Zs=0.0_wp; state%exc_pore_pres=0.0_wp
        IF (stoch_transport_flag) state%Zs=qp(idx_stoch)
        IF (pore_pressure_flag) state%exc_pore_pres=qp(idx_pore)

        IF (REAL(state%h,wp) <= EPSILON(1.0_wp)) THEN
            state%rho_m=rho_a_amb; state%inv_rhom=1.0_wp/rho_a_amb
            state%rho_c=rho_a_amb; state%inv_rho_c=1.0_wp/rho_a_amb
            state%cp_m=sp_heat_a; state%cp_c=sp_heat_a
            state%alphas=0.0_wp; state%alphag=0.0_wp; state%red_grav=0.0_wp
        ELSE
            CALL evaluate_mixture(state%T,state%xs,state%xg,state%rho_m,state%inv_rhom, &
                state%rho_c,state%inv_rho_c,state%cp_m,state%cp_c,state%alphas,state%alphag)
            state%red_grav=(state%rho_m-rho_a_amb)*state%inv_rhom*grav
        END IF
        state%Ri=richardson_complex(state%h,state%u,state%v,state%red_grav)
        geometry=0.0_wp
        IF (PRESENT(bx)) geometry(1)=bx; IF (PRESENT(by)) geometry(2)=by
        IF (PRESENT(bxx)) geometry(3)=bxx; IF (PRESENT(bxy)) geometry(4)=bxy; IF (PRESENT(byy)) geometry(5)=byy
        CALL update_kinematics_complex(state,geometry(1),geometry(2),geometry(3),geometry(4),geometry(5))
    END SUBROUTINE build_physical_state_from_qp_complex

    SUBROUTINE build_physical_state_from_q_complex(q, state, bx, by, bxx, bxy, byy)
        COMPLEX(wp), INTENT(IN) :: q(n_vars)
        TYPE(physical_state_complex), INTENT(INOUT) :: state
        REAL(wp), INTENT(IN), OPTIONAL :: bx, by, bxx, bxy, byy
        COMPLEX(wp) :: mass, T, cp_m, cp_c, rho_m, inv_rho_m, rho_c, inv_rho_c
        COMPLEX(wp) :: xs(n_solid), xg(n_add_gas), alphas(n_solid), alphag(n_add_gas)
        REAL(wp) :: geometry(5)

        CALL require_complex_storage(state)
        mass=q(1); state%xs=0.0_wp; state%xg=0.0_wp
        IF (REAL(mass,wp)>EPSILON(1.0_wp)) THEN
            xs=q(idx_solid_first:idx_solid_last)/mass; xg=q(idx_add_gas_first:idx_add_gas_last)/mass
            CALL evaluate_mixture((1.0_wp,0.0_wp),xs,xg,rho_m,inv_rho_m,rho_c,inv_rho_c,cp_m,cp_c,alphas,alphag)
            T=q(4)/(mass*cp_m)
            CALL evaluate_mixture(T,xs,xg,rho_m,inv_rho_m,rho_c,inv_rho_c,cp_m,cp_c,alphas,alphag)
            state%h=mass*inv_rho_m; state%u=q(2)/mass; state%v=q(3)/mass; state%T=T
            state%xs=xs; state%xg=xg; state%alphas=alphas; state%alphag=alphag
            state%rho_m=rho_m; state%inv_rhom=inv_rho_m; state%rho_c=rho_c; state%inv_rho_c=inv_rho_c
            state%cp_m=cp_m; state%cp_c=cp_c; state%red_grav=(rho_m-rho_a_amb)*inv_rho_m*grav
        ELSE
            state%h=0.0_wp; state%u=0.0_wp; state%v=0.0_wp; state%T=0.0_wp
            state%rho_m=rho_a_amb; state%inv_rhom=1.0_wp/rho_a_amb; state%rho_c=rho_a_amb
            state%inv_rho_c=1.0_wp/rho_a_amb; state%cp_m=sp_heat_a; state%cp_c=sp_heat_a
            state%alphas=0.0_wp; state%alphag=0.0_wp; state%red_grav=0.0_wp
        END IF
        state%Ri=richardson_complex(state%h,state%u,state%v,state%red_grav)
        geometry=0.0_wp
        IF (PRESENT(bx)) geometry(1)=bx; IF (PRESENT(by)) geometry(2)=by
        IF (PRESENT(bxx)) geometry(3)=bxx; IF (PRESENT(bxy)) geometry(4)=bxy; IF (PRESENT(byy)) geometry(5)=byy
        CALL update_kinematics_complex(state,geometry(1),geometry(2),geometry(3),geometry(4),geometry(5))
    END SUBROUTINE build_physical_state_from_q_complex

    SUBROUTINE update_kinematics_real(state,bx,by,bxx,bxy,byy)
        TYPE(physical_state_real), INTENT(INOUT) :: state
        REAL(wp), INTENT(IN) :: bx,by,bxx,bxy,byy
        INCLUDE 'kinematics.inc'
    END SUBROUTINE update_kinematics_real

    SUBROUTINE update_kinematics_complex(state,bx,by,bxx,bxy,byy)
        TYPE(physical_state_complex), INTENT(INOUT) :: state
        REAL(wp), INTENT(IN) :: bx,by,bxx,bxy,byy
        INCLUDE 'kinematics.inc'
    END SUBROUTINE update_kinematics_complex

    REAL(wp) FUNCTION richardson_real(h,u,v,red_grav)
        REAL(wp), INTENT(IN) :: h,u,v,red_grav
        REAL(wp) :: speed2
        speed2=u**2+v**2; richardson_real=0.0_wp
        IF (speed2>0.0_wp) richardson_real=MIN(1.0E15_wp,MAX(0.0_wp,red_grav*h/speed2))
    END FUNCTION richardson_real

    COMPLEX(wp) FUNCTION richardson_complex(h,u,v,red_grav)
        COMPLEX(wp), INTENT(IN) :: h,u,v,red_grav
        COMPLEX(wp) :: speed2, ri
        INCLUDE 'richardson.inc'
        richardson_complex=ri
    END FUNCTION richardson_complex

END MODULE state_new_physical
