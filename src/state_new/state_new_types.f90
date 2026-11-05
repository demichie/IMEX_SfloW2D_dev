!*******************************************************************************
!> Types for the clean-slate state architecture.
!
!> This module is intentionally parallel to state_conversion_2d.  It contains
!> only state data ownership and workspace lifetime; numerical consumers are
!> not migrated by this first step.
!*******************************************************************************
MODULE state_new_types

    USE parameters_2d, ONLY : wp, n_solid, n_add_gas

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: physical_state_real
    PUBLIC :: physical_state_complex
    PUBLIC :: state_workspace_real
    PUBLIC :: state_workspace_complex
    PUBLIC :: initialize_state_workspace_real
    PUBLIC :: initialize_state_workspace_complex
    PUBLIC :: reset_physical_state_real
    PUBLIC :: reset_physical_state_complex

    !> Complete cell-local physical state for real-valued evaluation.
    TYPE, PUBLIC :: physical_state_real
        LOGICAL :: layout_valid=.FALSE.
        REAL(wp) :: h=0.0_wp, u=0.0_wp, v=0.0_wp, T=0.0_wp
        REAL(wp) :: rho_m=0.0_wp, inv_rhom=0.0_wp
        REAL(wp) :: rho_c=0.0_wp, inv_rho_c=0.0_wp
        REAL(wp) :: cp_m=0.0_wp, cp_c=0.0_wp
        REAL(wp) :: red_grav=0.0_wp, Ri=0.0_wp
        REAL(wp) :: Zs=0.0_wp, exc_pore_pres=0.0_wp
        REAL(wp) :: w=0.0_wp, mod_vel_hor=0.0_wp, mod_vel2=0.0_wp
        REAL(wp) :: mod_vel=0.0_wp, curvature_acceleration=0.0_wp
        REAL(wp), ALLOCATABLE :: xs(:), xg(:), alphas(:), alphag(:)
    END TYPE physical_state_real

    !> Complete cell-local physical state for complex-step evaluation.
    TYPE, PUBLIC :: physical_state_complex
        LOGICAL :: layout_valid=.FALSE.
        COMPLEX(wp) :: h=(0.0_wp,0.0_wp), u=(0.0_wp,0.0_wp)
        COMPLEX(wp) :: v=(0.0_wp,0.0_wp), T=(0.0_wp,0.0_wp)
        COMPLEX(wp) :: rho_m=(0.0_wp,0.0_wp), inv_rhom=(0.0_wp,0.0_wp)
        COMPLEX(wp) :: rho_c=(0.0_wp,0.0_wp), inv_rho_c=(0.0_wp,0.0_wp)
        COMPLEX(wp) :: cp_m=(0.0_wp,0.0_wp), cp_c=(0.0_wp,0.0_wp)
        COMPLEX(wp) :: red_grav=(0.0_wp,0.0_wp), Ri=(0.0_wp,0.0_wp)
        COMPLEX(wp) :: Zs=(0.0_wp,0.0_wp), exc_pore_pres=(0.0_wp,0.0_wp)
        COMPLEX(wp) :: w=(0.0_wp,0.0_wp), mod_vel_hor=(0.0_wp,0.0_wp)
        COMPLEX(wp) :: mod_vel2=(0.0_wp,0.0_wp), mod_vel=(0.0_wp,0.0_wp)
        COMPLEX(wp) :: curvature_acceleration=(0.0_wp,0.0_wp)
        COMPLEX(wp), ALLOCATABLE :: xs(:), xg(:), alphas(:), alphag(:)
    END TYPE physical_state_complex

    !> Reusable per-worker scratch for real state construction.
    TYPE, PUBLIC :: state_workspace_real
        TYPE(physical_state_real) :: state
    END TYPE state_workspace_real

    !> Reusable per-worker scratch for complex-step state construction.
    TYPE, PUBLIC :: state_workspace_complex
        TYPE(physical_state_complex) :: state
    END TYPE state_workspace_complex

CONTAINS

    SUBROUTINE reset_physical_state_real(state)
        TYPE(physical_state_real), INTENT(INOUT) :: state
        LOGICAL :: layout_valid
        layout_valid=state%layout_valid
        state%h=0.0_wp; state%u=0.0_wp; state%v=0.0_wp; state%T=0.0_wp
        state%rho_m=0.0_wp; state%inv_rhom=0.0_wp; state%rho_c=0.0_wp; state%inv_rho_c=0.0_wp
        state%cp_m=0.0_wp; state%cp_c=0.0_wp; state%red_grav=0.0_wp; state%Ri=0.0_wp
        state%Zs=0.0_wp; state%exc_pore_pres=0.0_wp; state%w=0.0_wp
        state%mod_vel_hor=0.0_wp; state%mod_vel2=0.0_wp; state%mod_vel=0.0_wp
        state%curvature_acceleration=0.0_wp; state%layout_valid=layout_valid
        IF (ALLOCATED(state%xs)) state%xs=0.0_wp
        IF (ALLOCATED(state%xg)) state%xg=0.0_wp
        IF (ALLOCATED(state%alphas)) state%alphas=0.0_wp
        IF (ALLOCATED(state%alphag)) state%alphag=0.0_wp
    END SUBROUTINE reset_physical_state_real

    SUBROUTINE reset_physical_state_complex(state)
        TYPE(physical_state_complex), INTENT(INOUT) :: state
        LOGICAL :: layout_valid
        layout_valid=state%layout_valid
        state%h=0.0_wp; state%u=0.0_wp; state%v=0.0_wp; state%T=0.0_wp
        state%rho_m=0.0_wp; state%inv_rhom=0.0_wp; state%rho_c=0.0_wp; state%inv_rho_c=0.0_wp
        state%cp_m=0.0_wp; state%cp_c=0.0_wp; state%red_grav=0.0_wp; state%Ri=0.0_wp
        state%Zs=0.0_wp; state%exc_pore_pres=0.0_wp; state%w=0.0_wp
        state%mod_vel_hor=0.0_wp; state%mod_vel2=0.0_wp; state%mod_vel=0.0_wp
        state%curvature_acceleration=0.0_wp; state%layout_valid=layout_valid
        IF (ALLOCATED(state%xs)) state%xs=0.0_wp
        IF (ALLOCATED(state%xg)) state%xg=0.0_wp
        IF (ALLOCATED(state%alphas)) state%alphas=0.0_wp
        IF (ALLOCATED(state%alphag)) state%alphag=0.0_wp
    END SUBROUTINE reset_physical_state_complex

    !> Allocate and clear the arrays owned by a reusable REAL state workspace.
    !> This routine belongs to setup/configuration, never to a cell-evaluation loop.
    SUBROUTINE initialize_state_workspace_real(workspace)
        TYPE(state_workspace_real), INTENT(INOUT) :: workspace

        IF (.NOT.ALLOCATED(workspace%state%xs) .OR. .NOT.ALLOCATED(workspace%state%xg) .OR. &
            .NOT.ALLOCATED(workspace%state%alphas) .OR. .NOT.ALLOCATED(workspace%state%alphag)) THEN
            IF (ALLOCATED(workspace%state%xs)) DEALLOCATE(workspace%state%xs)
            IF (ALLOCATED(workspace%state%xg)) DEALLOCATE(workspace%state%xg)
            IF (ALLOCATED(workspace%state%alphas)) DEALLOCATE(workspace%state%alphas)
            IF (ALLOCATED(workspace%state%alphag)) DEALLOCATE(workspace%state%alphag)
            ALLOCATE(workspace%state%xs(n_solid), workspace%state%xg(n_add_gas), &
                workspace%state%alphas(n_solid), workspace%state%alphag(n_add_gas))
        ELSE IF (SIZE(workspace%state%xs) /= n_solid .OR. &
                 SIZE(workspace%state%xg) /= n_add_gas) THEN
            ERROR STOP 'state_new REAL workspace has incompatible composition layout'
        END IF

        workspace%state%xs=0.0_wp
        workspace%state%xg=0.0_wp
        workspace%state%alphas=0.0_wp
        workspace%state%alphag=0.0_wp
        workspace%state%layout_valid=.TRUE.
    END SUBROUTINE initialize_state_workspace_real

    !> Allocate and clear the arrays owned by a reusable COMPLEX state workspace.
    !> This routine belongs to setup/configuration, never to a cell-evaluation loop.
    SUBROUTINE initialize_state_workspace_complex(workspace)
        TYPE(state_workspace_complex), INTENT(INOUT) :: workspace

        IF (.NOT.ALLOCATED(workspace%state%xs) .OR. .NOT.ALLOCATED(workspace%state%xg) .OR. &
            .NOT.ALLOCATED(workspace%state%alphas) .OR. .NOT.ALLOCATED(workspace%state%alphag)) THEN
            IF (ALLOCATED(workspace%state%xs)) DEALLOCATE(workspace%state%xs)
            IF (ALLOCATED(workspace%state%xg)) DEALLOCATE(workspace%state%xg)
            IF (ALLOCATED(workspace%state%alphas)) DEALLOCATE(workspace%state%alphas)
            IF (ALLOCATED(workspace%state%alphag)) DEALLOCATE(workspace%state%alphag)
            ALLOCATE(workspace%state%xs(n_solid), workspace%state%xg(n_add_gas), &
                workspace%state%alphas(n_solid), workspace%state%alphag(n_add_gas))
        ELSE IF (SIZE(workspace%state%xs) /= n_solid .OR. &
                 SIZE(workspace%state%xg) /= n_add_gas) THEN
            ERROR STOP 'state_new COMPLEX workspace has incompatible composition layout'
        END IF

        workspace%state%xs=(0.0_wp,0.0_wp)
        workspace%state%xg=(0.0_wp,0.0_wp)
        workspace%state%alphas=(0.0_wp,0.0_wp)
        workspace%state%alphag=(0.0_wp,0.0_wp)
        workspace%state%layout_valid=.TRUE.
    END SUBROUTINE initialize_state_workspace_complex

END MODULE state_new_types
