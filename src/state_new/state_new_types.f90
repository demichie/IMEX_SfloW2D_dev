!*******************************************************************************
!> Types for the clean-slate state architecture.
!
!> This module is intentionally parallel to state_conversion_2d.  It contains
!> only state data ownership and workspace lifetime; numerical consumers are
!> not migrated by this first step.
!*******************************************************************************
MODULE state_new_types

    USE parameters_2d, ONLY : wp

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: physical_state_real
    PUBLIC :: physical_state_complex
    PUBLIC :: state_workspace_real
    PUBLIC :: state_workspace_complex

    !> Complete cell-local physical state for real-valued evaluation.
    TYPE, PUBLIC :: physical_state_real
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
        REAL(wp), ALLOCATABLE :: xs(:), xg(:), alphas(:), alphag(:)
    END TYPE state_workspace_real

    !> Reusable per-worker scratch for complex-step state construction.
    TYPE, PUBLIC :: state_workspace_complex
        TYPE(physical_state_complex) :: state
        COMPLEX(wp), ALLOCATABLE :: xs(:), xg(:), alphas(:), alphag(:)
    END TYPE state_workspace_complex

END MODULE state_new_types
