!> @brief Reusable REAL and COMPLEX physical-state workspace types.
!*******************************************************************************
!> Types for the clean-slate state architecture.
!
!> Canonical REAL/COMPLEX cell-local physical states, initialized for reuse.
!> Grid-sized authoritative fields are owned by state_fields.
!*******************************************************************************
MODULE state_types

    USE parameters_2d, ONLY : wp, n_solid, n_add_gas

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: physical_state_real
    PUBLIC :: physical_state_complex
    PUBLIC :: initialize_physical_state_real, initialize_physical_state_complex
    PUBLIC :: reset_physical_state_real
    PUBLIC :: reset_physical_state_complex
    PUBLIC :: promote_physical_state

    !> Complete cell-local physical state for real-valued evaluation.
    TYPE, PUBLIC :: physical_state_real
        !> Workspace allocation/layout validated at setup; never allocated in the cell loop.
        LOGICAL :: layout_valid=.FALSE.
        !> Thickness (m), horizontal velocity (m/s) and temperature (K).
        REAL(wp) :: h=0.0_wp, u=0.0_wp, v=0.0_wp, T=0.0_wp
        !> Bulk density and its inverse (kg/m3 and m3/kg).
        REAL(wp) :: rho_m=0.0_wp, inv_rhom=0.0_wp
        !> Carrier density and its inverse (kg/m3 and m3/kg).
        REAL(wp) :: rho_c=0.0_wp, inv_rho_c=0.0_wp
        !> Mixture/carrier heat capacity (J/kg/K).
        REAL(wp) :: cp_m=0.0_wp, cp_c=0.0_wp
        !> Reduced gravity (m/s2) and Richardson number (dimensionless).
        REAL(wp) :: red_grav=0.0_wp, Ri=0.0_wp
        !> Stochastic transport and excess pore-pressure state components.
        REAL(wp) :: Zs=0.0_wp, exc_pore_pres=0.0_wp
        !> Vertical velocity (m/s), horizontal speed (m/s) and speed squared (m2/s2).
        REAL(wp) :: w=0.0_wp, mod_vel_hor=0.0_wp, mod_vel2=0.0_wp
        !> Total speed (m/s) and bed-curvature acceleration (m/s2).
        REAL(wp) :: mod_vel=0.0_wp, curvature_acceleration=0.0_wp
        !> Mass and volume fraction arrays, sized n_solid and n_add_gas at setup.
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

CONTAINS

    !> @brief Prepare and validate REAL composition workspace arrays at setup.
    !> @details Allocate only at setup. An already allocated incompatible layout is fatal.
    !> @param[inout] state Initialized cell-local physical-state workspace; values are updated in
    !> place.
    SUBROUTINE initialize_physical_state_real(state)
        TYPE(physical_state_real), INTENT(INOUT) :: state

        IF (.NOT.ALLOCATED(state%xs) .OR. .NOT.ALLOCATED(state%xg) .OR.                                                &
            .NOT.ALLOCATED(state%alphas) .OR. .NOT.ALLOCATED(state%alphag)) THEN
            IF (ALLOCATED(state%xs)) DEALLOCATE(state%xs)
            IF (ALLOCATED(state%xg)) DEALLOCATE(state%xg)
            IF (ALLOCATED(state%alphas)) DEALLOCATE(state%alphas)
            IF (ALLOCATED(state%alphag)) DEALLOCATE(state%alphag)
            ALLOCATE(state%xs(n_solid), state%xg(n_add_gas),                                                           &
                state%alphas(n_solid), state%alphag(n_add_gas))
        ELSE IF (SIZE(state%xs) /= n_solid .OR.                                                                        &
                SIZE(state%xg) /= n_add_gas .OR.                                                                       &
                SIZE(state%alphas) /= n_solid .OR.                                                                     &
                SIZE(state%alphag) /= n_add_gas) THEN
            ERROR STOP 'state REAL workspace has incompatible composition layout'
        END IF

        state%xs=0.0_wp
        state%xg=0.0_wp
        state%alphas=0.0_wp
        state%alphag=0.0_wp
        state%layout_valid=.TRUE.
    END SUBROUTINE initialize_physical_state_real

    !> @brief Prepare and validate COMPLEX composition workspace arrays at setup.
    !> @details Allocate only at setup. An already allocated incompatible layout is fatal.
    !> @param[inout] state Initialized cell-local physical-state workspace; values are updated in
    !> place.
    SUBROUTINE initialize_physical_state_complex(state)
        TYPE(physical_state_complex), INTENT(INOUT) :: state

        IF (.NOT.ALLOCATED(state%xs) .OR. .NOT.ALLOCATED(state%xg) .OR.                                                &
            .NOT.ALLOCATED(state%alphas) .OR. .NOT.ALLOCATED(state%alphag)) THEN
            IF (ALLOCATED(state%xs)) DEALLOCATE(state%xs)
            IF (ALLOCATED(state%xg)) DEALLOCATE(state%xg)
            IF (ALLOCATED(state%alphas)) DEALLOCATE(state%alphas)
            IF (ALLOCATED(state%alphag)) DEALLOCATE(state%alphag)
            ALLOCATE(state%xs(n_solid), state%xg(n_add_gas),                                                           &
                state%alphas(n_solid), state%alphag(n_add_gas))
        ELSE IF (SIZE(state%xs) /= n_solid .OR.                                                                        &
                SIZE(state%xg) /= n_add_gas .OR.                                                                       &
                SIZE(state%alphas) /= n_solid .OR.                                                                     &
                SIZE(state%alphag) /= n_add_gas) THEN
            ERROR STOP 'state COMPLEX workspace has incompatible composition layout'
        END IF

        state%xs=(0.0_wp,0.0_wp)
        state%xg=(0.0_wp,0.0_wp)
        state%alphas=(0.0_wp,0.0_wp)
        state%alphag=(0.0_wp,0.0_wp)
        state%layout_valid=.TRUE.
    END SUBROUTINE initialize_physical_state_complex

    !> @brief Clear values in an initialized REAL workspace without reallocating.
    !> @param[inout] state Initialized cell-local physical-state workspace; values are updated in
    !> place.
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

    !> @brief Clear values in an initialized COMPLEX workspace without reallocating.
    !> @param[inout] state Initialized cell-local physical-state workspace; values are updated in
    !> place.
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

    !> @brief Promote a REAL physical state into a COMPLEX workspace.
    !> @details Requires a validated source; allocatable array assignment can allocate target.
    !> @param[in] source Initialized REAL physical state to promote without changing its values.
    !> @param[inout] target COMPLEX physical-state destination, including allocatable phase arrays.
    SUBROUTINE promote_physical_state(source,target)
        TYPE(physical_state_real), INTENT(IN) :: source
        TYPE(physical_state_complex), INTENT(INOUT) :: target

        IF (.NOT.source%layout_valid) ERROR STOP 'Cannot promote an uninitialized physical state'
        ! An explicit law may promote into a local complex state. Allocatable
        ! arrays follow intrinsic assignment semantics, as in the old kernel.
        target%layout_valid=.TRUE.
        target%h=CMPLX(source%h,0.0_wp,wp)
        target%u=CMPLX(source%u,0.0_wp,wp)
        target%v=CMPLX(source%v,0.0_wp,wp)
        target%T=CMPLX(source%T,0.0_wp,wp)
        target%rho_m=CMPLX(source%rho_m,0.0_wp,wp)
        target%inv_rhom=CMPLX(source%inv_rhom,0.0_wp,wp)
        target%rho_c=CMPLX(source%rho_c,0.0_wp,wp)
        target%inv_rho_c=CMPLX(source%inv_rho_c,0.0_wp,wp)
        target%cp_m=CMPLX(source%cp_m,0.0_wp,wp)
        target%cp_c=CMPLX(source%cp_c,0.0_wp,wp)
        target%red_grav=CMPLX(source%red_grav,0.0_wp,wp)
        target%Ri=CMPLX(source%Ri,0.0_wp,wp)
        target%Zs=CMPLX(source%Zs,0.0_wp,wp)
        target%exc_pore_pres=CMPLX(source%exc_pore_pres,0.0_wp,wp)
        target%w=CMPLX(source%w,0.0_wp,wp)
        target%mod_vel_hor=CMPLX(source%mod_vel_hor,0.0_wp,wp)
        target%mod_vel2=CMPLX(source%mod_vel2,0.0_wp,wp)
        target%mod_vel=CMPLX(source%mod_vel,0.0_wp,wp)
        target%curvature_acceleration=CMPLX(source%curvature_acceleration,0.0_wp,wp)
        target%xs=CMPLX(source%xs,0.0_wp,wp)
        target%xg=CMPLX(source%xg,0.0_wp,wp)
        target%alphas=CMPLX(source%alphas,0.0_wp,wp)
        target%alphag=CMPLX(source%alphag,0.0_wp,wp)
    END SUBROUTINE promote_physical_state

END MODULE state_types
