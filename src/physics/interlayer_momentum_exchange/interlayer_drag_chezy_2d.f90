!> @brief Exact analytic solution of the quadratic Chezy interlayer momentum drag.
!> @details Integrates D=f*rho_2*|u_2-u_1|*(u_2-u_1) with constant
!>          layer masses and density during the drag-only Lie substep.
!>          The returned impulse is added to L1 and subtracted from L2.
MODULE interlayer_drag_chezy_2d
    USE parameters_2d, ONLY : wp
    USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: evaluate_interlayer_drag_chezy_analytic

CONTAINS

    !> @brief Calculate the exact paired momentum impulse for Chezy drag.
    !> @details No mass, composition, or thermal energy is exchanged.
    !>          Invalid physical inputs are rejected without a partial update.
    !> @param[in] mass_l1 L1 mixture mass per unit horizontal area [kg/m2].
    !> @param[in] mass_l2 L2 mixture mass per unit horizontal area [kg/m2].
    !> @param[in] rho_l2 L2 mixture density [kg/m3].
    !> @param[in] velocity_l1 L1 horizontal velocity [m/s].
    !> @param[in] velocity_l2 L2 horizontal velocity [m/s].
    !> @param[in] coefficient Nonnegative dimensionless Chezy interfacial coefficient.
    !> @param[in] dt Drag-only Lie substep [s].
    !> @param[out] impulse Momentum increment of L1, opposite for L2 [kg/(m s)].
    !> @param[out] accepted True for valid physical inputs and a finite impulse.
    PURE SUBROUTINE evaluate_interlayer_drag_chezy_analytic(mass_l1,mass_l2,rho_l2,velocity_l1,velocity_l2, &
        coefficient,dt,impulse,accepted)
        REAL(wp), INTENT(IN) :: mass_l1,mass_l2,rho_l2,velocity_l1(2),velocity_l2(2),coefficient,dt
        REAL(wp), INTENT(OUT) :: impulse(2)
        LOGICAL, INTENT(OUT) :: accepted
        REAL(wp) :: slip(2),slip_norm,decay,relaxation,reduced_mass

        impulse=0.0_wp
        accepted=.FALSE.
        IF (.NOT.ieee_is_finite(mass_l1) .OR. .NOT.ieee_is_finite(mass_l2) .OR. &
            .NOT.ieee_is_finite(rho_l2) .OR. .NOT.ieee_is_finite(coefficient) .OR. &
            .NOT.ieee_is_finite(dt) .OR. .NOT.ALL(ieee_is_finite(velocity_l1)) .OR. &
            .NOT.ALL(ieee_is_finite(velocity_l2))) RETURN
        IF (mass_l1<=0.0_wp .OR. mass_l2<=0.0_wp .OR. rho_l2<=0.0_wp) RETURN
        IF (coefficient<0.0_wp .OR. dt<0.0_wp) RETURN

        slip=velocity_l2-velocity_l1
        slip_norm=SQRT(SUM(slip**2))
        IF (coefficient==0.0_wp .OR. dt==0.0_wp .OR. slip_norm==0.0_wp) THEN
            accepted=.TRUE.
            RETURN
        END IF

        ! The reduced mass avoids subtracting nearly equal barycentric velocities.
        reduced_mass=mass_l1*(mass_l2/(mass_l1+mass_l2))
        relaxation=coefficient*rho_l2*(1.0_wp/mass_l1+1.0_wp/mass_l2)*slip_norm*dt
        IF (.NOT.ieee_is_finite(relaxation)) RETURN
        decay=relaxation/(1.0_wp+relaxation)
        impulse=reduced_mass*decay*slip
        accepted=ALL(ieee_is_finite(impulse))
        IF (.NOT.accepted) impulse=0.0_wp
    END SUBROUTINE evaluate_interlayer_drag_chezy_analytic
END MODULE interlayer_drag_chezy_2d
