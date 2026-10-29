    !> \file velocity_from_conservative.inc
    !> \brief Shared REAL/COMPLEX near-dry velocity recovery.
    !> \details Included inside the typed velocity wrappers in state_conversion_2d.
    !> \par Inputs
    !> qc (M,M*u,M*v,...) and mass regularization constants eps_sing, eps_sing4.
    !> \par Outputs
    !> u and v [m s^-1]; both wrappers use identical desingularization algebra.
    !--------------------------------------------------------------------------
    ! Shared executable core for the desingularized velocity recovered from a
    ! conservative state. It is included by the REAL and COMPLEX wrappers in
    ! state_conversion_2d so both paths use exactly the same near-dry formula.
    !--------------------------------------------------------------------------

    IF (REAL(qc(idx_totMassEqn), wp) .GT. eps_sing) THEN

       u = qc(idx_uEqn) / qc(idx_totMassEqn)
       v = qc(idx_vEqn) / qc(idx_totMassEqn)

    ELSE

       ! Smoothly suppress velocity as M tends to zero; do not clip a momentum
       ! component independently or discard its complex-step perturbation.
       u = SQRT(2.0_wp) * qc(idx_totMassEqn) * qc(idx_uEqn) / SQRT(qc(idx_totMassEqn)**4 + eps_sing4)
       v = SQRT(2.0_wp) * qc(idx_totMassEqn) * qc(idx_vEqn) / SQRT(qc(idx_totMassEqn)**4 + eps_sing4)

    END IF
