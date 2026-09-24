    !--------------------------------------------------------------------------
    ! Shared executable core for the desingularized velocity recovered from a
    ! conservative state. It is included by the REAL and COMPLEX wrappers in
    ! state_conversion_2d so both paths use exactly the same near-dry formula.
    !--------------------------------------------------------------------------

    IF (REAL(qc(1), wp) .GT. eps_sing) THEN

       u = qc(2) / qc(1)
       v = qc(3) / qc(1)

    ELSE

       u = SQRT(2.0_wp) * qc(1) * qc(2) / SQRT(qc(1)**4 + eps_sing4)
       v = SQRT(2.0_wp) * qc(1) * qc(3) / SQRT(qc(1)**4 + eps_sing4)

    END IF
