    IF (REAL(mass,wp)>eps_sing) THEN
        u=q(idx_uEqn)/mass
        v=q(idx_vEqn)/mass
    ELSE
        u=SQRT(2.0_wp)*mass*q(idx_uEqn)/SQRT(mass**4+eps_sing4)
        v=SQRT(2.0_wp)*mass*q(idx_vEqn)/SQRT(mass**4+eps_sing4)
    END IF
