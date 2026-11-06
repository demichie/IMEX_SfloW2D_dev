    CALL evaluate_mixture(state%T,state%xs,state%xg,state%rho_m,state%inv_rhom, &
        state%rho_c,state%inv_rho_c,state%cp_m,state%cp_c,state%alphas,state%alphag)
    IF (REAL(state%cp_m,wp)<=0.0_wp) RETURN
    state%red_grav=evaluate_reduced_gravity(state%rho_m,state%inv_rhom)
