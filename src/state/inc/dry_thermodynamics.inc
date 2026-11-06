    state%T=T_ambient
    state%rho_m=rho_a_amb
    state%inv_rhom=1.0_wp/rho_a_amb
    state%rho_c=rho_a_amb
    state%inv_rho_c=1.0_wp/rho_a_amb
    IF (liquid_flag) THEN
        state%cp_m=sp_heat_l
        state%cp_c=sp_heat_l
    ELSE
        state%cp_m=sp_heat_a
        state%cp_c=sp_heat_a
    END IF
    state%alphas=0.0_wp
    state%alphag=0.0_wp
    state%red_grav=0.0_wp
