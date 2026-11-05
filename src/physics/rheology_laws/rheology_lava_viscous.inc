! The frozen-temperature viscous contribution is shared between precisions.
rate_x=0.0_wp
rate_y=0.0_wp
IF (REAL(state%h,wp)>h_threshold) THEN
    gamma=3.0_wp*nu_ref/state%h*EXP(-visc_par*(state%T-T_ref))
ELSE
    gamma=3.0_wp*nu_ref/h_threshold*EXP(-visc_par*(state%T-T_ref))
END IF
IF (ABS(state%mod_vel)>0.0_wp) THEN
    rate_x=-state%rho_m*gamma*state%u
    rate_y=-state%rho_m*gamma*state%v
END IF
