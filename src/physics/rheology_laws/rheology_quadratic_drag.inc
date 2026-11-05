! Single expression for the quadratic implicit contribution in both kinds.
rate_x=0.0_wp
rate_y=0.0_wp
IF (ABS(state%mod_vel)==0.0_wp) RETURN
rate_x=-state%rho_m*(state%u/state%mod_vel_hor) &
    *friction_factor*state%mod_vel2
rate_y=-state%rho_m*(state%v/state%mod_vel_hor) &
    *friction_factor*state%mod_vel2
