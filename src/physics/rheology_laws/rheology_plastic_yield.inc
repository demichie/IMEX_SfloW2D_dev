! Single expression for the plastic implicit contribution in both kinds.
rate_x=0.0_wp
rate_y=0.0_wp
IF (ABS(state%mod_vel)==0.0_wp) RETURN
rate_x=-state%rho_m*config%tau*(state%u/state%mod_vel)
rate_y=-state%rho_m*config%tau*(state%v/state%mod_vel)
