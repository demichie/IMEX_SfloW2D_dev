! Empirical drag compiled in REAL for explicit assembly and COMPLEX for Newton.
rate_x=0.0_wp
rate_y=0.0_wp
IF (ABS(state%mod_vel)==0.0_wp) RETURN
friction_coefficient=1.0E-3_wp/(1.0_wp+10.0_wp*state%h)*state%mod_vel
rate_x=-state%rho_m*friction_coefficient*(state%u/state%mod_vel_hor)
rate_y=-state%rho_m*friction_coefficient*(state%v/state%mod_vel_hor)
