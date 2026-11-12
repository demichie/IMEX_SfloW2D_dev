! Viscous and turbulent lahar drag, with frozen state, in either arithmetic kind.
rate_x=0.0_wp
rate_y=0.0_wp
IF (wp==sp) THEN
    depth_threshold=1.0E-10_wp
ELSE
    depth_threshold=1.0E-20_wp
END IF

fluid_viscosity=eval_lahar_fluid_viscosity(config,state%T,SUM(state%alphas))
IF (REAL(state%h,wp)>depth_threshold) THEN
    inverse_depth=1.0_wp/state%h
    viscous_slope=config%Kappa*fluid_viscosity*state%mod_vel*0.125_wp &
        *state%inv_rhom*inv_grav*inverse_depth**2
    turbulent_slope=(config%n_td**2)*state%mod_vel**2*inverse_depth**four_thirds
ELSE
    viscous_slope=config%Kappa*fluid_viscosity*state%mod_vel &
        /(8.0_wp*state%rho_m*grav*depth_threshold**2)
    turbulent_slope=(config%n_td**2)*state%mod_vel**2*depth_threshold**neg_four_thirds
END IF

total_slope=viscous_slope+turbulent_slope
IF (ABS(state%mod_vel)==0.0_wp) RETURN
friction_term=grav*state%rho_m*state%h*total_slope/state%mod_vel_hor
rate_x=-state%u*friction_term
rate_y=-state%v*friction_term
