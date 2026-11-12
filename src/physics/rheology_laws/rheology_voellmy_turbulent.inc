! The same turbulence law is evaluated in REAL for explicit rates and
! in COMPLEX for Newton. Branch only on a zero velocity norm.
rate_x=0.0_wp
rate_y=0.0_wp
IF (ABS(state%mod_vel)==0.0_wp) RETURN

IF (use_stochastic_xi) THEN
    xi_effective=config%xi+stochastic_value
    IF (xi_effective<1.0_wp) xi_effective=1.0_wp
    IF (xi_effective>1.0E5_wp) xi_effective=1.0E5_wp
ELSE
    xi_effective=config%xi
END IF

turbulent_stress=state%rho_m*state%red_grav/xi_effective*state%mod_vel2
rate_x=-turbulent_stress*(state%u/state%mod_vel)
rate_y=-turbulent_stress*(state%v/state%mod_vel)
