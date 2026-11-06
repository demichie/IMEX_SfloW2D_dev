!> @file conservative_from_primitive.inc
!> @brief Compute conserved mass, momenta, thermal content and species from primitive variables.
    CALL evaluate_mixture(qp(idx_T),xs,xg,rho_m,inv_rho_m,rho_c,inv_rho_c,cp_m,cp_c,alphas,alphag)
    IF (REAL(cp_m,wp)<=0.0_wp) RETURN
    mass=rho_m*qp(idx_h)
    q(idx_totMassEqn)=mass
    q(idx_uEqn)=mass*qp(idx_u)
    q(idx_vEqn)=mass*qp(idx_v)
    q(idx_engyEqn)=mass*cp_m*qp(idx_T)
    q(idx_solid_first:idx_solid_last)=mass*xs
    q(idx_add_gas_first:idx_add_gas_last)=mass*xg
    IF (stoch_transport_flag) q(idx_stoch)=mass*qp(idx_stoch)
    IF (pore_pressure_flag) q(idx_pore)=mass*qp(idx_pore)
