    ! This executable fragment is instantiated by the REAL and COMPLEX
    ! phys_var cores. Keep type-specific declarations in the wrappers.

    alphal = 0.0_wp
    alphag = 0.0_wp
    Zs = 0.0_wp
    exc_pore_pres = 0.0_wp
    xl = 0.0_wp
    xg = 0.0_wp

    is_wet = REAL(qj(1), wp) .GT. EPSILON(1.0_wp)

    IF ( .NOT. is_wet ) THEN

       h = 0.0_wp
       u = 0.0_wp
       v = 0.0_wp
       alphas = 0.0_wp
       rho_m = rho_a_amb
       T = T_ambient
       inv_rhom = 1.0_wp / rho_m
       rho_c = rho_a_amb
       inv_rho_c = 1.0_wp / rho_c

       RETURN

    END IF

    inv_qj1 = 1.0_wp / qj(1)

    xs(1:n_solid) = qj(idx_alfas_first:idx_alfas_last) * inv_qj1

    IF ( n_add_gas .GT. 0 ) THEN

       xg(1:n_add_gas) = qj(idx_addGas_first:idx_addGas_last) * inv_qj1

    END IF

    IF ( stoch_transport_flag ) Zs = qj(idx_stoch) * inv_qj1
    IF ( pore_pressure_flag ) exc_pore_pres = qj(idx_pore) * inv_qj1

    IF ( gas_flag .AND. liquid_flag ) THEN

       xl = qj(n_vars) * inv_qj1

    END IF

    CALL eval_mixture_heat_capacity(xs, xg, xl, xc, mixture_sp_heat_c,          &
         sp_heat_mix)

    IF ( REAL(qj(1), wp) .GT. eps_sing ) THEN

       IF ( energy_flag ) THEN

          T = ( qj(4) - 0.5_wp * (qj(2)**2 + qj(3)**2) * inv_qj1 ) /           &
               ( qj(1) * sp_heat_mix )

       ELSE

          T = qj(4) / ( qj(1) * sp_heat_mix )

       END IF

       IF ( REAL(T, wp) .LE. 0.0_wp ) T = T_ambient

    ELSE

       T = T_ambient

    END IF

    CALL eval_mixture_from_mass_fractions(T, xs, xg, xl, rho_m, inv_rhom,       &
         rho_c, inv_rho_c, alphas, alphag, alphal)

    h = qj(1) * inv_rhom

    IF ( REAL(qj(1), wp) .GT. eps_sing ) THEN

       u = qj(2) * inv_qj1
       v = qj(3) * inv_qj1

    ELSE

       u = SQRT(2.0_wp) * qj(1) * qj(2) / SQRT(qj(1)**4 + eps_sing4)
       v = SQRT(2.0_wp) * qj(1) * qj(3) / SQRT(qj(1)**4 + eps_sing4)

    END IF
