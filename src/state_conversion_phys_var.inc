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
       xs_tot = 0.0_wp

       RETURN

    END IF

    inv_qj1 = 1.0_wp / qj(1)

    xs(1:n_solid) = qj(idx_alfas_first:idx_alfas_last) * inv_qj1

    IF ( n_add_gas .GT. 0 ) THEN

       xg(1:n_add_gas) = qj(idx_addGas_first:idx_addGas_last) * inv_qj1

    END IF

    IF ( stoch_transport_flag ) Zs = qj(idx_stoch) * inv_qj1
    IF ( pore_pressure_flag ) exc_pore_pres = qj(idx_pore) * inv_qj1

    xs_tot = SUM(xs)

    IF ( gas_flag .AND. liquid_flag ) THEN

       xl = qj(n_vars) * inv_qj1
       xc = 1.0_wp - xs_tot - xl

    ELSE

       xc = 1.0_wp - xs_tot

    END IF

    ! Do not use DOT_PRODUCT here or below: for a complex first argument it
    ! applies conjugation and therefore corrupts complex-step derivatives.
    IF ( gas_flag ) THEN

       carrier_sp_heat_mass =                                                   &
            ( xc - SUM(xg(1:n_add_gas)) ) * sp_heat_a +                        &
            SUM(xg(1:n_add_gas) * sp_heat_g(1:n_add_gas))

    ELSE

       carrier_sp_heat_mass = sp_heat_l * xc

    END IF

    sp_heat_mix = SUM(xs(1:n_solid) * sp_heat_s(1:n_solid))

    IF ( gas_flag .AND. liquid_flag ) THEN

       sp_heat_mix = sp_heat_mix + xl * sp_heat_l

    END IF

    ! carrier_sp_heat_mass is already weighted by the carrier-component mass
    ! fractions, so it enters the mixture heat capacity exactly once.
    sp_heat_mix = sp_heat_mix + carrier_sp_heat_mass

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

    IF ( gas_flag ) THEN

       ! Use the same carrier-fraction guard for REAL and COMPLEX arithmetic.
       ! Testing the real part keeps the selected branch unchanged by an
       ! infinitesimal complex-step perturbation.
       IF ( REAL(xc, wp) .GT. EPSILON(1.0_wp) ) THEN

          sp_gas_const_c =                                                      &
               ( (xc - SUM(xg(1:n_add_gas))) * sp_gas_const_a +                &
               SUM(xg(1:n_add_gas) *                                           &
               sp_gas_const_g(1:n_add_gas)) ) / xc

       ELSE

          sp_gas_const_c = sp_gas_const_a

       END IF

       inv_rho_c = sp_gas_const_c * T * inv_pres
       rho_c = 1.0_wp / inv_rho_c
       inv_rho_g(1:n_add_gas) = sp_gas_const_g(1:n_add_gas) * T * inv_pres

    ELSE

       rho_c = rho_l
       inv_rho_c = inv_rho_l

    END IF

    inv_rhom = SUM(xs(1:n_solid) * inv_rho_s(1:n_solid)) +                     &
         xc * inv_rho_c

    IF ( gas_flag .AND. liquid_flag ) inv_rhom = inv_rhom + xl * inv_rho_l

    rho_m = 1.0_wp / inv_rhom

    IF ( gas_flag .AND. liquid_flag ) alphal = xl * rho_m * inv_rho_l

    alphas(1:n_solid) = xs(1:n_solid) * rho_m * inv_rho_s(1:n_solid)
    alphag(1:n_add_gas) = xg(1:n_add_gas) * rho_m * inv_rho_g(1:n_add_gas)

    h = qj(1) * inv_rhom

    IF ( REAL(qj(1), wp) .GT. eps_sing ) THEN

       u = qj(2) * inv_qj1
       v = qj(3) * inv_qj1

    ELSE

       u = SQRT(2.0_wp) * qj(1) * qj(2) / SQRT(qj(1)**4 + eps_sing4)
       v = SQRT(2.0_wp) * qj(1) * qj(3) / SQRT(qj(1)**4 + eps_sing4)

    END IF
