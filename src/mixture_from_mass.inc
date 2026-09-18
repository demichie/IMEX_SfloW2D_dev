    xc = 1.0_wp - SUM(xs)
    IF ( gas_flag .AND. liquid_flag ) xc = xc - xl

    IF ( gas_flag ) THEN

       IF ( REAL(xc, wp) .GT. EPSILON(1.0_wp) ) THEN
          sp_gas_const_c =                                                   &
               ( (xc - SUM(xg)) * sp_gas_const_a +                           &
               SUM(xg * sp_gas_const_g) ) / xc
       ELSE
          sp_gas_const_c = sp_gas_const_a
       END IF

       inv_rho_c = sp_gas_const_c * T * inv_pres
       rho_c = 1.0_wp / inv_rho_c
       inv_rho_g = sp_gas_const_g * T * inv_pres

    ELSE

       rho_c = rho_l
       inv_rho_c = inv_rho_l
       inv_rho_g = 0.0_wp

    END IF

    inv_rhom = SUM(xs * inv_rho_s) + xc * inv_rho_c
    IF ( gas_flag .AND. liquid_flag ) inv_rhom = inv_rhom + xl * inv_rho_l
    rho_m = 1.0_wp / inv_rhom

    alphas = xs * rho_m * inv_rho_s
    alphag = xg * rho_m * inv_rho_g
    alphal = 0.0_wp
    IF ( gas_flag .AND. liquid_flag ) alphal = xl * rho_m * inv_rho_l
