    xc = 1.0_wp - SUM(xs)
    IF ( gas_flag .AND. liquid_flag ) xc = xc - xl

    IF ( gas_flag ) THEN

       carrier_sp_heat_mass =                                                &
            ( xc - SUM(xg) ) * sp_heat_a + SUM(xg * sp_heat_g)

       IF ( REAL(xc, wp) .GT. EPSILON(1.0_wp) ) THEN
          sp_heat_c_mix = carrier_sp_heat_mass / xc
       ELSE
          sp_heat_c_mix = sp_heat_a
       END IF

    ELSE

       carrier_sp_heat_mass = xc * sp_heat_l
       sp_heat_c_mix = sp_heat_l

    END IF

    sp_heat_mix = SUM(xs * sp_heat_s) + carrier_sp_heat_mass
    IF ( gas_flag .AND. liquid_flag ) sp_heat_mix = sp_heat_mix + xl * sp_heat_l
