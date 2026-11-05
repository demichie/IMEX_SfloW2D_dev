    xc = 1.0_wp - SUM(xs)

    IF ( gas_flag ) THEN
        carrier_sp_heat_mass = (xc - SUM(xg)) * sp_heat_a + SUM(xg * sp_heat_g)
        IF (REAL(xc,wp) > EPSILON(1.0_wp)) THEN
            cp_c = carrier_sp_heat_mass / xc
        ELSE
            cp_c = sp_heat_a
        END IF
    ELSE
        carrier_sp_heat_mass = xc * sp_heat_l
        cp_c = sp_heat_l
    END IF

    cp_m = SUM(xs * sp_heat_s) + carrier_sp_heat_mass
