    xc = 1.0_wp - SUM(xs)
    IF (gas_flag) THEN
        IF (REAL(xc,wp)>EPSILON(1.0_wp)) THEN
            gas_const_c=((xc-SUM(xg))*sp_gas_const_a+SUM(xg*sp_gas_const_g))/xc
        ELSE
            gas_const_c=sp_gas_const_a
        END IF
        inv_rho_c=gas_const_c*T*inv_pres
    ELSE
        inv_rho_c=inv_rho_l
    END IF
    inv_rho_m=SUM(xs*inv_rho_s)+xc*inv_rho_c
    rho_m=1.0_wp/inv_rho_m
