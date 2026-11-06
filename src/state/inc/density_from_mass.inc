!> @file density_from_mass.inc
!> @brief Evaluate carrier and mixture specific volumes from mass fractions and temperature.
    xc = 1.0_wp - SUM(xs)
    IF (gas_flag) THEN
        IF (REAL(xc,wp)>EPSILON(1.0_wp)) THEN
            gas_const_c=((xc-SUM(xg))*sp_gas_const_a+SUM(xg*sp_gas_const_g))/xc
        ELSE
            gas_const_c=sp_gas_const_a
        END IF
        inv_rho_c=evaluate_gas_specific_volume(T,gas_const_c)
    ELSE
        inv_rho_c=inv_rho_l
    END IF
    inv_rho_m=SUM(xs*inv_rho_s)+xc*inv_rho_c
    rho_m=1.0_wp/inv_rho_m
