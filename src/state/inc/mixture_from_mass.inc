!> @file mixture_from_mass.inc
!> @brief Compute mixture density and phase volume fractions from mass fractions.
    ! Single thermodynamic density kernel shared with evaluate_mixture_density.
    ! Keep the REAL/COMPLEX arithmetic and carrier fallback in one place.
    INCLUDE 'inc/density_from_mass.inc'

    ! Only phase-specific properties remain here. inv_rho_c, inv_rho_m and
    ! rho_m have already been computed by the shared density kernel.
    IF (gas_flag) THEN
        rho_c = 1.0_wp / inv_rho_c
        inv_rho_g = evaluate_gas_specific_volume(T,sp_gas_const_g)
    ELSE
        rho_c = rho_l
        inv_rho_g = 0.0_wp
    END IF

    alphas = xs * rho_m * inv_rho_s
    alphag = xg * rho_m * inv_rho_g
