    !> \file mixture_from_mass.inc
    !> \brief Shared REAL/COMPLEX closure for mass-fraction mixture properties.
    !> \details Included inside the typed wrappers in state_conversion_2d.
    !> \par Inputs
    !> T, xs, xg, xl and the configured component densities/gas constants.
    !> \par Outputs
    !> rho_m, inv_rhom, rho_c, inv_rho_c, alphas, alphag and alphal.
    !> \par Local scratch
    !> xc, sp_gas_const_c and inv_rho_g are declared by the including wrapper.
    !--------------------------------------------------------------------------
    ! Shared executable core for mixture properties derived from mass fractions:
    ! carrier/mixture densities and the corresponding volume fractions.
    ! It is included by both the REAL and COMPLEX wrapper procedures in
    ! state_conversion_2d; declarations and type selection belong there.
    ! Keep this fragment type-independent so both paths use identical formulas.
    !--------------------------------------------------------------------------

    ! Additional gases are part of xc, not a separate dispersed phase.
    ! Ambient air therefore has mass fraction xc-SUM(xg).
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

    ! Mass-fraction mixing is linear in specific volume (1/rho), not in rho.
    inv_rhom = SUM(xs * inv_rho_s) + xc * inv_rho_c
    IF ( gas_flag .AND. liquid_flag ) inv_rhom = inv_rhom + xl * inv_rho_l
    rho_m = 1.0_wp / inv_rhom

    ! Recover volume fractions only after the whole-mixture density is known.
    alphas = xs * rho_m * inv_rho_s
    alphag = xg * rho_m * inv_rho_g
    alphal = 0.0_wp
    IF ( gas_flag .AND. liquid_flag ) alphal = xl * rho_m * inv_rho_l
