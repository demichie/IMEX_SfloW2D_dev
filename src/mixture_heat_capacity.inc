    !> \file mixture_heat_capacity.inc
    !> \brief Shared REAL/COMPLEX mass-weighted heat-capacity closure.
    !> \details Included inside the typed wrappers in state_conversion_2d.
    !> \par Inputs
    !> xs, xg, xl, phase flags and configured component heat capacities.
    !> \par Outputs
    !> xc, sp_heat_c_mix and sp_heat_mix [J kg^-1 K^-1 for heat capacities].
    !> \par Local scratch
    !> carrier_sp_heat_mass is declared by the including wrapper.
    !--------------------------------------------------------------------------
    ! Shared executable core for carrier and whole-mixture heat capacities.
    ! It is included by both the REAL and COMPLEX wrapper procedures in
    ! state_conversion_2d; all declarations and type selection belong there.
    ! Keep this fragment type-independent so both paths use identical formulas.
    !--------------------------------------------------------------------------

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

    ! carrier_sp_heat_mass already includes its carrier mass fraction; do not
    ! multiply it by xc again when forming the mixture heat capacity.
    sp_heat_mix = SUM(xs * sp_heat_s) + carrier_sp_heat_mass
    IF ( gas_flag .AND. liquid_flag ) sp_heat_mix = sp_heat_mix + xl * sp_heat_l
