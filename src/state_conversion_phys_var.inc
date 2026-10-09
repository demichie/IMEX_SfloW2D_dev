    !> \file state_conversion_phys_var.inc
    !> \brief Shared REAL/COMPLEX conservative-to-physical state decoder.
    !> \details Included by r_phys_var and c_phys_var, which own typed declarations.
    !> \par Inputs
    !> qj, component indices, phase flags, ambient state and material parameters.
    !> \par Outputs
    !> h, u, v, T, rho_m, rho_c, their reciprocals, alphas, alphag, alphal,
    !> Zs and exc_pore_pres. Temporary fractions and heat capacities are local.
    !--------------------------------------------------------------------------
    ! Shared executable core that derives the physical state from conservative
    ! variables. It is included by both r_phys_var and c_phys_var; declarations
    ! and type selection belong to those REAL and COMPLEX wrapper procedures.
    ! Keep this fragment type-independent so both paths use identical formulas.
    !--------------------------------------------------------------------------

    alphal = 0.0_wp
    alphag = 0.0_wp
    Zs = 0.0_wp
    exc_pore_pres = 0.0_wp
    xl = 0.0_wp
    xg = 0.0_wp

    ! Branch on the real base state, but retain imaginary perturbations in the
    ! algebra below so the COMPLEX wrapper can supply complex-step derivatives.
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

    ! The dry return above must precede every division by column mass M=qj(1).
    inv_qj1 = 1.0_wp / qj(1)

    xs(1:n_solid) = qj(idx_solid_first:idx_solid_last) * inv_qj1

    IF ( n_add_gas .GT. 0 ) THEN

       xg(1:n_add_gas) = qj(idx_add_gas_first:idx_add_gas_last) * inv_qj1

    END IF

    IF ( stoch_transport_flag ) Zs = qj(idx_stoch) * inv_qj1
    IF ( pore_pressure_flag ) exc_pore_pres = qj(idx_pore) * inv_qj1

    IF ( gas_flag .AND. liquid_flag ) THEN

       xl = qj(n_vars) * inv_qj1

    END IF

    CALL eval_mixture_heat_capacity(                                           &
         ! IN
         xs, xg, xl,                                                           &
         ! OUT
         xc, mixture_sp_heat_c, sp_heat_mix)

    ! qj(4)=M*cp*T contains thermal energy only, with no kinetic-energy term.
    ! The near-dry temperature fallback avoids amplification by negligible mass.
    IF ( REAL(qj(1), wp) .GT. eps_sing ) THEN

       T = qj(4) / ( qj(1) * sp_heat_mix )

       IF ( REAL(T, wp) .LE. 0.0_wp ) T = T_ambient

    ELSE

       T = T_ambient

    END IF

    CALL eval_mixture_properties_from_mass_fractions(                          &
         ! IN
         T, xs, xg, xl,                                                        &
         ! OUT
         rho_m, inv_rhom, rho_c, inv_rho_c, alphas, alphag, alphal)

    h = qj(1) * inv_rhom

    CALL velocity_from_conservative(qj, u, v)
