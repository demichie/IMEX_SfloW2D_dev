!> \brief Assemble the atomic post-RK mass-exchange source rates.
!> \details Physical exchange and equation-rate assembly only; the STEP
!> executor retains ownership of inventory and bed-state updates.
MODULE mass_exchange_laws_2d
  USE diagnostics_2d, ONLY : fatal_error
  USE parameters_2d
  USE constitutive_parameters_2d
  USE state_conversion_2d
  USE mass_exchange_balance_2d, ONLY : assemble_mass_exchange_balance
  USE mass_exchange_lava_deposition_2d, ONLY: lava_deposition_flag, eval_mass_exchange_lava_deposition_explicit
  USE mass_exchange_rates_2d, ONLY: mass_exchange_rates_type
  USE mass_exchange_carrier_loss_2d, ONLY : eval_mass_exchange_carrier_loss_explicit
  USE mass_exchange_erosion_2d, ONLY : eval_mass_exchange_erosion_explicit
  USE mass_exchange_deposition_2d, ONLY : eval_mass_exchange_deposition_explicit
  USE mass_exchange_entrainment_2d, ONLY : eval_mass_exchange_entrainment_explicit
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_mass_exchange_terms
CONTAINS
  !******************************************************************************
  !> \brief Compute limited erosion, deposition, entrainment and associated equation sources.
  !> \details Evaluates per-class solid erosion and settling/deposition, carrier
  !>          erosion and loss, optional air entrainment, vaporization,
  !>          thermal and pore-pressure sources, and the corresponding
  !>          conservative equation source. Erosion is limited by the available
  !>          erodible inventory over dt; deposition and carrier loss are
  !>          limited by the mass available in the flow and by the maximum
  !>          packing constraint. The bed-elevation rate is
  !>          \f$\dot{B}=(\sum_iD_i-\sum_iE_i)/(1-\phi_e)\f$.
  !>          If the flow is dry or its total solid fraction is at/below the
  !>          configured cutoff, every output is zero.
  !> \date 20/01/2018
  !> @author Mattia de' Michieli Vitturi
  !> \param[in] qpj Local physical state of length n_vars+2 with component mass fractions.
  !> \param[in] B_zone Connected-water-region indicator used by the exchange laws.
  !> \param[in] B_prime_x Filtered bed derivative dB/dx [dimensionless].
  !> \param[in] B_prime_y Filtered bed derivative dB/dy [dimensionless].
  !> \param[in] erodible Available erodible solid inventory expressed as thickness per class [m].
  !> \param[in] dt Time increment [s].
  !> \param[out] erosion_term Per-solid erosion rates after inventory limiting [m s^-1].
  !> \param[out] deposition_term Per-solid deposition rates [m s^-1].
  !> \param[out] continuous_phase_erosion_term Carrier-phase volume rate entrained with erosion
  !>                                           from the substrate [m s^-1].
  !> \param[out] continuous_phase_loss_term Carrier-phase volume loss rate, including configured
  !>                                        loss and any active pore-gas outflow [m s^-1].
  !> \param[out] eqns_term Conservative source vector [equation] for the cell-local mass-exchange
  !>                        update; includes mass, momentum, energy and component sources.
  !> \param[out] topo_term Proposed cell bed-elevation rate, positive for bed growth [m s^-1].
  !> \note B_zone selects the connected-water zone used by vaporization. Substrate temperature
  !>       is selected locally so parallel calls never mutate shared material configuration.
  !******************************************************************************

  SUBROUTINE eval_mass_exchange_terms( qpj , B_zone , B_prime_x , B_prime_y ,   &
       erodible , dt , erosion_term , deposition_term ,                         &
       continuous_phase_erosion_term , continuous_phase_loss_term , eqns_term , &
       topo_term, lava_rate  )

    USE parameters_2d, ONLY : erodible_deposit_flag , liquid_vaporization_flag ,&
         pore_pressure_flag , gas_loss_flag

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(IN) :: B_prime_x
    REAL(wp), INTENT(IN) :: B_prime_y
    REAL(wp), INTENT(IN) :: erodible(n_solid)
    REAL(wp), INTENT(IN) :: dt

    REAL(wp), INTENT(OUT) :: erosion_term(n_solid)
    REAL(wp), INTENT(OUT) :: deposition_term(n_solid)
    REAL(wp), INTENT(OUT) :: continuous_phase_erosion_term
    REAL(wp), INTENT(OUT) :: continuous_phase_loss_term
    REAL(wp), INTENT(OUT) :: eqns_term(n_eqns)
    REAL(wp), INTENT(OUT) :: topo_term, lava_rate
    TYPE(mass_exchange_rates_type) :: lava_rates


    REAL(wp) :: mod_vel
    REAL(wp) :: mod_vel2

    REAL(wp) :: r_h          !< real-value flow thickness
    REAL(wp) :: r_u          !< real-value x-velocity
    REAL(wp) :: r_v          !< real-value y-velocity
    REAL(wp) :: r_W          !< real-value z-velocity
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_alphag(n_add_gas) !< real-value add.gas volume fractions
    REAL(wp) :: r_alphal          !< real-value liquid volume fraction
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_T          !< real-value mixture temperature [K]
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(wp) :: r_Zs         !< real-value stochastic variable
    REAL(wp) :: r_exc_pore_pres   !< real-value excess p pres

    REAL(wp) :: alphas_tot   !< total solid volume fraction

    REAL(wp) :: r_Ri
    REAL(wp) :: r_red_grav

    REAL(wp) :: r_sp_heat_c
    REAL(wp) :: r_sp_heat_mix

    LOGICAL :: solid_processes_enabled, liquid_loss_requested
    REAL(wp) :: air_entr

    REAL(wp) :: T_erodible_local !< cell-local substrate temperature [K]


    erosion_term(1:n_solid) = 0.0_wp
    deposition_term(1:n_solid) = 0.0_wp
    continuous_phase_erosion_term = 0.0_wp
    continuous_phase_loss_term = 0.0_wp
    eqns_term(1:n_eqns) = 0.0_wp
    topo_term = 0.0_wp
    lava_rate = 0.0_wp
    r_Zs = 0.0_wp
    r_exc_pore_pres = 0.0_wp

    IF ( qpj(1) .LE. epsilon(1.0_wp) ) THEN

       RETURN

    END IF

    IF (lava_deposition_flag) THEN
      IF (.NOT.lava_flag) CALL fatal_error('lava deposition requires LAVA_FLAG=T')
      CALL eval_mass_exchange_lava_deposition_explicit(qpj(1),qpj(4),dt,lava_rates)
      lava_rate = lava_rates%continuous_phase_deposition_rate
      eqns_term(1) = lava_rates%mixture_mass_rate
      eqns_term(2) = qpj(idx_u)*MIN(0.0_wp,lava_rates%mixture_mass_rate)
      eqns_term(3) = qpj(idx_v)*MIN(0.0_wp,lava_rates%mixture_mass_rate)
      eqns_term(4) = lava_rates%energy_rate
      topo_term = lava_rates%dB_dt
      RETURN
    END IF

    r_h = qpj(1)
    r_u = qpj(idx_u)
    r_v = qpj(idx_v)

    CALL primitive_to_volume_fractions(qpj, r_alphas, r_alphag, r_alphal)

    alphas_tot = SUM(r_alphas)

    ! The solid-fraction cutoff governs solid exchange, not liquid carrier loss.
    ! Preserve the legacy cutoff for all other processes when liquid loss
    ! is not requested, and never consume substrate below the cutoff.
    solid_processes_enabled = alphas_tot > alphastot_min
    liquid_loss_requested = .FALSE.
    IF (liquid_flag .AND. ALLOCATED(loss_rate)) &
         liquid_loss_requested = loss_rate > 0.0_wp
    IF (.NOT.solid_processes_enabled .AND. .NOT.liquid_loss_requested) RETURN

    IF ( stoch_transport_flag ) r_Zs = qpj(idx_stoch)

    IF ( pore_pressure_flag ) r_exc_pore_pres = qpj(idx_pore)

    IF ( slope_correction_flag ) THEN

       r_w = r_u * B_prime_x + r_v * B_prime_y

    ELSE

       r_w = 0.0_wp

    END IF

    mod_vel2 = r_u**2 + r_v**2 + r_w**2
    mod_vel = SQRT( mod_vel2 )

    IF (solid_processes_enabled) THEN
       CALL eval_mass_exchange_erosion_explicit(r_h, mod_vel, alphas_tot, erodible, dt, &
            erosion_term, continuous_phase_erosion_term)
    END IF

    r_T = qpj(4)

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,r_sp_heat_c,             &
         r_sp_heat_mix)

    IF (solid_processes_enabled) THEN
       CALL eval_mass_exchange_deposition_explicit(r_h, r_T, r_alphas, alphas_tot, &
            r_rho_m, r_rho_c, dt, deposition_term)
    END IF

    CALL eval_mass_exchange_carrier_loss_explicit(r_h, r_T, r_alphas, alphas_tot, &
         r_rho_c, r_exc_pore_pres, deposition_term, dt, &
         continuous_phase_loss_term)


    ! Re-erodible deposits use the current flow temperature by model policy,
    ! but the configured substrate temperature is shared read-only state.
    T_erodible_local = T_erodible
    IF ( erodible_deposit_flag ) T_erodible_local = r_T

    air_entr = 0.0_wp
    IF (solid_processes_enabled) &
         air_entr = eval_mass_exchange_entrainment_explicit(r_h, r_Ri, mod_vel2)

    CALL assemble_mass_exchange_balance(B_zone, r_u, r_v, r_T, r_rho_c, &
         r_sp_heat_c, T_erodible_local, air_entr, r_Zs, r_exc_pore_pres, &
         erosion_term, deposition_term, continuous_phase_erosion_term, &
         continuous_phase_loss_term, eqns_term, topo_term)

    RETURN

  END SUBROUTINE eval_mass_exchange_terms
END MODULE mass_exchange_laws_2d
