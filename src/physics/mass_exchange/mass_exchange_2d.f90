!> \brief Assemble the atomic post-RK mass-exchange source rates.
!> \details Physical exchange and equation-rate assembly only; the STEP
!> executor retains ownership of inventory and bed-state updates.
MODULE mass_exchange_2d
  USE diagnostics_2d, ONLY : fatal_error
  USE parameters_2d
  USE constitutive_parameters_2d
  USE state_conversion_2d
  USE mass_exchange_lava_deposition_2d, ONLY: lava_deposition_flag, eval_mass_exchange_lava_deposition_explicit
  USE mass_exchange_rates_2d, ONLY: mass_exchange_rates_type
  USE mass_exchange_carrier_loss_2d, ONLY : eval_mass_exchange_carrier_loss_explicit
  USE mass_exchange_erosion_2d, ONLY : eval_mass_exchange_erosion_explicit
  USE mass_exchange_deposition_2d, ONLY : eval_mass_exchange_deposition_explicit
  USE mass_exchange_entrainment_2d, ONLY : eval_mass_exchange_entrainment_explicit
  USE mass_exchange_vaporization_2d, ONLY: eval_mass_exchange_vaporization_rates_explicit
  USE pore_pressure_laws_2d, ONLY: eval_pore_pressure_mass_exchange_source
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
  SUBROUTINE assemble_mass_exchange_balance(B_zone, r_u, r_v, r_T, r_rho_c, &
       r_sp_heat_c, T_erodible_local, air_entr, r_Zs, r_exc_pore_pres, &
       erosion_term, deposition_term, continuous_phase_erosion_term, &
       continuous_phase_loss_term, eqns_term, topo_term)
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(IN) :: r_u, r_v, r_T, r_rho_c, r_sp_heat_c
    REAL(wp), INTENT(IN) :: T_erodible_local, air_entr, r_Zs, r_exc_pore_pres
    REAL(wp), INTENT(IN) :: erosion_term(n_solid), deposition_term(n_solid)
    REAL(wp), INTENT(IN) :: continuous_phase_erosion_term, continuous_phase_loss_term
    REAL(wp), INTENT(OUT) :: eqns_term(n_eqns), topo_term
    REAL(wp) :: dep_tot, ers_tot, rho_dep_tot, rho_ers_tot
    REAL(wp) :: vapor_mass_rate, vapor_component_rate, vapor_energy_rate

    eqns_term(1:n_eqns) = 0.0_wp

    ! solid total volume deposition rate [m s-1]
    dep_tot = SUM( deposition_term )
    ! solid total volume deposition rate [m s-1]
    ers_tot = SUM( erosion_term )

    ! solid total mass deposition rate [kg m-2 s-1]
    rho_dep_tot = DOT_PRODUCT( rho_s , deposition_term )
    ! solid total mass erosion rate [kg m-2 s-1]
    rho_ers_tot = DOT_PRODUCT( rho_s , erosion_term )

    ! total mass equation source term [kg m-2 s-1]:
    ! deposition, erosion and entrainment are considered
    eqns_term(1) = rho_a_amb * air_entr + rho_ers_tot - rho_dep_tot             &
         + rho_c_sub * continuous_phase_erosion_term                            &
         - r_rho_c * continuous_phase_loss_term

    ! x-momenutm equation source term [kg m-1 s-2]:
    ! only deposition contribute to change in momentum, erosion does not carry
    ! any momentum inside the flow
    eqns_term(2) = - r_u * ( rho_dep_tot + r_rho_c * continuous_phase_loss_term )

    ! y-momentum equation source term [kg m-1 s-2]:
    ! only deposition contribute to change in momentum, erosion does not carry
    ! any momentum inside the flow
    eqns_term(3) = - r_v * ( rho_dep_tot + r_rho_c * continuous_phase_loss_term )

    ! Thermal-energy equation source term [kg s-3]:
    ! deposition, erosion and entrainment are considered
    eqns_term(4) = -r_T * ( SUM(rho_s * sp_heat_s * deposition_term) +       &
         r_rho_c * r_sp_heat_c * continuous_phase_loss_term ) +               &
         T_erodible_local * ( SUM(rho_s * sp_heat_s * erosion_term) +         &
         rho_c_sub * r_sp_heat_c * continuous_phase_erosion_term ) +           &
         T_ambient * sp_heat_a * rho_a_amb * air_entr

    ! solid phase mass equation source term [kg m-2 s-1]:
    ! due to solid erosion and deposition
    eqns_term(idx_solidEqn_first:idx_solidEqn_last) = rho_s(1:n_solid) *        &
         ( erosion_term(1:n_solid) - deposition_term(1:n_solid) )

    CALL eval_mass_exchange_vaporization_rates_explicit(B_zone, r_T, deposition_term, &
         vapor_mass_rate, vapor_component_rate, vapor_energy_rate)
    eqns_term(1) = eqns_term(1) + vapor_mass_rate
    eqns_term(4) = eqns_term(4) + vapor_energy_rate
    IF (gas_flag .AND. liquid_vaporization_flag .AND. B_zone /= 0) &
         eqns_term(5+n_solid) = eqns_term(5+n_solid) + vapor_component_rate

    IF ( stoch_transport_flag) THEN

       ! Evaluate the equation term related to the noise transport equation
       ! (if transport flag is false nothing is done)
       eqns_term(idx_stochEqn) = eqns_term(1) * r_Zs

    END IF

    IF ( pore_pressure_flag ) THEN

       ! Equation for q1*exc_pore_press
       eqns_term(idx_poreEqn) = eval_pore_pressure_mass_exchange_source(       &
            eqns_term(1), r_exc_pore_pres)

    END IF

    ! erodible layer thickness source terms [m s-1]:
    ! due to erosion and deposition of solid+continuous phase
    topo_term = ( dep_tot - ers_tot ) / ( 1.0_wp - erodible_porosity )

  END SUBROUTINE assemble_mass_exchange_balance
END MODULE mass_exchange_2d
