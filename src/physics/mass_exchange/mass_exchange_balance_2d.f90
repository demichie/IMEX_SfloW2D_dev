!> Assemble conservative and topographic sources from one atomic mass-exchange batch.
!> Rate limiting remains owned by the process kernels and occurs before this call.
MODULE mass_exchange_balance_2d
  USE parameters_2d
  USE constitutive_parameters_2d
  USE pore_pressure_laws_2d, ONLY : eval_pore_pressure_mass_exchange_source
  USE mass_exchange_vaporization_2d, ONLY : add_mass_exchange_vaporization
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: assemble_mass_exchange_balance
CONTAINS
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

    CALL add_mass_exchange_vaporization(B_zone, r_T, deposition_term, eqns_term)

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
END MODULE mass_exchange_balance_2d
