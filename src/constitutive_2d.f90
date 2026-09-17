!********************************************************************************
!> \brief Compatibility facade for constitutive model modules
!********************************************************************************
MODULE constitutive_2d

  USE constitutive_parameters_2d
  USE state_conversion_2d
  USE equation_terms_2d

  USE equation_metadata_2d, ONLY : equation_partition_type

  USE parameters_2d, ONLY : wp, sp, tolh
  USE parameters_2d, ONLY : n_eqns, n_vars, n_solid, n_add_gas, n_quad,         &
       n_stoch_vars, n_pore_vars
  USE parameters_2d, ONLY : rheology_flag, rheology_model, energy_flag,         &
       liquid_flag, gas_flag, alpha_flag, slope_correction_flag,                &
       curvature_term_flag, stochastic_flag, mean_field_flag,                  &
       stoch_transport_flag, pore_pressure_flag, sutherland_flag

  USE parameters_2d, ONLY : idx_h, idx_hu, idx_hv, idx_T, idx_alfas_first,      &
       idx_alfas_last, idx_addGas_first, idx_addGas_last, idx_stoch, idx_pore,  &
       idx_u, idx_v

  USE parameters_2d, ONLY : idx_totMassEqn, idx_uEqn, idx_vEqn, idx_engyEqn,    &
       idx_solidEqn_first, idx_solidEqn_last, idx_addGasEqn_first,              &
       idx_addGasEqn_last, idx_stochEqn, idx_poreEqn

  IMPLICIT NONE

END MODULE constitutive_2d
