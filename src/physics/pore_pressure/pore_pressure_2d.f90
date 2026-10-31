!> Coordinator for the optional excess pore-pressure source family.
!> Physical Darcy closures are implemented by pore_pressure_darcy_2d.
MODULE pore_pressure_2d
 USE pore_pressure_darcy_2d, ONLY: &
  solve_frozen_darcy_mass_pressure, &
  eval_pore_pressure_dissipation_explicit, &
  eval_pore_pressure_dissipation_semi_implicit, &
  eval_pore_pressure_dissipation_implicit, &
  solve_pore_pressure_dissipation_analytic, &
  advance_pore_pressure_diffusion, &
  eval_pore_pressure_gas_loss, &
  eval_pore_pressure_inlet_source, &
  eval_pore_pressure_mass_exchange_source
 IMPLICIT NONE
 PRIVATE
 PUBLIC :: solve_frozen_darcy_mass_pressure
 PUBLIC :: eval_pore_pressure_dissipation_explicit
 PUBLIC :: eval_pore_pressure_dissipation_semi_implicit
 PUBLIC :: eval_pore_pressure_dissipation_implicit
 PUBLIC :: solve_pore_pressure_dissipation_analytic
 PUBLIC :: advance_pore_pressure_diffusion
 PUBLIC :: eval_pore_pressure_gas_loss
 PUBLIC :: eval_pore_pressure_inlet_source
 PUBLIC :: eval_pore_pressure_mass_exchange_source
END MODULE pore_pressure_2d
