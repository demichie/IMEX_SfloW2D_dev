!> Thermal source coordinator and compatibility functions for dependent physics.
MODULE thermal_laws_2d
  USE thermal_cooling_2d, ONLY : eval_thermal_cooling_explicit, &
       eval_thermal_cooling_semi_implicit, eval_thermal_cooling_implicit, &
       solve_thermal_cooling_analytic, advance_lava_cooling_temperature, &
       eval_cooling_state_explicit, eval_cooling_state_semi_implicit, &
       eval_cooling_state_implicit, solve_cooling_state_analytic
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_thermal_cooling_explicit, eval_thermal_cooling_semi_implicit
  PUBLIC :: eval_thermal_cooling_implicit, solve_thermal_cooling_analytic
  PUBLIC :: advance_lava_cooling_temperature
  PUBLIC :: eval_cooling_state_explicit, eval_cooling_state_semi_implicit
  PUBLIC :: eval_cooling_state_implicit, solve_cooling_state_analytic
END MODULE thermal_laws_2d
