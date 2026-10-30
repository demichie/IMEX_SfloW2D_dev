!> Thermal source coordinator and compatibility functions for dependent physics.
MODULE thermal_laws_2d
  USE parameters_2d, ONLY : wp
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
  PUBLIC :: eval_inlet_thermal_energy_source, eval_rheology_cooling_energy_rate
  PUBLIC :: eval_vaporization_mass_rate, eval_vaporization_energy_rate
  REAL(wp), PARAMETER :: vaporization_liquid_temperature = 290.0_wp
  REAL(wp), PARAMETER :: vaporization_boiling_temperature = 373.15_wp
  REAL(wp), PARAMETER :: vaporization_latent_heat = 2264705.0_wp
  REAL(wp), PARAMETER :: vaporization_liquid_heat_capacity = 4184.0_wp
CONTAINS
PURE FUNCTION eval_inlet_thermal_energy_source(time_coefficient, injection_rate, &
       mixture_density, mixture_heat_capacity, source_temperature) RESULT(source)

    REAL(wp), INTENT(IN) :: time_coefficient
    REAL(wp), INTENT(IN) :: injection_rate
    REAL(wp), INTENT(IN) :: mixture_density
    REAL(wp), INTENT(IN) :: mixture_heat_capacity
    REAL(wp), INTENT(IN) :: source_temperature
    REAL(wp) :: source

    source = time_coefficient * injection_rate * mixture_density *             &
         mixture_heat_capacity * source_temperature

  END FUNCTION eval_inlet_thermal_energy_source

PURE FUNCTION eval_rheology_cooling_energy_rate(model_name, temperature,     &
       radiative_coefficient, convective_coefficient, ambient_temperature)   &
       RESULT(rate)

    CHARACTER(*), INTENT(IN) :: model_name
    COMPLEX(wp), INTENT(IN) :: temperature
    REAL(wp), INTENT(IN) :: radiative_coefficient, convective_coefficient
    REAL(wp), INTENT(IN) :: ambient_temperature
    COMPLEX(wp) :: rate

    SELECT CASE (TRIM(model_name))
    CASE ('lava','rheology_lava')
       rate = eval_thermal_cooling_implicit(temperature, radiative_coefficient, &
            convective_coefficient, ambient_temperature)
    CASE DEFAULT
       rate = CMPLX(0.0_wp, 0.0_wp, wp)
    END SELECT

  END FUNCTION eval_rheology_cooling_energy_rate

PURE FUNCTION eval_vaporization_mass_rate(flow_temperature, solid_density,   &
       solid_heat_capacity, deposition_rate, steam_fraction) RESULT(rate)

    REAL(wp), INTENT(IN) :: flow_temperature
    REAL(wp), INTENT(IN) :: solid_density(:)
    REAL(wp), INTENT(IN) :: solid_heat_capacity(:)
    REAL(wp), INTENT(IN) :: deposition_rate(:)
    REAL(wp), INTENT(IN) :: steam_fraction
    REAL(wp) :: rate

    rate = steam_fraction * MAX(0.0_wp, flow_temperature -                    &
         vaporization_boiling_temperature) *                                 &
         SUM(solid_density * solid_heat_capacity * deposition_rate) /         &
         (vaporization_liquid_heat_capacity *                                 &
         (vaporization_boiling_temperature -                                  &
         vaporization_liquid_temperature) + vaporization_latent_heat)

  END FUNCTION eval_vaporization_mass_rate

PURE FUNCTION eval_vaporization_energy_rate(vapor_mass_rate,               &
       gas_heat_capacity) RESULT(rate)

    REAL(wp), INTENT(IN) :: vapor_mass_rate, gas_heat_capacity
    REAL(wp) :: rate

    rate = vapor_mass_rate * gas_heat_capacity *                             &
         vaporization_boiling_temperature

  END FUNCTION eval_vaporization_energy_rate

END MODULE thermal_laws_2d
