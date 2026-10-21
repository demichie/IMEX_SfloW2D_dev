!********************************************************************************
!> \brief General thermal source terms for the depth-averaged equations.
!********************************************************************************

MODULE thermal_sources_2d

  USE parameters_2d, ONLY : wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_inlet_thermal_energy_source
  PUBLIC :: eval_mass_exchange_thermal_energy
  PUBLIC :: eval_vaporization_sources

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

  PURE FUNCTION eval_mass_exchange_thermal_energy(flow_temperature,            &
       solid_density, solid_heat_capacity, deposition_rate, erosion_rate,      &
       carrier_density, carrier_heat_capacity, carrier_loss_rate,              &
       substrate_temperature, substrate_carrier_density,                       &
       carrier_erosion_rate, ambient_temperature, air_heat_capacity,            &
       ambient_air_density, air_entrainment_rate) RESULT(source)

    REAL(wp), INTENT(IN) :: flow_temperature
    REAL(wp), INTENT(IN) :: solid_density(:)
    REAL(wp), INTENT(IN) :: solid_heat_capacity(:)
    REAL(wp), INTENT(IN) :: deposition_rate(:)
    REAL(wp), INTENT(IN) :: erosion_rate(:)
    REAL(wp), INTENT(IN) :: carrier_density
    REAL(wp), INTENT(IN) :: carrier_heat_capacity
    REAL(wp), INTENT(IN) :: carrier_loss_rate
    REAL(wp), INTENT(IN) :: substrate_temperature
    REAL(wp), INTENT(IN) :: substrate_carrier_density
    REAL(wp), INTENT(IN) :: carrier_erosion_rate
    REAL(wp), INTENT(IN) :: ambient_temperature
    REAL(wp), INTENT(IN) :: air_heat_capacity
    REAL(wp), INTENT(IN) :: ambient_air_density
    REAL(wp), INTENT(IN) :: air_entrainment_rate
    REAL(wp) :: source

    source = -flow_temperature * ( SUM( solid_density * solid_heat_capacity *  &
         deposition_rate ) + carrier_density * carrier_heat_capacity *         &
         carrier_loss_rate ) + substrate_temperature * ( SUM( solid_density *  &
         solid_heat_capacity * erosion_rate ) + substrate_carrier_density *    &
         carrier_heat_capacity * carrier_erosion_rate ) + ambient_temperature &
         * air_heat_capacity * ambient_air_density * air_entrainment_rate

  END FUNCTION eval_mass_exchange_thermal_energy

  PURE SUBROUTINE eval_vaporization_sources(flow_temperature, solid_density,    &
       solid_heat_capacity, deposition_rate, steam_fraction, gas_heat_capacity,&
       vapor_mass_rate, vapor_energy_source)

    REAL(wp), INTENT(IN) :: flow_temperature
    REAL(wp), INTENT(IN) :: solid_density(:)
    REAL(wp), INTENT(IN) :: solid_heat_capacity(:)
    REAL(wp), INTENT(IN) :: deposition_rate(:)
    REAL(wp), INTENT(IN) :: steam_fraction
    REAL(wp), INTENT(IN) :: gas_heat_capacity
    REAL(wp), INTENT(OUT) :: vapor_mass_rate
    REAL(wp), INTENT(OUT) :: vapor_energy_source

    REAL(wp), PARAMETER :: liquid_temperature = 290.0_wp
    REAL(wp), PARAMETER :: boiling_temperature = 373.15_wp
    REAL(wp), PARAMETER :: latent_heat = 2264705.0_wp
    REAL(wp), PARAMETER :: liquid_heat_capacity = 4184.0_wp

    vapor_mass_rate = steam_fraction * MAX(0.0_wp, flow_temperature -         &
         boiling_temperature) * SUM(solid_density * solid_heat_capacity *      &
         deposition_rate) / (liquid_heat_capacity * (boiling_temperature -    &
         liquid_temperature) + latent_heat)
    vapor_energy_source = vapor_mass_rate * gas_heat_capacity *                &
         boiling_temperature

  END SUBROUTINE eval_vaporization_sources

END MODULE thermal_sources_2d
