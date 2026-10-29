MODULE mass_exchange_rates_2d
  USE parameters_2d, ONLY: wp, n_solid, n_add_gas
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: mass_exchange_rates_type, reset_mass_exchange_rates
  TYPE mass_exchange_rates_type
    REAL(wp) :: mixture_mass_rate=0.0_wp, liquid_mass_rate=0.0_wp
    REAL(wp) :: energy_rate=0.0_wp, dB_dt=0.0_wp
    REAL(wp) :: continuous_phase_erosion_rate=0.0_wp
    REAL(wp) :: continuous_phase_deposition_rate=0.0_wp
    REAL(wp), ALLOCATABLE :: solid_mass_rate(:), add_gas_mass_rate(:)
    REAL(wp), ALLOCATABLE :: solid_erosion_rate(:), solid_deposition_rate(:)
  END TYPE
CONTAINS
  SUBROUTINE reset_mass_exchange_rates(r)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    IF (ALLOCATED(r%solid_mass_rate)) THEN
      IF (SIZE(r%solid_mass_rate)/=n_solid) DEALLOCATE(r%solid_mass_rate)
    END IF
    IF (ALLOCATED(r%add_gas_mass_rate)) THEN
      IF (SIZE(r%add_gas_mass_rate)/=n_add_gas) DEALLOCATE(r%add_gas_mass_rate)
    END IF
    IF (ALLOCATED(r%solid_erosion_rate)) THEN
      IF (SIZE(r%solid_erosion_rate)/=n_solid) DEALLOCATE(r%solid_erosion_rate)
    END IF
    IF (ALLOCATED(r%solid_deposition_rate)) THEN
      IF (SIZE(r%solid_deposition_rate)/=n_solid) DEALLOCATE(r%solid_deposition_rate)
    END IF
    IF (.NOT.ALLOCATED(r%solid_mass_rate)) ALLOCATE(r%solid_mass_rate(n_solid))
    IF (.NOT.ALLOCATED(r%add_gas_mass_rate)) ALLOCATE(r%add_gas_mass_rate(n_add_gas))
    IF (.NOT.ALLOCATED(r%solid_erosion_rate)) ALLOCATE(r%solid_erosion_rate(n_solid))
    IF (.NOT.ALLOCATED(r%solid_deposition_rate)) ALLOCATE(r%solid_deposition_rate(n_solid))
    r%mixture_mass_rate=0.0_wp
    r%liquid_mass_rate=0.0_wp
    r%energy_rate=0.0_wp
    r%dB_dt=0.0_wp
    r%continuous_phase_erosion_rate=0.0_wp
    r%continuous_phase_deposition_rate=0.0_wp
    r%solid_mass_rate=0.0_wp
    r%add_gas_mass_rate=0.0_wp
    r%solid_erosion_rate=0.0_wp
    r%solid_deposition_rate=0.0_wp
  END SUBROUTINE
END MODULE mass_exchange_rates_2d
