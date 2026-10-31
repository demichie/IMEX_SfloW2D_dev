!> Analytic frozen-coefficient Darcy subsystem.
!> Only total mass and excess pore pressure are integrated as independent ODEs.
!> The caller is responsible for common carrier-inventory/packing limiters,
!> and for consistent conservative momentum, energy, and species updates.
MODULE pore_pressure_darcy_analytic_2d
  USE parameters_2d, ONLY : wp
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: solve_frozen_darcy_mass_pressure
CONTAINS
  !> Exact map for dM/dt=-A*p, dp/dt=-lambda*p, Qp=M*p.
  !> A is set to zero when gas_loss_flag is false.
  !> Feasibility checks are delegated to the coordinator: the map must not
  !> be applied to a state with insufficient carrier-gas inventory.
  SUBROUTINE solve_frozen_darcy_mass_pressure(mass_before, pressure_before, &
       interval, dissipation_rate, loss_coefficient, gas_loss_enabled, &
       mass_after, pressure_after, conservative_pressure_after, &
       gas_mass_lost, admissible)
    REAL(wp), INTENT(IN) :: mass_before, pressure_before, interval
    REAL(wp), INTENT(IN) :: dissipation_rate, loss_coefficient
    LOGICAL, INTENT(IN) :: gas_loss_enabled
    REAL(wp), INTENT(OUT) :: mass_after, pressure_after
    REAL(wp), INTENT(OUT) :: conservative_pressure_after, gas_mass_lost
    LOGICAL, INTENT(OUT) :: admissible
    REAL(wp) :: exponential_factor, pressure_integral, coefficient

    mass_after=mass_before
    pressure_after=pressure_before
    conservative_pressure_after=mass_before*pressure_before
    gas_mass_lost=0.0_wp
    admissible=.FALSE.

    IF (.NOT.ieee_is_finite(mass_before) .OR. &
        .NOT.ieee_is_finite(pressure_before) .OR. &
        .NOT.ieee_is_finite(interval) .OR. &
        .NOT.ieee_is_finite(dissipation_rate) .OR. &
        .NOT.ieee_is_finite(loss_coefficient)) RETURN
    IF (mass_before<=0.0_wp .OR. interval<0.0_wp .OR. &
        dissipation_rate<0.0_wp .OR. loss_coefficient<0.0_wp) RETURN

    ! The existing v2/dev diffusion closure does not dissipate
    ! negative excess pore pressure; preserve that guard.
    IF (pressure_before<=0.0_wp .OR. interval==0.0_wp) THEN
      admissible=.TRUE.
      RETURN
    END IF

    exponential_factor=EXP(-dissipation_rate*interval)
    pressure_after=pressure_before*exponential_factor

    coefficient=0.0_wp
    IF (gas_loss_enabled) coefficient=loss_coefficient
    IF (dissipation_rate>0.0_wp) THEN
      ! Stable evaluation of (1-exp(-lambda*dt))/lambda,
      ! including small lambda*dt without cancellation.
      pressure_integral=-EXPM1(-dissipation_rate*interval)/dissipation_rate
    ELSE
      pressure_integral=interval
    END IF

    gas_mass_lost=coefficient*pressure_before*pressure_integral
    mass_after=mass_before-gas_mass_lost
    admissible=(mass_after>0.0_wp .AND. ieee_is_finite(mass_after))
    IF (.NOT.admissible) RETURN
    conservative_pressure_after=mass_after*pressure_after
  END SUBROUTINE solve_frozen_darcy_mass_pressure
END MODULE pore_pressure_darcy_analytic_2d
