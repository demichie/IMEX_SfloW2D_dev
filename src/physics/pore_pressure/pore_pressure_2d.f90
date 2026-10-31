!> Dispatch pore-pressure analytic models through the conservative rate interface.
!> The model list contains canonical names, allowing future physical models
!> to be added without changing any caller. Darcy is the only model today.
MODULE pore_pressure_2d
  USE parameters_2d, ONLY : wp, n_solid, n_eqns, pore_pressure_flag
  USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
  USE diagnostics_2d, ONLY : fatal_error
  USE pore_pressure_darcy_2d, ONLY : &
       darcy_analytic => solve_pore_pressure_darcy_analytic
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: solve_pore_pressure_analytic
CONTAINS

  !> Complete conservative average rate for the selected analytic operator.
  !> For sequential operators, reconstruct the physical state before
  !> extending this dispatch; do not reuse stale primitive coefficients.
  SUBROUTINE solve_pore_pressure_analytic(models,state,q_before,interval,rate,admissible, &
       include_gas_loss)
    CHARACTER(*), INTENT(IN) :: models(:)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: q_before(n_eqns),interval
    REAL(wp), INTENT(OUT) :: rate(n_eqns)
    LOGICAL, INTENT(OUT) :: admissible
    LOGICAL, INTENT(IN), OPTIONAL :: include_gas_loss
    INTEGER :: i
    rate=0.0_wp
    admissible=.TRUE.
    IF (.NOT.pore_pressure_flag) RETURN
    IF (SIZE(models)>1) THEN
      CALL fatal_error('Coupled pore-pressure analytic models require state reconstruction')
    END IF
    DO i=1,SIZE(models)
      SELECT CASE(TRIM(models(i)))
      CASE('pore_pressure_darcy')
        CALL darcy_analytic(state,q_before,interval,rate,admissible,include_gas_loss)
      CASE DEFAULT
        CALL fatal_error('Unknown analytic pore-pressure model: '//TRIM(models(i)))
      END SELECT
      IF (.NOT.admissible) RETURN
    END DO
  END SUBROUTINE solve_pore_pressure_analytic

END MODULE pore_pressure_2d
