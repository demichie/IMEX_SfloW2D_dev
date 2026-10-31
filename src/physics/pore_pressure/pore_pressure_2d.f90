!> Dispatch pore-pressure models through the common four-treatment interface.
!> The model list contains canonical names, allowing future physical models
!> to be added without changing any caller. Darcy is the only model today.
MODULE pore_pressure_2d
  USE parameters_2d, ONLY : wp, n_solid, n_eqns, pore_pressure_flag
  USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
  USE diagnostics_2d, ONLY : fatal_error
  USE pore_pressure_darcy_2d, ONLY : &
       darcy_explicit => eval_pore_pressure_darcy_explicit, &
       darcy_semi_implicit => eval_pore_pressure_darcy_semi_implicit, &
       darcy_implicit => eval_pore_pressure_darcy_implicit, &
       darcy_analytic => solve_pore_pressure_darcy_analytic, &
       solve_frozen_darcy_mass_pressure, apply_frozen_darcy_withdrawal, &
       eval_pore_pressure_gas_loss, &
       eval_pore_pressure_inlet_source, eval_pore_pressure_mass_exchange_source
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_pore_pressure_explicit, eval_pore_pressure_semi_implicit
  PUBLIC :: eval_pore_pressure_implicit, solve_pore_pressure_analytic
  ! Compatibility entry points while the full conservative source contract
  ! and the shared carrier-inventory limiters are migrated.
  PUBLIC :: solve_frozen_darcy_mass_pressure, apply_frozen_darcy_withdrawal
  PUBLIC :: eval_pore_pressure_gas_loss
  PUBLIC :: eval_pore_pressure_inlet_source, eval_pore_pressure_mass_exchange_source
CONTAINS

  SUBROUTINE eval_pore_pressure_explicit(models,state,rate,include_gas_loss)
    CHARACTER(*), INTENT(IN) :: models(:)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(OUT) :: rate(n_eqns)
    LOGICAL, INTENT(IN), OPTIONAL :: include_gas_loss
    REAL(wp) :: model_rate(n_eqns)
    INTEGER :: i
    rate=0.0_wp
    IF (.NOT.pore_pressure_flag) RETURN
    DO i=1,SIZE(models)
      SELECT CASE(TRIM(models(i)))
      CASE('pore_pressure_darcy')
        CALL darcy_explicit(state,model_rate,include_gas_loss)
      CASE DEFAULT
        CALL fatal_error('Unknown explicit pore-pressure model: '//TRIM(models(i)))
      END SELECT
      rate=rate+model_rate
    END DO
  END SUBROUTINE eval_pore_pressure_explicit

  SUBROUTINE eval_pore_pressure_semi_implicit(models,state,rate)
    CHARACTER(*), INTENT(IN) :: models(:)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(OUT) :: rate(n_eqns)
    REAL(wp) :: model_rate(n_eqns)
    INTEGER :: i
    rate=0.0_wp
    IF (.NOT.pore_pressure_flag) RETURN
    DO i=1,SIZE(models)
      SELECT CASE(TRIM(models(i)))
      CASE('pore_pressure_darcy')
        CALL darcy_semi_implicit(state,model_rate)
      CASE DEFAULT
        CALL fatal_error('Unknown semi-implicit pore-pressure model: '//TRIM(models(i)))
      END SELECT
      rate=rate+model_rate
    END DO
  END SUBROUTINE eval_pore_pressure_semi_implicit

  SUBROUTINE eval_pore_pressure_implicit(models,state,rate)
    CHARACTER(*), INTENT(IN) :: models(:)
    TYPE(complex_physical_state_type), INTENT(IN) :: state
    COMPLEX(wp), INTENT(OUT) :: rate(n_eqns)
    COMPLEX(wp) :: model_rate(n_eqns)
    INTEGER :: i
    rate=CMPLX(0.0_wp,0.0_wp,wp)
    IF (.NOT.pore_pressure_flag) RETURN
    DO i=1,SIZE(models)
      SELECT CASE(TRIM(models(i)))
      CASE('pore_pressure_darcy')
        CALL darcy_implicit(state,model_rate)
      CASE DEFAULT
        CALL fatal_error('Unknown implicit pore-pressure model: '//TRIM(models(i)))
      END SELECT
      rate=rate+model_rate
    END DO
  END SUBROUTINE eval_pore_pressure_implicit

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

  ! Transitional adapters: existing tests and source assembly still call
  ! the old scalar-only interface until the full process rates are wired.








END MODULE pore_pressure_2d
