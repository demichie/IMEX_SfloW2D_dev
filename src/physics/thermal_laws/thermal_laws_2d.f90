!> Dispatch thermal source models through a common four-treatment interface.
!> Add future models in the SELECT CASE statements below.
MODULE thermal_laws_2d
  USE parameters_2d, ONLY : wp
  USE diagnostics_2d, ONLY : fatal_error
  USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
  USE thermal_cooling_2d, ONLY : eval_thermal_cooling_explicit, &
       eval_thermal_cooling_semi_implicit, eval_thermal_cooling_implicit, &
       solve_thermal_cooling_analytic, advance_lava_cooling_temperature, &
       eval_cooling_state_explicit, eval_cooling_state_semi_implicit, &
       eval_cooling_state_implicit, solve_cooling_state_analytic
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_thermal_explicit, eval_thermal_semi_implicit
  PUBLIC :: eval_thermal_explicit_processes, eval_thermal_semi_implicit_processes
  PUBLIC :: eval_thermal_implicit_processes, solve_thermal_analytic_processes
  PUBLIC :: eval_thermal_implicit, solve_thermal_analytic
  ! Retain scalar entry points for existing callers and regression tests.
  PUBLIC :: eval_thermal_cooling_explicit, eval_thermal_cooling_semi_implicit
  PUBLIC :: eval_thermal_cooling_implicit, solve_thermal_cooling_analytic
  PUBLIC :: advance_lava_cooling_temperature
CONTAINS
  SUBROUTINE eval_thermal_explicit(model_name,state,radiative,convective,ambient,rate)
    CHARACTER(*), INTENT(IN) :: model_name
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: radiative,convective,ambient
    REAL(wp), INTENT(OUT) :: rate
    SELECT CASE (TRIM(model_name))
    CASE ('thermal_cooling')
      CALL eval_cooling_state_explicit(state,radiative,convective,ambient,rate)
    CASE DEFAULT
      CALL fatal_error('Unknown explicit thermal model: '//TRIM(model_name))
    END SELECT
  END SUBROUTINE eval_thermal_explicit

  SUBROUTINE eval_thermal_semi_implicit(model_name,state,rate)
    CHARACTER(*), INTENT(IN) :: model_name
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(OUT) :: rate
    SELECT CASE (TRIM(model_name))
    CASE ('thermal_cooling')
      CALL eval_cooling_state_semi_implicit(state,rate)
    CASE DEFAULT
      CALL fatal_error('Unknown semi-implicit thermal model: '//TRIM(model_name))
    END SELECT
  END SUBROUTINE eval_thermal_semi_implicit

  SUBROUTINE eval_thermal_implicit(model_name,state,radiative,convective,ambient,rate)
    CHARACTER(*), INTENT(IN) :: model_name
    TYPE(complex_physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: radiative,convective,ambient
    COMPLEX(wp), INTENT(OUT) :: rate
    SELECT CASE (TRIM(model_name))
    CASE ('thermal_cooling')
      CALL eval_cooling_state_implicit(state,radiative,convective,ambient,rate)
    CASE DEFAULT
      CALL fatal_error('Unknown implicit thermal model: '//TRIM(model_name))
    END SELECT
  END SUBROUTINE eval_thermal_implicit

  SUBROUTINE solve_thermal_analytic(model_name,state,interval,heat_capacity, &
       radiative,convective,ambient,temperature_after)
    CHARACTER(*), INTENT(IN) :: model_name
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: interval,heat_capacity,radiative,convective,ambient
    REAL(wp), INTENT(OUT) :: temperature_after
    SELECT CASE (TRIM(model_name))
    CASE ('thermal_cooling')
      CALL solve_cooling_state_analytic(state,interval,heat_capacity, &
           radiative,convective,ambient,temperature_after)
    CASE DEFAULT
      CALL fatal_error('Unknown analytic thermal model: '//TRIM(model_name))
    END SELECT
  END SUBROUTINE solve_thermal_analytic

  ! All instantaneous thermal contributions are evaluated at one unchanged state.
  SUBROUTINE eval_thermal_explicit_processes(models,state,radiative,convective,ambient,rate)
    CHARACTER(*), INTENT(IN) :: models(:)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: radiative,convective,ambient
    REAL(wp), INTENT(OUT) :: rate
    REAL(wp) :: process_rate
    INTEGER :: i
    rate=0.0_wp
    DO i=1,SIZE(models)
      CALL eval_thermal_explicit(TRIM(models(i)),state,radiative,convective,ambient,process_rate)
      rate=rate+process_rate
    END DO
  END SUBROUTINE eval_thermal_explicit_processes

  SUBROUTINE eval_thermal_semi_implicit_processes(models,state,rate)
    CHARACTER(*), INTENT(IN) :: models(:)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(OUT) :: rate
    REAL(wp) :: process_rate
    INTEGER :: i
    rate=0.0_wp
    DO i=1,SIZE(models)
      CALL eval_thermal_semi_implicit(TRIM(models(i)),state,process_rate)
      rate=rate+process_rate
    END DO
  END SUBROUTINE eval_thermal_semi_implicit_processes

  SUBROUTINE eval_thermal_implicit_processes(models,state,radiative,convective,ambient,rate)
    CHARACTER(*), INTENT(IN) :: models(:)
    TYPE(complex_physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: radiative,convective,ambient
    COMPLEX(wp), INTENT(OUT) :: rate
    COMPLEX(wp) :: process_rate
    INTEGER :: i
    rate=CMPLX(0.0_wp,0.0_wp,wp)
    DO i=1,SIZE(models)
      CALL eval_thermal_implicit(TRIM(models(i)),state,radiative,convective,ambient,process_rate)
      rate=rate+process_rate
    END DO
  END SUBROUTINE eval_thermal_implicit_processes

  ! Exact operators are composed on the updated temperature in declaration order.
  SUBROUTINE solve_thermal_analytic_processes(models,state,interval,heat_capacity, &
       radiative,convective,ambient,temperature_after)
    CHARACTER(*), INTENT(IN) :: models(:)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: interval,heat_capacity,radiative,convective,ambient
    REAL(wp), INTENT(OUT) :: temperature_after
    TYPE(physical_state_type) :: current_state
    INTEGER :: i
    current_state=state
    temperature_after=state%T
    DO i=1,SIZE(models)
      CALL solve_thermal_analytic(TRIM(models(i)),current_state,interval,heat_capacity, &
           radiative,convective,ambient,temperature_after)
      current_state%T=temperature_after
    END DO
  END SUBROUTINE solve_thermal_analytic_processes

END MODULE thermal_laws_2d
