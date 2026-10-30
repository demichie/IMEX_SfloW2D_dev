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
END MODULE thermal_laws_2d
