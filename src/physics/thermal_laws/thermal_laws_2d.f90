!> Dispatch thermal source models through a common four-treatment interface.
!> Add future models in the SELECT CASE statements below.
MODULE thermal_laws_2d
  USE parameters_2d, ONLY : wp
  USE diagnostics_2d, ONLY : fatal_error
  USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
  USE thermal_cooling_2d, ONLY : eval_thermal_cooling_explicit, &
       eval_thermal_cooling_semi_implicit, eval_thermal_cooling_implicit, &
       solve_thermal_cooling_analytic
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_thermal_explicit, eval_thermal_semi_implicit
  PUBLIC :: eval_thermal_implicit, solve_thermal_analytic
CONTAINS

  SUBROUTINE eval_thermal_explicit(models,state,radiative,convective,ambient,rate)
    CHARACTER(*), INTENT(IN) :: models(:)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: radiative,convective,ambient
    REAL(wp), INTENT(OUT) :: rate
    REAL(wp) :: process_rate
    INTEGER :: i
    rate=0.0_wp
    DO i=1,SIZE(models)
      SELECT CASE (TRIM(models(i)))
      CASE ('thermal_cooling')
        CALL eval_thermal_cooling_explicit(state,radiative,convective,ambient,process_rate)
      CASE DEFAULT
        CALL fatal_error('Unknown explicit thermal model: '//TRIM(models(i)))
      END SELECT
      rate=rate+process_rate
    END DO
  END SUBROUTINE eval_thermal_explicit


  SUBROUTINE eval_thermal_semi_implicit(models,state,rate)
    CHARACTER(*), INTENT(IN) :: models(:)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(OUT) :: rate
    REAL(wp) :: process_rate
    INTEGER :: i
    rate=0.0_wp
    DO i=1,SIZE(models)
      SELECT CASE (TRIM(models(i)))
      CASE ('thermal_cooling')
        CALL eval_thermal_cooling_semi_implicit(state,process_rate)
      CASE DEFAULT
        CALL fatal_error('Unknown semi-implicit thermal model: '//TRIM(models(i)))
      END SELECT
      rate=rate+process_rate
    END DO
  END SUBROUTINE eval_thermal_semi_implicit


  SUBROUTINE eval_thermal_implicit(models,state,radiative,convective,ambient,rate)
    CHARACTER(*), INTENT(IN) :: models(:)
    TYPE(complex_physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: radiative,convective,ambient
    COMPLEX(wp), INTENT(OUT) :: rate
    COMPLEX(wp) :: process_rate
    INTEGER :: i
    rate=CMPLX(0.0_wp,0.0_wp,wp)
    DO i=1,SIZE(models)
      SELECT CASE (TRIM(models(i)))
      CASE ('thermal_cooling')
        CALL eval_thermal_cooling_implicit(state,radiative,convective,ambient,process_rate)
      CASE DEFAULT
        CALL fatal_error('Unknown implicit thermal model: '//TRIM(models(i)))
      END SELECT
      rate=rate+process_rate
    END DO
  END SUBROUTINE eval_thermal_implicit


  SUBROUTINE solve_thermal_analytic(models,state,interval,heat_capacity, &
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
      SELECT CASE (TRIM(models(i)))
      CASE ('thermal_cooling')
        CALL solve_thermal_cooling_analytic(current_state,interval,heat_capacity, &
           radiative,convective,ambient,temperature_after)
      CASE DEFAULT
        CALL fatal_error('Unknown analytic thermal model: '//TRIM(models(i)))
      END SELECT
      current_state%T=temperature_after
    END DO
  END SUBROUTINE solve_thermal_analytic

END MODULE thermal_laws_2d
