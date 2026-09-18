PROGRAM test_diagnostics

  USE diagnostics_2d, ONLY : interactive_debug_flag, debug_pause, fatal_error

  IMPLICIT NONE

  CHARACTER(LEN=16) :: mode

  CALL get_command_argument(1, mode)

  SELECT CASE (TRIM(mode))

  CASE ('noninteractive')
     interactive_debug_flag = .FALSE.
     CALL debug_pause('disabled pause must return immediately')

  CASE ('interactive')
     interactive_debug_flag = .TRUE.
     CALL debug_pause('enabled serial pause')

  CASE ('parallel')
     interactive_debug_flag = .TRUE.
     !$OMP PARALLEL
     CALL debug_pause('parallel pause must be skipped')
     !$OMP END PARALLEL

  CASE ('fatal')
     CALL fatal_error('expected test failure')

  CASE DEFAULT
     ERROR STOP 'unknown diagnostics test mode'

  END SELECT

  WRITE(*,*) 'PASS: diagnostics mode ', TRIM(mode)

END PROGRAM test_diagnostics
