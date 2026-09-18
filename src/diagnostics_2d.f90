!********************************************************************************
!> \brief Runtime diagnostics and controlled debug pauses
!>
!> Verbosity only controls diagnostic output elsewhere in the solver. Optional
!> interactive pauses are enabled independently and are never attempted from
!> an OpenMP parallel region. Fatal conditions always terminate with a nonzero
!> status so batch executions cannot remain blocked on standard input.
!********************************************************************************
MODULE diagnostics_2d

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit
  USE omp_lib, ONLY : omp_in_parallel

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PUBLIC :: interactive_debug_flag = .FALSE.

  PUBLIC :: debug_pause
  PUBLIC :: fatal_error

CONTAINS

  !******************************************************************************
  !> \brief Pause for interactive inspection when explicitly enabled.
  !>
  !> Pauses requested from an OpenMP parallel region are skipped because stdin
  !> access from worker threads is unsafe and can deadlock a batch execution.
  !******************************************************************************
  SUBROUTINE debug_pause(message)

    CHARACTER(LEN=*), INTENT(IN) :: message

    CHARACTER(LEN=1) :: response
    INTEGER :: io_status

    IF (.NOT. interactive_debug_flag) RETURN

    IF (omp_in_parallel()) THEN
       !$OMP CRITICAL(diagnostics_output)
       WRITE(error_unit,*) 'DEBUG PAUSE SKIPPED inside OpenMP region: ',     &
            TRIM(message)
       CALL FLUSH(error_unit)
       !$OMP END CRITICAL(diagnostics_output)
       RETURN
    END IF

    WRITE(output_unit,*) 'DEBUG PAUSE: ', TRIM(message)
    WRITE(output_unit,*) 'Press ENTER to continue.'
    CALL FLUSH(output_unit)
    READ(input_unit,'(A)',IOSTAT=io_status) response

    IF (io_status /= 0) THEN
       WRITE(error_unit,*) 'WARNING: Unable to read interactive debug input.'
       CALL FLUSH(error_unit)
    END IF

  END SUBROUTINE debug_pause

  !******************************************************************************
  !> \brief Report an unrecoverable condition and terminate with failure status.
  !******************************************************************************
  SUBROUTINE fatal_error(message)

    CHARACTER(LEN=*), INTENT(IN) :: message

    !$OMP CRITICAL(diagnostics_output)
    CALL FLUSH(output_unit)
    WRITE(error_unit,*) 'FATAL ERROR: ', TRIM(message)
    CALL FLUSH(error_unit)
    !$OMP END CRITICAL(diagnostics_output)

    IF (interactive_debug_flag .AND. (.NOT. omp_in_parallel()))             &
         CALL debug_pause('fatal condition: '//TRIM(message))

    ERROR STOP 1

  END SUBROUTINE fatal_error

END MODULE diagnostics_2d
