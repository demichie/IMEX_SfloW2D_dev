!> \brief Exercise stage-count validation through real input and workspace entry points.
!> \details Arguments select input parsing or direct workspace initialization and N_RK.
!>          Invalid workspace requests must fail before allocation: dimensions are unset.
PROGRAM test_stage_count
  USE parameters_2d, ONLY : n_RK, n_vars, n_eqns
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE time_integration_2d, ONLY : time_integration_workspace_type
  USE inpout_2d, ONLY : init_param, read_param
  IMPLICIT NONE
  TYPE(time_integration_workspace_type) :: integration
  CHARACTER(LEN=32) :: mode, argument
  INTEGER :: expected

  CALL get_command_argument(1,mode)
  CALL get_command_argument(2,argument)
  READ(argument,*) expected
  SELECT CASE (TRIM(mode))
  CASE ('input')
     CALL init_param
     CALL read_param
     IF (n_RK /= expected) ERROR STOP 'input stage count changed unexpectedly'
  CASE ('workspace')
     n_RK = expected
     ! Valid requests need sizes; invalid requests must be rejected before using them.
     IF (expected >= 2 .AND. expected <= 4) THEN
        n_vars=5; n_eqns=5; comp_cells_x=2; comp_cells_y=2
     END IF
     CALL integration%initialize
     CALL integration%finalize
  CASE DEFAULT
     ERROR STOP 'unknown stage-count test mode'
  END SELECT
  WRITE(*,*) 'PASS: stage count accepted:', n_RK
END PROGRAM test_stage_count
