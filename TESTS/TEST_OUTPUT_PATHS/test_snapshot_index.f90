PROGRAM test_snapshot_index

  USE output_paths, ONLY: snapshot_index_label
  USE output_paths, ONLY: initialize_output_paths, begin_snapshot, snapshot_file_path
  USE parameters_2d, ONLY: wp

  IMPLICIT NONE

  CHARACTER(LEN=1024) :: base_directory
  CHARACTER(LEN=1024) :: snapshot_paths(4)
  INTEGER :: path_count

  CALL GET_COMMAND_ARGUMENT(1, base_directory)
  IF (LEN_TRIM(base_directory) == 0) ERROR STOP 'Missing output base directory'
  CALL initialize_output_paths(TRIM(base_directory), 'snapshot-index-test', .FALSE., .FALSE.)

  path_count = 0
  CALL expect_snapshot(0, '0000')
  CALL expect_snapshot(9999, '9999')
  CALL expect_snapshot(10000, '10000')
  CALL expect_snapshot(10001, '10001')

CONTAINS

  SUBROUTINE expect_snapshot(index, expected)
    INTEGER, INTENT(IN) :: index
    CHARACTER(LEN=*), INTENT(IN) :: expected
    CHARACTER(LEN=:), ALLOCATABLE :: actual
    CHARACTER(LEN=:), ALLOCATABLE :: marker_path
    INTEGER :: i, io_status, unit

    actual = snapshot_index_label(index)
    IF (actual /= expected) THEN
      WRITE (*, '(A,I0,A,A,A,A)') 'FAIL: index ', index, &
           ' produced "', actual, '", expected "', expected, '"'
      ERROR STOP 1
    END IF

    CALL begin_snapshot(index, REAL(index, wp))
    marker_path = snapshot_file_path('marker')
    DO i = 1, path_count
      IF (TRIM(snapshot_paths(i)) == marker_path) &
           ERROR STOP 'Snapshot indices produced duplicate directories'
    END DO
    path_count = path_count + 1
    snapshot_paths(path_count) = marker_path
    OPEN(NEWUNIT=unit, FILE=marker_path, STATUS='NEW', ACTION='WRITE', IOSTAT=io_status)
    IF (io_status /= 0) ERROR STOP 'Snapshot directory was not created uniquely'
    CLOSE(unit)
  END SUBROUTINE expect_snapshot

END PROGRAM test_snapshot_index
