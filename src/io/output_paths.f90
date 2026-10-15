!> \brief Centralized per-run output path construction and safe directory setup.
!
!> The module keeps all solver-generated files below the configured run directory.
MODULE output_paths

  USE, INTRINSIC :: iso_c_binding, ONLY: c_char, c_int, c_loc, c_null_char, c_ptr, c_size_t
  USE parameters_2d, ONLY: wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: initialize_output_paths, output_file_path, begin_snapshot, snapshot_file_path
  PUBLIC :: ensure_output_directory

  CHARACTER(LEN=1024), SAVE :: run_output_directory = ''
  CHARACTER(LEN=:), ALLOCATABLE, SAVE :: active_snapshot_directory

  INTERFACE
    FUNCTION c_prepare_run_directory(base, run, overwrite, allow_existing, &
                                     error_message, error_size)             &
        BIND(C, NAME='imex_prepare_run_directory') RESULT(status)
      IMPORT :: c_char, c_int, c_ptr, c_size_t
      TYPE(c_ptr), VALUE :: base, run, error_message
      INTEGER(c_int), VALUE :: overwrite, allow_existing
      INTEGER(c_size_t), VALUE :: error_size
      INTEGER(c_int) :: status
    END FUNCTION c_prepare_run_directory

    FUNCTION c_ensure_directory(path, error_message, error_size)            &
        BIND(C, NAME='imex_ensure_directory') RESULT(status)
      IMPORT :: c_int, c_ptr, c_size_t
      TYPE(c_ptr), VALUE :: path, error_message
      INTEGER(c_size_t), VALUE :: error_size
      INTEGER(c_int) :: status
    END FUNCTION c_ensure_directory
  END INTERFACE

CONTAINS

  !> \brief Create and validate the run output directory tree.
  !> \param[in] base_directory Parent output directory.
  !> \param[in] run_name Validated run identifier.
  !> \param[in] overwrite Whether an existing run directory may be reused.
  !> \param[in] allow_existing Whether a binary restart may resume in place.
  SUBROUTINE initialize_output_paths(base_directory, run_name, overwrite, allow_existing)
    CHARACTER(LEN=*), INTENT(IN) :: base_directory, run_name
    LOGICAL, INTENT(IN) :: overwrite, allow_existing
    CHARACTER(KIND=c_char), ALLOCATABLE, TARGET :: c_base(:), c_run(:), c_error(:)
    CHARACTER(LEN=1024) :: error_message
    INTEGER(c_int) :: status
    INTEGER :: i

    IF (LEN_TRIM(base_directory) == 0 .OR. LEN_TRIM(run_name) == 0) &
         CALL fail('OUTPUT_BASE_DIRECTORY and RUN_NAME must not be empty')

    IF (LEN_TRIM(base_directory) + LEN_TRIM(run_name) + 1 > &
        LEN(run_output_directory)) CALL fail('Output run path is too long')

    IF (INDEX(base_directory, ACHAR(0)) /= 0 .OR. INDEX(run_name, ACHAR(0)) /= 0) &
         CALL fail('NUL characters are not allowed in output paths')

    ALLOCATE(c_base(LEN_TRIM(base_directory) + 1), &
             c_run(LEN_TRIM(run_name) + 1), c_error(LEN(error_message)))
    DO i = 1, LEN_TRIM(base_directory)
      c_base(i) = base_directory(i:i)
    END DO
    c_base(SIZE(c_base)) = c_null_char
    DO i = 1, LEN_TRIM(run_name)
      c_run(i) = run_name(i:i)
    END DO
    c_run(SIZE(c_run)) = c_null_char
    c_error = c_null_char

    status = c_prepare_run_directory(c_loc(c_base(1)), c_loc(c_run(1)), &
         MERGE(1_c_int, 0_c_int, overwrite), &
         MERGE(1_c_int, 0_c_int, allow_existing), c_loc(c_error(1)), &
         INT(SIZE(c_error), c_size_t))
    IF (status /= 0_c_int) THEN
      error_message = ''
      DO i = 1, SIZE(c_error)
        IF (c_error(i) == c_null_char) EXIT
        error_message(i:i) = c_error(i)
      END DO
      CALL fail(TRIM(error_message))
    END IF

    run_output_directory = TRIM(base_directory)
    IF (run_output_directory(LEN_TRIM(run_output_directory):LEN_TRIM(run_output_directory)) &
        /= '/') run_output_directory = TRIM(run_output_directory)//'/'
    run_output_directory = TRIM(run_output_directory)//TRIM(run_name)
  END SUBROUTINE initialize_output_paths

  !> \brief Construct a path inside a named run-output category.
  !> \param[in] category One of the predefined output subdirectories.
  !> \param[in] filename Basename of the output file.
  !> \return Path to the requested output file.
  FUNCTION output_file_path(category, filename) RESULT(path)
    CHARACTER(LEN=*), INTENT(IN) :: category, filename
    CHARACTER(LEN=:), ALLOCATABLE :: path

    IF (LEN_TRIM(run_output_directory) == 0) CALL fail('Output paths are not initialized')
    IF (LEN_TRIM(run_output_directory) + LEN_TRIM(category) + LEN_TRIM(filename) + 2 &
        > LEN(run_output_directory)) CALL fail('Output file path is too long')
    path = TRIM(run_output_directory)//'/'//TRIM(category)//'/'//TRIM(filename)
  END FUNCTION output_file_path

  !> \brief Create the directory for a scheduled snapshot.
  !> \param[in] index Zero-based snapshot index.
  !> \param[in] time Exact scheduled snapshot time [s].
  SUBROUTINE begin_snapshot(index, time)
    INTEGER, INTENT(IN) :: index
    REAL(wp), INTENT(IN) :: time
    CHARACTER(LEN=4) :: index_string
    CHARACTER(LEN=32) :: time_label
    INTEGER :: i

    WRITE(index_string, '(I4.4)') index
    WRITE(time_label, '(ES24.16E3)') time
    time_label = ADJUSTL(time_label)
    DO i = 1, LEN_TRIM(time_label)
      SELECT CASE (time_label(i:i))
      CASE (' ')
        time_label(i:i) = '_'
      CASE ('+')
        time_label(i:i) = 'p'
      CASE ('-')
        time_label(i:i) = 'm'
      CASE ('.')
        time_label(i:i) = 'd'
      END SELECT
    END DO
    IF (LEN_TRIM(run_output_directory) + LEN('/snapshots/') + LEN_TRIM(index_string) + &
        LEN('_t_') + LEN_TRIM(time_label) > LEN(run_output_directory)) &
         CALL fail('Snapshot directory path is too long')
    active_snapshot_directory = TRIM(run_output_directory)//'/snapshots/'// &
         TRIM(index_string)//'_t_'//TRIM(time_label)
    CALL ensure_output_directory(active_snapshot_directory)
  END SUBROUTINE begin_snapshot

  !> \brief Construct a path inside the active snapshot directory.
  !> \param[in] filename Basename of the snapshot file.
  !> \return Path to the requested snapshot file.
  FUNCTION snapshot_file_path(filename) RESULT(path)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=:), ALLOCATABLE :: path

    IF (.NOT. ALLOCATED(active_snapshot_directory)) &
         CALL fail('Snapshot output path is not initialized')
    IF (LEN_TRIM(active_snapshot_directory) + 1 + LEN_TRIM(filename) > &
        LEN(run_output_directory)) CALL fail('Snapshot file path is too long')
    path = TRIM(active_snapshot_directory)//'/'//TRIM(filename)
  END FUNCTION snapshot_file_path

  !> \brief Ensure a directory exists without following symlink components.
  !> \param[in] directory Directory path to create or validate.
  SUBROUTINE ensure_output_directory(directory)
    CHARACTER(LEN=*), INTENT(IN) :: directory
    CHARACTER(KIND=c_char), ALLOCATABLE, TARGET :: c_path(:), c_error(:)
    CHARACTER(LEN=1024) :: error_message
    INTEGER(c_int) :: status
    INTEGER :: i

    IF (LEN_TRIM(directory) == 0 .OR. LEN_TRIM(directory) > LEN(run_output_directory) .OR. &
        INDEX(directory, ACHAR(0)) /= 0) &
         CALL fail('Invalid output directory path')

    ALLOCATE(c_path(LEN_TRIM(directory) + 1), c_error(LEN(error_message)))
    DO i = 1, LEN_TRIM(directory)
      c_path(i) = directory(i:i)
    END DO
    c_path(SIZE(c_path)) = c_null_char
    c_error = c_null_char
    status = c_ensure_directory(c_loc(c_path(1)), c_loc(c_error(1)), &
         INT(SIZE(c_error), c_size_t))
    IF (status /= 0_c_int) THEN
      error_message = ''
      DO i = 1, SIZE(c_error)
        IF (c_error(i) == c_null_char) EXIT
        error_message(i:i) = c_error(i)
      END DO
      CALL fail(TRIM(error_message))
    END IF
  END SUBROUTINE ensure_output_directory

  !> \brief Report an output path error using the solver's standard fatal diagnostic.
  !> \param[in] message Description of the invalid path or failed operation.
  SUBROUTINE fail(message)
    USE diagnostics_2d, ONLY: fatal_error
    CHARACTER(LEN=*), INTENT(IN) :: message
    CALL fatal_error('Output path error: '//TRIM(message))
  END SUBROUTINE fail

END MODULE output_paths
