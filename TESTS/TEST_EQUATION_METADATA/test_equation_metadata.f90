PROGRAM test_equation_metadata

  USE equation_metadata_2d, ONLY : equation_partition_type

  IMPLICIT NONE

  TYPE(equation_partition_type) :: partition

  LOGICAL :: implicit_mask(7)

  implicit_mask = [ .FALSE., .TRUE., .TRUE., .FALSE., .FALSE., .FALSE.,    &
       .TRUE. ]

  CALL partition%initialize( implicit_mask )

  CALL assert_true( 'partition initialized', partition%is_initialized() )
  CALL assert_true( 'equation count', partition%n_equations .EQ. 7 )
  CALL assert_true( 'implicit count', partition%n_implicit .EQ. 3 )
  CALL assert_true( 'explicit count', partition%n_explicit .EQ. 4 )
  CALL assert_true( 'implicit mask', ALL(partition%implicit .EQV. implicit_mask) )
  CALL assert_true( 'implicit map', ALL(partition%implicit_map .EQ. [ 2, 3, 7 ]) )
  CALL assert_true( 'explicit map', ALL(partition%explicit_map .EQ. [ 1, 4, 5, 6 ]) )

  ! Reconfiguration must replace the previous allocation and support an empty
  ! implicit partition.
  implicit_mask = .FALSE.
  CALL partition%initialize( implicit_mask )

  CALL assert_true( 'all-explicit partition initialized',                   &
       partition%is_initialized() )
  CALL assert_true( 'zero implicit equations', partition%n_implicit .EQ. 0 )
  CALL assert_true( 'all equations explicit', partition%n_explicit .EQ. 7 )
  CALL assert_true( 'zero-sized implicit map', SIZE(partition%implicit_map) .EQ. 0 )
  CALL assert_true( 'identity explicit map',                                &
       ALL(partition%explicit_map .EQ. [ 1, 2, 3, 4, 5, 6, 7 ]) )

  ! The opposite degenerate partition must also be valid.
  implicit_mask = .TRUE.
  CALL partition%initialize( implicit_mask )

  CALL assert_true( 'all-implicit partition initialized',                   &
       partition%is_initialized() )
  CALL assert_true( 'all equations implicit', partition%n_implicit .EQ. 7 )
  CALL assert_true( 'zero explicit equations', partition%n_explicit .EQ. 0 )
  CALL assert_true( 'identity implicit map',                                &
       ALL(partition%implicit_map .EQ. [ 1, 2, 3, 4, 5, 6, 7 ]) )
  CALL assert_true( 'zero-sized explicit map', SIZE(partition%explicit_map) .EQ. 0 )

  CALL partition%finalize

  CALL assert_true( 'partition finalized', .NOT. partition%is_initialized() )
  CALL assert_true( 'finalized counts', partition%n_equations .EQ. 0        &
       .AND. partition%n_implicit .EQ. 0 .AND. partition%n_explicit .EQ. 0 )

  WRITE(*,*) 'PASS: equation partition metadata verified'

CONTAINS

  SUBROUTINE assert_true( label, condition )

    CHARACTER(LEN=*), INTENT(IN) :: label
    LOGICAL, INTENT(IN) :: condition

    IF ( .NOT. condition ) THEN
       WRITE(*,*) 'FAIL: ', label
       ERROR STOP 1
    END IF

  END SUBROUTINE assert_true

END PROGRAM test_equation_metadata
