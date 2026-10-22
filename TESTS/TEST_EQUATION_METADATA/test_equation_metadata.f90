PROGRAM test_equation_metadata

  USE equation_metadata_2d, ONLY : equation_partition_type,                  &
       IM_GROUP_METHOD_ANALYTIC, IM_GROUP_METHOD_NEWTON

  IMPLICIT NONE

  TYPE(equation_partition_type) :: partition

  LOGICAL :: implicit_mask(7)
  LOGICAL :: implicit_mask_five(5)

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
  CALL assert_true( 'default group count', partition%n_IM_groups .EQ. 1 )
  CALL assert_true( 'default group membership',                            &
       ALL(partition%IM_group .EQ. [ 0, 1, 1, 0, 0, 0, 1 ]) )
  CALL assert_true( 'default group method',                                &
       partition%IM_group_method(1) .EQ. IM_GROUP_METHOD_NEWTON )
  CALL assert_true( 'default group size',                                  &
       ALL(partition%n_IM_group_equations .EQ. [ 3 ]) )
  CALL assert_true( 'default compact group map',                           &
       ALL(partition%IM_group_map(:,1) .EQ. [ 2, 3, 7 ]) )

  ! Reconfiguration must replace the previous allocation and support an empty
  ! implicit partition.
  implicit_mask = .FALSE.
  CALL partition%initialize( implicit_mask )

  CALL assert_true( 'all-explicit partition initialized',                   &
       partition%is_initialized() )
  CALL assert_true( 'zero implicit equations', partition%n_implicit .EQ. 0 )
  CALL assert_true( 'all equations explicit', partition%n_explicit .EQ. 7 )
  CALL assert_true( 'no groups for all-explicit partition',                 &
       partition%n_IM_groups .EQ. 0 .AND. SIZE(partition%IM_group_method) .EQ. 0 )
  CALL assert_true( 'empty compact group map',                              &
       SIZE(partition%IM_group_map,1) .EQ. 0 .AND.                         &
       SIZE(partition%IM_group_map,2) .EQ. 0 )
  CALL assert_true( 'zero group for all-explicit equations',                &
       ALL(partition%IM_group .EQ. 0) )
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
  CALL assert_true( 'all-implicit default group',                           &
       partition%n_IM_groups .EQ. 1 .AND. ALL(partition%IM_group .EQ. 1) )
  CALL assert_true( 'all-implicit group map',                              &
       ALL(partition%IM_group_map(:,1) .EQ. [ 1, 2, 3, 4, 5, 6, 7 ]) )

  ! The documented D11 layout keeps energy, pore pressure and the coupled
  ! momenta in ordered groups with independent solve methods.
  implicit_mask_five = [ .FALSE., .TRUE., .TRUE., .TRUE., .TRUE. ]
  CALL partition%initialize( implicit_mask_five, [ 0, 3, 3, 1, 2 ],        &
       [ IM_GROUP_METHOD_ANALYTIC, IM_GROUP_METHOD_ANALYTIC,                &
       IM_GROUP_METHOD_NEWTON ] )

  CALL assert_true( 'multi-group partition initialized', partition%is_initialized() )
  CALL assert_true( 'multi-group count', partition%n_IM_groups .EQ. 3 )
  CALL assert_true( 'documented group membership',                          &
       ALL(partition%IM_group .EQ. [ 0, 3, 3, 1, 2 ]) )
  CALL assert_true( 'documented group methods',                              &
       ALL(partition%IM_group_method .EQ. [ IM_GROUP_METHOD_ANALYTIC,       &
       IM_GROUP_METHOD_ANALYTIC, IM_GROUP_METHOD_NEWTON ]) )
  CALL assert_true( 'documented group sizes',                               &
       ALL(partition%n_IM_group_equations .EQ. [ 1, 1, 2 ]) )
  CALL assert_true( 'documented compact group maps',                        &
       ALL(partition%IM_group_map(:,1) .EQ. [ 4, 0, 0, 0 ]) .AND.           &
       ALL(partition%IM_group_map(:,2) .EQ. [ 5, 0, 0, 0 ]) .AND.           &
       ALL(partition%IM_group_map(:,3) .EQ. [ 2, 3, 0, 0 ]) )

  ! The public consistency check rejects metadata mutations that violate the
  ! partition, dense group numbering or supported method identifiers.
  partition%IM_group(1) = 1
  CALL assert_true( 'explicit equation cannot have a group',                &
       .NOT. partition%is_initialized() )
  partition%IM_group(1) = 0
  partition%implicit(1) = .TRUE.
  CALL assert_true( 'implicit count must match partition mask',             &
       .NOT. partition%is_initialized() )
  partition%implicit(1) = .FALSE.

  partition%IM_group(2:3) = 2
  CALL assert_true( 'group numbering must be dense',                        &
       .NOT. partition%is_initialized() )
  partition%IM_group = [ 0, 3, 3, 1, 2 ]

  partition%IM_group_method(1) = 99
  CALL assert_true( 'unsupported group method rejected',                    &
       .NOT. partition%is_initialized() )
  partition%IM_group_method(1) = IM_GROUP_METHOD_ANALYTIC
  partition%IM_group_map(1,3) = 7
  CALL assert_true( 'compact map must match group membership',              &
       .NOT. partition%is_initialized() )
  partition%IM_group_map(1,3) = 2
  CALL assert_true( 'restored multi-group partition valid', partition%is_initialized() )

  CALL partition%finalize

  CALL assert_true( 'partition finalized', .NOT. partition%is_initialized() )
  CALL assert_true( 'finalized counts', partition%n_equations .EQ. 0        &
       .AND. partition%n_implicit .EQ. 0 .AND. partition%n_explicit .EQ. 0 &
       .AND. partition%n_IM_groups .EQ. 0 )

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
