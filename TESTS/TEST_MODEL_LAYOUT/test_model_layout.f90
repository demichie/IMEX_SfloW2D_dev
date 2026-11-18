PROGRAM test_model_layout

  USE model_layout_2d, ONLY : model_layout_type, supported_layer_count

  IMPLICIT NONE

  TYPE(model_layout_type) :: layout

  CALL assert_true('single layer supported', supported_layer_count(1))
  CALL assert_true('zero layers rejected', .NOT. supported_layer_count(0))
  CALL assert_true('two layers rejected for now', .NOT. supported_layer_count(2))

  CALL layout%initialize(1, 6, 6, 1, 1, 0, 0)

  CALL assert_true('basic layout initialized', layout%is_initialized())
  CALL assert_true('basic counts', layout%n_layers == 1 .AND.               &
       layout%n_vars == 6 .AND. layout%n_physical_variables == 6      &
       .AND. layout%n_eqns == 6)
  CALL assert_true('solid range', layout%idx_solid_first == 5 .AND.         &
       layout%idx_solid_last == 5)
  CALL assert_true('gas range', layout%idx_add_gas_first == 6        &
       .AND. layout%idx_add_gas_last == 6)
  CALL assert_true('physical velocity indices', layout%idx_u == 2 &
       .AND. layout%idx_v == 3)

  ! Reconfigure the same descriptor with every current optional component.
  CALL layout%initialize(1, 9, 9, 2, 1, 1, 1)

  CALL assert_true('full layout initialized', layout%is_initialized())
  CALL assert_true('full component counts', layout%n_solid == 2 .AND.       &
       layout%n_add_gas == 1 .AND. layout%n_stoch_vars == 1         &
       .AND. layout%n_pore_vars == 1 )
  CALL assert_true('full solid range', layout%idx_solid_first == 5 .AND.    &
       layout%idx_solid_last == 6)
  CALL assert_true('full optional indices',                                &
       layout%idx_add_gas_first == 7 .AND.                          &
       layout%idx_add_gas_last == 7 .AND.                           &
       layout%idx_stoch == 8 .AND.                                    &
       layout%idx_pore == 9)
  CALL assert_true('full velocity indices', layout%idx_u == 2    &
       .AND. layout%idx_v == 3)

  CALL layout%finalize
  CALL assert_true('layout finalized', .NOT. layout%is_initialized())

  WRITE(*,*) 'PASS: runtime model layout metadata verified'

CONTAINS

  SUBROUTINE assert_true(label, condition)

    CHARACTER(LEN=*), INTENT(IN) :: label
    LOGICAL, INTENT(IN) :: condition

    IF (.NOT. condition) THEN
       WRITE(*,*) 'FAIL: ', label
       ERROR STOP 1
    END IF

  END SUBROUTINE assert_true

END PROGRAM test_model_layout
