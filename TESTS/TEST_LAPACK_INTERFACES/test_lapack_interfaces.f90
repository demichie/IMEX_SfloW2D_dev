PROGRAM test_lapack_interfaces

  USE lapack_interfaces_2d, ONLY : solve_dense_system
  USE parameters_2d, ONLY : sp, dp

  IMPLICIT NONE

  REAL(sp) :: matrix_sp(2,2), rhs_sp(2)
  REAL(dp) :: matrix_dp(2,2), rhs_dp(2)
  INTEGER :: pivot(2), info

  matrix_sp = RESHAPE([ 3.0_sp, 1.0_sp, 1.0_sp, 2.0_sp ], [ 2, 2 ])
  rhs_sp = [ 9.0_sp, 8.0_sp ]
  CALL solve_dense_system(matrix_sp, rhs_sp, pivot, info)
  CALL assert_true('SGESV status', info == 0)
  CALL assert_true('SGESV solution',                                &
       MAXVAL(ABS(rhs_sp - [ 2.0_sp, 3.0_sp ])) < 1.0E-5_sp)

  matrix_dp = RESHAPE([ 3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp ], [ 2, 2 ])
  rhs_dp = [ 9.0_dp, 8.0_dp ]
  CALL solve_dense_system(matrix_dp, rhs_dp, pivot, info)
  CALL assert_true('DGESV status', info == 0)
  CALL assert_true('DGESV solution',                                &
       MAXVAL(ABS(rhs_dp - [ 2.0_dp, 3.0_dp ])) < 1.0E-12_dp)

  WRITE(*,*) 'PASS: single- and double-precision LAPACK interfaces verified'

CONTAINS

  SUBROUTINE assert_true(label, condition)

    CHARACTER(LEN=*), INTENT(IN) :: label
    LOGICAL, INTENT(IN) :: condition

    IF (.NOT. condition) THEN
       WRITE(*,*) 'FAIL: ', label
       ERROR STOP 1
    END IF

  END SUBROUTINE assert_true

END PROGRAM test_lapack_interfaces
