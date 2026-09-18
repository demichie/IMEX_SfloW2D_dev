!********************************************************************************
!> \brief Type-safe interfaces to the LAPACK dense linear solvers
!>
!> The generic solve_dense_system interface selects SGESV or DGESV from the
!> kind of the matrix. The wrappers adapt the solver's rank-one right-hand side
!> to LAPACK's rank-two B argument and keep raw LAPACK interfaces out of the
!> nonlinear solver.
!********************************************************************************
MODULE lapack_interfaces_2d

  USE parameters_2d, ONLY : sp, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: solve_dense_system

  INTERFACE solve_dense_system
     MODULE PROCEDURE solve_dense_system_sp
     MODULE PROCEDURE solve_dense_system_dp
  END INTERFACE solve_dense_system

  INTERFACE

     SUBROUTINE SGESV(n, nrhs, a, lda, ipiv, b, ldb, info)
       IMPORT :: sp
       INTEGER, INTENT(IN) :: n, nrhs, lda, ldb
       INTEGER, INTENT(OUT) :: ipiv(*)
       REAL(sp), INTENT(INOUT) :: a(lda,*), b(ldb,*)
       INTEGER, INTENT(OUT) :: info
     END SUBROUTINE SGESV

     SUBROUTINE DGESV(n, nrhs, a, lda, ipiv, b, ldb, info)
       IMPORT :: dp
       INTEGER, INTENT(IN) :: n, nrhs, lda, ldb
       INTEGER, INTENT(OUT) :: ipiv(*)
       REAL(dp), INTENT(INOUT) :: a(lda,*), b(ldb,*)
       INTEGER, INTENT(OUT) :: info
     END SUBROUTINE DGESV

  END INTERFACE

CONTAINS

  !******************************************************************************
  !> \brief Solve A*x=b in single precision using SGESV.
  !******************************************************************************
  SUBROUTINE solve_dense_system_sp(a, b, ipiv, info)

    REAL(sp), CONTIGUOUS, INTENT(INOUT) :: a(:,:)
    REAL(sp), CONTIGUOUS, INTENT(INOUT) :: b(:)
    INTEGER, INTENT(OUT) :: ipiv(:)
    INTEGER, INTENT(OUT) :: info

    REAL(sp) :: rhs(SIZE(b),1)
    INTEGER :: n

    n = SIZE(b)
    IF ((SIZE(a,1) /= n) .OR. (SIZE(a,2) /= n) .OR. (SIZE(ipiv) < n)) THEN
       info = -1
       RETURN
    END IF

    rhs(:,1) = b
    CALL SGESV(n, 1, a, n, ipiv, rhs, n, info)
    b = rhs(:,1)

  END SUBROUTINE solve_dense_system_sp

  !******************************************************************************
  !> \brief Solve A*x=b in double precision using DGESV.
  !******************************************************************************
  SUBROUTINE solve_dense_system_dp(a, b, ipiv, info)

    REAL(dp), CONTIGUOUS, INTENT(INOUT) :: a(:,:)
    REAL(dp), CONTIGUOUS, INTENT(INOUT) :: b(:)
    INTEGER, INTENT(OUT) :: ipiv(:)
    INTEGER, INTENT(OUT) :: info

    REAL(dp) :: rhs(SIZE(b),1)
    INTEGER :: n

    n = SIZE(b)
    IF ((SIZE(a,1) /= n) .OR. (SIZE(a,2) /= n) .OR. (SIZE(ipiv) < n)) THEN
       info = -1
       RETURN
    END IF

    rhs(:,1) = b
    CALL DGESV(n, 1, a, n, ipiv, rhs, n, info)
    b = rhs(:,1)

  END SUBROUTINE solve_dense_system_dp

END MODULE lapack_interfaces_2d
