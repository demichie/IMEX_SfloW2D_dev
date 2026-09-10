PROGRAM test_nonconservative_path

  USE parameters_2d, ONLY : wp
  USE nonconservative_2d, ONLY : eval_path_contribution
  USE nonconservative_2d, ONLY : PATH_DIR_X, PATH_DIR_Y

  IMPLICIT NONE

  REAL(wp) :: contribution(2)
  REAL(wp) :: contribution_a(2), contribution_b(2)
  REAL(wp) :: contribution_reverse(2)
  REAL(wp) :: error_quad(3)
  REAL(wp) :: exact
  REAL(wp) :: stateL(2), stateR(2)
  REAL(wp) :: tolerance

  INTEGER :: polynomial_degree

  tolerance = 4096.0_wp * EPSILON(1.0_wp)

  ! Equal endpoints must produce a zero path contribution.
  stateL = [ 0.37_wp, -0.21_wp ]
  stateR = stateL
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, exponential_provider, &
       contribution )
  CALL assert_close_vector( 'equal-state path', contribution,               &
       [ 0.0_wp, 0.0_wp ], tolerance )

  ! Reversing a straight path must reverse its integral.
  stateL = [ -0.4_wp, 0.2_wp ]
  stateR = [ 0.9_wp, 1.1_wp ]
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, exponential_provider, &
       contribution )
  CALL eval_path_contribution( PATH_DIR_X, stateR, stateL, exponential_provider, &
       contribution_reverse )
  CALL assert_close_vector( 'path reversal', contribution_reverse,          &
       -contribution, tolerance * MAX(1.0_wp, MAXVAL(ABS(contribution))) )

  ! Gauss-Legendre rules integrate polynomials through degree 2*n-1.
  stateL = [ 0.0_wp, 0.0_wp ]
  stateR = [ 1.0_wp, 0.0_wp ]

  polynomial_degree = 1
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, polynomial_provider, &
       contribution, n_quad=1 )
  CALL assert_close_scalar( 'one-point polynomial quadrature', contribution(1), &
       0.5_wp, tolerance )

  polynomial_degree = 3
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, polynomial_provider, &
       contribution, n_quad=2 )
  CALL assert_close_scalar( 'two-point polynomial quadrature', contribution(1), &
       0.25_wp, tolerance )

  polynomial_degree = 5
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, polynomial_provider, &
       contribution, n_quad=3 )
  CALL assert_close_scalar( 'three-point polynomial quadrature', contribution(1), &
       1.0_wp / 6.0_wp, tolerance )

  ! A nonlinear manufactured product must converge as quadrature order grows.
  exact = EXP(1.0_wp) - 1.0_wp
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, exponential_provider, &
       contribution, n_quad=1 )
  error_quad(1) = ABS(contribution(1) - exact)
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, exponential_provider, &
       contribution, n_quad=2 )
  error_quad(2) = ABS(contribution(1) - exact)
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, exponential_provider, &
       contribution, n_quad=3 )
  error_quad(3) = ABS(contribution(1) - exact)
  IF ( .NOT. ( error_quad(3) .LT. error_quad(2) .AND.                       &
       error_quad(2) .LT. error_quad(1) ) ) THEN
     WRITE(*,*) 'FAIL: nonlinear quadrature errors = ', error_quad
     ERROR STOP 1
  END IF

  ! The sum returned by one provider must equal independently integrated terms.
  stateL = [ -0.3_wp, 0.4_wp ]
  stateR = [ 0.8_wp, 1.2_wp ]
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, product_a_provider, &
       contribution_a )
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, product_b_provider, &
       contribution_b )
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, product_sum_provider, &
       contribution )
  CALL assert_close_vector( 'sum of registered products', contribution,     &
       contribution_a + contribution_b, tolerance                           &
       * MAX(1.0_wp, MAXVAL(ABS(contribution))) )

  ! Matrix-form and column-wise factorized products must be equivalent.
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, matrix_provider,  &
       contribution_a )
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, factorized_provider, &
       contribution_b )
  CALL assert_close_vector( 'matrix/factorized equivalence', contribution_a, &
       contribution_b, tolerance * MAX(1.0_wp, MAXVAL(ABS(contribution_a))) )

  ! A composite scalar G=x*y must match its expanded chain rule.
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, composite_provider, &
       contribution_a )
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, chain_rule_provider, &
       contribution_b )
  CALL assert_close_vector( 'composite/chain-rule equivalence', contribution_a, &
       contribution_b, tolerance * MAX(1.0_wp, MAXVAL(ABS(contribution_a))) )

  ! The direction is dispatched to the model without changing path orientation.
  stateL = [ 0.2_wp, 0.0_wp ]
  stateR = [ 1.3_wp, 0.0_wp ]
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, directional_provider, &
       contribution_a )
  CALL eval_path_contribution( PATH_DIR_Y, stateL, stateR, directional_provider, &
       contribution_b )
  CALL assert_close_vector( 'x-direction dispatch', contribution_a,         &
       [ 1.1_wp, 0.0_wp ], tolerance )
  CALL assert_close_vector( 'y-direction dispatch', contribution_b,         &
       [ 0.0_wp, -1.1_wp ], tolerance )

  ! A model-defined curved path must override the default straight segment.
  stateL = [ 0.0_wp, 0.0_wp ]
  stateR = [ 1.0_wp, 1.0_wp ]
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, path_sensitive_provider, &
       contribution_a )
  CALL eval_path_contribution( PATH_DIR_X, stateL, stateR, path_sensitive_provider, &
       contribution_b, path_provider=quadratic_path )
  CALL assert_close_scalar( 'straight-path callback default', contribution_a(1), &
       0.5_wp, tolerance )
  CALL assert_close_scalar( 'custom quadratic path', contribution_b(1),     &
       1.0_wp / 3.0_wp, tolerance )

  WRITE(*,*) 'PASS: generic nonconservative path engine verified'

CONTAINS

  SUBROUTINE exponential_provider( direction, path_state, dstate_ds, integrand )
    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(IN) :: path_state(:), dstate_ds(:)
    REAL(wp), INTENT(OUT) :: integrand(:)
    integrand = 0.0_wp
    integrand(1) = EXP(path_state(1)) * dstate_ds(1)
  END SUBROUTINE exponential_provider


  SUBROUTINE polynomial_provider( direction, path_state, dstate_ds, integrand )
    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(IN) :: path_state(:), dstate_ds(:)
    REAL(wp), INTENT(OUT) :: integrand(:)
    integrand = 0.0_wp
    integrand(1) = path_state(1)**polynomial_degree * dstate_ds(1)
  END SUBROUTINE polynomial_provider


  SUBROUTINE product_a_provider( direction, path_state, dstate_ds, integrand )
    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(IN) :: path_state(:), dstate_ds(:)
    REAL(wp), INTENT(OUT) :: integrand(:)
    integrand = 0.0_wp
    integrand(1) = 2.0_wp * path_state(1) * dstate_ds(1)
    integrand(2) = path_state(2) * dstate_ds(1)
  END SUBROUTINE product_a_provider


  SUBROUTINE product_b_provider( direction, path_state, dstate_ds, integrand )
    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(IN) :: path_state(:), dstate_ds(:)
    REAL(wp), INTENT(OUT) :: integrand(:)
    integrand = 0.0_wp
    integrand(1) = 3.0_wp * path_state(2)**2 * dstate_ds(2)
    integrand(2) = path_state(1) * dstate_ds(2)
  END SUBROUTINE product_b_provider


  SUBROUTINE product_sum_provider( direction, path_state, dstate_ds, integrand )
    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(IN) :: path_state(:), dstate_ds(:)
    REAL(wp), INTENT(OUT) :: integrand(:)
    REAL(wp) :: term(2)
    CALL product_a_provider( direction, path_state, dstate_ds, integrand )
    CALL product_b_provider( direction, path_state, dstate_ds, term )
    integrand = integrand + term
  END SUBROUTINE product_sum_provider


  SUBROUTINE matrix_provider( direction, path_state, dstate_ds, integrand )
    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(IN) :: path_state(:), dstate_ds(:)
    REAL(wp), INTENT(OUT) :: integrand(:)
    REAL(wp) :: matrix(2,2)
    matrix(:,1) = [ path_state(1), path_state(2) ]
    matrix(:,2) = [ 2.0_wp * path_state(2), -path_state(1) ]
    integrand = MATMUL(matrix, dstate_ds)
  END SUBROUTINE matrix_provider


  SUBROUTINE factorized_provider( direction, path_state, dstate_ds, integrand )
    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(IN) :: path_state(:), dstate_ds(:)
    REAL(wp), INTENT(OUT) :: integrand(:)
    integrand = [ path_state(1), path_state(2) ] * dstate_ds(1)              &
         + [ 2.0_wp * path_state(2), -path_state(1) ] * dstate_ds(2)
  END SUBROUTINE factorized_provider


  SUBROUTINE composite_provider( direction, path_state, dstate_ds, integrand )
    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(IN) :: path_state(:), dstate_ds(:)
    REAL(wp), INTENT(OUT) :: integrand(:)
    REAL(wp) :: dG_ds
    dG_ds = path_state(2) * dstate_ds(1) + path_state(1) * dstate_ds(2)
    integrand = [ 1.0_wp + path_state(1), path_state(2) ] * dG_ds
  END SUBROUTINE composite_provider


  SUBROUTINE chain_rule_provider( direction, path_state, dstate_ds, integrand )
    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(IN) :: path_state(:), dstate_ds(:)
    REAL(wp), INTENT(OUT) :: integrand(:)
    REAL(wp) :: coefficient(2)
    coefficient = [ 1.0_wp + path_state(1), path_state(2) ]
    integrand = coefficient * path_state(2) * dstate_ds(1)                  &
         + coefficient * path_state(1) * dstate_ds(2)
  END SUBROUTINE chain_rule_provider


  SUBROUTINE directional_provider( direction, path_state, dstate_ds, integrand )
    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(IN) :: path_state(:), dstate_ds(:)
    REAL(wp), INTENT(OUT) :: integrand(:)
    integrand = 0.0_wp
    SELECT CASE ( direction )
    CASE ( PATH_DIR_X )
       integrand(1) = dstate_ds(1)
    CASE ( PATH_DIR_Y )
       integrand(2) = -dstate_ds(1)
    END SELECT
  END SUBROUTINE directional_provider


  SUBROUTINE path_sensitive_provider( direction, path_state, dstate_ds, integrand )
    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(IN) :: path_state(:), dstate_ds(:)
    REAL(wp), INTENT(OUT) :: integrand(:)
    integrand = 0.0_wp
    integrand(1) = path_state(2) * dstate_ds(1)
  END SUBROUTINE path_sensitive_provider


  SUBROUTINE quadratic_path( path_left, path_right, s, path_state, dstate_ds )
    REAL(wp), INTENT(IN) :: path_left(:), path_right(:), s
    REAL(wp), INTENT(OUT) :: path_state(:), dstate_ds(:)
    path_state(1) = path_left(1) + s * (path_right(1) - path_left(1))
    path_state(2) = path_left(2) + s**2 * (path_right(2) - path_left(2))
    dstate_ds(1) = path_right(1) - path_left(1)
    dstate_ds(2) = 2.0_wp * s * (path_right(2) - path_left(2))
  END SUBROUTINE quadratic_path


  SUBROUTINE assert_close_scalar( label, value, expected, atol )
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: value, expected, atol
    IF ( ABS(value - expected) .GT. atol ) THEN
       WRITE(*,*) 'FAIL: ', label
       WRITE(*,*) ' value, expected, tolerance = ', value, expected, atol
       ERROR STOP 1
    END IF
  END SUBROUTINE assert_close_scalar


  SUBROUTINE assert_close_vector( label, value, expected, atol )
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: value(:), expected(:), atol
    IF ( MAXVAL(ABS(value - expected)) .GT. atol ) THEN
       WRITE(*,*) 'FAIL: ', label
       WRITE(*,*) ' value = ', value
       WRITE(*,*) ' expected = ', expected
       WRITE(*,*) ' tolerance = ', atol
       ERROR STOP 1
    END IF
  END SUBROUTINE assert_close_vector

END PROGRAM test_nonconservative_path
