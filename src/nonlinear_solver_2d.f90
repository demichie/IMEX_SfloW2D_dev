!********************************************************************************
!> \brief Local nonlinear solver for the implicit IMEX stage
!>
!> Newton iteration, line search, residual/Jacobian evaluation and small dense
!> linear solves.  This module is intentionally independent of the spatial
!> domain and operates on one cell state at a time.
!********************************************************************************
MODULE nonlinear_solver_2d

  USE parameters_2d, ONLY : wp, sp, n_eqns, n_vars, n_nh, verbose_level
  USE constitutive_2d, ONLY : implicit_flag, implicit_map

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: initialize_nonlinear_solver
  PUBLIC :: finalize_nonlinear_solver
  PUBLIC :: solve_rk_step

  !> Map from explicit variables to the full system.
  INTEGER, ALLOCATABLE :: explicit_map(:)

  !> Complex-step perturbation and its reciprocal.
  REAL(wp) :: h, one_by_h

  !> Relative threshold used to reject numerically singular small systems.
  REAL(wp), PARAMETER :: pivot_tol_factor = 64.0_wp

CONTAINS

  SUBROUTINE initialize_nonlinear_solver

    INTEGER :: i, j

    h = n_vars * EPSILON(1.0_wp)
    one_by_h = 1.0_wp / h

    IF ( ALLOCATED(explicit_map) ) DEALLOCATE(explicit_map)
    ALLOCATE( explicit_map(n_eqns-n_nh) )

    j = 0
    DO i = 1,n_eqns
       IF ( .NOT.implicit_flag(i) ) THEN
          j = j + 1
          explicit_map(j) = i
       END IF
    END DO

  END SUBROUTINE initialize_nonlinear_solver

  SUBROUTINE finalize_nonlinear_solver

    IF ( ALLOCATED(explicit_map) ) DEALLOCATE(explicit_map)

  END SUBROUTINE finalize_nonlinear_solver

  !******************************************************************************

  SUBROUTINE solve_rk_step( qj, qj_old, dt_step, a_diag, Rj_not_impl,          &
       Bprimej_x, Bprimej_y, Zij, fric_val,                                    &
       iterations_used, converged, linear_info, line_search_failed )

    USE parameters_2d, ONLY : max_nl_iter , tol_rel , tol_abs

    USE constitutive_2d, ONLY : rheology_model

    USE constitutive_2d, ONLY : integrate_friction_term

    IMPLICIT NONE

    REAL(wp), INTENT(INOUT) :: qj(n_vars)
    REAL(wp), INTENT(IN) :: qj_old(n_vars)
    REAL(wp), INTENT(IN) :: dt_step
    REAL(wp), INTENT(IN) :: a_diag
    REAL(wp), INTENT(IN) :: Rj_not_impl(n_eqns)
    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y
    REAL(wp), INTENT(IN):: Zij ! value stochastic process
    REAL(wp), INTENT(OUT) :: fric_val ! to save the value of the friction
    INTEGER, INTENT(OUT) :: iterations_used
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: linear_info
    LOGICAL, INTENT(OUT) :: line_search_failed

    REAL(wp) :: qj_init(n_vars)

    REAL(wp) :: qj_org(n_vars) , qj_rel(n_vars)

    REAL(wp) :: left_matrix(n_eqns,n_vars)
    REAL(wp) :: right_term(n_eqns)

    REAL(wp) :: scal_f

    REAL(wp) :: coeff_f(n_eqns)
    REAL(wp) :: residual_ref(n_eqns)
    REAL(wp) :: residual_tol(n_eqns)

    REAL(wp) :: qj_rel_NR_old(n_vars)
    REAL(wp) :: scal_f_old
    REAL(wp) :: desc_dir(n_vars)
    REAL(wp) :: grad_f(n_vars)

    INTEGER :: pivot(n_vars)

    REAL(wp) :: left_matrix_small22(n_nh,n_nh)

    REAL(wp) :: desc_dir_small2(n_nh)
    INTEGER :: pivot_small2(n_nh)

    REAL(wp) :: desc_dir_small1(n_vars-n_nh)

    INTEGER :: ok

    INTEGER :: i,j
    INTEGER :: idx
    INTEGER :: nl_iter

    REAL(wp), PARAMETER :: STPMX=100.0_wp
    REAL(wp) :: stpmax
    LOGICAL :: check

    ! REAL(wp) :: qpj(n_vars+2) , p_dyn

    REAL(wp) :: desc_dir_temp(n_vars)

    iterations_used = 0
    converged = .FALSE.
    linear_info = 0
    line_search_failed = .FALSE.

    IF ( rheology_model .EQ. 8 ) THEN

       CALL integrate_friction_term( qj , dt_step )
       converged = .TRUE.
       RETURN

    END IF

    coeff_f(1:n_eqns) = 1.0_wp

    grad_f(1:n_eqns) = 0.0_wp

    qj_init = qj

    !---- normalize the conservative variables ------

    qj_org = MAX( ABS(qj) , 1.0E-3_wp )

    qj_rel = qj / qj_org
    check = .FALSE.

    ! -----------------------------------------------
    newton_raphson_loop:DO nl_iter=1,max_nl_iter

       iterations_used = nl_iter

       IF ( verbose_level .GE. 2 ) WRITE(*,*) 'solve_rk_step: nl_iter',nl_iter

       CALL eval_f( qj , qj_old , dt_step , a_diag , coeff_f ,             &
            Rj_not_impl , Bprimej_x , Bprimej_y , right_term , scal_f , Zij ,   &
            fric_val )

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*) 'solve_rk_step: right_term',right_term

       END IF

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*) 'before_lnsrch: scal_f',scal_f

       END IF

       ! check the residual of the system

       IF ( nl_iter .EQ. 1 ) residual_ref = ABS(right_term)
       residual_tol = tol_abs + tol_rel * residual_ref

       IF ( ALL( ABS(right_term) .LE. residual_tol ) ) THEN

          IF ( verbose_level .GE. 3 ) WRITE(*,*) '1: check',check
          converged = .TRUE.
          EXIT newton_raphson_loop

       END IF

       ! ---- evaluate the descent direction ------------------------------------

       CALL eval_jacobian( qj_rel , qj_org , dt_step , a_diag , coeff_f ,   &
            Bprimej_x , Bprimej_y , left_matrix, Zij, fric_val )

       ! DGESV/SGESV overwrite the Jacobian with its LU factors in the fully
       ! implicit case. Form the line-search gradient before the linear solve.
       IF ( nl_iter .GT. 1 ) grad_f = MATMUL( right_term , left_matrix )

       IF ( n_nh .EQ. n_eqns ) THEN

          desc_dir_temp = - right_term

          IF ( wp .EQ. sp ) THEN

             CALL SGESV(n_eqns,1, left_matrix , n_eqns, pivot, desc_dir_temp ,  &
                  n_eqns, ok)

          ELSE

             CALL DGESV(n_eqns,1, left_matrix , n_eqns, pivot, desc_dir_temp ,  &
                  n_eqns, ok)

          END IF

          IF ( ok .NE. 0 ) THEN
             linear_info = ok
             qj = qj_init
             RETURN
          END IF

          desc_dir = desc_dir_temp

       ELSE

          DO i=1,n_nh

             desc_dir_small2(i) = right_term(implicit_map(i))

             DO j=1,n_nh
                left_matrix_small22(i,j) =                                     &
                     left_matrix(implicit_map(i),implicit_map(j))
             END DO

          END DO

          ! Non-implicit columns are diagonal by construction in
          ! eval_jacobian; therefore the A21 block is identically zero.
          DO i=1,n_vars-n_nh

             idx = explicit_map(i)
             desc_dir_small1(i) = right_term(idx)

             IF ( ABS(left_matrix(idx,idx)) .LE. TINY(1.0_wp) ) THEN
                linear_info = i
                qj = qj_init
                RETURN
             END IF

             desc_dir_small1(i) = desc_dir_small1(i) / left_matrix(idx,idx)

          END DO

          IF ( n_nh .EQ. 2 ) THEN

             CALL solve_2x2_pivoted( left_matrix_small22,                     &
                  desc_dir_small2, ok )

             IF ( ok .NE. 0 ) THEN
                linear_info = ok
                qj = qj_init
                RETURN
             END IF

          ELSEIF ( n_nh .EQ. 3 ) THEN

             CALL solve_3x3_pivoted( left_matrix_small22,                      &
                  desc_dir_small2, ok )

             IF ( ok .NE. 0 ) THEN
                linear_info = ok
                qj = qj_init
                RETURN
             END IF

          ELSE

             IF ( wp .EQ. sp ) THEN

                CALL SGESV(n_nh,1, left_matrix_small22 , n_nh , pivot_small2 ,  &
                     desc_dir_small2 , n_nh, ok)

             ELSE

                CALL DGESV(n_nh,1, left_matrix_small22 , n_nh , pivot_small2 ,  &
                     desc_dir_small2 , n_nh, ok)

             END IF

             IF ( ok .NE. 0 ) THEN
                linear_info = ok
                qj = qj_init
                RETURN
             END IF

          END IF

          desc_dir = 0.0_wp

          DO i=1,n_nh
             desc_dir(implicit_map(i)) = -desc_dir_small2(i)
          END DO

          DO i=1,n_vars-n_nh
             desc_dir(explicit_map(i)) = -desc_dir_small1(i)
          END DO

       END IF

       IF ( verbose_level .GE. 3 ) WRITE(*,*) 'desc_dir',desc_dir

       qj_rel_NR_old = qj_rel
       scal_f_old = scal_f

       IF ( nl_iter .GT. 1 ) THEN
          ! Search for the step lambda giving a suffic. decrease in the solution

          stpmax = STPMX * MAX( SQRT( DOT_PRODUCT(qj_rel,qj_rel) ) ,            &
               DBLE( SIZE(qj_rel) ) )

          CALL lnsrch( qj_rel_NR_old , qj_org , qj_old , scal_f_old , grad_f ,  &
               desc_dir , coeff_f , qj_rel , scal_f , right_term , stpmax ,     &
               check , dt_step , a_diag , Rj_not_impl , Bprimej_x , Bprimej_y,  &
               Zij, fric_val )

          IF ( check ) THEN
             qj = qj_rel * qj_org
             line_search_failed = .TRUE.
             IF ( verbose_level .GE. 2 )                                      &
                  WRITE(*,*) 'solve_rk_step: line search failed'
             RETURN
          END IF

       ELSE

          qj_rel = qj_rel_NR_old + desc_dir

          qj = qj_rel * qj_org

          CALL eval_f( qj , qj_old , dt_step , a_diag , coeff_f ,          &
               Rj_not_impl , Bprimej_x , Bprimej_y , right_term , scal_f, Zij,  &
               fric_val )

       END IF

       IF ( verbose_level .GE. 2 ) WRITE(*,*) 'after_lnsrch: scal_f',scal_f

       qj = qj_rel * qj_org

       IF ( verbose_level .GE. 3 ) THEN

          WRITE(*,*) 'qj',qj

       END IF

       IF ( ALL( ABS(right_term) .LE. residual_tol ) ) THEN

          IF ( verbose_level .GE. 3 ) WRITE(*,*) '1: check',check
          check= .FALSE.
          converged = .TRUE.
          EXIT newton_raphson_loop

       END IF

       IF ( MAXVAL( ABS( qj_rel(:) - qj_rel_NR_old(:) ) / MAX( ABS( qj_rel(:)) ,&
            1.0_wp ) ) < EPSILON(1.0_wp) ) THEN

          IF ( verbose_level .GE. 2 )                                         &
               WRITE(*,*) 'solve_rk_step: stagnation before convergence'
          RETURN

       END IF

    END DO newton_raphson_loop

    RETURN

  END SUBROUTINE solve_rk_step

  !******************************************************************************
  !> \brief Solve a 2x2 linear system with partial pivoting
  !
  !> The matrix and right-hand side are overwritten with the elimination
  !> factors and the solution, respectively. A nonzero info value identifies
  !> the first numerically singular pivot.
  !******************************************************************************

  SUBROUTINE solve_2x2_pivoted( matrix, rhs, info )

    IMPLICIT NONE

    REAL(wp), INTENT(INOUT) :: matrix(2,2)
    REAL(wp), INTENT(INOUT) :: rhs(2)
    INTEGER, INTENT(OUT) :: info

    REAL(wp) :: factor
    REAL(wp) :: matrix_scale
    REAL(wp) :: pivot_tol
    REAL(wp) :: swap_value

    info = 0
    matrix_scale = MAXVAL(ABS(matrix))

    IF ( matrix_scale .LE. TINY(1.0_wp) ) THEN
       info = 1
       RETURN
    END IF

    pivot_tol = pivot_tol_factor * EPSILON(1.0_wp) * matrix_scale

    IF ( ABS(matrix(2,1)) .GT. ABS(matrix(1,1)) ) THEN
       swap_value = matrix(1,1)
       matrix(1,1) = matrix(2,1)
       matrix(2,1) = swap_value

       swap_value = matrix(1,2)
       matrix(1,2) = matrix(2,2)
       matrix(2,2) = swap_value

       swap_value = rhs(1)
       rhs(1) = rhs(2)
       rhs(2) = swap_value
    END IF

    IF ( ABS(matrix(1,1)) .LE. pivot_tol ) THEN
       info = 1
       RETURN
    END IF

    factor = matrix(2,1) / matrix(1,1)
    matrix(2,1) = factor
    matrix(2,2) = matrix(2,2) - factor * matrix(1,2)
    rhs(2) = rhs(2) - factor * rhs(1)

    IF ( ABS(matrix(2,2)) .LE. pivot_tol ) THEN
       info = 2
       RETURN
    END IF

    rhs(2) = rhs(2) / matrix(2,2)
    rhs(1) = ( rhs(1) - matrix(1,2) * rhs(2) ) / matrix(1,1)

    RETURN

  END SUBROUTINE solve_2x2_pivoted


  !******************************************************************************
  !> \brief Solve a 3x3 linear system with partial pivoting
  !
  !> The matrix and right-hand side are overwritten with the elimination
  !> factors and the solution, respectively. A nonzero info value identifies
  !> the first singular pivot.
  !******************************************************************************

  SUBROUTINE solve_3x3_pivoted( matrix, rhs, info )

    IMPLICIT NONE

    REAL(wp), INTENT(INOUT) :: matrix(3,3)
    REAL(wp), INTENT(INOUT) :: rhs(3)
    INTEGER, INTENT(OUT) :: info

    INTEGER :: j
    INTEGER :: pivot_row
    REAL(wp) :: factor
    REAL(wp) :: matrix_scale
    REAL(wp) :: pivot_abs
    REAL(wp) :: pivot_tol
    REAL(wp) :: swap_value

    info = 0
    matrix_scale = MAXVAL(ABS(matrix))

    IF ( matrix_scale .LE. TINY(1.0_wp) ) THEN
       info = 1
       RETURN
    END IF

    pivot_tol = pivot_tol_factor * EPSILON(1.0_wp) * matrix_scale

    ! First elimination column.
    pivot_row = 1
    pivot_abs = ABS(matrix(1,1))

    IF ( ABS(matrix(2,1)) .GT. pivot_abs ) THEN
       pivot_row = 2
       pivot_abs = ABS(matrix(2,1))
    END IF

    IF ( ABS(matrix(3,1)) .GT. pivot_abs ) THEN
       pivot_row = 3
       pivot_abs = ABS(matrix(3,1))
    END IF

    IF ( pivot_abs .LE. pivot_tol ) THEN
       info = 1
       RETURN
    END IF

    IF ( pivot_row .NE. 1 ) THEN

       DO j=1,3
          swap_value = matrix(1,j)
          matrix(1,j) = matrix(pivot_row,j)
          matrix(pivot_row,j) = swap_value
       END DO

       swap_value = rhs(1)
       rhs(1) = rhs(pivot_row)
       rhs(pivot_row) = swap_value

    END IF

    factor = matrix(2,1) / matrix(1,1)
    matrix(2,1) = factor
    matrix(2,2) = matrix(2,2) - factor * matrix(1,2)
    matrix(2,3) = matrix(2,3) - factor * matrix(1,3)
    rhs(2) = rhs(2) - factor * rhs(1)

    factor = matrix(3,1) / matrix(1,1)
    matrix(3,1) = factor
    matrix(3,2) = matrix(3,2) - factor * matrix(1,2)
    matrix(3,3) = matrix(3,3) - factor * matrix(1,3)
    rhs(3) = rhs(3) - factor * rhs(1)

    ! Second elimination column.
    pivot_row = 2
    pivot_abs = ABS(matrix(2,2))

    IF ( ABS(matrix(3,2)) .GT. pivot_abs ) THEN
       pivot_row = 3
       pivot_abs = ABS(matrix(3,2))
    END IF

    IF ( pivot_abs .LE. pivot_tol ) THEN
       info = 2
       RETURN
    END IF

    IF ( pivot_row .NE. 2 ) THEN

       DO j=1,3
          swap_value = matrix(2,j)
          matrix(2,j) = matrix(pivot_row,j)
          matrix(pivot_row,j) = swap_value
       END DO

       swap_value = rhs(2)
       rhs(2) = rhs(pivot_row)
       rhs(pivot_row) = swap_value

    END IF

    factor = matrix(3,2) / matrix(2,2)
    matrix(3,2) = factor
    matrix(3,3) = matrix(3,3) - factor * matrix(2,3)
    rhs(3) = rhs(3) - factor * rhs(2)

    IF ( ABS(matrix(3,3)) .LE. pivot_tol ) THEN
       info = 3
       RETURN
    END IF

    ! Back substitution.
    rhs(3) = rhs(3) / matrix(3,3)
    rhs(2) = ( rhs(2) - matrix(2,3) * rhs(3) ) / matrix(2,2)
    rhs(1) = ( rhs(1) - matrix(1,2) * rhs(2)                         &
         - matrix(1,3) * rhs(3) ) / matrix(1,1)

    RETURN

  END SUBROUTINE solve_3x3_pivoted

  !******************************************************************************
  !> \brief Search the descent stepsize
  !
  !> This subroutine search for the lenght of the descent step in order to have
  !> a decrease in the nonlinear function.
  !> \param[in]     qj_rel_NR_old
  !> \param[in]     qj_org
  !> \param[in]     qj_old
  !> \param[in]     scal_f_old
  !> \param[in]     grad_f
  !> \param[in,out] desc_dir
  !> \param[in]     coeff_f
  !> \param[out]    qj_rel
  !> \param[out]    scal_f
  !> \param[out]    right_term
  !> \param[in]     stpmax
  !> \param[out]    check
  !> \param[in]     RJ_not_impl
  !> @author
  !> Mattia de' Michieli Vitturi
  !> \date 2019/12/16
  !******************************************************************************

  SUBROUTINE lnsrch( qj_rel_NR_old , qj_org , qj_old , scal_f_old , grad_f ,    &
       desc_dir , coeff_f , qj_rel , scal_f , right_term , stpmax , check ,     &
       dt_step , a_diag , Rj_not_impl , Bprimej_x , Bprimej_y, Zij, fric_val )

    IMPLICIT NONE

    !> Initial point
    REAL(wp), DIMENSION(:), INTENT(IN) :: qj_rel_NR_old

    !> Initial point
    REAL(wp), DIMENSION(:), INTENT(IN) :: qj_org

    !> Initial point
    REAL(wp), DIMENSION(:), INTENT(IN) :: qj_old

    !> Gradient at xold
    REAL(wp), DIMENSION(:), INTENT(IN) :: grad_f

    !> Value of the function at xold
    REAL(wp), INTENT(IN) :: scal_f_old

    !> Descent direction (usually Newton direction)
    REAL(wp), DIMENSION(:), INTENT(INOUT) :: desc_dir

    REAL(wp), INTENT(IN) :: stpmax

    !> Coefficients to rescale the nonlinear function
    REAL(wp), DIMENSION(:), INTENT(IN) :: coeff_f

    !> Updated solution
    REAL(wp), DIMENSION(:), INTENT(OUT) :: qj_rel

    !> Value of the scalar function at x
    REAL(wp), INTENT(OUT) :: scal_f

    !> Residual at the initial point, updated with the accepted step
    REAL(wp), INTENT(INOUT) :: right_term(n_eqns)

    !> Output quantity check is false on a normal exit
    LOGICAL, INTENT(OUT) :: check

    REAL(wp), INTENT(IN) :: dt_step
    REAL(wp), INTENT(IN) :: a_diag
    REAL(wp), INTENT(IN) :: Rj_not_impl(n_eqns)

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y

    ! vars for stochastic variable
    REAL(wp), INTENT(IN):: Zij ! value stochastic process
    REAL(wp), INTENT(INOUT) :: fric_val ! to save the value of the friction
    REAL(wp), PARAMETER :: TOLX=epsilon(qj_rel)

    INTEGER, DIMENSION(1) :: ndum
    REAL(wp) :: ALF , a,alam,alam2,alamin,b,disc
    REAL(wp) :: scal_f2
    REAL(wp) :: desc_dir_abs
    REAL(wp) :: rhs1 , rhs2 , slope, tmplam

    REAL(wp) :: qj(n_vars)
    REAL(wp) :: right_term_old(n_eqns)
    REAL(wp) :: fric_val_old

    ALF = 1.0e-4_wp

    IF ( size(grad_f) == size(desc_dir) .AND. size(grad_f) == size(qj_rel)      &
         .AND. size(qj_rel) == size(qj_rel_NR_old) ) THEN

       ndum = size(grad_f)

    ELSE

       WRITE(*,*) 'nrerror: an assert_eq failed with this tag:', 'lnsrch'
       STOP 'program terminated by assert_eq4'

    END IF

    check = .FALSE.
    right_term_old = right_term
    fric_val_old = fric_val

    desc_dir_abs = NORM2(desc_dir)

    IF ( desc_dir_abs > stpmax ) desc_dir(:) = desc_dir(:) * stpmax/desc_dir_abs

    slope = DOT_PRODUCT(grad_f,desc_dir)

    IF ( slope .GE. 0.0_wp ) THEN
       qj_rel = qj_rel_NR_old
       scal_f = scal_f_old
       right_term = right_term_old
       fric_val = fric_val_old
       check = .TRUE.
       RETURN
    END IF

    alamin = TOLX / MAXVAL(ABS( desc_dir(:))/MAX( ABS(qj_rel_NR_old(:)),1.0_wp ))

    IF ( alamin .EQ. 0.0_wp ) THEN

       qj_rel(:) = qj_rel_NR_old(:)
       scal_f = scal_f_old
       right_term = right_term_old
       fric_val = fric_val_old
       check = .TRUE.

       RETURN

    END IF

    alam = 1.0_wp
    alam2 = alam
    scal_f2 = scal_f_old

    optimal_step_search: DO

       IF ( verbose_level .GE. 4 ) THEN

          WRITE(*,*) 'alam',alam

       END IF

       qj_rel = qj_rel_NR_old + alam * desc_dir

       qj = qj_rel * qj_org

       CALL eval_f( qj , qj_old , dt_step , a_diag , coeff_f ,             &
            Rj_not_impl , Bprimej_x , Bprimej_y, right_term , scal_f, Zij,      &
            fric_val )

       IF ( verbose_level .GE. 4 ) THEN

          WRITE(*,*) 'lnsrch: effe_old,effe',scal_f_old,scal_f
          READ(*,*)

       END IF

       IF ( scal_f .LE. scal_f_old + ALF * alam * slope ) THEN
          ! Sufficient decrease according to the Armijo condition.

          IF ( verbose_level .GE. 4 ) THEN

             WRITE(*,*) 'sufficient function decrease'

          END IF

          EXIT optimal_step_search

       ELSE IF ( alam < alamin ) THEN
          ! convergence on Delta_x

          IF ( verbose_level .GE. 4 ) THEN

             WRITE(*,*) ' convergence on Delta_x',alam,alamin

          END IF

          qj_rel(:) = qj_rel_NR_old(:)
          scal_f = scal_f_old
          right_term = right_term_old
          fric_val = fric_val_old
          check = .TRUE.

          EXIT optimal_step_search

       ELSE

          IF ( alam .EQ. 1.0_wp ) THEN

             tmplam = - slope / ( 2.0_wp * ( scal_f - scal_f_old - slope ) )

          ELSE

             rhs1 = scal_f - scal_f_old - alam*slope
             rhs2 = scal_f2 - scal_f_old - alam2*slope

             a = ( rhs1/alam**2 - rhs2/alam2**2 ) / ( alam - alam2 )
             b = ( -alam2*rhs1/alam**2 + alam*rhs2/alam2**2 ) / ( alam - alam2 )

             IF ( a .EQ. 0.0_wp ) THEN

                tmplam = - slope / ( 2.0_wp * b )

             ELSE

                disc = b*b - 3.0_wp*a*slope

                IF ( disc .LT. 0.0_wp ) THEN

                   tmplam = 0.5_wp * alam

                ELSE IF ( b .LE. 0.0_wp ) THEN

                   tmplam = ( - b + SQRT(disc) ) / ( 3.0_wp * a )

                ELSE

                   tmplam = - slope / ( b + SQRT(disc) )

                ENDIF

             END IF

             IF ( tmplam .GT. 0.5_wp * alam ) tmplam = 0.5_wp * alam

          END IF

       END IF

       alam2 = alam
       scal_f2 = scal_f
       ! Keep the interpolated step.  A lower bound of 0.5 made every
       ! accepted update an exact halving after the upper bound above.
       alam = MAX( tmplam , 0.1_wp * alam )

    END DO optimal_step_search

    RETURN

  END SUBROUTINE lnsrch

  !******************************************************************************
  !> \brief Evaluate the nonlinear system
  !
  !> This subroutine evaluate the value of the nonlinear system in the state
  !> defined by the variables qj.
  !> \param[in]    qj          conservative variables
  !> \param[in]    qj_old      conservative variables at the old time step
  !> \param[in]    a_diag      implicit coefficient for the non-hyperbolic term
  !> \param[in]    coeff_f     coefficient to rescale the nonlinear functions
  !> \param[in]    Rj_not_impl explicit terms
  !> \param[out]   f_nl        values of the nonlinear functions
  !> \param[out]   scal_f      value of the scalar function f=0.5*<F,F>
  !> \date 2019/12/16
  !> @author
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_f( qj , qj_old , dt_step, a_diag , coeff_f , Rj_not_impl ,    &
       Bprimej_x ,                                                              &
       Bprimej_y , f_nl , scal_f, Zij, fric_val )

    USE constitutive_2d, ONLY : eval_implicit_terms

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: qj(n_vars)
    REAL(wp), INTENT(IN) :: qj_old(n_vars)
    REAL(wp), INTENT(IN) :: dt_step
    REAL(wp), INTENT(IN) :: a_diag
    REAL(wp), INTENT(IN) :: coeff_f(n_eqns)
    REAL(wp), INTENT(IN) :: Rj_not_impl(n_eqns)

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y


    REAL(wp), INTENT(OUT) :: f_nl(n_eqns)
    REAL(wp), INTENT(OUT) :: scal_f

    REAL(wp), INTENT(IN):: Zij ! value stochastic process
    REAL(wp), INTENT(OUT) :: fric_val ! to save the value of the friction

    REAL(wp) :: nh_term_impl(n_eqns)
    REAL(wp) :: Rj(n_eqns)

    CALL eval_implicit_terms( Bprimej_x , Bprimej_y, Zij, fric_val, r_qj = qj , &
         r_nh_term_impl=nh_term_impl )

    Rj = Rj_not_impl - a_diag * nh_term_impl

    f_nl = qj - qj_old + dt_step * Rj

    f_nl = coeff_f * f_nl

    scal_f = 0.5_wp * DOT_PRODUCT( f_nl , f_nl )

    RETURN

  END SUBROUTINE eval_f

  !******************************************************************************
  !> \brief Evaluate the jacobian
  !
  !> This subroutine evaluate the jacobian of the non-linear system
  !> with respect to the conservative variables.
  !
  !> \param[in]    qj_rel        relative variation (qj=qj_rel*qj_org)
  !> \param[in]    qj_org        conservative variables at the old time step
  !> \param[in]    coeff_f       coefficient to rescale the nonlinear functions
  !> \param[out]   left_matrix   matrix from the linearization of the system
  !
  !> \date 07/10/2016
  !> @author
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_jacobian( qj_rel , qj_org , dt_step, a_diag , coeff_f,   &
       Bprimej_x , Bprimej_y , left_matrix, Zij, fric_val)

    USE constitutive_2d, ONLY : eval_implicit_terms

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: qj_rel(n_vars)
    REAL(wp), INTENT(IN) :: qj_org(n_vars)
    REAL(wp), INTENT(IN) :: dt_step
    REAL(wp), INTENT(IN) :: a_diag
    REAL(wp), INTENT(IN) :: coeff_f(n_eqns)

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y

    REAL(wp), INTENT(OUT) :: left_matrix(n_eqns,n_vars)

    REAL(wp), INTENT(IN):: Zij ! value stochastic process
    REAL(wp), INTENT(OUT) :: fric_val ! to save the value of the friction

    REAL(wp) :: Jacob_relax(n_eqns,n_vars)
    COMPLEX(wp) :: nh_terms_cmplx_impl(n_eqns)
    COMPLEX(wp) :: qj_cmplx(n_vars) , qj_rel_cmplx(n_vars)
    COMPLEX(wp) :: qj_rel_cmplx_init(n_vars)

    INTEGER :: i

    ! initialize the matrix of the linearized system and the Jacobian

    left_matrix(1:n_eqns,1:n_vars) = 0.0_wp
    Jacob_relax(1:n_eqns,1:n_vars) = 0.0_wp

    ! evaluate the jacobian of the non-hyperbolic terms

    DO i=1,n_vars

       qj_rel_cmplx_init(i) = CMPLX(qj_rel(i),0.0_wp,wp)

    END DO

    DO i=1,n_vars

       left_matrix(i,i) = coeff_f(i) * qj_org(i)

       IF ( implicit_flag(i) ) THEN

          qj_rel_cmplx(1:n_vars) = qj_rel_cmplx_init(1:n_vars)
          qj_rel_cmplx(i) = CMPLX(qj_rel(i), h,wp)

          qj_cmplx = qj_rel_cmplx * qj_org

          CALL eval_implicit_terms( Bprimej_x , Bprimej_y, Zij, fric_val,       &
               c_qj = qj_cmplx , c_nh_term_impl = nh_terms_cmplx_impl )

          Jacob_relax(1:n_eqns,i) = coeff_f(1:n_eqns) *                        &
               AIMAG(nh_terms_cmplx_impl) * one_by_h

          left_matrix(1:n_eqns,i) = left_matrix(1:n_eqns,i) - dt_step * a_diag       &
               * Jacob_relax(1:n_eqns,i)

       END IF

    END DO

    RETURN

  END SUBROUTINE eval_jacobian


END MODULE nonlinear_solver_2d
