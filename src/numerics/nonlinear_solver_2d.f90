!********************************************************************************
!> \brief Cell-local backward-Euler source group nonlinear solver
!>
!> Newton iteration, line search, residual/Jacobian evaluation and small dense
!> linear solves.  This module is intentionally independent of the spatial
!> domain and operates on one cell state at a time.
!********************************************************************************

MODULE nonlinear_solver_2d

  USE diagnostics_2d, ONLY : debug_pause

  USE lapack_interfaces_2d, ONLY : solve_dense_system

  USE parameters_2d, ONLY : wp, n_eqns, n_vars, verbose_level


  USE source_registry_2d, ONLY : source_group_type
  USE source_catalog_2d, ONLY : SOURCE_METHOD_EXPLICIT, SOURCE_METHOD_NEWTON, &
       SOURCE_METHOD_ANALYTIC
  USE source_terms_2d, ONLY : analytic_term_context_type, &
       eval_source_analytic, eval_source_semi_implicit
  USE rheology_2d, ONLY : limit_rheology_momentum_reversal

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: initialize_nonlinear_solver
  PUBLIC :: finalize_nonlinear_solver
  PUBLIC :: solve_backward_euler_group

  !> Complex-step perturbation and its reciprocal.
  REAL(wp) :: h, one_by_h

  !> Relative threshold used to reject numerically singular small systems.
  REAL(wp), PARAMETER :: pivot_tol_factor = 64.0_wp

CONTAINS

  !> Solve one full backward-Euler Lie group without RK coefficients.
  SUBROUTINE solve_backward_euler_group(source_group, &
       q_cell,dt,Bprime_x,Bprime_y,stochastic_value, &
       iterations_used,converged,linear_info,line_search_failed)
    TYPE(source_group_type), INTENT(IN) :: source_group
    REAL(wp), INTENT(INOUT) :: q_cell(n_vars)
    REAL(wp), INTENT(IN) :: dt,Bprime_x,Bprime_y,stochastic_value
    INTEGER, INTENT(OUT) :: iterations_used,linear_info
    LOGICAL, INTENT(OUT) :: converged,line_search_failed
    REAL(wp) :: q_in(n_vars)
    q_in=q_cell
    CALL solve_implicit_groups([source_group],q_cell,q_in, &
         dt,Bprime_x,Bprime_y,0.0_wp,0.0_wp,0.0_wp,0.0_wp, &
         stochastic_value,iterations_used,converged,linear_info,line_search_failed)
  END SUBROUTINE solve_backward_euler_group

  !> \brief Set the complex-step perturbation and its reciprocal.
  !>
  !> \note Writes h and one_by_h, used to differentiate the implicit source with a complex
  !>       perturbation.

  SUBROUTINE initialize_nonlinear_solver

    h = n_vars * EPSILON(1.0_wp)
    one_by_h = 1.0_wp / h

  END SUBROUTINE initialize_nonlinear_solver

  !> \brief Provide the matching lifecycle hook for the cell-local nonlinear solver.
  !>
  !> \note Currently a no-op: the nonlinear solver owns scalar constants, not allocated arrays.

  SUBROUTINE finalize_nonlinear_solver

  END SUBROUTINE finalize_nonlinear_solver

  !> \brief Apply ordered cell-local source groups using Newton or analytic flows.
  !>
  !> \param[in] lie_groups Validated ordered Lie source groups; explicit groups are skipped.
  !> \param[in,out] qj Local conservative state vector; mass and momenta are per unit horizontal
  !>                   area.
  !> \param[in] q_group_input Conservative input state for the backward-Euler group residual.
  !> \param[in] dt_step Full Lie group timestep [s].
  !> \param[in] Bprimej_x Filtered cell bed derivative dB/dx [dimensionless].
  !> \param[in] Bprimej_y Filtered cell bed derivative dB/dy [dimensionless].
  !> \param[in] Bsecondj_xx Filtered bed curvature in xx [m^-1].
  !> \param[in] Bsecondj_xy Filtered bed curvature in xy [m^-1].
  !> \param[in] Bsecondj_yy Filtered bed curvature in yy [m^-1].
  !> \param[in] grav_coeff Effective gravitational acceleration [m s^-2].
  !> \param[in] Zij Local effective stochastic fluctuation supplied to the rheology.
  !> \param[in] apply_semi_implicit Apply the semi-implicit momentum contribution before the sweep.
  !> \param[out] iterations_used Number of Newton iterations attempted for this source-group cell solve.
  !> \param[out] converged True when every Newton group residual meets its bound.
  !> \param[out] linear_info Linear-solver status: zero for success, nonzero for a failed solve.
  !> \param[out] line_search_failed True when backtracking fails to accept a step.
  !******************************************************************************

  SUBROUTINE solve_implicit_groups(lie_groups, &
       qj, q_group_input, dt_step,                                                  &
Bprimej_x, Bprimej_y, Bsecondj_xx, Bsecondj_xy, &
       Bsecondj_yy, grav_coeff, Zij,                                         &
       iterations_used, converged, linear_info, line_search_failed,          &
       apply_semi_implicit )

    USE parameters_2d, ONLY : tol_rel, tol_abs
    USE state_conversion_2d, ONLY : qc_to_qp

    IMPLICIT NONE

    TYPE(source_group_type), INTENT(IN) :: lie_groups(:)
    REAL(wp), INTENT(INOUT) :: qj(n_vars)
    REAL(wp), INTENT(IN) :: q_group_input(n_vars)
    REAL(wp), INTENT(IN) :: dt_step
    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y
    REAL(wp), INTENT(IN) :: Bsecondj_xx, Bsecondj_xy, Bsecondj_yy
    REAL(wp), INTENT(IN) :: grav_coeff
    REAL(wp), INTENT(IN):: Zij ! value stochastic process
    LOGICAL, INTENT(IN), OPTIONAL :: apply_semi_implicit
    INTEGER, INTENT(OUT) :: iterations_used
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: linear_info
    LOGICAL, INTENT(OUT) :: line_search_failed

    REAL(wp) :: qj_init(n_vars), q_known(n_vars), q_group_before(n_vars)
    REAL(wp) :: coeff_f(n_eqns)
    REAL(wp) :: tol_abs_local, tol_rel_local
    REAL(wp) :: analytic_rate(n_eqns)
    TYPE(analytic_term_context_type) :: analytic_context
    REAL(wp) :: qp_current(n_vars+2)
    REAL(wp) :: dynamic_pressure, si_term(n_eqns)
    INTEGER :: group, group_size, equation, i
    INTEGER :: analytic_status
    INTEGER :: group_equations(n_eqns)
    LOGICAL :: group_has_momentum, apply_si_local

    iterations_used = 0
    converged = .FALSE.
    linear_info = 0
    line_search_failed = .FALSE.

    IF (n_eqns /= n_vars .OR. SIZE(lie_groups)==0) THEN
       linear_info = -1
       RETURN
    END IF

    coeff_f = 1.0_wp

    tol_abs_local = tol_abs
    tol_rel_local = tol_rel

    qj_init = qj
    q_known = q_group_input
    qj = q_known
    apply_si_local = .FALSE.
    IF (PRESENT(apply_semi_implicit)) apply_si_local = apply_semi_implicit
    analytic_context%Bprime_x=Bprimej_x
    analytic_context%Bprime_y=Bprimej_y
    analytic_context%Bsecond_xx=Bsecondj_xx
    analytic_context%Bsecond_xy=Bsecondj_xy
    analytic_context%Bsecond_yy=Bsecondj_yy
    analytic_context%grav_coeff=grav_coeff
    analytic_context%Zij=Zij

    ! All non-group equations remain at the incoming group state.

    DO group = 1, SIZE(lie_groups)
       IF (lie_groups(group)%method==SOURCE_METHOD_EXPLICIT) CYCLE
       IF (.NOT.ALLOCATED(lie_groups(group)%equation_ids)) THEN
          linear_info=-1
          qj=qj_init
          RETURN
       END IF
       group_size=SIZE(lie_groups(group)%equation_ids)
       IF (group_size<1 .OR. group_size>SIZE(group_equations)) THEN
          linear_info=-1
          qj=qj_init
          RETURN
       END IF
       group_equations(1:group_size)=lie_groups(group)%equation_ids

       IF (lie_groups(group)%method==SOURCE_METHOD_ANALYTIC) THEN
          ! Analytic operators act on the state produced by preceding groups.
          group_has_momentum=ANY(group_equations(1:group_size)==2) .OR. &
               ANY(group_equations(1:group_size)==3)
          IF (group_has_momentum .AND. &
               .NOT.(ANY(group_equations(1:group_size)==2) .AND. &
               ANY(group_equations(1:group_size)==3))) THEN
             linear_info=-6
             qj=qj_init
             RETURN
          END IF
          CALL eval_source_analytic(lie_groups(group),qj, &
               analytic_context,dt_step,analytic_rate,             &
               status=analytic_status)
          IF (analytic_status/=0) THEN
             linear_info=analytic_status
             qj=qj_init
             RETURN
          END IF
          qj=qj+dt_step*analytic_rate(1:n_vars)

       ELSEIF (lie_groups(group)%method==SOURCE_METHOD_NEWTON) THEN
          ! A group solves against the state left by preceding groups.
          q_group_before=qj
             si_term=0.0_wp
             group_has_momentum = ANY(group_equations(1:group_size)==2) .OR. &
                  ANY(group_equations(1:group_size)==3)
             IF (group_has_momentum .AND. apply_si_local) THEN
                CALL qc_to_qp(qj, qp_current, dynamic_pressure)
                CALL eval_source_semi_implicit(Bprimej_x, Bprimej_y,          &
                     Bsecondj_xx,Bsecondj_xy,Bsecondj_yy,grav_coeff,        &
                     qp_current,si_term,Zij,group=lie_groups(group))
             END IF
             DO i=1,group_size
                equation=group_equations(i)
                qj(equation)=q_group_before(equation)+dt_step*si_term(equation)
             END DO
             IF (group_has_momentum .AND. apply_si_local) THEN
                CALL limit_rheology_momentum_reversal(q_group_before(2:3),qj(2:3))
                DO i=1,group_size
                   equation=group_equations(i)
                   si_term(equation)=                                        &
                        (qj(equation)-q_group_before(equation))/(dt_step)
                END DO
             END IF
             ! Semi-implicit information seeds the Newton initial guess only.
             ! The backward-Euler residual contains the full implicit source once.
             CALL solve_newton_group(lie_groups(group),&
                  group_equations(1:group_size),qj,q_group_before,dt_step,&
                  Bprimej_x,Bprimej_y,Zij,coeff_f,      &
                  tol_abs_local,tol_rel_local,iterations_used,converged,   &
                  linear_info,line_search_failed)
             IF (.NOT.converged) THEN
                qj=qj_init
                RETURN
             END IF
       ELSE
          linear_info = -5
          qj = qj_init
          RETURN
       END IF
    END DO

    converged = .TRUE.

  END SUBROUTINE solve_implicit_groups

  !> \brief Solve one backward-Euler Newton group while holding other equations fixed.
  !> \param[in] group_equations Indices of the conservative equations in this group.
  !> \param[in,out] qj Local conservative state; only group equations are updated.
  !> \param[in] q_group_input Conservative reference state at the start of the step.
  !> \param[in] dt_step Full Lie group time interval [s].
  !> \param[in] Bprimej_x Filtered bed slope in x [dimensionless].
  !> \param[in] Bprimej_y Filtered bed slope in y [dimensionless].
  !> \param[in] Zij Local stochastic-process value supplied to the rheology.
  !> \param[in] coeff_f Per-equation residual scaling factors.
  !> \param[in] source_group Validated Newton source group being solved.
  !> \param[in] source_solver Validated solver descriptor selecting rates and target equations.
  !> \param[in] tol_abs_local Absolute residual tolerance.
  !> \param[in] tol_rel_local Relative residual tolerance.
  !> \param[in,out] iterations_used Cumulative Newton iterations attempted.
  !> \param[out] converged True when the group residual meets its tolerance.
  !> \param[out] linear_info Linear solver status; zero denotes success.
  !> \param[out] line_search_failed True if backtracking cannot accept a step.
  SUBROUTINE solve_newton_group(source_group,group_equations,   &
       qj,q_group_input,dt_step,Bprimej_x,Bprimej_y,Zij,coeff_f, &
       tol_abs_local,                                                         &
       tol_rel_local, iterations_used, converged, linear_info,               &
       line_search_failed)

    USE parameters_2d, ONLY : max_nl_iter

    IMPLICIT NONE

    TYPE(source_group_type), INTENT(IN) :: source_group
    INTEGER, INTENT(IN) :: group_equations(:)
    REAL(wp), INTENT(INOUT) :: qj(n_vars)
    REAL(wp), INTENT(IN) :: q_group_input(n_vars), dt_step
    REAL(wp), INTENT(IN) :: Bprimej_x, Bprimej_y, Zij
    REAL(wp), INTENT(IN) :: coeff_f(n_eqns), tol_abs_local, tol_rel_local
    INTEGER, INTENT(INOUT) :: iterations_used
    LOGICAL, INTENT(OUT) :: converged, line_search_failed
    INTEGER, INTENT(OUT) :: linear_info

    REAL(wp) :: qj_org(n_vars), qj_rel(n_vars), qj_rel_old(n_vars)
    REAL(wp) :: qj_trial(n_vars)
    REAL(wp) :: left_matrix(n_eqns,n_vars), right_term(n_eqns)
    REAL(wp) :: residual_ref(n_eqns), residual_tol(n_eqns)
    REAL(wp) :: desc_dir(n_vars), grad_f(n_vars), scal_f, scal_f_old
    REAL(wp) :: stpmax
    REAL(wp) :: group_matrix(SIZE(group_equations),SIZE(group_equations))
    REAL(wp) :: group_rhs(SIZE(group_equations))
    INTEGER :: group_pivot(SIZE(group_equations))
    INTEGER :: nl_iter, i, j, ok
    LOGICAL :: check
    REAL(wp), PARAMETER :: STPMX = 100.0_wp

    converged = .FALSE.
    line_search_failed = .FALSE.
    linear_info = 0
    grad_f = 0.0_wp
    qj_org = MAX(ABS(qj),1.0E-3_wp)
    qj_rel = qj / qj_org

    DO nl_iter = 1, max_nl_iter
       iterations_used = iterations_used + 1
       CALL eval_f(qj, q_group_input, dt_step, coeff_f,       &
            Bprimej_x, Bprimej_y, right_term, scal_f, Zij, group_equations, &
            source_group)

       IF (nl_iter .EQ. 1) residual_ref = ABS(right_term)
       residual_tol = tol_abs_local + tol_rel_local*residual_ref
       IF (ALL(ABS(right_term) .LE. residual_tol)) THEN
          converged = .TRUE.
          EXIT
       END IF

       CALL eval_jacobian(qj_rel, qj_org, dt_step, coeff_f,         &
            Bprimej_x, Bprimej_y, left_matrix, Zij, group_equations,       &
            source_group)
       IF (nl_iter .GT. 1) grad_f = MATMUL(right_term,left_matrix)

       DO i = 1, SIZE(group_equations)
          group_rhs(i) = -right_term(group_equations(i))
          DO j = 1, SIZE(group_equations)
             group_matrix(i,j) = left_matrix(group_equations(i),            &
                  group_equations(j))
          END DO
       END DO
       IF (SIZE(group_equations) .EQ. 2) THEN
          CALL solve_2x2_pivoted(group_matrix, group_rhs, ok)
       ELSEIF (SIZE(group_equations) .EQ. 3) THEN
          CALL solve_3x3_pivoted(group_matrix, group_rhs, ok)
       ELSE
          CALL solve_dense_system(group_matrix, group_rhs, group_pivot, ok)
       END IF
       IF (ok .NE. 0) THEN
          linear_info = ok
          RETURN
       END IF

       desc_dir = 0.0_wp
       DO i = 1, SIZE(group_equations)
          desc_dir(group_equations(i)) = group_rhs(i)
       END DO
       qj_rel_old = qj_rel
       scal_f_old = scal_f

       IF (nl_iter .GT. 1) THEN
          stpmax = STPMX * MAX(NORM2(qj_rel),REAL(SIZE(qj_rel),wp))
          CALL lnsrch(qj_rel_old, qj_org, qj, q_group_input, scal_f_old, grad_f,   &
               desc_dir, coeff_f, qj_rel, scal_f, right_term, stpmax,      &
               check, dt_step, Bprimej_x, Bprimej_y,  &
               Zij, group_equations,source_group)
          IF (check) THEN
             line_search_failed = .TRUE.
             RETURN
          END IF
       ELSE
          qj_rel = qj_rel_old + desc_dir
          qj_trial = qj
          DO i = 1, SIZE(group_equations)
             qj_trial(group_equations(i)) =                                &
                  qj_rel(group_equations(i))*qj_org(group_equations(i))
          END DO
          CALL eval_f(qj_trial, q_group_input, dt_step, coeff_f,           &
Bprimej_x, Bprimej_y, right_term, scal_f, Zij, &
               group_equations,source_group)
       END IF

       DO i = 1, SIZE(group_equations)
          qj(group_equations(i)) = qj_rel(group_equations(i))*             &
               qj_org(group_equations(i))
       END DO
       IF (ALL(ABS(right_term) .LE. residual_tol)) THEN
          converged = .TRUE.
          EXIT
       END IF
       IF (MAXVAL(ABS(qj_rel-qj_rel_old)/MAX(ABS(qj_rel),1.0_wp)) .LT.      &
            EPSILON(1.0_wp)) EXIT
    END DO

  END SUBROUTINE solve_newton_group

  !******************************************************************************
  !> \brief Solve a two-equation dense system with guarded partial pivoting.
  !
  !> The matrix and right-hand side are overwritten with the elimination
  !> factors and the solution, respectively. Rows are equilibrated before
  !> elimination so that dimensional differences do not affect pivot tests.
  !> A nonzero info value identifies the first numerically singular pivot.
  !>
  !> \param[in,out] matrix Coefficient matrix; overwritten by row scaling and elimination.
  !> \param[in,out] rhs Right-hand side on entry; solution vector on successful exit.
  !> \param[out] info Zero on success; positive for a singular pivot, negative for invalid
  !>                  arguments/shapes.
  !******************************************************************************

  SUBROUTINE solve_2x2_pivoted( matrix, rhs, info )

    IMPLICIT NONE

    REAL(wp), INTENT(INOUT) :: matrix(2,2)
    REAL(wp), INTENT(INOUT) :: rhs(2)
    INTEGER, INTENT(OUT) :: info

    INTEGER :: i
    REAL(wp) :: factor
    REAL(wp) :: pivot_tol
    REAL(wp) :: row_scale
    REAL(wp) :: swap_value

    info = 0

    DO i = 1,2
       row_scale = MAXVAL(ABS(matrix(i,:)))

       IF ( row_scale .LE. TINY(1.0_wp) ) THEN
          info = i
          RETURN
       END IF

       matrix(i,:) = matrix(i,:) / row_scale
       rhs(i) = rhs(i) / row_scale
    END DO

    pivot_tol = pivot_tol_factor * EPSILON(1.0_wp)

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
  !> \brief Solve a three-equation dense system with guarded partial pivoting.
  !
  !> The matrix and right-hand side are overwritten with the elimination
  !> factors and the solution, respectively. Rows are equilibrated before
  !> elimination so that dimensional differences do not affect pivot tests.
  !> A nonzero info value identifies the first numerically singular pivot.
  !>
  !> \param[in,out] matrix Coefficient matrix; overwritten by row scaling and elimination.
  !> \param[in,out] rhs Right-hand side on entry; solution vector on successful exit.
  !> \param[out] info Zero on success; positive for a singular pivot, negative for invalid
  !>                  arguments/shapes.
  !******************************************************************************

  SUBROUTINE solve_3x3_pivoted( matrix, rhs, info )

    IMPLICIT NONE

    REAL(wp), INTENT(INOUT) :: matrix(3,3)
    REAL(wp), INTENT(INOUT) :: rhs(3)
    INTEGER, INTENT(OUT) :: info

    INTEGER :: i
    INTEGER :: j
    INTEGER :: pivot_row
    REAL(wp) :: factor
    REAL(wp) :: pivot_abs
    REAL(wp) :: pivot_tol
    REAL(wp) :: row_scale
    REAL(wp) :: swap_value

    info = 0

    DO i = 1,3
       row_scale = MAXVAL(ABS(matrix(i,:)))

       IF ( row_scale .LE. TINY(1.0_wp) ) THEN
          info = i
          RETURN
       END IF

       matrix(i,:) = matrix(i,:) / row_scale
       rhs(i) = rhs(i) / row_scale
    END DO

    pivot_tol = pivot_tol_factor * EPSILON(1.0_wp)

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
  !> \brief Backtrack a Newton direction until the scaled residual merit function decreases.
  !
  !> This subroutine search for the lenght of the descent step in order to have
  !> a decrease in the nonlinear function.
  !> @author
  !> Mattia de' Michieli Vitturi
  !> \date 2019/12/16
  !>
  !> \param[in] qj_rel_NR_old Normalized state at the start of the current line search.
  !> \param[in] qj_org Per-variable scaling factors for normalized Newton variables.
  !> \param[in] q_group_input Conservative input state for the backward-Euler group residual.
  !> \param[in] qj_fixed Non-group conservative components held bitwise fixed during the line search.
  !> \param[in] scal_f_old Initial value of the half-squared scaled residual norm.
  !> \param[in] grad_f Gradient of the scalar residual merit function at the starting point.
  !> \param[in,out] desc_dir Newton descent direction, limited to the maximum permitted step length.
  !> \param[in] coeff_f Per-equation scaling factors applied to the group residual.
  !> \param[out] qj_rel Dimensionless normalized conservative state, qj=qj_rel*qj_org.
  !> \param[out] scal_f Half the squared Euclidean norm of the scaled nonlinear residual.
  !> \param[in,out] right_term Initial residual on entry; residual of the accepted line-search state
  !>                           on exit.
  !> \param[in] stpmax Maximum Euclidean length of the normalized trial step.
  !> \param[out] check True when the line search cannot find a usable decrease.
  !> \param[in] dt_step Full Lie group timestep [s].
  !> \param[in] Bprimej_x Filtered cell bed derivative dB/dx [dimensionless].
  !> \param[in] Bprimej_y Filtered cell bed derivative dB/dy [dimensionless].
  !> \param[in] Zij Local effective stochastic fluctuation supplied to the rheology.
  !> \param[in] active_equations Equation indices whose residual rows are active.
  !> \param[in] source_group Validated Newton source group for the selected rate evaluation.
    !******************************************************************************

  SUBROUTINE lnsrch( qj_rel_NR_old, qj_org, qj_fixed, q_group_input, scal_f_old,      &
       grad_f, desc_dir, coeff_f, qj_rel, scal_f, right_term, stpmax, check,  &
       dt_step  , Bprimej_x , Bprimej_y, Zij,           &
       active_equations,source_group )

    IMPLICIT NONE

    !> Initial point
    REAL(wp), DIMENSION(:), INTENT(IN) :: qj_rel_NR_old

    !> Initial point
    REAL(wp), DIMENSION(:), INTENT(IN) :: qj_org

    ! Non-group variables must remain bitwise unchanged during line search.
    REAL(wp), INTENT(IN) :: qj_fixed(n_vars)

    !> Initial point
    REAL(wp), DIMENSION(:), INTENT(IN) :: q_group_input

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
    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y

    ! vars for stochastic variable
    REAL(wp), INTENT(IN):: Zij ! value stochastic process
    INTEGER, INTENT(IN) :: active_equations(:)
    TYPE(source_group_type), INTENT(IN) :: source_group
    REAL(wp), PARAMETER :: TOLX=epsilon(qj_rel)

    INTEGER, DIMENSION(1) :: ndum
    INTEGER :: i, equation
    REAL(wp) :: ALF , a,alam,alam2,alamin,b,disc
    REAL(wp) :: scal_f2
    REAL(wp) :: desc_dir_abs
    REAL(wp) :: rhs1 , rhs2 , slope, tmplam

    REAL(wp) :: qj(n_vars)
    REAL(wp) :: right_term_old(n_eqns)

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

    desc_dir_abs = NORM2(desc_dir)

    IF ( desc_dir_abs > stpmax ) desc_dir(:) = desc_dir(:) * stpmax/desc_dir_abs

    slope = DOT_PRODUCT(grad_f,desc_dir)

    IF ( slope .GE. 0.0_wp ) THEN
       qj_rel = qj_rel_NR_old
       scal_f = scal_f_old
       right_term = right_term_old
       check = .TRUE.
       RETURN
    END IF

    alamin = TOLX / MAXVAL(ABS( desc_dir(:))/MAX( ABS(qj_rel_NR_old(:)),1.0_wp ))

    IF ( alamin .EQ. 0.0_wp ) THEN

       qj_rel(:) = qj_rel_NR_old(:)
       scal_f = scal_f_old
       right_term = right_term_old
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

       qj = qj_fixed
       DO i = 1, SIZE(active_equations)
          equation = active_equations(i)
          qj(equation) = qj_rel(equation)*qj_org(equation)
       END DO

       CALL eval_f( qj , q_group_input , dt_step  , coeff_f ,             &
Bprimej_x , Bprimej_y, right_term , scal_f, Zij,  &
            active_equations,source_group )

       IF ( verbose_level .GE. 4 ) THEN

          WRITE(*,*) 'lnsrch: effe_old,effe',scal_f_old,scal_f
          CALL debug_pause('nonlinear line-search iteration')

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
  !> \brief Evaluate the scaled backward-Euler group residual and its scalar merit function.
  !
  !> This subroutine evaluate the value of the nonlinear system in the state
  !> defined by the variables qj.
  !> \date 2019/12/16
  !> @author
  !> Mattia de' Michieli Vitturi
  !>
  !> \param[in] qj Local conservative state vector; mass and momenta are per unit horizontal area.
  !> \param[in] q_group_input Conservative input state for the backward-Euler group residual.
  !> \param[in] dt_step Full Lie group timestep [s].
  !> \param[in] coeff_f Per-equation scaling factors applied to the group residual.
  !> \param[in] Bprimej_x Filtered cell bed derivative dB/dx [dimensionless].
  !> \param[in] Bprimej_y Filtered cell bed derivative dB/dy [dimensionless].
  !> \param[out] f_nl Scaled nonlinear group residual, one entry per equation.
  !> \param[out] scal_f Half the squared Euclidean norm of the scaled nonlinear residual.
  !> \param[in] Zij Local effective stochastic fluctuation supplied to the rheology.
  !> \param[in] active_equations Equation indices whose residual rows and Jacobian columns are active.
  !> \param[in] source_group Validated Newton source group for the selected rate evaluation.
    !******************************************************************************

  SUBROUTINE eval_f( qj , q_group_input , dt_step , coeff_f ,    &
       Bprimej_x ,                                                              &
       Bprimej_y , f_nl , scal_f, Zij, active_equations,source_group )

    USE source_terms_2d, ONLY : eval_source_implicit

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: qj(n_vars)
    REAL(wp), INTENT(IN) :: q_group_input(n_vars)
    REAL(wp), INTENT(IN) :: dt_step
    REAL(wp), INTENT(IN) :: coeff_f(n_eqns)
    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y


    REAL(wp), INTENT(OUT) :: f_nl(n_eqns)
    REAL(wp), INTENT(OUT) :: scal_f

    REAL(wp), INTENT(IN):: Zij ! value stochastic process
    INTEGER, INTENT(IN) :: active_equations(:)
    TYPE(source_group_type), INTENT(IN) :: source_group

    REAL(wp) :: nh_term_impl(n_eqns)
    REAL(wp) :: Rj(n_eqns)
    ! Keep the conservative difference separate from the forcing even under
    ! -Ofast reassociation. Otherwise q + (dt*R-q_old) can lose a small forcing
    ! before subtracting equal large states and falsely report zero residual.
    REAL(wp), VOLATILE :: delta_q(n_vars)

    CALL eval_source_implicit( Bprimej_x , Bprimej_y, Zij, r_qj = qj ,           &
         r_nh_term_impl=nh_term_impl,group=source_group )

    Rj = -nh_term_impl

    delta_q = qj - q_group_input
    f_nl = delta_q + dt_step * Rj
    f_nl = coeff_f * f_nl
    WHERE (.NOT. equation_is_active(active_equations))
       f_nl = 0.0_wp
    END WHERE

    scal_f = 0.5_wp * DOT_PRODUCT( f_nl , f_nl )

    RETURN

  END SUBROUTINE eval_f

  !******************************************************************************
  !> \brief Build the scaled backward-Euler group Jacobian using complex-step differentiation.
  !
  !> This subroutine evaluate the jacobian of the non-linear system
  !> with respect to the conservative variables.
  !
  !
  !> \date 07/10/2016
  !> @author
  !> Mattia de' Michieli Vitturi
  !>
  !> \param[in] qj_rel Dimensionless normalized conservative state, qj=qj_rel*qj_org.
  !> \param[in] qj_org Per-variable scaling factors for normalized Newton variables.
  !> \param[in] dt_step Full Lie group timestep [s].
  !> \param[in] coeff_f Per-equation scaling factors applied to the group residual.
  !> \param[in] Bprimej_x Filtered cell bed derivative dB/dx [dimensionless].
  !> \param[in] Bprimej_y Filtered cell bed derivative dB/dy [dimensionless].
  !> \param[out] left_matrix Jacobian of the scaled residual with respect to the normalized state.
  !> \param[in] Zij Local effective stochastic fluctuation supplied to the rheology.
  !> \param[in] active_equations Equation indices whose Jacobian columns are differentiated.
  !> \param[in] source_group Validated Newton source group for the selected rate evaluation.
    !> \details Uses complex-step perturbations of normalized conservative variables. The resulting
  !>          imaginary derivatives form the active Jacobian columns without finite-difference
  !>          subtraction.
  !******************************************************************************

  SUBROUTINE eval_jacobian( qj_rel, qj_org, dt_step, coeff_f,        &
       Bprimej_x, Bprimej_y, left_matrix, Zij, active_equations,             &
       source_group )

    USE source_terms_2d, ONLY : eval_source_implicit

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: qj_rel(n_vars)
    REAL(wp), INTENT(IN) :: qj_org(n_vars)
    REAL(wp), INTENT(IN) :: dt_step
    REAL(wp), INTENT(IN) :: coeff_f(n_eqns)

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y

    REAL(wp), INTENT(OUT) :: left_matrix(n_eqns,n_vars)

    REAL(wp), INTENT(IN):: Zij ! value stochastic process
    INTEGER, INTENT(IN) :: active_equations(:)
    TYPE(source_group_type), INTENT(IN) :: source_group

    REAL(wp) :: Jacob_relax(n_eqns,n_vars)
    COMPLEX(wp) :: nh_terms_cmplx_impl(n_eqns)
    COMPLEX(wp) :: qj_cmplx(n_vars) , qj_rel_cmplx(n_vars)
    COMPLEX(wp) :: qj_rel_cmplx_init(n_vars)

    INTEGER :: i, j, equation

    ! initialize the matrix of the linearized system and the Jacobian

    left_matrix(1:n_eqns,1:n_vars) = 0.0_wp
    Jacob_relax(1:n_eqns,1:n_vars) = 0.0_wp

    ! Differentiate only implicit columns. A tiny imaginary perturbation of
    ! one normalized conservative variable produces a derivative without
    ! subtracting nearly equal residuals; the typed closures must retain it.

    DO i=1,n_vars

       qj_rel_cmplx_init(i) = CMPLX(qj_rel(i),0.0_wp,wp)

    END DO

    DO i = 1, SIZE(active_equations)
       equation = active_equations(i)
       left_matrix(equation,equation) = coeff_f(equation) * qj_org(equation)

       qj_rel_cmplx = qj_rel_cmplx_init
       qj_rel_cmplx(equation) = CMPLX(qj_rel(equation), h, wp)
       qj_cmplx = qj_rel_cmplx * qj_org

       CALL eval_source_implicit(Bprimej_x, Bprimej_y, Zij,                   &
            c_qj=qj_cmplx, c_nh_term_impl=nh_terms_cmplx_impl,               &
            group=source_group)

       DO j = 1, SIZE(active_equations)
          Jacob_relax(active_equations(j),equation) =                         &
               coeff_f(active_equations(j)) *                                &
               AIMAG(nh_terms_cmplx_impl(active_equations(j))) * one_by_h
       END DO
       DO j = 1, n_eqns
          left_matrix(j,equation) = left_matrix(j,equation) -                 &
               dt_step*Jacob_relax(j,equation)
       END DO
    END DO

    RETURN

  END SUBROUTINE eval_jacobian

  !> \brief Return a mask selecting residual rows in one implicit group.
  !> \param[in] active_equations Equation indices in the group.
  !> \return Logical mask over all equations, true for active indices.
  PURE FUNCTION equation_is_active(active_equations) RESULT(mask)
    INTEGER, INTENT(IN) :: active_equations(:)
    LOGICAL :: mask(n_eqns)
    INTEGER :: i

    mask = .FALSE.
    DO i = 1, SIZE(active_equations)
       mask(active_equations(i)) = .TRUE.
    END DO
  END FUNCTION equation_is_active


END MODULE nonlinear_solver_2d
