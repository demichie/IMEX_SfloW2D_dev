!********************************************************************************
!> \brief Local nonlinear solver for the implicit IMEX stage
!>
!> Newton iteration, line search, residual/Jacobian evaluation and small dense
!> linear solves.  This module is intentionally independent of the spatial
!> domain and operates on one cell state at a time.
!********************************************************************************

MODULE nonlinear_solver_2d

  USE diagnostics_2d, ONLY : debug_pause

  USE lapack_interfaces_2d, ONLY : solve_dense_system

  USE parameters_2d, ONLY : wp, n_eqns, n_vars, verbose_level

  USE equation_metadata_2d, ONLY : equation_partition_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: initialize_nonlinear_solver
  PUBLIC :: finalize_nonlinear_solver
  PUBLIC :: solve_imex_stage
  PUBLIC :: solve_implicit_stage
  PUBLIC :: solve_rk_step

  !> Complex-step perturbation and its reciprocal.
  REAL(wp) :: h, one_by_h

  !> Relative threshold used to reject numerically singular small systems.
  REAL(wp), PARAMETER :: pivot_tol_factor = 64.0_wp

CONTAINS

  !******************************************************************************
  !> \brief Coordinate one IMEX stage and return its effective IM contribution.
  !> \details Apply the semi-implicit contribution and its corrections before
  !>          the ordered implicit-group sweep. The returned IM term is the
  !>          net conservative increment over the stage interval.
  !******************************************************************************
  SUBROUTINE solve_imex_stage(equation_partition, q_known, q_old, qp_known,  &
       dt_step, a_diag, known_stage_residual, Bprimej_x, Bprimej_y,            &
       Bsecondj_xx, Bsecondj_xy, Bsecondj_yy, grav_coeff, Zij,                &
       zero_x_momentum, zero_y_momentum, q_stage, effective_im_term,          &
       iterations_used, converged, linear_info, line_search_failed)

    USE equation_terms_2d, ONLY : eval_nh_semi_impl_terms
    USE rheology_2d, ONLY : correct_friction_momentum

    IMPLICIT NONE

    TYPE(equation_partition_type), INTENT(IN) :: equation_partition
    REAL(wp), INTENT(IN) :: q_known(n_vars), q_old(n_vars)
    REAL(wp), INTENT(IN) :: qp_known(n_vars+2)
    REAL(wp), INTENT(IN) :: dt_step, a_diag
    REAL(wp), INTENT(IN) :: known_stage_residual(n_eqns)
    REAL(wp), INTENT(IN) :: Bprimej_x, Bprimej_y
    REAL(wp), INTENT(IN) :: Bsecondj_xx, Bsecondj_xy, Bsecondj_yy
    REAL(wp), INTENT(IN) :: grav_coeff, Zij
    LOGICAL, INTENT(IN) :: zero_x_momentum, zero_y_momentum
    REAL(wp), INTENT(OUT) :: q_stage(n_vars), effective_im_term(n_eqns)
    INTEGER, INTENT(OUT) :: iterations_used, linear_info
    LOGICAL, INTENT(OUT) :: converged, line_search_failed

    REAL(wp) :: si_term(n_eqns)
    REAL(wp) :: q_si(n_vars), residual_not_impl(n_eqns)

    CALL eval_nh_semi_impl_terms(Bprimej_x, Bprimej_y, Bsecondj_xx,          &
         Bsecondj_xy, Bsecondj_yy, grav_coeff, qp_known, si_term, Zij)

    q_si = q_known + dt_step * a_diag * si_term(1:n_vars)
    CALL correct_friction_momentum(q_known(2:3), q_si(2:3))
    si_term(1:n_vars) = (q_si - q_known) / (dt_step * a_diag)

    q_stage = q_si
    residual_not_impl = known_stage_residual - a_diag * si_term
    CALL solve_implicit_groups(equation_partition, q_stage, q_old, dt_step,   &
         a_diag, residual_not_impl, Bprimej_x, Bprimej_y, Zij,                &
         iterations_used, converged, linear_info, line_search_failed)

    IF (zero_y_momentum) q_stage(3) = 0.0_wp
    IF (zero_x_momentum) q_stage(2) = 0.0_wp

    CALL correct_friction_momentum(q_si(2:3), q_stage(2:3))

    effective_im_term = 0.0_wp
    effective_im_term(1:n_vars) = (q_stage - q_known) / (dt_step * a_diag)

  END SUBROUTINE solve_imex_stage

  !******************************************************************************
  !> \brief Backward-compatible name for the local IMEX stage coordinator.
  !******************************************************************************
  SUBROUTINE solve_implicit_stage(equation_partition, q_known, q_old,         &
       qp_known, dt_step, a_diag, known_stage_residual, Bprimej_x, Bprimej_y, &
       Bsecondj_xx, Bsecondj_xy, Bsecondj_yy, grav_coeff, Zij,               &
       zero_x_momentum, zero_y_momentum, q_stage, effective_im_term,         &
       iterations_used, converged, linear_info, line_search_failed)

    TYPE(equation_partition_type), INTENT(IN) :: equation_partition
    REAL(wp), INTENT(IN) :: q_known(n_vars), q_old(n_vars)
    REAL(wp), INTENT(IN) :: qp_known(n_vars+2)
    REAL(wp), INTENT(IN) :: dt_step, a_diag
    REAL(wp), INTENT(IN) :: known_stage_residual(n_eqns)
    REAL(wp), INTENT(IN) :: Bprimej_x, Bprimej_y
    REAL(wp), INTENT(IN) :: Bsecondj_xx, Bsecondj_xy, Bsecondj_yy
    REAL(wp), INTENT(IN) :: grav_coeff, Zij
    LOGICAL, INTENT(IN) :: zero_x_momentum, zero_y_momentum
    REAL(wp), INTENT(OUT) :: q_stage(n_vars), effective_im_term(n_eqns)
    INTEGER, INTENT(OUT) :: iterations_used, linear_info
    LOGICAL, INTENT(OUT) :: converged, line_search_failed

    CALL solve_imex_stage(equation_partition, q_known, q_old, qp_known,      &
         dt_step, a_diag, known_stage_residual, Bprimej_x, Bprimej_y,         &
         Bsecondj_xx, Bsecondj_xy, Bsecondj_yy, grav_coeff, Zij,              &
         zero_x_momentum, zero_y_momentum, q_stage, effective_im_term,        &
         iterations_used, converged, linear_info, line_search_failed)

  END SUBROUTINE solve_implicit_stage

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

  !> \brief Backward-compatible entry point for the implicit-group sweep.
  !******************************************************************************
  SUBROUTINE solve_rk_step(equation_partition, qj, qj_old, dt_step, a_diag,  &
       Rj_not_impl, Bprimej_x, Bprimej_y, Zij, iterations_used, converged,  &
       linear_info, line_search_failed)

    TYPE(equation_partition_type), INTENT(IN) :: equation_partition
    REAL(wp), INTENT(INOUT) :: qj(n_vars)
    REAL(wp), INTENT(IN) :: qj_old(n_vars)
    REAL(wp), INTENT(IN) :: dt_step, a_diag
    REAL(wp), INTENT(IN) :: Rj_not_impl(n_eqns)
    REAL(wp), INTENT(IN) :: Bprimej_x, Bprimej_y, Zij
    INTEGER, INTENT(OUT) :: iterations_used, linear_info
    LOGICAL, INTENT(OUT) :: converged, line_search_failed

    CALL solve_implicit_groups(equation_partition, qj, qj_old, dt_step,      &
         a_diag, Rj_not_impl, Bprimej_x, Bprimej_y, Zij, iterations_used,   &
         converged, linear_info, line_search_failed)

  END SUBROUTINE solve_rk_step

  !> \brief Apply the ordered implicit-group sweep using Newton or closed-form DIRK solves.
  !>
  !> \param[in] equation_partition Explicit/implicit equation mask and compact index maps.
  !> \param[in,out] qj Local conservative state vector; mass and momenta are per unit horizontal
  !>                   area.
  !> \param[in] qj_old Cell conservative reference state used in the implicit-stage residual.
  !> \param[in] dt_step Timestep entering the implicit-stage equation [s].
  !> \param[in] a_diag Diagonal coefficient of the current implicit Runge-Kutta stage.
  !> \param[in] Rj_not_impl Already known stage contribution excluding the current implicit source.
  !> \param[in] Bprimej_x Filtered cell bed derivative dB/dx [dimensionless].
  !> \param[in] Bprimej_y Filtered cell bed derivative dB/dy [dimensionless].
  !> \param[in] Zij Local effective stochastic fluctuation supplied to the rheology.
  !> \param[out] iterations_used Number of Newton iterations attempted for this cell/stage.
  !> \param[out] converged True when every Newton stage residual meets its bound.
  !> \param[out] linear_info Linear-solver status: zero for success, nonzero for a failed solve.
  !> \param[out] line_search_failed True when backtracking fails to accept a step.
  !******************************************************************************

  SUBROUTINE solve_implicit_groups(equation_partition, qj, qj_old, dt_step,  &
       a_diag, Rj_not_impl,                                                   &
       Bprimej_x, Bprimej_y, Zij,                                              &
       iterations_used, converged, linear_info, line_search_failed )

    USE parameters_2d, ONLY : max_nl_iter, tol_rel, tol_abs, n_RK,            &
         n_solid, n_add_gas, idx_pore
    USE equation_metadata_2d, ONLY : IM_GROUP_METHOD_ANALYTIC,                &
         IM_GROUP_METHOD_NEWTON
    USE parameters_2d, ONLY : rheology_model_id, rheology_id_lava, rheology_flag
    USE constitutive_parameters_2d, ONLY : T_env, radiative_term_coeff,       &
         convective_term_coeff
    USE state_conversion_2d, ONLY : r_phys_var
    USE thermal_sources_2d, ONLY : solve_lava_analytic_temperature
    USE pore_pressure_2d, ONLY : solve_analytic_pore_pressure

    IMPLICIT NONE

    TYPE(equation_partition_type), INTENT(IN) :: equation_partition
    REAL(wp), INTENT(INOUT) :: qj(n_vars)
    REAL(wp), INTENT(IN) :: qj_old(n_vars)
    REAL(wp), INTENT(IN) :: dt_step
    REAL(wp), INTENT(IN) :: a_diag
    REAL(wp), INTENT(IN) :: Rj_not_impl(n_eqns)
    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y
    REAL(wp), INTENT(IN):: Zij ! value stochastic process
    INTEGER, INTENT(OUT) :: iterations_used
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: linear_info
    LOGICAL, INTENT(OUT) :: line_search_failed

    REAL(wp) :: qj_init(n_vars), q_known(n_vars)
    REAL(wp) :: coeff_f(n_eqns)
    REAL(wp) :: tol_abs_local, tol_rel_local
    REAL(wp) :: lava_h, lava_u, lava_v, lava_temperature, lava_density
    REAL(wp) :: lava_alphas(n_solid), lava_alphag(n_add_gas), lava_alphal
    REAL(wp) :: lava_red_grav, lava_p_dyn, lava_Zs, lava_pore_pressure
    REAL(wp) :: pore_h, pore_u, pore_v, pore_temperature, pore_density
    REAL(wp) :: pore_alphas(n_solid), pore_alphag(n_add_gas), pore_alphal
    REAL(wp) :: pore_red_grav, pore_p_dyn, pore_Zs, pore_pressure
    REAL(wp) :: temperature_stage, heat_capacity
    REAL(wp) :: pressure_known, pressure_stage
    INTEGER :: group, group_size, equation, i
    INTEGER :: group_equations(equation_partition%n_implicit)

    iterations_used = 0
    converged = .FALSE.
    linear_info = 0
    line_search_failed = .FALSE.

    IF (.NOT. equation_partition%is_initialized() .OR.                     &
         equation_partition%n_equations .NE. n_eqns .OR.                   &
         n_eqns .NE. n_vars .OR. equation_partition%n_implicit .GT. n_vars) THEN
       linear_info = -1
       RETURN
    END IF

    coeff_f = 1.0_wp

    ! The four-stage third-order method needs a more accurate solution of
    ! each implicit stage than the default tolerances provide. Keep the
    ! configured tolerances for lower-order methods; tighten only this solve.
    tol_abs_local = tol_abs
    tol_rel_local = tol_rel
    IF ( n_RK .EQ. 4 ) THEN
       tol_abs_local = 1.0E-8_wp
       tol_rel_local = 1.0E-12_wp
    END IF

    qj_init = qj
    q_known = qj_old - dt_step * Rj_not_impl(1:n_vars)
    qj = q_known

    ! Explicit equations have no current-stage implicit source and are solved
    ! directly before the ordered implicit-group sweep.
    DO i = 1, equation_partition%n_explicit
       equation = equation_partition%explicit_map(i)
       qj(equation) = q_known(equation)
    END DO

    IF (equation_partition%n_IM_groups .EQ. 0) THEN
       converged = .TRUE.
       RETURN
    END IF

    DO group = 1, equation_partition%n_IM_groups
       group_size = equation_partition%n_IM_group_equations(group)
       group_equations(1:group_size) =                                      &
            equation_partition%IM_group_map(1:group_size,group)

       IF (equation_partition%IM_group_method(group) .EQ.                   &
            IM_GROUP_METHOD_ANALYTIC) THEN
          IF (group_size .NE. 1) THEN
             linear_info = -6
             qj = qj_init
             RETURN
          END IF
          equation = group_equations(1)
          SELECT CASE (equation)
          CASE (4)
             IF (rheology_model_id .NE. rheology_id_lava) THEN
                linear_info = -2
                qj = qj_init
                RETURN
             END IF
             IF (.NOT. rheology_flag) THEN
                CYCLE
             END IF
             CALL r_phys_var(qj, lava_h, lava_u, lava_v, lava_alphas,           &
                  lava_density, lava_temperature, lava_alphal, lava_alphag,   &
                  lava_red_grav, lava_p_dyn, lava_Zs, lava_pore_pressure)
             IF (lava_temperature .LE. 0.0_wp) THEN
                linear_info = -3
                qj = qj_init
                RETURN
             END IF
             heat_capacity = qj(4) / lava_temperature
             CALL solve_lava_analytic_temperature(lava_temperature,             &
                  dt_step*a_diag, heat_capacity, radiative_term_coeff,          &
                  convective_term_coeff, T_env, temperature_stage)
             qj(4) = heat_capacity * temperature_stage
          CASE DEFAULT
             IF (equation .NE. idx_pore) THEN
                linear_info = -4
                qj = qj_init
                RETURN
             END IF
             IF (qj(1) .GT. 0.0_wp) THEN
                CALL r_phys_var(qj, pore_h, pore_u, pore_v, pore_alphas,        &
                     pore_density, pore_temperature, pore_alphal, pore_alphag, &
                     pore_red_grav, pore_p_dyn, pore_Zs, pore_pressure)
                pressure_known = qj(idx_pore) / qj(1)
                CALL solve_analytic_pore_pressure(pressure_known,              &
                     dt_step*a_diag, pore_h, pore_temperature, pore_alphas,    &
                     pressure_stage)
                qj(idx_pore) = qj(1) * pressure_stage
             END IF
             END SELECT

       ELSEIF (equation_partition%IM_group_method(group) .EQ.              &
               IM_GROUP_METHOD_NEWTON) THEN
          DO i = 1, group_size
             equation = group_equations(i)
             qj(equation) = qj_init(equation)
          END DO
          CALL solve_newton_group(group_equations(1:group_size), qj,        &
               qj_old, dt_step, a_diag, Rj_not_impl, Bprimej_x, Bprimej_y, &
               Zij, coeff_f, tol_abs_local, tol_rel_local, iterations_used,&
               converged, linear_info, line_search_failed)
          IF (.NOT. converged) THEN
             qj = qj_init
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

  !> \brief Solve one implicit Newton group while holding other equations fixed.
  SUBROUTINE solve_newton_group(group_equations, qj, qj_old, dt_step, a_diag,&
       Rj_not_impl, Bprimej_x, Bprimej_y, Zij, coeff_f, tol_abs_local,       &
       tol_rel_local, iterations_used, converged, linear_info,               &
       line_search_failed)

    USE parameters_2d, ONLY : max_nl_iter

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: group_equations(:)
    REAL(wp), INTENT(INOUT) :: qj(n_vars)
    REAL(wp), INTENT(IN) :: qj_old(n_vars), dt_step, a_diag
    REAL(wp), INTENT(IN) :: Rj_not_impl(n_eqns), Bprimej_x, Bprimej_y, Zij
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
       CALL eval_f(qj, qj_old, dt_step, a_diag, coeff_f, Rj_not_impl,       &
            Bprimej_x, Bprimej_y, right_term, scal_f, Zij, group_equations)

       IF (nl_iter .EQ. 1) residual_ref = ABS(right_term)
       residual_tol = tol_abs_local + tol_rel_local*residual_ref
       IF (ALL(ABS(right_term) .LE. residual_tol)) THEN
          converged = .TRUE.
          EXIT
       END IF

       CALL eval_jacobian(qj_rel, qj_org, dt_step, a_diag, coeff_f,         &
            Bprimej_x, Bprimej_y, left_matrix, Zij, group_equations)
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
          CALL lnsrch(qj_rel_old, qj_org, qj, qj_old, scal_f_old, grad_f,   &
               desc_dir, coeff_f, qj_rel, scal_f, right_term, stpmax,      &
               check, dt_step, a_diag, Rj_not_impl, Bprimej_x, Bprimej_y,  &
               Zij, group_equations)
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
          CALL eval_f(qj_trial, qj_old, dt_step, a_diag, coeff_f,           &
               Rj_not_impl, Bprimej_x, Bprimej_y, right_term, scal_f, Zij, &
               group_equations)
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
  !> \param[in] qj_old Cell conservative reference state used in the implicit-stage residual.
  !> \param[in] scal_f_old Initial value of the half-squared scaled residual norm.
  !> \param[in] grad_f Gradient of the scalar residual merit function at the starting point.
  !> \param[in,out] desc_dir Newton descent direction, limited to the maximum permitted step length.
  !> \param[in] coeff_f Per-equation scaling factors applied to the stage residual.
  !> \param[out] qj_rel Dimensionless normalized conservative state, qj=qj_rel*qj_org.
  !> \param[out] scal_f Half the squared Euclidean norm of the scaled nonlinear residual.
  !> \param[in,out] right_term Initial residual on entry; residual of the accepted line-search state
  !>                           on exit.
  !> \param[in] stpmax Maximum Euclidean length of the normalized trial step.
  !> \param[out] check True when the line search cannot find a usable decrease.
  !> \param[in] dt_step Timestep entering the implicit-stage equation [s].
  !> \param[in] a_diag Diagonal coefficient of the current implicit Runge-Kutta stage.
  !> \param[in] Rj_not_impl Already known stage contribution excluding the current implicit source.
  !> \param[in] Bprimej_x Filtered cell bed derivative dB/dx [dimensionless].
  !> \param[in] Bprimej_y Filtered cell bed derivative dB/dy [dimensionless].
  !> \param[in] Zij Local effective stochastic fluctuation supplied to the rheology.
  !******************************************************************************

  SUBROUTINE lnsrch( qj_rel_NR_old, qj_org, qj_fixed, qj_old, scal_f_old,      &
       grad_f, desc_dir, coeff_f, qj_rel, scal_f, right_term, stpmax, check,  &
       dt_step , a_diag , Rj_not_impl , Bprimej_x , Bprimej_y, Zij,           &
       active_equations )

    IMPLICIT NONE

    !> Initial point
    REAL(wp), DIMENSION(:), INTENT(IN) :: qj_rel_NR_old

    !> Initial point
    REAL(wp), DIMENSION(:), INTENT(IN) :: qj_org

    ! Non-group variables must remain bitwise unchanged during line search.
    REAL(wp), INTENT(IN) :: qj_fixed(n_vars)

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
    INTEGER, INTENT(IN) :: active_equations(:)
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

       CALL eval_f( qj , qj_old , dt_step , a_diag , coeff_f ,             &
            Rj_not_impl , Bprimej_x , Bprimej_y, right_term , scal_f, Zij,  &
            active_equations )

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
  !> \brief Evaluate the scaled implicit-stage residual and its scalar merit function.
  !
  !> This subroutine evaluate the value of the nonlinear system in the state
  !> defined by the variables qj.
  !> \date 2019/12/16
  !> @author
  !> Mattia de' Michieli Vitturi
  !>
  !> \param[in] qj Local conservative state vector; mass and momenta are per unit horizontal area.
  !> \param[in] qj_old Cell conservative reference state used in the implicit-stage residual.
  !> \param[in] dt_step Timestep entering the implicit-stage equation [s].
  !> \param[in] a_diag Diagonal coefficient of the current implicit Runge-Kutta stage.
  !> \param[in] coeff_f Per-equation scaling factors applied to the stage residual.
  !> \param[in] Rj_not_impl Already known stage contribution excluding the current implicit source.
  !> \param[in] Bprimej_x Filtered cell bed derivative dB/dx [dimensionless].
  !> \param[in] Bprimej_y Filtered cell bed derivative dB/dy [dimensionless].
  !> \param[out] f_nl Scaled nonlinear stage residual, one entry per equation.
  !> \param[out] scal_f Half the squared Euclidean norm of the scaled nonlinear residual.
  !> \param[in] Zij Local effective stochastic fluctuation supplied to the rheology.
  !******************************************************************************

  SUBROUTINE eval_f( qj , qj_old , dt_step, a_diag , coeff_f , Rj_not_impl ,    &
       Bprimej_x ,                                                              &
       Bprimej_y , f_nl , scal_f, Zij, active_equations )

    USE equation_terms_2d, ONLY : eval_implicit_terms

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
    INTEGER, INTENT(IN) :: active_equations(:)

    REAL(wp) :: nh_term_impl(n_eqns)
    REAL(wp) :: Rj(n_eqns)
    ! Keep the conservative difference separate from the forcing even under
    ! -Ofast reassociation. Otherwise q + (dt*R-q_old) can lose a small forcing
    ! before subtracting equal large states and falsely report zero residual.
    REAL(wp), VOLATILE :: delta_q(n_vars)

    CALL eval_implicit_terms( Bprimej_x , Bprimej_y, Zij, r_qj = qj ,           &
         r_nh_term_impl=nh_term_impl )

    Rj = Rj_not_impl - a_diag * nh_term_impl

    delta_q = qj - qj_old
    f_nl = delta_q + dt_step * Rj
    f_nl = coeff_f * f_nl
    WHERE (.NOT. equation_is_active(active_equations))
       f_nl = 0.0_wp
    END WHERE

    scal_f = 0.5_wp * DOT_PRODUCT( f_nl , f_nl )

    RETURN

  END SUBROUTINE eval_f

  !******************************************************************************
  !> \brief Build the scaled implicit-stage Jacobian using complex-step differentiation.
  !
  !> This subroutine evaluate the jacobian of the non-linear system
  !> with respect to the conservative variables.
  !
  !
  !> \date 07/10/2016
  !> @author
  !> Mattia de' Michieli Vitturi
  !>
  !> \param[in] equation_partition Explicit/implicit equation mask and compact index maps.
  !> \param[in] qj_rel Dimensionless normalized conservative state, qj=qj_rel*qj_org.
  !> \param[in] qj_org Per-variable scaling factors for normalized Newton variables.
  !> \param[in] dt_step Timestep entering the implicit-stage equation [s].
  !> \param[in] a_diag Diagonal coefficient of the current implicit Runge-Kutta stage.
  !> \param[in] coeff_f Per-equation scaling factors applied to the stage residual.
  !> \param[in] Bprimej_x Filtered cell bed derivative dB/dx [dimensionless].
  !> \param[in] Bprimej_y Filtered cell bed derivative dB/dy [dimensionless].
  !> \param[out] left_matrix Jacobian of the scaled residual with respect to the normalized state.
  !> \param[in] Zij Local effective stochastic fluctuation supplied to the rheology.
  !******************************************************************************

  SUBROUTINE eval_jacobian( qj_rel, qj_org, dt_step, a_diag, coeff_f,        &
       Bprimej_x, Bprimej_y, left_matrix, Zij, active_equations )

    USE equation_terms_2d, ONLY : eval_implicit_terms

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
    INTEGER, INTENT(IN) :: active_equations(:)

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

       CALL eval_implicit_terms(Bprimej_x, Bprimej_y, Zij,                   &
            c_qj=qj_cmplx, c_nh_term_impl=nh_terms_cmplx_impl)

       DO j = 1, SIZE(active_equations)
          Jacob_relax(active_equations(j),equation) =                         &
               coeff_f(active_equations(j)) *                                &
               AIMAG(nh_terms_cmplx_impl(active_equations(j))) * one_by_h
       END DO
       DO j = 1, n_eqns
          left_matrix(j,equation) = left_matrix(j,equation) -                 &
               dt_step*a_diag*Jacob_relax(j,equation)
       END DO
    END DO

    RETURN

  END SUBROUTINE eval_jacobian

  !> \brief Return a mask selecting residual rows in one implicit group.
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
