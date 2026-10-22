!********************************************************************************
!> \brief Exact-flow helpers for frozen-coefficient rheology sources.
!********************************************************************************

MODULE rheology_analytic_utils_2d

  USE parameters_2d, ONLY : wp
  USE diagnostics_2d, ONLY : fatal_error
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: solve_quadratic_in_interval
  PUBLIC :: solve_cubic_in_interval
  PUBLIC :: lambert_w0
  PUBLIC :: solve_exact_velocity_flow
  PUBLIC :: EXACT_FLOW_CONSTANT, EXACT_FLOW_LINEAR, EXACT_FLOW_AFFINE
  PUBLIC :: EXACT_FLOW_QUADRATIC, EXACT_FLOW_QUADRATIC_POLYNOMIAL
  PUBLIC :: EXACT_FLOW_FROUDE, EXACT_FLOW_LUCAS
  PUBLIC :: EXACT_FLOW_MUI, EXACT_FLOW_MUI_VOELLMY
  PUBLIC :: validate_analytic_velocity_inputs
  PUBLIC :: fatal_analytic_input_error

  INTEGER, PARAMETER :: EXACT_FLOW_CONSTANT = 1
  INTEGER, PARAMETER :: EXACT_FLOW_LINEAR = 2
  INTEGER, PARAMETER :: EXACT_FLOW_AFFINE = 3
  INTEGER, PARAMETER :: EXACT_FLOW_QUADRATIC = 4
  INTEGER, PARAMETER :: EXACT_FLOW_QUADRATIC_POLYNOMIAL = 5
  INTEGER, PARAMETER :: EXACT_FLOW_FROUDE = 6
  INTEGER, PARAMETER :: EXACT_FLOW_LUCAS = 7
  INTEGER, PARAMETER :: EXACT_FLOW_MUI = 8
  INTEGER, PARAMETER :: EXACT_FLOW_MUI_VOELLMY = 9

CONTAINS

  !******************************************************************************
  !> \brief Report a fatal error raised by an exact rheology flow routine.
  !******************************************************************************
  SUBROUTINE fatal_analytic_input_error(message)

    CHARACTER(LEN=*), INTENT(IN) :: message

    CALL fatal_error(message)

  END SUBROUTINE fatal_analytic_input_error

  !******************************************************************************
  !> \brief Validate inputs to a module-level exact velocity solve.
  !******************************************************************************
  SUBROUTINE validate_analytic_velocity_inputs(integration_step, flow_depth, &
       speed_ratio, speed_before, rheology_name)

    REAL(wp), INTENT(IN) :: integration_step, flow_depth, speed_ratio
    REAL(wp), INTENT(IN) :: speed_before
    CHARACTER(LEN=*), INTENT(IN) :: rheology_name

    IF (integration_step .LT. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Negative integration step in '//TRIM(rheology_name)//' exact solver')
    END IF
    IF (integration_step .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN

    IF (flow_depth .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive flow depth in '//TRIM(rheology_name)//' analytic solver')
    END IF
    IF (speed_ratio .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive speed ratio in '//TRIM(rheology_name)//' analytic solver')
    END IF
    IF (.NOT. ieee_is_finite(integration_step) .OR.                         &
         .NOT. ieee_is_finite(flow_depth) .OR.                              &
         .NOT. ieee_is_finite(speed_ratio) .OR.                             &
         .NOT. ieee_is_finite(speed_before)) THEN
       CALL fatal_analytic_input_error('Nonfinite input in '//TRIM(rheology_name)//' exact solver')
    END IF

  END SUBROUTINE validate_analytic_velocity_inputs

  !******************************************************************************
  !> \brief Integrate a frozen-coefficient, direction-preserving velocity source.
  !> \details Simple constant, linear, affine, and quadratic laws use explicit
  !>          exact-flow formulas. Other supported laws are advanced by
  !>          inverting their separated scalar time integral, without taking
  !>          time-discretization substeps.
  !******************************************************************************
  SUBROUTINE solve_exact_velocity_flow(flow_model, coefficients, time_step, &
       speed_before, speed_after)

    INTEGER, INTENT(IN) :: flow_model
    REAL(wp), INTENT(IN) :: coefficients(:), time_step, speed_before
    REAL(wp), INTENT(OUT) :: speed_after

    REAL(wp) :: constant_drag, linear_drag, quadratic_drag, exponent
    REAL(wp) :: exponential_offset, rate_constant, value

    speed_after = speed_before
    IF (time_step .LT. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Negative time step in exact velocity integrator')
    END IF
    IF (.NOT. ieee_is_finite(time_step) .OR.                               &
         .NOT. ieee_is_finite(speed_before) .OR.                           &
         .NOT. ALL(ieee_is_finite(coefficients))) THEN
       CALL fatal_analytic_input_error('Nonfinite input to exact velocity integrator')
    END IF
    SELECT CASE (flow_model)
    CASE (EXACT_FLOW_CONSTANT, EXACT_FLOW_LINEAR, EXACT_FLOW_QUADRATIC)
       IF (SIZE(coefficients) .NE. 1)                                      &
            CALL fatal_analytic_input_error('Invalid coefficient count for exact velocity flow')
    CASE (EXACT_FLOW_AFFINE)
       IF (SIZE(coefficients) .NE. 2)                                      &
            CALL fatal_analytic_input_error('Invalid coefficient count for exact velocity flow')
    CASE (EXACT_FLOW_QUADRATIC_POLYNOMIAL)
       IF (SIZE(coefficients) .NE. 3)                                      &
            CALL fatal_analytic_input_error('Invalid coefficient count for exact velocity flow')
    CASE (EXACT_FLOW_FROUDE)
       IF (SIZE(coefficients) .NE. 5)                                      &
            CALL fatal_analytic_input_error('Invalid coefficient count for exact velocity flow')
    CASE (EXACT_FLOW_LUCAS)
       IF (SIZE(coefficients) .NE. 7)                                      &
            CALL fatal_analytic_input_error('Invalid coefficient count for exact velocity flow')
    CASE (EXACT_FLOW_MUI)
       IF (SIZE(coefficients) .NE. 5)                                      &
            CALL fatal_analytic_input_error('Invalid coefficient count for exact velocity flow')
    CASE (EXACT_FLOW_MUI_VOELLMY)
       IF (SIZE(coefficients) .NE. 6)                                      &
            CALL fatal_analytic_input_error('Invalid coefficient count for exact velocity flow')
    CASE DEFAULT
       CALL fatal_analytic_input_error('Unknown exact velocity flow model')
    END SELECT
    IF (time_step .EQ. 0.0_wp .OR. speed_before .LE. 0.0_wp) RETURN

    SELECT CASE (flow_model)
    CASE (EXACT_FLOW_CONSTANT)
       speed_after = MAX(0.0_wp, speed_before - coefficients(1)*time_step)
    CASE (EXACT_FLOW_LINEAR)
       speed_after = speed_before * EXP(-coefficients(1)*time_step)
    CASE (EXACT_FLOW_AFFINE)
       constant_drag = coefficients(1)
       linear_drag = coefficients(2)
       IF (linear_drag .GT. 0.0_wp) THEN
          exponent = EXP(-linear_drag*time_step)
          speed_after = MAX(0.0_wp, (speed_before + constant_drag/linear_drag)* &
               exponent - constant_drag/linear_drag)
       ELSE
          speed_after = MAX(0.0_wp, speed_before-constant_drag*time_step)
       END IF
    CASE (EXACT_FLOW_QUADRATIC)
       quadratic_drag = coefficients(1)
       speed_after = speed_before /                                     &
            (1.0_wp + quadratic_drag*speed_before*time_step)
    CASE (EXACT_FLOW_QUADRATIC_POLYNOMIAL)
       CALL solve_quadratic_polynomial_flow(coefficients, time_step,        &
            speed_before, speed_after)
    CASE (EXACT_FLOW_FROUDE)
       IF (coefficients(3) .GE. 0.0_wp .AND.                               &
            coefficients(4) .GE. coefficients(5) .AND.                     &
            coefficients(1) .GT. 0.0_wp) THEN
          exponential_offset = (coefficients(4)-coefficients(5))*          &
               EXP(-coefficients(3))
          rate_constant = coefficients(1)
          IF (coefficients(2) .LE. 0.0_wp .OR. exponential_offset .EQ. 0.0_wp) THEN
             value = coefficients(5) + exponential_offset
             speed_after = MAX(0.0_wp, speed_before-rate_constant*value*time_step)
          ELSEIF (coefficients(5) .GT. 0.0_wp) THEN
             value = coefficients(5)*EXP(coefficients(2)*speed_before) +   &
                  exponential_offset
             value = value*EXP(-rate_constant*coefficients(2)*             &
                  coefficients(5)*time_step) - exponential_offset
             IF (value .LE. coefficients(5)) THEN
                speed_after = 0.0_wp
             ELSE
                speed_after = LOG(value/coefficients(5))/coefficients(2)
             END IF
          ELSE
             value = EXP(coefficients(2)*speed_before) -                   &
                  rate_constant*exponential_offset*coefficients(2)*time_step
             IF (value .LE. 1.0_wp) THEN
                speed_after = 0.0_wp
             ELSE
                speed_after = LOG(value)/coefficients(2)
             END IF
          END IF
       ELSE
          CALL solve_by_exact_time_integral(flow_model, coefficients,       &
               time_step, speed_before, speed_after)
       END IF
    CASE (EXACT_FLOW_LUCAS, EXACT_FLOW_MUI, EXACT_FLOW_MUI_VOELLMY)
       CALL solve_by_exact_time_integral(flow_model, coefficients,          &
            time_step, speed_before, speed_after)
    END SELECT

  END SUBROUTINE solve_exact_velocity_flow

  !******************************************************************************
  !> \brief Advance a polynomial drag law c + a*r + b*r**2.
  !******************************************************************************
  SUBROUTINE solve_quadratic_polynomial_flow(coefficients, time_step,       &
       speed_before, speed_after)

    REAL(wp), INTENT(IN) :: coefficients(:), time_step, speed_before
    REAL(wp), INTENT(OUT) :: speed_after

    INTEGER, PARAMETER :: max_iterations = 80
    REAL(wp) :: constant_drag, linear_drag, quadratic_drag
    REAL(wp) :: lower, upper, midpoint, target_primitive, primitive
    REAL(wp) :: exponential_factor, tolerance
    INTEGER :: iteration

    constant_drag = coefficients(1)
    linear_drag = coefficients(2)
    quadratic_drag = coefficients(3)
    speed_after = speed_before

    IF (MIN(constant_drag, linear_drag, quadratic_drag) .LT. 0.0_wp) THEN
       CALL solve_by_exact_time_integral(EXACT_FLOW_QUADRATIC_POLYNOMIAL,   &
            coefficients, time_step, speed_before, speed_after)
       RETURN
    END IF
    IF (constant_drag .EQ. 0.0_wp) THEN
       IF (linear_drag .GT. 0.0_wp) THEN
          exponential_factor = EXP(-linear_drag*time_step)
          speed_after = speed_before*exponential_factor /                   &
               (1.0_wp + (quadratic_drag*speed_before/linear_drag)*         &
               (1.0_wp-exponential_factor))
       ELSEIF (quadratic_drag .GT. 0.0_wp) THEN
          speed_after = speed_before /                                      &
               (1.0_wp+quadratic_drag*speed_before*time_step)
       END IF
       RETURN
    END IF
    IF (quadratic_drag .EQ. 0.0_wp) THEN
       IF (linear_drag .GT. 0.0_wp) THEN
          exponential_factor = EXP(-linear_drag*time_step)
          speed_after = MAX(0.0_wp, (speed_before+constant_drag/linear_drag)*&
               exponential_factor-constant_drag/linear_drag)
       ELSE
          speed_after = MAX(0.0_wp, speed_before-constant_drag*time_step)
       END IF
       RETURN
    END IF

    primitive = quadratic_decay_time(speed_before, constant_drag,           &
         linear_drag, quadratic_drag)
    target_primitive = primitive-time_step
    IF (target_primitive .LE. 0.0_wp) THEN
       speed_after = 0.0_wp
       RETURN
    END IF

    lower = 0.0_wp
    upper = speed_before
    tolerance = 32.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp, speed_before)
    DO iteration = 1, max_iterations
       midpoint = lower+0.5_wp*(upper-lower)
       IF (midpoint .EQ. lower .OR. midpoint .EQ. upper) EXIT
       primitive = quadratic_decay_time(midpoint, constant_drag,            &
            linear_drag, quadratic_drag)
       IF (primitive .LT. target_primitive) THEN
          lower = midpoint
       ELSE
          upper = midpoint
       END IF
       IF (upper-lower .LE. tolerance) EXIT
    END DO
    IF (iteration .GT. max_iterations) THEN
       CALL fatal_analytic_input_error('Polynomial exact velocity inversion failed to converge')
    END IF
    speed_after = lower+0.5_wp*(upper-lower)

  END SUBROUTINE solve_quadratic_polynomial_flow

  PURE FUNCTION quadratic_decay_time(speed, constant_drag, linear_drag,     &
       quadratic_drag) RESULT(time)

    REAL(wp), INTENT(IN) :: speed, constant_drag, linear_drag, quadratic_drag
    REAL(wp) :: time, discriminant, root_discriminant, numerator, denominator

    discriminant = 4.0_wp*quadratic_drag*constant_drag-linear_drag**2
    IF (ABS(discriminant) .LE. 64.0_wp*EPSILON(1.0_wp)*                     &
         MAX(1.0_wp, linear_drag**2,                                       &
         4.0_wp*quadratic_drag*constant_drag)) THEN
       time = 2.0_wp/linear_drag - 2.0_wp/                                  &
            (2.0_wp*quadratic_drag*speed+linear_drag)
    ELSEIF (discriminant .GT. 0.0_wp) THEN
       root_discriminant = SQRT(discriminant)
       time = 2.0_wp/root_discriminant *                                    &
            (ATAN((2.0_wp*quadratic_drag*speed+linear_drag)/                &
            root_discriminant)-ATAN(linear_drag/root_discriminant))
    ELSE
       root_discriminant = SQRT(-discriminant)
       numerator = 2.0_wp*quadratic_drag*speed + linear_drag
       denominator = (numerator-root_discriminant) /                        &
            (numerator+root_discriminant)
       time = (LOG(denominator)-LOG((linear_drag-root_discriminant)/        &
            (linear_drag+root_discriminant)))/root_discriminant
    END IF

  END FUNCTION quadratic_decay_time

  !******************************************************************************
  !> \brief Invert the exact separated time integral for nonlinear speed laws.
  !******************************************************************************
  SUBROUTINE solve_by_exact_time_integral(flow_model, coefficients, time_step, &
       speed_before, speed_after)

    INTEGER, INTENT(IN) :: flow_model
    REAL(wp), INTENT(IN) :: coefficients(:), time_step, speed_before
    REAL(wp), INTENT(OUT) :: speed_after

    INTEGER, PARAMETER :: max_iterations = 80
    REAL(wp), PARAMETER :: relative_tolerance = 32.0_wp*EPSILON(1.0_wp)
    REAL(wp) :: lower, upper, midpoint, elapsed, stop_time, tolerance
    INTEGER :: iteration

    speed_after = speed_before
    IF (velocity_decay_rate(flow_model, coefficients, speed_before) .LE.   &
         0.0_wp) RETURN

    IF (velocity_decay_rate(flow_model, coefficients, 0.0_wp) .GT. 0.0_wp) THEN
       stop_time = integrate_inverse_decay_rate(flow_model, coefficients,   &
            0.0_wp, speed_before)
       IF (time_step .GE. stop_time) THEN
          speed_after = 0.0_wp
          RETURN
       END IF
    END IF

    lower = 0.0_wp
    upper = speed_before
    DO iteration = 1, max_iterations
       midpoint = lower + 0.5_wp*(upper-lower)
       IF (midpoint .EQ. lower .OR. midpoint .EQ. upper) EXIT
       elapsed = integrate_inverse_decay_rate(flow_model, coefficients,     &
            midpoint, speed_before)
       IF (elapsed .GT. time_step) THEN
          lower = midpoint
       ELSE
          upper = midpoint
       END IF
       tolerance = relative_tolerance*MAX(1.0_wp, speed_before)
       IF (upper-lower .LE. tolerance) EXIT
    END DO
    IF (iteration .GT. max_iterations) THEN
       CALL fatal_analytic_input_error('Exact velocity time-integral inversion failed to converge')
    END IF
    speed_after = lower + 0.5_wp*(upper-lower)

  END SUBROUTINE solve_by_exact_time_integral

  !******************************************************************************
  !> \brief Evaluate the frozen-coefficient scalar speed-decay rate.
  !******************************************************************************
  PURE FUNCTION velocity_decay_rate(flow_model, coefficients, speed)        &
       RESULT(rate)

    INTEGER, INTENT(IN) :: flow_model
    REAL(wp), INTENT(IN) :: coefficients(:), speed
    REAL(wp) :: rate, friction_speed

    SELECT CASE (flow_model)
    CASE (EXACT_FLOW_QUADRATIC_POLYNOMIAL)
       rate = coefficients(1) + coefficients(2)*speed +                    &
            coefficients(3)*speed**2
    CASE (EXACT_FLOW_FROUDE)
       IF (coefficients(2)*speed+coefficients(3) .LE. 0.0_wp) THEN
          rate = coefficients(1)*coefficients(4)
       ELSE
          rate = coefficients(1) * (coefficients(5) +                      &
               (coefficients(4)-coefficients(5))*EXP(                      &
               -coefficients(2)*speed-coefficients(3)))
       END IF
    CASE (EXACT_FLOW_LUCAS)
       friction_speed = coefficients(2)*speed + coefficients(3)
       IF (coefficients(7) .GT. 0.5_wp .OR.                               &
            friction_speed .LT. coefficients(6)) THEN
          rate = coefficients(1)*coefficients(4)
       ELSEIF (friction_speed .GT. 0.0_wp) THEN
          rate = coefficients(1)*(coefficients(5) +                        &
               (coefficients(4)-coefficients(5))*coefficients(6)/friction_speed)
       ELSE
          rate = coefficients(1)*coefficients(4)
       END IF
    CASE (EXACT_FLOW_MUI)
       IF (coefficients(4)+coefficients(5)*speed .GT. TINY(1.0_wp)) THEN
          rate = (coefficients(1) + coefficients(2)*speed +                 &
               coefficients(3)*speed**2) /                                &
               (coefficients(4) + coefficients(5)*speed)
       ELSEIF (coefficients(5) .GT. 0.0_wp) THEN
          rate = coefficients(2)/coefficients(5)
       ELSE
          rate = 0.0_wp
       END IF
    CASE (EXACT_FLOW_MUI_VOELLMY)
       IF (coefficients(4)+coefficients(5)*speed .GT. TINY(1.0_wp)) THEN
          rate = (coefficients(1) + coefficients(2)*speed +                 &
               coefficients(3)*speed**2) /                                &
               (coefficients(4) + coefficients(5)*speed) +                &
               coefficients(6)*speed**2
       ELSEIF (coefficients(5) .GT. 0.0_wp) THEN
          rate = coefficients(2)/coefficients(5) + coefficients(6)*speed**2
       ELSE
          rate = 0.0_wp
       END IF
    CASE DEFAULT
       rate = 0.0_wp
    END SELECT

  END FUNCTION velocity_decay_rate

  !******************************************************************************
  !> \brief Integrate reciprocal decay rate with adaptive Simpson quadrature.
  !******************************************************************************
  FUNCTION integrate_inverse_decay_rate(flow_model, coefficients, lower,   &
       upper) RESULT(integral)

    INTEGER, INTENT(IN) :: flow_model
    REAL(wp), INTENT(IN) :: coefficients(:), lower, upper
    REAL(wp) :: integral, branch_speed

    IF (upper .LE. lower) THEN
       integral = 0.0_wp
       RETURN
    END IF

    branch_speed = -HUGE(1.0_wp)
    SELECT CASE (flow_model)
    CASE (EXACT_FLOW_FROUDE)
       IF (coefficients(2) .GT. 0.0_wp)                                &
            branch_speed = -coefficients(3)/coefficients(2)
    CASE (EXACT_FLOW_LUCAS)
       IF (coefficients(2) .GT. 0.0_wp .AND. coefficients(7) .LE. 0.5_wp) &
            branch_speed = (coefficients(6)-coefficients(3))/coefficients(2)
    END SELECT
    IF (branch_speed .GT. lower .AND. branch_speed .LT. upper) THEN
       integral = integrate_inverse_interval(flow_model, coefficients,     &
            lower, branch_speed) + integrate_inverse_interval(flow_model,  &
            coefficients, branch_speed, upper)
    ELSE
       integral = integrate_inverse_interval(flow_model, coefficients,    &
            lower, upper)
    END IF

  END FUNCTION integrate_inverse_decay_rate

  FUNCTION integrate_inverse_interval(flow_model, coefficients, lower, upper) &
       RESULT(integral)

    INTEGER, INTENT(IN) :: flow_model
    REAL(wp), INTENT(IN) :: coefficients(:), lower, upper
    REAL(wp) :: integral, f_lower, f_middle, f_upper, whole, tolerance
    REAL(wp) :: middle

    f_lower = inverse_rate(flow_model, coefficients, lower)
    f_upper = inverse_rate(flow_model, coefficients, upper)
    middle = lower + 0.5_wp*(upper-lower)
    f_middle = inverse_rate(flow_model, coefficients, middle)
    whole = simpson_estimate(lower, upper, f_lower, f_middle, f_upper)
    tolerance = 2.0E-11_wp*MAX(1.0_wp, ABS(whole))
    integral = adaptive_simpson(flow_model, coefficients, lower, upper,    &
         f_lower, f_middle, f_upper, whole, tolerance, 24)

  END FUNCTION integrate_inverse_interval

  RECURSIVE FUNCTION adaptive_simpson(flow_model, coefficients, lower, upper, &
       f_lower, f_middle, f_upper, whole, tolerance, depth) RESULT(integral)

    INTEGER, INTENT(IN) :: flow_model, depth
    REAL(wp), INTENT(IN) :: coefficients(:), lower, upper
    REAL(wp), INTENT(IN) :: f_lower, f_middle, f_upper, whole, tolerance
    REAL(wp) :: integral, left_middle, right_middle, f_left, f_right
    REAL(wp) :: left_estimate, right_estimate, delta, midpoint

    midpoint = lower + 0.5_wp*(upper-lower)
    left_middle = lower + 0.5_wp*(midpoint-lower)
    right_middle = midpoint + 0.5_wp*(upper-midpoint)
    f_left = inverse_rate(flow_model, coefficients, left_middle)
    f_right = inverse_rate(flow_model, coefficients, right_middle)
    left_estimate = simpson_estimate(lower, midpoint, f_lower, f_left,     &
         f_middle)
    right_estimate = simpson_estimate(midpoint, upper, f_middle, f_right,  &
         f_upper)
    delta = left_estimate + right_estimate - whole

    IF (depth .LE. 0 .AND. ABS(delta) .GT. 15.0_wp*tolerance) THEN
       CALL fatal_analytic_input_error('Exact velocity time-integral quadrature failed to converge')
    END IF
    IF (depth .LE. 0 .OR. ABS(delta) .LE. 15.0_wp*tolerance) THEN
       integral = left_estimate + right_estimate + delta/15.0_wp
       RETURN
    END IF

    integral = adaptive_simpson(flow_model, coefficients, lower, midpoint, &
         f_lower, f_left, f_middle, left_estimate, 0.5_wp*tolerance, depth-1) &
         + adaptive_simpson(flow_model, coefficients, midpoint, upper,     &
         f_middle, f_right, f_upper, right_estimate, 0.5_wp*tolerance, depth-1)

  END FUNCTION adaptive_simpson

  FUNCTION inverse_rate(flow_model, coefficients, speed) RESULT(value)

    INTEGER, INTENT(IN) :: flow_model
    REAL(wp), INTENT(IN) :: coefficients(:), speed
    REAL(wp) :: value, rate

    rate = velocity_decay_rate(flow_model, coefficients, speed)
    IF (.NOT. ieee_is_finite(rate)) THEN
       CALL fatal_analytic_input_error('Nonfinite source rate in exact velocity integrator')
    END IF
    IF (rate .LT. -TINY(1.0_wp)) THEN
       CALL fatal_analytic_input_error('Negative decay rate in exact velocity integrator')
    END IF
    IF (rate .LE. TINY(1.0_wp)) THEN
       value = HUGE(1.0_wp)
    ELSE
       value = 1.0_wp/rate
    END IF

  END FUNCTION inverse_rate

  PURE FUNCTION simpson_estimate(lower, upper, f_lower, f_middle, f_upper)  &
       RESULT(value)

    REAL(wp), INTENT(IN) :: lower, upper, f_lower, f_middle, f_upper
    REAL(wp) :: value

    value = (upper-lower)*(f_lower+4.0_wp*f_middle+f_upper)/6.0_wp

  END FUNCTION simpson_estimate

  !******************************************************************************
  !> \brief Select the largest real root of a quadratic in a closed interval.
  !> \param[in] a,b,c Polynomial coefficients for a*x**2+b*x+c=0.
  !> \param[in] lower,upper Closed interval for admissible roots.
  !> \param[out] root Largest admissible root.
  !> \param[out] found True when at least one admissible real root exists.
  !******************************************************************************
  PURE SUBROUTINE solve_quadratic_in_interval(a, b, c, lower, upper, root, found)

    REAL(wp), INTENT(IN) :: a, b, c, lower, upper
    REAL(wp), INTENT(OUT) :: root
    LOGICAL, INTENT(OUT) :: found

    REAL(wp) :: discriminant, discriminant_tolerance, sqrt_discriminant
    REAL(wp) :: q, roots(2), interval_tolerance
    INTEGER :: i

    root = 0.0_wp
    found = .FALSE.
    interval_tolerance = 32.0_wp * EPSILON(1.0_wp) *                       &
         MAX(1.0_wp, ABS(lower), ABS(upper))

    IF (a .EQ. 0.0_wp) THEN
       IF (b .EQ. 0.0_wp) THEN
          IF (c .EQ. 0.0_wp) THEN
             root = upper
             found = .TRUE.
          END IF
          RETURN
       END IF
       roots(1) = -c / b
       roots(2) = roots(1)
    ELSE
       discriminant = b*b - 4.0_wp*a*c
       discriminant_tolerance = 64.0_wp * EPSILON(1.0_wp) *                &
            MAX(TINY(1.0_wp), ABS(b*b), ABS(4.0_wp*a*c))
       IF (discriminant .LT. -discriminant_tolerance) RETURN
       discriminant = MAX(0.0_wp, discriminant)

       sqrt_discriminant = SQRT(discriminant)
       q = -0.5_wp * (b + SIGN(sqrt_discriminant, b))
       IF (q .EQ. 0.0_wp) THEN
          roots = -b / (2.0_wp*a)
       ELSE
          roots(1) = q / a
          roots(2) = c / q
       END IF
    END IF

    DO i = 1, 2
       IF (roots(i) .LT. lower - interval_tolerance) CYCLE
       IF (roots(i) .GT. upper + interval_tolerance) CYCLE
       IF (.NOT. found .OR. roots(i) .GT. root) THEN
          root = MIN(upper, MAX(lower, roots(i)))
          found = .TRUE.
       END IF
    END DO

  END SUBROUTINE solve_quadratic_in_interval

  !******************************************************************************
  !> \brief Select the largest real root of a cubic in a closed interval.
  !> \param[in] a,b,c,d Polynomial coefficients for a*x**3+b*x**2+c*x+d=0.
  !> \param[in] lower,upper Closed interval for admissible roots.
  !> \param[out] root Largest admissible root.
  !> \param[out] found True when at least one admissible real root exists.
  !******************************************************************************
  PURE SUBROUTINE solve_cubic_in_interval(a, b, c, d, lower, upper, root, found)

    REAL(wp), INTENT(IN) :: a, b, c, d, lower, upper
    REAL(wp), INTENT(OUT) :: root
    LOGICAL, INTENT(OUT) :: found

    REAL(wp) :: quadratic_root, discriminant, discriminant_tolerance
    REAL(wp) :: p, q, shift, angle, radius, u, roots(3)
    REAL(wp) :: interval_tolerance
    LOGICAL :: quadratic_found
    INTEGER :: i

    root = 0.0_wp
    found = .FALSE.
    IF (a .EQ. 0.0_wp) THEN
       CALL solve_quadratic_in_interval(b, c, d, lower, upper,               &
            quadratic_root, quadratic_found)
       root = quadratic_root
       found = quadratic_found
       RETURN
    END IF

    shift = b / (3.0_wp*a)
    p = c / a - 3.0_wp*shift**2
    q = 2.0_wp*shift**3 - shift*(c/a) + d/a
    discriminant = 0.25_wp*q*q + (p/3.0_wp)**3
    discriminant_tolerance = 64.0_wp * EPSILON(1.0_wp) *                   &
         MAX(TINY(1.0_wp), ABS(0.25_wp*q*q), ABS((p/3.0_wp)**3))
    interval_tolerance = 32.0_wp * EPSILON(1.0_wp) *                       &
         MAX(1.0_wp, ABS(lower), ABS(upper))

    IF (discriminant .GT. discriminant_tolerance) THEN
       roots(1) = SIGN(ABS(-0.5_wp*q + SQRT(discriminant))**(1.0_wp/3.0_wp), &
            -0.5_wp*q + SQRT(discriminant)) +                              &
            SIGN(ABS(-0.5_wp*q - SQRT(discriminant))**(1.0_wp/3.0_wp),      &
            -0.5_wp*q - SQRT(discriminant)) - shift
       roots(2:3) = roots(1)
    ELSEIF (discriminant .GE. -discriminant_tolerance) THEN
       u = SIGN(ABS(-0.5_wp*q)**(1.0_wp/3.0_wp), -0.5_wp*q)
       roots(1) = 2.0_wp*u - shift
       roots(2) = -u - shift
       roots(3) = roots(2)
    ELSE
       radius = 2.0_wp * SQRT(-p/3.0_wp)
       angle = ACOS(MAX(-1.0_wp, MIN(1.0_wp,                              &
            (3.0_wp*q/(2.0_wp*p))*SQRT(-3.0_wp/p))))
       DO i = 1, 3
          roots(i) = radius * COS((angle - 2.0_wp*ACOS(-1.0_wp)*           &
               REAL(i-1,wp))/3.0_wp) - shift
       END DO
    END IF

    DO i = 1, 3
       IF (roots(i) .LT. lower - interval_tolerance) CYCLE
       IF (roots(i) .GT. upper + interval_tolerance) CYCLE
       IF (.NOT. found .OR. roots(i) .GT. root) THEN
          root = MIN(upper, MAX(lower, roots(i)))
          found = .TRUE.
       END IF
    END DO

  END SUBROUTINE solve_cubic_in_interval

  !******************************************************************************
  !> \brief Evaluate the principal real branch of Lambert W.
  !> \details Uses Halley's iteration to evaluate the special function defined
  !> by \f$W\exp(W)=x\f$ for \f$x\ge-1/e\f$.
  !> \param[in] argument Real argument.
  !> \param[out] value Principal real value \f$W_0(x)\f$.
  !> \param[out] success True when the real branch was evaluated.
  !******************************************************************************
  PURE SUBROUTINE lambert_w0(argument, value, success)

    REAL(wp), INTENT(IN) :: argument
    REAL(wp), INTENT(OUT) :: value
    LOGICAL, INTENT(OUT) :: success

    REAL(wp) :: branch_point, residual, exponential_value, denominator
    REAL(wp) :: correction, tolerance, p
    INTEGER :: iteration

    branch_point = -1.0_wp / EXP(1.0_wp)
    tolerance = 8.0_wp * EPSILON(1.0_wp)
    value = 0.0_wp
    success = .FALSE.

    IF (argument .LT. branch_point - tolerance) RETURN
    IF (argument .EQ. 0.0_wp) THEN
       success = .TRUE.
       RETURN
    END IF
    IF (argument .LE. branch_point + tolerance) THEN
       value = -1.0_wp
       success = .TRUE.
       RETURN
    END IF

    IF (argument .LT. -0.1_wp) THEN
       p = SQRT(2.0_wp * (EXP(1.0_wp)*argument + 1.0_wp))
       value = -1.0_wp + p - p*p/3.0_wp + 11.0_wp*p**3/72.0_wp
    ELSEIF (argument .LT. 1.0_wp) THEN
       value = LOG(1.0_wp + argument)
    ELSEIF (argument .LT. 3.0_wp) THEN
       value = 0.5_wp
    ELSE
       value = LOG(argument) - LOG(LOG(argument))
    END IF

    DO iteration = 1, 100
       exponential_value = EXP(value)
       residual = value*exponential_value - argument
       denominator = exponential_value*(value + 1.0_wp) -                &
            (value + 2.0_wp)*residual/(2.0_wp*value + 2.0_wp)
       IF (denominator .EQ. 0.0_wp) RETURN

       correction = residual / denominator
       value = value - correction
       IF (ABS(correction) .LE. tolerance * MAX(1.0_wp, ABS(value))) THEN
          success = .TRUE.
          RETURN
       END IF
    END DO

  END SUBROUTINE lambert_w0

END MODULE rheology_analytic_utils_2d
