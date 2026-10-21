!********************************************************************************
!> \brief Scalar closed-form helpers for analytic rheology stage solutions.
!> \details Provides real roots of quadratic and cubic stage equations over a
!> bounded physical interval, plus evaluation of the principal real Lambert W
!> function required by the Froude-weakening law.
!********************************************************************************

MODULE rheology_analytic_utils_2d

  USE parameters_2d, ONLY : wp
  USE diagnostics_2d, ONLY : fatal_error

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: solve_quadratic_in_interval
  PUBLIC :: solve_cubic_in_interval
  PUBLIC :: lambert_w0
  PUBLIC :: mui_velocity_stage_polynomial
  PUBLIC :: validate_analytic_velocity_inputs
  PUBLIC :: fatal_analytic_input_error

CONTAINS

  !******************************************************************************
  !> \brief Report a fatal error raised by an analytic rheology stage routine.
  !******************************************************************************
  SUBROUTINE fatal_analytic_input_error(message)

    CHARACTER(LEN=*), INTENT(IN) :: message

    CALL fatal_error(message)

  END SUBROUTINE fatal_analytic_input_error

  !******************************************************************************
  !> \brief Validate common inputs to a module-level analytic velocity stage.
  !******************************************************************************
  SUBROUTINE validate_analytic_velocity_inputs(integration_step, flow_depth, &
       speed_ratio, speed_before, rheology_name)

    REAL(wp), INTENT(IN) :: integration_step, flow_depth, speed_ratio
    REAL(wp), INTENT(IN) :: speed_before
    CHARACTER(LEN=*), INTENT(IN) :: rheology_name

    IF (integration_step .LT. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Negative DIRK stage coefficient in '//TRIM(rheology_name)//' analytic solver')
    END IF
    IF (integration_step .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN

    IF (flow_depth .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive flow depth in '//TRIM(rheology_name)//' analytic solver')
    END IF
    IF (speed_ratio .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive speed ratio in '//TRIM(rheology_name)//' analytic solver')
    END IF

  END SUBROUTINE validate_analytic_velocity_inputs

  !******************************************************************************
  !> \brief Form the polynomial for a regularized \f$\mu(I)\f$ stage.
  !> \details The explicit stresses, particle properties, depth, and speed
  !> ratio are frozen. The output coefficients describe the equation
  !> \f$a_2 r^2+a_1 r+a_0=0\f$ after clearing the positive
  !> \f$I_0+I(r)\f$ denominator.
  !> \param integration_step Effective DIRK integration step (dt * a_ii).
  !> \param speed_before Input horizontal speed.
  !> \param flow_depth Actual depth in the conservative momentum.
  !> \param effective_depth Regularized depth used by the constitutive law.
  !> \param speed_ratio Explicit total-to-horizontal speed ratio.
  !> \param mixture_density Explicit mixture density.
  !> \param effective_normal_stress Explicit regularized normal stress.
  !> \param gravity_coefficient Gravity multiplier.
  !> \param particle_diameter Explicit representative particle diameter.
  !> \param particle_density Explicit representative particle density.
  !> \param I_0,mu_s,mu_2,muI_inf Regularized \f$\mu(I)\f$ parameters.
  !> \param quadratic,linear,constant Returned polynomial coefficients.
  !> \param inertial_slope Returned coefficient in \f$I(r)=c_Ir\f$.
  !> \param basal_factor Returned stress-to-momentum acceleration factor.
  !******************************************************************************
  PURE SUBROUTINE mui_velocity_stage_polynomial(integration_step,          &
       speed_before, flow_depth, effective_depth, speed_ratio, mixture_density, &
       effective_normal_stress, gravity_coefficient, particle_diameter,    &
       particle_density, I_0, mu_s, mu_2, muI_inf, quadratic, linear,      &
       constant, inertial_slope, basal_factor)

    REAL(wp), INTENT(IN) :: integration_step, speed_before, flow_depth
    REAL(wp), INTENT(IN) :: effective_depth, speed_ratio, mixture_density
    REAL(wp), INTENT(IN) :: effective_normal_stress, gravity_coefficient
    REAL(wp), INTENT(IN) :: particle_diameter, particle_density
    REAL(wp), INTENT(IN) :: I_0, mu_s, mu_2, muI_inf
    REAL(wp), INTENT(OUT) :: quadratic, linear, constant
    REAL(wp), INTENT(OUT) :: inertial_slope, basal_factor

    REAL(wp) :: sqrt_gravity, sqrt_stress_ratio

    sqrt_gravity = SQRT(gravity_coefficient)
    sqrt_stress_ratio = SQRT(effective_normal_stress / particle_density)
    inertial_slope = 2.5_wp * particle_diameter * speed_ratio /             &
         (effective_depth * sqrt_gravity * sqrt_stress_ratio)
    basal_factor = effective_normal_stress /                                 &
         (mixture_density * flow_depth * speed_ratio)

    quadratic = inertial_slope + integration_step * basal_factor *           &
         muI_inf * inertial_slope**2
    linear = I_0 - inertial_slope * speed_before + integration_step *       &
         basal_factor * mu_2 * inertial_slope
    constant = -I_0 * speed_before + integration_step * basal_factor *      &
         mu_s * I_0

  END SUBROUTINE mui_velocity_stage_polynomial

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
  !> by \f$W\exp(W)=x\f$ for \f$x\ge-1/e\f$. The iteration evaluates W only;
  !> it does not solve a rheology stage equation.
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
