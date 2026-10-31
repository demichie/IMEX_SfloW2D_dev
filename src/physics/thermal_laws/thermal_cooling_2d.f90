!> Radiative and convective cooling source model and exact integrator.
MODULE thermal_cooling_2d
  USE parameters_2d, ONLY : wp
  USE diagnostics_2d, ONLY : fatal_error
  USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_thermal_cooling_explicit, eval_thermal_cooling_semi_implicit
  PUBLIC :: eval_thermal_cooling_implicit, eval_thermal_cooling_analytic
CONTAINS

  SUBROUTINE eval_thermal_cooling_explicit(state,radiative,convective,ambient,rate)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: radiative,convective,ambient
    REAL(wp), INTENT(OUT) :: rate
    rate=-radiative*(state%T**4-ambient**4)-convective*(state%T-ambient)
  END SUBROUTINE eval_thermal_cooling_explicit


  SUBROUTINE eval_thermal_cooling_semi_implicit(state,rate)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(OUT) :: rate
    rate=0.0_wp
  END SUBROUTINE eval_thermal_cooling_semi_implicit


  SUBROUTINE eval_thermal_cooling_implicit(state,radiative,convective,ambient,rate)
    TYPE(complex_physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: radiative,convective,ambient
    COMPLEX(wp), INTENT(OUT) :: rate
    rate=-radiative*(state%T**4-ambient**4)-convective*(state%T-ambient)
  END SUBROUTINE eval_thermal_cooling_implicit


  SUBROUTINE eval_thermal_cooling_analytic(state,interval,heat_capacity, &
       radiative,convective,ambient,temperature_after)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: interval,heat_capacity,radiative,convective,ambient
    REAL(wp), INTENT(OUT) :: temperature_after
    REAL(wp) :: temperature
    REAL(wp) :: beta, cubic_root, lower, upper, midpoint
    REAL(wp) :: initial_primitive, target_primitive, midpoint_primitive
    REAL(wp) :: elapsed_scale, tolerance
    INTEGER :: iteration


    temperature = state%T

    IF (.NOT. ieee_is_finite(temperature) .OR.                              &
         .NOT. ieee_is_finite(interval) .OR.                        &
         .NOT. ieee_is_finite(heat_capacity) .OR.                    &
         .NOT. ieee_is_finite(radiative) .OR.                  &
         .NOT. ieee_is_finite(convective) .OR.                 &
         .NOT. ieee_is_finite(ambient)) THEN
       CALL fatal_error('Nonfinite input to exact lava cooling integrator')
    END IF
    IF (temperature .LE. 0.0_wp .OR. ambient .LE. 0.0_wp) THEN
       CALL fatal_error('Nonpositive temperature in exact lava cooling integrator')
    END IF
    IF (interval .LT. 0.0_wp .OR. heat_capacity .LE. 0.0_wp) THEN
       CALL fatal_error('Invalid time or heat capacity in exact lava cooling integrator')
    END IF
    IF (radiative .LT. 0.0_wp .OR.                            &
         convective .LT. 0.0_wp) THEN
       CALL fatal_error('Negative cooling coefficient in exact lava cooling integrator')
    END IF

    temperature_after = temperature
    IF (interval .EQ. 0.0_wp .OR.                                   &
         (radiative .EQ. 0.0_wp .AND.                          &
         convective .EQ. 0.0_wp) .OR.                          &
         temperature .EQ. ambient) RETURN

    IF (radiative .EQ. 0.0_wp) THEN
       temperature_after = ambient + (temperature -             &
            ambient) * EXP(-convective *          &
            interval / heat_capacity)
       RETURN
    END IF

    beta = convective / (radiative *               &
         ambient**3)
    IF (.NOT. ieee_is_finite(beta)) THEN
       CALL fatal_error('Cooling coefficient ratio overflow in lava integrator')
    END IF
    CALL lava_cooling_cubic_root(beta, cubic_root)

    initial_primitive = lava_cooling_time_primitive(                        &
         temperature / ambient, beta, cubic_root)
    elapsed_scale = radiative * ambient**3 *        &
         interval / heat_capacity
    IF (.NOT. ieee_is_finite(elapsed_scale)) THEN
       temperature_after = ambient
       RETURN
    END IF
    target_primitive = initial_primitive - elapsed_scale

    IF (temperature .GT. ambient) THEN
       lower = 1.0_wp
       upper = temperature / ambient
    ELSE
       lower = temperature / ambient
       upper = 1.0_wp
    END IF

    tolerance = 8.0_wp * EPSILON(1.0_wp) *                                  &
         MAX(1.0_wp, ABS(lower), ABS(upper))
    DO iteration = 1, 256
       midpoint = lower + 0.5_wp * (upper - lower)
       IF (midpoint .EQ. lower .OR. midpoint .EQ. upper .OR.                &
            ABS(upper - lower) .LE. tolerance) EXIT
       IF (midpoint .EQ. 1.0_wp) THEN
          IF (temperature .GT. ambient) THEN
             lower = midpoint
          ELSE
             upper = midpoint
          END IF
          CYCLE
       END IF

       midpoint_primitive = lava_cooling_time_primitive(midpoint, beta,       &
            cubic_root)
       IF (temperature .GT. ambient) THEN
          IF (midpoint_primitive .LT. target_primitive) THEN
             lower = midpoint
          ELSE
             upper = midpoint
          END IF
       ELSE
          IF (midpoint_primitive .GT. target_primitive) THEN
             lower = midpoint
          ELSE
             upper = midpoint
          END IF
       END IF
    END DO
    IF (iteration .GT. 256) THEN
       CALL fatal_error('Exact lava cooling inversion failed to converge')
    END IF

    midpoint = lower + 0.5_wp * (upper - lower)
    temperature_after = ambient * midpoint
    IF (temperature .GT. ambient) THEN
       temperature_after = MAX(ambient,                        &
            MIN(temperature, temperature_after))
    ELSE
       temperature_after = MIN(ambient,                        &
            MAX(temperature, temperature_after))
    END IF

    END SUBROUTINE eval_thermal_cooling_analytic





SUBROUTINE lava_cooling_cubic_root(beta, root)

    REAL(wp), INTENT(IN) :: beta
    REAL(wp), INTENT(OUT) :: root

    REAL(wp) :: lower, upper, midpoint, scale, polynomial
    INTEGER :: iteration

    IF (beta .EQ. 0.0_wp) THEN
       root = -1.0_wp
       RETURN
    END IF

    scale = MAX(1.0_wp, (1.0_wp + beta)**(1.0_wp / 3.0_wp))
    lower = -2.0_wp
    upper = -1.0_wp / scale
    DO iteration = 1, 256
       midpoint = lower + 0.5_wp * (upper - lower)
       IF (midpoint .EQ. lower .OR. midpoint .EQ. upper) EXIT
       polynomial = midpoint**3 + midpoint**2/scale + midpoint/scale**2 +  &
            1.0_wp
       IF (polynomial .LT. 0.0_wp) THEN
          lower = midpoint
       ELSE
          upper = midpoint
       END IF
    END DO
    IF (iteration .GT. 256) THEN
       CALL fatal_error('Cooling cubic root calculation failed to converge')
    END IF
    root = scale * (lower + 0.5_wp * (upper - lower))

  END SUBROUTINE lava_cooling_cubic_root

PURE FUNCTION lava_cooling_time_primitive(x, beta, cubic_root) RESULT(value)

    REAL(wp), INTENT(IN) :: x, beta, cubic_root
    REAL(wp) :: value

    REAL(wp) :: quadratic_linear, quadratic_constant, discriminant
    REAL(wp) :: partial_a, partial_b, partial_c, partial_d, residual

    quadratic_linear = 1.0_wp + cubic_root
    quadratic_constant = 1.0_wp + cubic_root + cubic_root**2
    discriminant = 4.0_wp * quadratic_constant - quadratic_linear**2

    partial_a = 1.0_wp / (4.0_wp + beta)
    partial_b = 1.0_wp / ((cubic_root - 1.0_wp) *                        &
         (3.0_wp * cubic_root**2 + 2.0_wp * cubic_root + 1.0_wp))
    partial_c = -partial_a - partial_b
    partial_d = -partial_a * (cubic_root + 2.0_wp) -                     &
         partial_b * (2.0_wp * cubic_root + 1.0_wp)
    residual = partial_d - 0.5_wp * partial_c * quadratic_linear

    value = partial_a * LOG(ABS(x - 1.0_wp)) +                            &
         partial_b * LOG(ABS(x - cubic_root)) +                           &
         0.5_wp * partial_c * LOG(x**2 + quadratic_linear*x +             &
         quadratic_constant) + 2.0_wp * residual / SQRT(discriminant) *   &
         ATAN((2.0_wp*x + quadratic_linear) / SQRT(discriminant))

  END FUNCTION lava_cooling_time_primitive

  ! All public model treatments take a physical state; scalar calculations are private.



END MODULE thermal_cooling_2d
