!> \brief Per-step source-window coefficients shared by injection laws.
!> \details Updated before the RK stages by the existing window initializer.
MODULE mass_injection_coefficients_2d
  USE parameters_2d, ONLY : wp
  IMPLICIT NONE
  REAL(wp) :: bottom_coeff_step = -1.0_wp
  REAL(wp) :: bottom_vel_step = -1.0_wp
  REAL(wp) :: lateral_coeff_step = -1.0_wp
  REAL(wp), ALLOCATABLE :: fissure_coeff_step(:)
CONTAINS
  SUBROUTINE set_source_time_coefficients( time , dt )

    USE parameters_2d, ONLY : time_param , bottom_radial_source_flag ,          &
         radial_source_flag , n_intervals , bottom_fissural_source_flag ,      &
         n_fissures , time_param_fissures

    REAL(wp), INTENT(IN) :: time
    REAL(wp), INTENT(IN) :: dt

    INTEGER :: i_fissure

    IF ( dt .LE. 0.0_wp ) RETURN

    IF ( bottom_radial_source_flag ) THEN

       IF ( n_intervals .GT. 0 ) THEN

          bottom_vel_step = interval_mean_velocity( time , dt )

       ELSE

          bottom_coeff_step = window_mean( time + time_param(4) , dt ,         &
               time_param(1) , time_param(2) , time_param(3) , .FALSE. ,      &
               HUGE(1.0_wp) )

       END IF

    END IF

    IF ( bottom_fissural_source_flag ) THEN

       IF ( .NOT. ALLOCATED( fissure_coeff_step ) )                            &
            ALLOCATE( fissure_coeff_step(n_fissures) )

       DO i_fissure = 1, n_fissures

          fissure_coeff_step(i_fissure) = window_mean( time +                 &
               time_param_fissures(4,i_fissure) , dt ,                        &
               time_param_fissures(1,i_fissure) ,                             &
               time_param_fissures(2,i_fissure) ,                             &
               time_param_fissures(3,i_fissure) , .FALSE. , HUGE(1.0_wp) )

       END DO

    END IF

    ! The lateral radial source has cosine ramps and stops at TIME_PARAM(4)
    IF ( radial_source_flag ) lateral_coeff_step = window_mean( time , dt ,     &
         time_param(1) , time_param(2) , time_param(3) , .TRUE. ,              &
         time_param(4) )

  END SUBROUTINE set_source_time_coefficients

  !******************************************************************************
  !> \brief Mean of a periodic source window over [v0, v0+dt].
  !>
  !> In each period the coefficient rises over the ramp (linearly, or as a half
  !> cosine), is 1 until act-ramp, falls to 0 at act and stays 0 to the end of
  !> the period; it is 0 from cutoff on. With ramp = 0 the window is on for
  !> 0 <= u <= act, as in the point evaluation.
  !>
  !> \param[in] v0 Start of the step, phase included [s].
  !> \param[in] dt Length of the step [s].
  !> \param[in] per Period [s].
  !> \param[in] act Active time in each period [s].
  !> \param[in] ramp Ramp time [s].
  !> \param[in] cosine Cosine ramps when true, linear ramps otherwise.
  !> \param[in] cutoff Time from which the window is 0 [s].
  !> \return Mean coefficient over the step.
  !> \note Assumes per > 0 and dt > 0; the result is clipped to [0,1].
  !******************************************************************************

  REAL(wp) FUNCTION window_mean( v0 , dt , per , act , ramp , cosine , cutoff )

    REAL(wp), INTENT(IN) :: v0 , dt , per , act , ramp , cutoff
    LOGICAL, INTENT(IN) :: cosine

    REAL(wp) :: v1 , n0 , n1 , u0 , u1 , total

    v1 = MIN( v0 + dt , cutoff )

    IF ( v1 .LE. v0 ) THEN

       window_mean = 0.0_wp
       RETURN

    END IF

    ! A window with no ramp that is active over the whole period is always on
    IF ( ( ramp .LE. 0.0_wp ) .AND. ( act .GE. per ) .AND.                    &
         ( v0 + dt .LE. cutoff ) ) THEN

       window_mean = 1.0_wp
       RETURN

    END IF

    n0 = AINT( v0 / per )
    n1 = AINT( v1 / per )
    u0 = v0 - n0 * per
    u1 = v1 - n1 * per

    IF ( n1 .EQ. n0 ) THEN

       ! A step that stays on the plateau gets exactly 1
       IF ( ( u0 .GE. ramp ) .AND. ( u1 .LE. act - ramp ) .AND.                &
            ( v0 + dt .LE. cutoff ) ) THEN

          window_mean = 1.0_wp
          RETURN

       END IF

       total = window_integral( u0 , u1 , act , ramp , cosine )

    ELSE

       total = window_integral( u0 , per , act , ramp , cosine )               &
            + ( n1 - n0 - 1.0_wp ) * ( act - ramp )                            &
            + window_integral( 0.0_wp , u1 , act , ramp , cosine )

    END IF

    window_mean = MIN( MAX( total / dt , 0.0_wp ) , 1.0_wp )

  END FUNCTION window_mean

  !******************************************************************************
  !> \brief Integral of one period of a source window between a and b.
  !>
  !> Each segment is integrated in closed form, with the cosine ramps written as
  !> sin/cos products so that short intervals lose no precision.
  !>
  !> \param[in] a Lower limit, 0 <= a <= b <= period [s].
  !> \param[in] b Upper limit [s].
  !> \param[in] act Active time [s].
  !> \param[in] ramp Ramp time [s].
  !> \param[in] cosine Cosine ramps when true, linear ramps otherwise.
  !> \return Integral of the coefficient [s].
  !> \note The limits are expected to satisfy 0 <= a <= b. A nonpositive ramp selects
  !>       a rectangular window.
  !******************************************************************************

  REAL(wp) FUNCTION window_integral( a , b , act , ramp , cosine )

    USE geometry_2d, ONLY : pi_g

    REAL(wp), INTENT(IN) :: a , b , act , ramp
    LOGICAL, INTENT(IN) :: cosine

    REAL(wp) :: x0 , x1

    window_integral = 0.0_wp

    IF ( ramp .LE. 0.0_wp ) THEN

       x1 = MIN( b , act )
       IF ( x1 .GT. a ) window_integral = x1 - a
       RETURN

    END IF

    ! ramp up
    x0 = MAX( a , 0.0_wp )
    x1 = MIN( b , ramp )

    IF ( x1 .GT. x0 ) THEN

       IF ( cosine ) THEN
          window_integral = window_integral + 0.5_wp * ( x1 - x0 )             &
               - ramp / pi_g * COS( 0.5_wp * pi_g * ( x0 + x1 ) / ramp )      &
               * SIN( 0.5_wp * pi_g * ( x1 - x0 ) / ramp )
       ELSE
          window_integral = window_integral                                    &
               + 0.5_wp * ( x1 - x0 ) * ( x1 + x0 ) / ramp
       END IF

    END IF

    ! plateau
    x0 = MAX( a , ramp )
    x1 = MIN( b , act - ramp )
    IF ( x1 .GT. x0 ) window_integral = window_integral + ( x1 - x0 )

    ! ramp down
    x0 = MAX( a , act - ramp )
    x1 = MIN( b , act )

    IF ( x1 .GT. x0 ) THEN

       IF ( cosine ) THEN
          window_integral = window_integral + 0.5_wp * ( x1 - x0 )             &
               + ramp / pi_g * COS( 0.5_wp * pi_g * ( x0 + x1 - 2.0_wp *       &
               ( act - ramp ) ) / ramp ) * SIN( 0.5_wp * pi_g * ( x1 - x0 )    &
               / ramp )
       ELSE
          window_integral = window_integral                                    &
               + 0.5_wp * ( x1 - x0 ) * ( 2.0_wp * act - x0 - x1 ) / ramp
       END IF

    END IF

  END FUNCTION window_integral

  !******************************************************************************
  !> \brief Mean velocity of the bottom radial source intervals over a step.
  !>
  !> Interval i holds VEL_INTERVALS(i) from T_INTERVALS(i) to T_INTERVALS(i+1),
  !> the last one to the end of the run; before T_INTERVALS(1) the source is
  !> off. A step inside one interval gets its velocity exactly.
  !>
  !> \param[in] t0 Start of the step [s].
  !> \param[in] dt Length of the step [s].
  !> \return Mean velocity over the step [m s^-1].
  !> \note Uses the configured interval start times and velocities; a positive dt is required.
  !******************************************************************************

  REAL(wp) FUNCTION interval_mean_velocity( t0 , dt )

    USE parameters_2d, ONLY : n_intervals , t_intervals , vel_intervals

    REAL(wp), INTENT(IN) :: t0 , dt

    REAL(wp) :: t1 , x0 , x1 , total
    INTEGER :: i_int

    t1 = t0 + dt
    total = 0.0_wp

    DO i_int = 1, n_intervals

       x0 = MAX( t0 , t_intervals(i_int) )

       IF ( i_int .LT. n_intervals ) THEN
          x1 = MIN( t1 , t_intervals(i_int+1) )
       ELSE
          x1 = t1
       END IF

       IF ( x1 .GT. x0 ) THEN

          IF ( ( x0 .EQ. t0 ) .AND. ( x1 .EQ. t1 ) ) THEN
             interval_mean_velocity = vel_intervals(i_int)
             RETURN
          END IF

          total = total + vel_intervals(i_int) * ( x1 - x0 )

       END IF

    END DO

    interval_mean_velocity = total / dt

  END FUNCTION interval_mean_velocity
END MODULE mass_injection_coefficients_2d
