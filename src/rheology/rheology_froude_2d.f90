!********************************************************************************
!> \brief Froude-dependent Coulomb basal-friction model (model 9).
!> \details The friction coefficient weakens exponentially with the Froude
!> number, \f$\mathrm{Fr}=|\boldsymbol{u}|/\sqrt{g' h}\f$:
!> \f[
!> \mu(\mathrm{Fr})=\mu_\infty+(\mu_0-\mu_\infty)
!>                 \exp(-\mathrm{Fr}/\mathrm{Fr}_0).
!> \f]
!> The horizontal source is
!> \f$\boldsymbol{S}=-\rho_m\mu h G(g'+a_c)
!> \boldsymbol{u}_h/|\boldsymbol{u}_h|\f$, where \f$G\f$ is the supplied
!> gravity coefficient and \f$a_c\f$ is the curvature acceleration. Optional
!> stochastic perturbations shift the Froude number for sufficiently deep
!> flow and clamp it at zero. The module provides
!> \ref eval_froude_source to evaluate the two horizontal source components.
!> \ref solve_froude_analytic_velocity solves the velocity DIRK stage on the
!> exponential branch with the principal real Lambert W function.
!********************************************************************************

MODULE rheology_froude_2d

  USE parameters_2d, ONLY : wp
  USE rheology_parameters_2d, ONLY : mu_0, mu_inf, Fr_0
  USE rheology_analytic_utils_2d, ONLY : lambert_w0,                        &
       validate_analytic_velocity_inputs, fatal_analytic_input_error

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_froude_source
  PUBLIC :: solve_froude_analytic_velocity

CONTAINS

  !******************************************************************************
  !> \brief Solve the Froude-dependent DIRK stage analytically for velocity.
  !> \details Holds depth, gravity, curvature, and the speed ratio fixed. On
  !> the weakening branch, the scalar stage equation is solved with the
  !> principal real Lambert W function; the zero-Froude clamp is checked
  !> separately. All inputs other than the two horizontal velocity components
  !> are explicit and frozen over the stage.
  !> \param integration_step Effective DIRK integration step (dt * a_ii).
  !> \param flow_depth Explicit flow depth.
  !> \param speed_ratio Explicit total-to-horizontal speed ratio.
  !> \param velocity_before Input horizontal velocity.
  !> \param reduced_gravity Explicit reduced gravity.
  !> \param gravity_coefficient Explicit gravity multiplier.
  !> \param curvature_term Explicit curvature acceleration.
  !> \param stochastic Enable the explicit Froude shift.
  !> \param stochastic_value Explicit stochastic increment.
  !> \param velocity_after Returned horizontal velocity.
  !******************************************************************************
  SUBROUTINE solve_froude_analytic_velocity(integration_step, flow_depth,  &
       speed_ratio, velocity_before, reduced_gravity, gravity_coefficient,   &
       curvature_term, stochastic, stochastic_value, velocity_after)

    REAL(wp), INTENT(IN) :: integration_step, flow_depth, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2), reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient, curvature_term
    LOGICAL, INTENT(IN) :: stochastic
    REAL(wp), INTENT(IN) :: stochastic_value
    REAL(wp), INTENT(OUT) :: velocity_after(2)

    REAL(wp) :: speed_before, speed_after, shift, base_drag, scale
    REAL(wp) :: froude_slope, froude_value_slope, constant_branch_speed
    REAL(wp) :: remaining_speed, weakening_amplitude, argument, w_value
    REAL(wp) :: exponent, lower_branch_limit
    LOGICAL :: w_found, valid_root

    velocity_after = velocity_before
    speed_before = SQRT(SUM(velocity_before**2))
    CALL validate_analytic_velocity_inputs(integration_step, flow_depth,    &
         speed_ratio, speed_before, 'Froude')
    IF (integration_step .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
    IF (flow_depth .LE. 0.0_wp .OR. reduced_gravity .LE. 0.0_wp .OR.       &
         speed_ratio .LE. 0.0_wp) RETURN
    IF (Fr_0 .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive Fr_0 in Froude analytic solver')
    END IF

    base_drag = gravity_coefficient * (reduced_gravity + curvature_term)
    scale = integration_step * base_drag
    IF (scale .LE. 0.0_wp) RETURN

    IF (mu_0 .EQ. mu_inf) THEN
       speed_after = MAX(0.0_wp, speed_before - scale*mu_0)
       velocity_after = (speed_after / speed_before) * velocity_before
       RETURN
    END IF

    shift = 0.0_wp
    IF (stochastic .AND. flow_depth .GT. 0.01_wp) shift = stochastic_value
    froude_value_slope = speed_ratio / SQRT(reduced_gravity*flow_depth)
    froude_slope = froude_value_slope / Fr_0
    constant_branch_speed = speed_before - scale*mu_0
    valid_root = constant_branch_speed .GE. 0.0_wp .AND.                   &
         (froude_value_slope*constant_branch_speed + shift) .LE. 0.0_wp

    IF (valid_root) THEN
       speed_after = constant_branch_speed
    ELSE
       remaining_speed = speed_before - scale*mu_inf
       weakening_amplitude = scale * (mu_0 - mu_inf) * EXP(-shift/Fr_0)
       IF (weakening_amplitude .EQ. 0.0_wp) THEN
          speed_after = remaining_speed
          valid_root = speed_after .GE. 0.0_wp .AND.                       &
               (froude_value_slope*speed_after + shift) .GE. 0.0_wp
       ELSE
          exponent = -froude_slope * remaining_speed
          lower_branch_limit = LOG(HUGE(1.0_wp)) -                         &
               LOG(ABS(froude_slope*weakening_amplitude))
          IF (exponent .GT. lower_branch_limit) THEN
             IF (weakening_amplitude .LT. 0.0_wp) THEN
                CALL fatal_analytic_input_error(                            &
                     'Lambert W argument overflow in Froude analytic solver')
             END IF
             valid_root = .FALSE.
          ELSE
             argument = -froude_slope * weakening_amplitude * EXP(exponent)
             IF (argument .LT. -1.0_wp/EXP(1.0_wp) -                       &
                  8.0_wp*EPSILON(1.0_wp)) THEN
                valid_root = .FALSE.
             ELSE
                CALL lambert_w0(argument, w_value, w_found)
                IF (.NOT. w_found) THEN
                   CALL fatal_analytic_input_error(                        &
                        'Lambert W evaluation failed in Froude analytic solver')
                END IF
                speed_after = remaining_speed + w_value/froude_slope
                valid_root = speed_after .GE. 0.0_wp .AND.                 &
                     speed_after .LE. speed_before .AND.                   &
                     (froude_value_slope*speed_after + shift) .GE. 0.0_wp
             END IF
          END IF
       END IF

       IF (.NOT. valid_root) speed_after = 0.0_wp
    END IF

    speed_after = MIN(speed_before, MAX(0.0_wp, speed_after))
    velocity_after = (speed_after / speed_before) * velocity_before

  END SUBROUTINE solve_froude_analytic_velocity

  !******************************************************************************
  !> \brief Evaluate model 9's Froude-dependent Coulomb source.
  !> \details Computes the Froude number and the exponentially weakening
  !> friction coefficient, then returns the x/y components of the basal
  !> friction source. A stochastic perturbation is applied only when enabled;
  !> the source is zero unless speed and reduced gravity are positive.
  !> \param flow_depth Flow depth \f$h\f$.
  !> \param velocity_x x component of horizontal velocity.
  !> \param velocity_y y component of horizontal velocity.
  !> \param speed Magnitude used in the Froude number.
  !> \param horizontal_speed Magnitude of horizontal velocity.
  !> \param mixture_density Mixture density \f$\rho_m\f$.
  !> \param reduced_gravity Reduced gravity \f$g'\f$.
  !> \param gravity_coefficient Gravity multiplier \f$G\f$.
  !> \param curvature_term Curvature acceleration \f$a_c\f$.
  !> \param stochastic Enable the stochastic Froude perturbation.
  !> \param stochastic_value Stochastic increment to the Froude number.
  !> \param source_x Returned x momentum source.
  !> \param source_y Returned y momentum source.
  !******************************************************************************

  SUBROUTINE eval_froude_source( flow_depth, velocity_x, velocity_y, speed,  &
       horizontal_speed, mixture_density, reduced_gravity, gravity_coefficient, &
       curvature_term, stochastic, stochastic_value, source_x, source_y )

    REAL(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y, speed
    REAL(wp), INTENT(IN) :: horizontal_speed, mixture_density, reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient, curvature_term
    LOGICAL, INTENT(IN) :: stochastic
    REAL(wp), INTENT(IN) :: stochastic_value
    REAL(wp), INTENT(OUT) :: source_x, source_y

    REAL(wp) :: froude_number
    REAL(wp) :: friction_coefficient
    REAL(wp) :: friction_stress

    source_x = 0.0_wp
    source_y = 0.0_wp

    IF ( speed .GT. 0.0_wp .AND. reduced_gravity .GT. 0.0_wp ) THEN
       froude_number = speed / SQRT(reduced_gravity * flow_depth)

       IF ( stochastic ) THEN
          IF ( flow_depth .GT. 0.01_wp ) froude_number = froude_number + stochastic_value
          IF ( froude_number .LT. 0.0_wp ) froude_number = 0.0_wp
       END IF

       friction_coefficient = mu_inf + (mu_0 - mu_inf) * EXP(-froude_number/Fr_0)
       friction_stress = mixture_density * ( friction_coefficient * flow_depth * &
            gravity_coefficient * ( reduced_gravity + curvature_term ) )
       source_x = source_x - friction_stress * velocity_x / horizontal_speed
       source_y = source_y - friction_stress * velocity_y / horizontal_speed
    END IF

  END SUBROUTINE eval_froude_source

END MODULE rheology_froude_2d
