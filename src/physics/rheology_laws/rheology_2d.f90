!********************************************************************************
!> \brief Rheology dispatch, model naming and friction-only momentum correction.
!> \details Dispatches common momentum-rate and analytic-velocity interfaces
!> to the selected constitutive law, resolves canonical model names to IDs,
!> and provides a projection for momentum updates caused exclusively by
!> friction. The projection constrains the candidate update to
!> \f$\boldsymbol{m}_{after}=\alpha\boldsymbol{m}_{before}\f$,
!> \f$0\le\alpha\le1\f$, so it cannot rotate, reverse or increase the
!> momentum. Callers must pass the momentum immediately before and after the
!> isolated friction substep, not the full stage state.
!********************************************************************************

MODULE rheology_2d

  USE parameters_2d, ONLY : wp, rheology_id_none, rheology_id_voellmy,       &
       rheology_id_plastic, rheology_id_lava, rheology_id_lahar,             &
       rheology_id_empirical, rheology_id_quadratic, rheology_id_froude,     &
       rheology_id_lucas, rheology_id_mui, rheology_id_mui_voellmy
  USE diagnostics_2d, ONLY : fatal_error
  USE rheology_voellmy_2d, ONLY : eval_voellmy_implicit_momentum_rate,       &
       eval_voellmy_semi_implicit_momentum_rate,                            &
       solve_voellmy_analytic_velocity
  USE rheology_plastic_2d, ONLY : eval_plastic_implicit_momentum_rate,       &
       solve_plastic_analytic_velocity
  USE rheology_lava_2d, ONLY : eval_lava_implicit_momentum_rate,             &
       eval_lava_semi_implicit_momentum_rate,                               &
       solve_lava_analytic_velocity
  USE rheology_lahar_2d, ONLY : eval_lahar_implicit_momentum_rate,            &
       eval_lahar_semi_implicit_momentum_rate,                              &
       solve_lahar_analytic_velocity
  USE rheology_empirical_2d, ONLY : eval_empirical_implicit_momentum_rate,   &
       solve_empirical_analytic_velocity
  USE rheology_quadratic_2d, ONLY : eval_quadratic_implicit_momentum_rate,   &
       solve_quadratic_analytic_velocity
  USE rheology_froude_2d, ONLY : eval_froude_momentum_rate,                  &
       solve_froude_analytic_velocity
  USE rheology_lucas_2d, ONLY : eval_lucas_momentum_rate,                    &
       solve_lucas_analytic_velocity
  USE rheology_mui_2d, ONLY : eval_mui_momentum_rate,                        &
       solve_mui_analytic_velocity
  USE rheology_mui_voellmy_2d, ONLY : eval_mui_voellmy_momentum_rate,        &
       eval_mui_voellmy_implicit_momentum_rate,                             &
       solve_mui_voellmy_analytic_velocity

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: resolve_rheology_model, rheology_name_from_id
  PUBLIC :: correct_friction_momentum
  PUBLIC :: curvature_acceleration, eval_implicit_rheology_rate
  PUBLIC :: eval_semi_implicit_rheology_rate
  PUBLIC :: solve_analytic_rheology_velocity

CONTAINS

  !******************************************************************************
  !> \brief Contract the filtered bed Hessian with the Cartesian velocity.
  !> \param[in] Bsecond_xx Filtered bed Hessian xx entry [m^-1].
  !> \param[in] Bsecond_xy Filtered bed Hessian xy entry [m^-1].
  !> \param[in] Bsecond_yy Filtered bed Hessian yy entry [m^-1].
  !> \param[in] velocity_x Cartesian depth-averaged x velocity [m s^-1].
  !> \param[in] velocity_y Cartesian depth-averaged y velocity [m s^-1].
  !> \return Signed bed-curvature acceleration [m s^-2].
  !******************************************************************************
  PURE FUNCTION curvature_acceleration(Bsecond_xx, Bsecond_xy, Bsecond_yy,  &
       velocity_x, velocity_y) RESULT(acceleration)

    REAL(wp), INTENT(IN) :: Bsecond_xx, Bsecond_xy, Bsecond_yy
    REAL(wp), INTENT(IN) :: velocity_x, velocity_y
    REAL(wp) :: acceleration

    acceleration = Bsecond_xx*velocity_x**2 +                              &
         2.0_wp*Bsecond_xy*velocity_x*velocity_y + Bsecond_yy*velocity_y**2

  END FUNCTION curvature_acceleration

  !******************************************************************************
  !> \brief Evaluate the selected implicit rheology momentum rate.
  !> \details All state values are explicit arguments so the dispatch has no
  !> per-cell workspace allocation or hidden cell-state conversion.
  !> \param[in] model_id Selected rheology identifier.
  !> \param[in] flow_depth Flow depth [m].
  !> \param[in] temperature Mixture temperature [K].
  !> \param[in] mixture_density Mixture density [kg m^-3].
  !> \param[in] inverse_density Reciprocal mixture density [m^3 kg^-1].
  !> \param[in] velocity_x x velocity [m s^-1].
  !> \param[in] velocity_y y velocity [m s^-1].
  !> \param[in] speed Total flow-speed magnitude [m s^-1].
  !> \param[in] horizontal_speed Horizontal flow-speed magnitude [m s^-1].
  !> \param[in] speed_squared Squared flow-speed magnitude [m^2 s^-2].
  !> \param[in] reduced_gravity Reduced gravitational acceleration [m s^-2].
  !> \param[in] solid_fractions Solid volume fractions, one per class.
  !> \param[in] moving True when the flow is moving.
  !> \param[in] use_stochastic Enable the stochastic rheology correction.
  !> \param[in] stochastic_value Local stochastic-process value.
  !> \param[out] rate x/y conservative momentum source rates [kg m^-1 s^-2].
  !> \note State arguments are complex-valued to preserve complex-step derivatives.
  !******************************************************************************
  SUBROUTINE eval_implicit_rheology_rate(model_id, flow_depth, temperature,   &
       mixture_density, inverse_density, velocity_x, velocity_y, speed,      &
       horizontal_speed, speed_squared, reduced_gravity, solid_fractions,    &
       moving, use_stochastic, stochastic_value, rate)

    INTEGER, INTENT(IN) :: model_id
    COMPLEX(wp), INTENT(IN) :: flow_depth, temperature, mixture_density
    COMPLEX(wp), INTENT(IN) :: inverse_density, velocity_x, velocity_y
    COMPLEX(wp), INTENT(IN) :: speed, horizontal_speed, speed_squared
    COMPLEX(wp), INTENT(IN) :: reduced_gravity, solid_fractions(:)
    LOGICAL, INTENT(IN) :: moving, use_stochastic
    REAL(wp), INTENT(IN) :: stochastic_value
    COMPLEX(wp), INTENT(OUT) :: rate(2)

    rate = CMPLX(0.0_wp, 0.0_wp, wp)
    SELECT CASE(model_id)
    CASE(rheology_id_voellmy)
       CALL eval_voellmy_implicit_momentum_rate(velocity_x, velocity_y,      &
            mixture_density, reduced_gravity, speed, speed_squared, moving,  &
            use_stochastic, stochastic_value, rate(1), rate(2))
    CASE(rheology_id_mui_voellmy)
       CALL eval_mui_voellmy_implicit_momentum_rate(velocity_x, velocity_y,  &
            mixture_density, reduced_gravity, speed, speed_squared, moving,  &
            use_stochastic, stochastic_value, rate(1), rate(2))
    CASE(rheology_id_plastic)
       CALL eval_plastic_implicit_momentum_rate(velocity_x, velocity_y,      &
            speed, mixture_density, moving, rate(1), rate(2))
    CASE(rheology_id_lava)
       CALL eval_lava_implicit_momentum_rate(flow_depth, temperature,        &
            mixture_density, velocity_x, velocity_y, speed, moving,          &
            rate(1), rate(2))
    CASE(rheology_id_lahar)
       CALL eval_lahar_implicit_momentum_rate(flow_depth, velocity_x,        &
            velocity_y, speed, horizontal_speed, temperature, mixture_density, &
            inverse_density, solid_fractions, moving, rate(1), rate(2))
    CASE(rheology_id_empirical)
       CALL eval_empirical_implicit_momentum_rate(flow_depth, velocity_x,    &
            velocity_y, horizontal_speed, speed, mixture_density, moving,   &
            rate(1), rate(2))
    CASE(rheology_id_quadratic)
       CALL eval_quadratic_implicit_momentum_rate(velocity_x, velocity_y,   &
            speed, horizontal_speed, speed_squared, mixture_density, moving, &
            rate(1), rate(2))
    CASE(rheology_id_none, rheology_id_froude, rheology_id_lucas, rheology_id_mui)
       CONTINUE
    CASE DEFAULT
       CALL fatal_error('Invalid rheology ID in implicit momentum-rate dispatch')
    END SELECT

  END SUBROUTINE eval_implicit_rheology_rate

  !******************************************************************************
  !> \brief Evaluate the selected semi-implicit rheology momentum rate.
  !> \details Curvature acceleration is used by the MUI-family laws only,
  !> matching their existing equation-term assembly semantics.
  !> \param[in] model_id Selected rheology identifier.
  !> \param[in] flow_depth Flow depth [m].
  !> \param[in] velocity_x x velocity [m s^-1].
  !> \param[in] velocity_y y velocity [m s^-1].
  !> \param[in] speed Total flow-speed magnitude [m s^-1].
  !> \param[in] horizontal_speed Horizontal flow-speed magnitude [m s^-1].
  !> \param[in] mixture_density Mixture density [kg m^-3].
  !> \param[in] reduced_gravity Reduced gravitational acceleration [m s^-2].
  !> \param[in] gravity_coefficient Effective gravity multiplier [dimensionless].
  !> \param[in] solid_fractions Solid volume fractions, one per class.
  !> \param[in] curvature_term Enable curvature acceleration for applicable laws.
  !> \param[in] curvature_xx Bed curvature in xx [m^-1].
  !> \param[in] curvature_xy Bed curvature in xy [m^-1].
  !> \param[in] curvature_yy Bed curvature in yy [m^-1].
  !> \param[in] use_pore_pressure Enable the pore-pressure correction.
  !> \param[in] excess_pore_pressure Excess pore pressure [Pa].
  !> \param[in] use_stochastic Enable the stochastic rheology correction.
  !> \param[in] stochastic_value Local stochastic-process value.
  !> \param[out] rate x/y conservative momentum source rates [kg m^-1 s^-2].
  !******************************************************************************
  SUBROUTINE eval_semi_implicit_rheology_rate(model_id, flow_depth,          &
       velocity_x, velocity_y, speed, horizontal_speed, mixture_density,    &
       reduced_gravity, gravity_coefficient, solid_fractions, curvature_term, &
       curvature_xx, curvature_xy, curvature_yy, use_pore_pressure,          &
       excess_pore_pressure, use_stochastic, stochastic_value, rate)

    INTEGER, INTENT(IN) :: model_id
    REAL(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y, speed
    REAL(wp), INTENT(IN) :: horizontal_speed, mixture_density, reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient, solid_fractions(:)
    LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure, use_stochastic
    REAL(wp), INTENT(IN) :: curvature_xx, curvature_xy, curvature_yy
    REAL(wp), INTENT(IN) :: excess_pore_pressure, stochastic_value
    REAL(wp), INTENT(OUT) :: rate(2)

    REAL(wp) :: curvature_term_value

    rate = 0.0_wp
    SELECT CASE(model_id)
    CASE(rheology_id_voellmy)
       CALL eval_voellmy_semi_implicit_momentum_rate(flow_depth, velocity_x, &
            velocity_y, speed, mixture_density, reduced_gravity,             &
            gravity_coefficient, curvature_term, curvature_xx, curvature_xy, &
            curvature_yy, use_pore_pressure, excess_pore_pressure,           &
            rate(1), rate(2))
    CASE(rheology_id_lava)
       CALL eval_lava_semi_implicit_momentum_rate(velocity_x, velocity_y,   &
            speed, rate(1), rate(2))
    CASE(rheology_id_lahar)
       CALL eval_lahar_semi_implicit_momentum_rate(flow_depth, velocity_x,  &
            velocity_y, horizontal_speed, mixture_density,                  &
            SUM(solid_fractions), speed .GT. 0.0_wp, rate(1), rate(2))
    CASE(rheology_id_froude)
       CALL eval_froude_momentum_rate(flow_depth, velocity_x, velocity_y,   &
            speed, horizontal_speed, mixture_density, reduced_gravity,      &
            gravity_coefficient, 0.0_wp, use_stochastic, stochastic_value,  &
            rate(1), rate(2))
    CASE(rheology_id_lucas)
       CALL eval_lucas_momentum_rate(flow_depth, velocity_x, velocity_y,    &
            speed, horizontal_speed, mixture_density, reduced_gravity,      &
            gravity_coefficient, 0.0_wp, use_stochastic, stochastic_value,  &
            rate(1), rate(2))
    CASE(rheology_id_mui)
       curvature_term_value = 0.0_wp
       IF (curvature_term) curvature_term_value = curvature_acceleration(   &
            curvature_xx, curvature_xy, curvature_yy, velocity_x, velocity_y)
       CALL eval_mui_momentum_rate(flow_depth, velocity_x, velocity_y,     &
            speed, mixture_density, reduced_gravity, gravity_coefficient,  &
            solid_fractions, curvature_term_value, use_pore_pressure,      &
            excess_pore_pressure, rate(1), rate(2))
    CASE(rheology_id_mui_voellmy)
       curvature_term_value = 0.0_wp
       IF (curvature_term) curvature_term_value = curvature_acceleration(   &
            curvature_xx, curvature_xy, curvature_yy, velocity_x, velocity_y)
       CALL eval_mui_voellmy_momentum_rate(flow_depth, velocity_x,          &
            velocity_y, speed, mixture_density, reduced_gravity,            &
            gravity_coefficient, solid_fractions, curvature_term_value,    &
            use_pore_pressure, excess_pore_pressure, rate(1), rate(2))
    CASE(rheology_id_none, rheology_id_plastic, rheology_id_empirical,       &
         rheology_id_quadratic)
       CONTINUE
    CASE DEFAULT
       CALL fatal_error('Invalid rheology ID in semi-implicit momentum-rate dispatch')
    END SELECT

  END SUBROUTINE eval_semi_implicit_rheology_rate

  !******************************************************************************
  !> \brief Integrate the selected rheology velocity law analytically.
  !> \details Dispatches to the model-specific exact-flow solver. All state
  !>          values are passed explicitly; unsupported model IDs leave the
  !>          output equal to the input and return model_supported as false.
  !> \param[in] model_id Rheology model ID.
  !> \param[in] integration_step Exact-flow interval [s].
  !> \param[in] flow_depth Flow depth [m].
  !> \param[in] temperature Flow temperature [K].
  !> \param[in] mixture_density Mixture density [kg m^-3].
  !> \param[in] speed_ratio Total-to-horizontal speed ratio.
  !> \param[in] velocity_before Horizontal velocity before integration [m s^-1].
  !> \param[in] reduced_gravity Reduced gravity [m s^-2].
  !> \param[in] gravity_coefficient Gravity coefficient used by rheology laws.
  !> \param[in] solid_fractions Per-class solid volume fractions.
  !> \param[in] curvature_acceleration Signed bed-curvature acceleration [m s^-2].
  !> \param[in] use_pore_pressure Whether excess pore pressure is active.
  !> \param[in] excess_pore_pressure Excess pore pressure.
  !> \param[in] use_stochastic Whether stochastic forcing is active.
  !> \param[in] stochastic_value Local stochastic rheology value.
  !> \param[out] velocity_after Horizontal velocity after integration [m s^-1].
  !> \param[out] model_supported True when model_id has an analytic velocity law.
  !******************************************************************************
  SUBROUTINE solve_analytic_rheology_velocity(model_id, integration_step,    &
       flow_depth, temperature, mixture_density, speed_ratio, velocity_before, &
       reduced_gravity, gravity_coefficient, solid_fractions,                &
       curvature_acceleration, use_pore_pressure, excess_pore_pressure,      &
       use_stochastic, stochastic_value, velocity_after, model_supported)

    INTEGER, INTENT(IN) :: model_id
    REAL(wp), INTENT(IN) :: integration_step, flow_depth, temperature
    REAL(wp), INTENT(IN) :: mixture_density, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2), reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient, solid_fractions(:)
    REAL(wp), INTENT(IN) :: curvature_acceleration, excess_pore_pressure
    REAL(wp), INTENT(IN) :: stochastic_value
    LOGICAL, INTENT(IN) :: use_pore_pressure, use_stochastic
    REAL(wp), INTENT(OUT) :: velocity_after(2)
    LOGICAL, INTENT(OUT) :: model_supported

    velocity_after = velocity_before
    model_supported = .TRUE.

    SELECT CASE (model_id)
    CASE (rheology_id_none)
       CONTINUE
    CASE (rheology_id_voellmy)
       CALL solve_voellmy_analytic_velocity(integration_step, flow_depth,    &
            speed_ratio, velocity_before, mixture_density, reduced_gravity,  &
            gravity_coefficient, curvature_acceleration, use_pore_pressure, &
            excess_pore_pressure, use_stochastic, stochastic_value,         &
            velocity_after)
    CASE (rheology_id_plastic)
       CALL solve_plastic_analytic_velocity(integration_step, flow_depth,    &
            speed_ratio, velocity_before, velocity_after)
    CASE (rheology_id_lava)
       CALL solve_lava_analytic_velocity(integration_step, flow_depth,      &
            temperature, mixture_density, speed_ratio, velocity_before,     &
            velocity_after)
    CASE (rheology_id_lahar)
       CALL solve_lahar_analytic_velocity(integration_step, flow_depth,     &
            temperature, SUM(solid_fractions), mixture_density, speed_ratio,&
            velocity_before, velocity_after)
    CASE (rheology_id_empirical)
       CALL solve_empirical_analytic_velocity(integration_step, flow_depth, &
            speed_ratio, velocity_before, velocity_after)
    CASE (rheology_id_quadratic)
       CALL solve_quadratic_analytic_velocity(integration_step, flow_depth, &
            speed_ratio, velocity_before, velocity_after)
    CASE (rheology_id_froude)
       CALL solve_froude_analytic_velocity(integration_step, flow_depth,    &
            speed_ratio, velocity_before, reduced_gravity,                  &
            gravity_coefficient, curvature_acceleration, use_stochastic,   &
            stochastic_value, velocity_after)
    CASE (rheology_id_lucas)
       CALL solve_lucas_analytic_velocity(integration_step, flow_depth,     &
            speed_ratio, velocity_before, reduced_gravity,                  &
            gravity_coefficient, curvature_acceleration, use_stochastic,   &
            stochastic_value, velocity_after)
    CASE (rheology_id_mui)
       CALL solve_mui_analytic_velocity(integration_step, flow_depth,       &
            speed_ratio, velocity_before, mixture_density, reduced_gravity, &
            gravity_coefficient, solid_fractions, curvature_acceleration,   &
            use_pore_pressure, excess_pore_pressure, velocity_after)
    CASE (rheology_id_mui_voellmy)
       CALL solve_mui_voellmy_analytic_velocity(integration_step,           &
            flow_depth, speed_ratio, velocity_before, mixture_density,      &
            reduced_gravity, gravity_coefficient, solid_fractions,           &
            curvature_acceleration, use_pore_pressure, excess_pore_pressure, &
            use_stochastic, stochastic_value, velocity_after)
    CASE DEFAULT
       model_supported = .FALSE.
    END SELECT

  END SUBROUTINE solve_analytic_rheology_velocity

  !******************************************************************************
  !> \brief Constrain a momentum update from an isolated friction substep.
  !> \details Projects the candidate post-friction momentum onto the ray
  !> defined by the pre-friction momentum and clamps its scale factor to
  !> \f$[0,1]\f$. This preserves direction, prevents reversal, and ensures
  !> friction cannot increase momentum magnitude. If the pre-friction
  !> momentum is zero, the post-friction momentum is set to zero.
  !> \param[in] momentum_before_friction Momentum immediately before friction.
  !> \param[in,out] momentum_after_friction Candidate momentum immediately after friction.
  !******************************************************************************
  PURE SUBROUTINE correct_friction_momentum(momentum_before_friction,         &
       momentum_after_friction)

    USE parameters_2d, ONLY : wp
    REAL(wp), INTENT(IN) :: momentum_before_friction(2)
    REAL(wp), INTENT(INOUT) :: momentum_after_friction(2)

    REAL(wp) :: momentum_before_squared
    REAL(wp) :: friction_scale

    momentum_before_squared = DOT_PRODUCT(momentum_before_friction,           &
         momentum_before_friction)
    IF (momentum_before_squared .EQ. 0.0_wp) THEN
       momentum_after_friction = 0.0_wp
       RETURN
    END IF

    friction_scale = DOT_PRODUCT(momentum_before_friction,                   &
         momentum_after_friction) / momentum_before_squared
    friction_scale = MAX(0.0_wp, MIN(1.0_wp, friction_scale))
    momentum_after_friction = friction_scale * momentum_before_friction

  END SUBROUTINE correct_friction_momentum

  !******************************************************************************
  !> \brief Resolve a user-supplied rheology name to its internal model ID.
  !>
  !> \param[in,out] name Input model name; trimmed and canonicalized to lowercase.
  !> \param[out] model_id Internal model ID, or -1 for an unsupported name.
  !******************************************************************************

  SUBROUTINE resolve_rheology_model( name, model_id )

    CHARACTER(LEN=*), INTENT(INOUT) :: name
    INTEGER, INTENT(OUT) :: model_id

    CHARACTER(LEN=LEN(name)) :: canonical_name
    INTEGER :: character_index, character_code

    canonical_name = ADJUSTL(name)
    DO character_index = 1, LEN_TRIM(canonical_name)
       character_code = IACHAR(canonical_name(character_index:character_index))
       IF ( character_code .GE. IACHAR('A') .AND. &
            character_code .LE. IACHAR('Z') ) THEN
          canonical_name(character_index:character_index) =                  &
               ACHAR(character_code + IACHAR('a') - IACHAR('A'))
       END IF
    END DO
    name = TRIM(canonical_name)

    SELECT CASE (TRIM(name))
    CASE ('none')
       model_id = rheology_id_none
    CASE ('voellmy')
       model_id = rheology_id_voellmy
    CASE ('plastic')
       model_id = rheology_id_plastic
    CASE ('lava')
       model_id = rheology_id_lava
    CASE ('lahar')
       model_id = rheology_id_lahar
    CASE ('empirical')
       model_id = rheology_id_empirical
    CASE ('quadratic')
       model_id = rheology_id_quadratic
    CASE ('froude')
       model_id = rheology_id_froude
    CASE ('lucas')
       model_id = rheology_id_lucas
    CASE ('mui')
       model_id = rheology_id_mui
    CASE ('mui_voellmy')
       model_id = rheology_id_mui_voellmy
    CASE DEFAULT
       model_id = -1
    END SELECT

  END SUBROUTINE resolve_rheology_model


  !******************************************************************************
  !> \brief Return the canonical input name corresponding to an internal ID.
  !>
  !> \param[in] model_id Internal rheology model ID.
  !> \return Canonical name, or 'unknown' for an unsupported ID.
  !******************************************************************************

  PURE FUNCTION rheology_name_from_id( model_id ) RESULT(name)

    INTEGER, INTENT(IN) :: model_id
    CHARACTER(LEN=32) :: name

    SELECT CASE (model_id)
    CASE (rheology_id_none)
       name = 'none'
    CASE (rheology_id_voellmy)
       name = 'voellmy'
    CASE (rheology_id_plastic)
       name = 'plastic'
    CASE (rheology_id_lava)
       name = 'lava'
    CASE (rheology_id_lahar)
       name = 'lahar'
    CASE (rheology_id_empirical)
       name = 'empirical'
    CASE (rheology_id_quadratic)
       name = 'quadratic'
    CASE (rheology_id_froude)
       name = 'froude'
    CASE (rheology_id_lucas)
       name = 'lucas'
    CASE (rheology_id_mui)
       name = 'mui'
    CASE (rheology_id_mui_voellmy)
       name = 'mui_voellmy'
    CASE DEFAULT
       name = 'unknown'
    END SELECT

  END FUNCTION rheology_name_from_id

END MODULE rheology_2d
