!********************************************************************************
!> \brief Rheology dispatch, model naming and friction-only momentum correction.
!> \details Dispatches common momentum-rate and analytic-velocity interfaces
!> to the selected constitutive law using canonical source names,
!> and provides a projection for momentum updates caused exclusively by
!> friction. The projection constrains the candidate update to
!> \f$\boldsymbol{m}_{after}=\alpha\boldsymbol{m}_{before}\f$,
!> \f$0\le\alpha\le1\f$, so it cannot rotate, reverse or increase the
!> momentum. Callers must pass the momentum immediately before and after the
!> isolated friction substep, not the full stage state.
!********************************************************************************

MODULE rheology_2d

  USE parameters_2d, ONLY : wp
  USE diagnostics_2d, ONLY : fatal_error
  USE rheology_voellmy_2d, ONLY : eval_rheology_voellmy_implicit,       &
       eval_rheology_voellmy_semi_implicit,                            &
       solve_rheology_voellmy_analytic
  USE rheology_plastic_2d, ONLY : eval_rheology_plastic_implicit,       &
       solve_rheology_plastic_analytic
  USE rheology_lava_2d, ONLY : eval_rheology_lava_implicit,             &
       eval_rheology_lava_semi_implicit,                               &
       solve_rheology_lava_analytic
  USE rheology_lahar_2d, ONLY : eval_rheology_lahar_implicit,            &
       eval_rheology_lahar_semi_implicit,                              &
       solve_rheology_lahar_analytic
  USE rheology_empirical_2d, ONLY : eval_rheology_empirical_implicit,   &
       solve_rheology_empirical_analytic
  USE rheology_quadratic_2d, ONLY : eval_rheology_quadratic_implicit,   &
       solve_rheology_quadratic_analytic
  USE rheology_newmodel_2d, ONLY : eval_rheology_newmodel_implicit
  USE rheology_froude_2d, ONLY : eval_rheology_froude_semi_implicit,                  &
       solve_rheology_froude_analytic
  USE rheology_lucas_2d, ONLY : eval_rheology_lucas_semi_implicit,                    &
       solve_rheology_lucas_analytic
  USE rheology_mui_2d, ONLY : eval_rheology_mui_semi_implicit,                        &
       solve_rheology_mui_analytic
  USE rheology_mui_voellmy_2d, ONLY : eval_rheology_mui_voellmy_semi_implicit,        &
       eval_rheology_mui_voellmy_implicit,                             &
       solve_rheology_mui_voellmy_analytic

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_rheology_explicit
  PUBLIC :: eval_rheology_semi_implicit
  PUBLIC :: eval_rheology_implicit
  PUBLIC :: solve_rheology_analytic
  PUBLIC :: curvature_acceleration

  ! Uniform treatment surface for this source family.

CONTAINS

  SUBROUTINE eval_rheology_explicit()
    CALL fatal_error('rheology: explicit treatment is unsupported')
  END SUBROUTINE eval_rheology_explicit

  SUBROUTINE eval_rheology_semi_implicit(model_name, flow_depth,          &
       velocity_x, velocity_y, speed, horizontal_speed, mixture_density,    &
       reduced_gravity, gravity_coefficient, solid_fractions, curvature_term, &
       curvature_xx, curvature_xy, curvature_yy, use_pore_pressure,          &
       excess_pore_pressure, use_stochastic, stochastic_value, rate)

    CHARACTER(*), INTENT(IN) :: model_name
    REAL(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y, speed
    REAL(wp), INTENT(IN) :: horizontal_speed, mixture_density, reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient, solid_fractions(:)
    LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure, use_stochastic
    REAL(wp), INTENT(IN) :: curvature_xx, curvature_xy, curvature_yy
    REAL(wp), INTENT(IN) :: excess_pore_pressure, stochastic_value
    REAL(wp), INTENT(OUT) :: rate(2)

    REAL(wp) :: curvature_term_value

    rate = 0.0_wp
    SELECT CASE(TRIM(model_name))
    CASE('rheology_voellmy')
       CALL eval_rheology_voellmy_semi_implicit(flow_depth, velocity_x, &
            velocity_y, speed, mixture_density, reduced_gravity,             &
            gravity_coefficient, curvature_term, curvature_xx, curvature_xy, &
            curvature_yy, use_pore_pressure, excess_pore_pressure,           &
            rate(1), rate(2))
    CASE('rheology_lava')
       CALL eval_rheology_lava_semi_implicit(velocity_x, velocity_y,   &
            speed, rate(1), rate(2))
    CASE('rheology_lahar')
       CALL eval_rheology_lahar_semi_implicit(flow_depth, velocity_x,  &
            velocity_y, horizontal_speed, mixture_density,                  &
            SUM(solid_fractions), speed .GT. 0.0_wp, rate(1), rate(2))
    CASE('rheology_froude')
       CALL eval_rheology_froude_semi_implicit(flow_depth, velocity_x, velocity_y,   &
            speed, horizontal_speed, mixture_density, reduced_gravity,      &
            gravity_coefficient, 0.0_wp, use_stochastic, stochastic_value,  &
            rate(1), rate(2))
    CASE('rheology_lucas')
       CALL eval_rheology_lucas_semi_implicit(flow_depth, velocity_x, velocity_y,    &
            speed, horizontal_speed, mixture_density, reduced_gravity,      &
            gravity_coefficient, 0.0_wp, use_stochastic, stochastic_value,  &
            rate(1), rate(2))
    CASE('rheology_mui')
       curvature_term_value = 0.0_wp
       IF (curvature_term) curvature_term_value = curvature_acceleration(   &
            curvature_xx, curvature_xy, curvature_yy, velocity_x, velocity_y)
       CALL eval_rheology_mui_semi_implicit(flow_depth, velocity_x, velocity_y,     &
            speed, mixture_density, reduced_gravity, gravity_coefficient,  &
            solid_fractions, curvature_term_value, use_pore_pressure,      &
            excess_pore_pressure, rate(1), rate(2))
    CASE('rheology_mui_voellmy')
       curvature_term_value = 0.0_wp
       IF (curvature_term) curvature_term_value = curvature_acceleration(   &
            curvature_xx, curvature_xy, curvature_yy, velocity_x, velocity_y)
       CALL eval_rheology_mui_voellmy_semi_implicit(flow_depth, velocity_x,          &
            velocity_y, speed, mixture_density, reduced_gravity,            &
            gravity_coefficient, solid_fractions, curvature_term_value,    &
            use_pore_pressure, excess_pore_pressure, rate(1), rate(2))
    CASE('none', 'rheology_plastic', 'rheology_empirical',       &
         'rheology_quadratic', 'rheology_newmodel')
       CONTINUE
    CASE DEFAULT
       CALL fatal_error('Invalid rheology name in semi-implicit momentum-rate dispatch')
    END SELECT

  END SUBROUTINE eval_rheology_semi_implicit

  SUBROUTINE eval_rheology_implicit(model_name, flow_depth, temperature,   &
       mixture_density, inverse_density, velocity_x, velocity_y, speed,      &
       horizontal_speed, speed_squared, reduced_gravity, solid_fractions,    &
       moving, use_stochastic, stochastic_value, rate)

    CHARACTER(*), INTENT(IN) :: model_name
    COMPLEX(wp), INTENT(IN) :: flow_depth, temperature, mixture_density
    COMPLEX(wp), INTENT(IN) :: inverse_density, velocity_x, velocity_y
    COMPLEX(wp), INTENT(IN) :: speed, horizontal_speed, speed_squared
    COMPLEX(wp), INTENT(IN) :: reduced_gravity, solid_fractions(:)
    LOGICAL, INTENT(IN) :: moving, use_stochastic
    REAL(wp), INTENT(IN) :: stochastic_value
    COMPLEX(wp), INTENT(OUT) :: rate(2)

    rate = CMPLX(0.0_wp, 0.0_wp, wp)
    SELECT CASE(TRIM(model_name))
    CASE('rheology_voellmy')
       CALL eval_rheology_voellmy_implicit(velocity_x, velocity_y,      &
            mixture_density, reduced_gravity, speed, speed_squared, moving,  &
            use_stochastic, stochastic_value, rate(1), rate(2))
    CASE('rheology_mui_voellmy')
       CALL eval_rheology_mui_voellmy_implicit(velocity_x, velocity_y,  &
            mixture_density, reduced_gravity, speed, speed_squared, moving,  &
            use_stochastic, stochastic_value, rate(1), rate(2))
    CASE('rheology_plastic')
       CALL eval_rheology_plastic_implicit(velocity_x, velocity_y,      &
            speed, mixture_density, moving, rate(1), rate(2))
    CASE('rheology_lava')
       CALL eval_rheology_lava_implicit(flow_depth, temperature,        &
            mixture_density, velocity_x, velocity_y, moving,                &
            rate(1), rate(2))
    CASE('rheology_lahar')
       CALL eval_rheology_lahar_implicit(flow_depth, velocity_x,        &
            velocity_y, speed, horizontal_speed, temperature, mixture_density, &
            inverse_density, solid_fractions, moving, rate(1), rate(2))
    CASE('rheology_empirical')
       CALL eval_rheology_empirical_implicit(flow_depth, velocity_x,    &
            velocity_y, horizontal_speed, speed, mixture_density, moving,   &
            rate(1), rate(2))
    CASE('rheology_newmodel')
       CALL eval_rheology_newmodel_implicit(velocity_x,velocity_y, &
            mixture_density,rate(1),rate(2))
    CASE('rheology_quadratic')
       CALL eval_rheology_quadratic_implicit(velocity_x, velocity_y,   &
            horizontal_speed, speed_squared, mixture_density, moving,       &
            rate(1), rate(2))
    CASE('none', 'rheology_froude', 'rheology_lucas', 'rheology_mui')
       CONTINUE
    CASE DEFAULT
       CALL fatal_error('Invalid rheology name in implicit momentum-rate dispatch')
    END SELECT

  END SUBROUTINE eval_rheology_implicit

  SUBROUTINE solve_rheology_analytic(model_name, integration_step,    &
       flow_depth, temperature, mixture_density, speed_ratio, velocity_before, &
       reduced_gravity, gravity_coefficient, solid_fractions,                &
       curvature_acceleration, use_pore_pressure, excess_pore_pressure,      &
       use_stochastic, stochastic_value, velocity_after, model_supported)

    CHARACTER(*), INTENT(IN) :: model_name
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

    SELECT CASE(TRIM(model_name))
    CASE ('none')
       CONTINUE
    CASE ('rheology_voellmy')
       CALL solve_rheology_voellmy_analytic(integration_step, flow_depth,    &
            speed_ratio, velocity_before, mixture_density, reduced_gravity,  &
            gravity_coefficient, curvature_acceleration, use_pore_pressure, &
            excess_pore_pressure, use_stochastic, stochastic_value,         &
            velocity_after)
    CASE ('rheology_plastic')
       CALL solve_rheology_plastic_analytic(integration_step, flow_depth,    &
            speed_ratio, velocity_before, velocity_after)
    CASE ('rheology_lava')
       CALL solve_rheology_lava_analytic(integration_step, flow_depth,      &
            temperature, mixture_density, speed_ratio, velocity_before,     &
            velocity_after)
    CASE ('rheology_lahar')
       CALL solve_rheology_lahar_analytic(integration_step, flow_depth,     &
            temperature, SUM(solid_fractions), mixture_density, speed_ratio,&
            velocity_before, velocity_after)
    CASE ('rheology_empirical')
       CALL solve_rheology_empirical_analytic(integration_step, flow_depth, &
            speed_ratio, velocity_before, velocity_after)
    CASE ('rheology_quadratic')
       CALL solve_rheology_quadratic_analytic(integration_step, flow_depth, &
            speed_ratio, velocity_before, velocity_after)
    CASE ('rheology_froude')
       CALL solve_rheology_froude_analytic(integration_step, flow_depth,    &
            speed_ratio, velocity_before, reduced_gravity,                  &
            gravity_coefficient, curvature_acceleration, use_stochastic,   &
            stochastic_value, velocity_after)
    CASE ('rheology_lucas')
       CALL solve_rheology_lucas_analytic(integration_step, flow_depth,     &
            speed_ratio, velocity_before, reduced_gravity,                  &
            gravity_coefficient, curvature_acceleration, use_stochastic,   &
            stochastic_value, velocity_after)
    CASE ('rheology_mui')
       CALL solve_rheology_mui_analytic(integration_step, flow_depth,       &
            speed_ratio, velocity_before, mixture_density, reduced_gravity, &
            gravity_coefficient, solid_fractions, curvature_acceleration,   &
            use_pore_pressure, excess_pore_pressure, velocity_after)
    CASE ('rheology_mui_voellmy')
       CALL solve_rheology_mui_voellmy_analytic(integration_step,           &
            flow_depth, speed_ratio, velocity_before, mixture_density,      &
            reduced_gravity, gravity_coefficient, solid_fractions,           &
            curvature_acceleration, use_pore_pressure, excess_pore_pressure, &
            use_stochastic, stochastic_value, velocity_after)
    CASE DEFAULT
       model_supported = .FALSE.
    END SELECT

  END SUBROUTINE solve_rheology_analytic


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
  !> \param[in] model_name Selected rheology identifier.
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

  !******************************************************************************
  !> \brief Evaluate the selected semi-implicit rheology momentum rate.
  !> \details Curvature acceleration is used by the MUI-family laws only,
  !> matching their existing equation-term assembly semantics.
  !> \param[in] model_name Selected rheology identifier.
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

  !******************************************************************************
  !> \brief Integrate the selected rheology velocity law analytically.
  !> \details Dispatches to the model-specific exact-flow solver. All state
  !>          values are passed explicitly; unsupported model IDs leave the
  !>          output equal to the input and return model_supported as false.
  !> \param[in] model_name Rheology model ID.
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
  !> \param[out] model_supported True when model_name has an analytic velocity law.
  !******************************************************************************




END MODULE rheology_2d
