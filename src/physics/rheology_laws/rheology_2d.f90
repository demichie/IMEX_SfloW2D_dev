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
  USE rheology_voellmy_2d, ONLY : eval_rheology_voellmy_explicit, &
        eval_rheology_voellmy_implicit,       &
       eval_rheology_voellmy_semi_implicit,                            &
       solve_rheology_voellmy_analytic
  USE rheology_plastic_2d, ONLY : eval_rheology_plastic_explicit, &
        eval_rheology_plastic_implicit,       &
       solve_rheology_plastic_analytic
  USE rheology_lava_2d, ONLY : eval_rheology_lava_explicit, &
        eval_rheology_lava_implicit,             &
       eval_rheology_lava_semi_implicit,                               &
       solve_rheology_lava_analytic
  USE rheology_lahar_2d, ONLY : eval_rheology_lahar_explicit, &
        eval_rheology_lahar_implicit,            &
       eval_rheology_lahar_semi_implicit,                              &
       solve_rheology_lahar_analytic
  USE rheology_empirical_2d, ONLY : eval_rheology_empirical_explicit, &
        eval_rheology_empirical_implicit,   &
       solve_rheology_empirical_analytic
  USE rheology_quadratic_2d, ONLY : eval_rheology_quadratic_explicit, &
        eval_rheology_quadratic_implicit,   &
       solve_rheology_quadratic_analytic
  USE rheology_newmodel_2d, ONLY : eval_rheology_newmodel_explicit, &
        eval_rheology_newmodel_implicit
  USE rheology_froude_2d, ONLY : eval_rheology_froude_explicit, &
        eval_rheology_froude_semi_implicit,                  &
       solve_rheology_froude_analytic
  USE rheology_lucas_2d, ONLY : eval_rheology_lucas_explicit, &
        eval_rheology_lucas_semi_implicit,                    &
       solve_rheology_lucas_analytic
  USE rheology_mui_2d, ONLY : eval_rheology_mui_explicit, &
        eval_rheology_mui_semi_implicit,                        &
       solve_rheology_mui_analytic
  USE rheology_mui_voellmy_2d, ONLY : eval_rheology_mui_voellmy_explicit, &
        eval_rheology_mui_voellmy_semi_implicit,        &
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

  ! Both Newton components are evaluated on the same, unchanged physical state.
  SUBROUTINE eval_rheology_explicit(model_name, h, T, &
       u, v, mod_vel, mod_hor_vel, mod_vel2, &
       rhom, inv_rhom, red_grav, grav_coeff, &
       alphas, curvature_term, B_second_xx, B_second_xy, B_second_yy, &
       use_pore_pressure, exc_pore_pres, moving, use_stochastic, &
       stochastic_value, rate)
    CHARACTER(*), INTENT(IN) :: model_name
    REAL(wp), INTENT(IN) :: h, T, u, v
    REAL(wp), INTENT(IN) :: mod_vel, mod_hor_vel, mod_vel2
    REAL(wp), INTENT(IN) :: rhom, inv_rhom
    REAL(wp), INTENT(IN) :: red_grav, grav_coeff
    REAL(wp), INTENT(IN) :: alphas(:)
    REAL(wp), INTENT(IN) :: B_second_xx, B_second_xy, B_second_yy
    REAL(wp), INTENT(IN) :: exc_pore_pres, stochastic_value
    LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure, moving, use_stochastic
    REAL(wp), INTENT(OUT) :: rate(2)
    SELECT CASE(TRIM(model_name))
    CASE('rheology_voellmy')
       CALL eval_rheology_voellmy_explicit( &
            h, T, u, v, mod_vel, mod_hor_vel, &
            mod_vel2, rhom, inv_rhom, red_grav, &
            grav_coeff, alphas, curvature_term, B_second_xx, &
            B_second_xy, B_second_yy, use_pore_pressure, exc_pore_pres, moving, &
            use_stochastic, stochastic_value, rate )
    CASE('rheology_plastic')
       CALL eval_rheology_plastic_explicit( &
            h, T, u, v, mod_vel, mod_hor_vel, &
            mod_vel2, rhom, inv_rhom, red_grav, &
            grav_coeff, alphas, curvature_term, B_second_xx, &
            B_second_xy, B_second_yy, use_pore_pressure, exc_pore_pres, moving, &
            use_stochastic, stochastic_value, rate )
    CASE('rheology_lava')
       CALL eval_rheology_lava_explicit( &
            h, T, u, v, mod_vel, mod_hor_vel, &
            mod_vel2, rhom, inv_rhom, red_grav, &
            grav_coeff, alphas, curvature_term, B_second_xx, &
            B_second_xy, B_second_yy, use_pore_pressure, exc_pore_pres, moving, &
            use_stochastic, stochastic_value, rate )
    CASE('rheology_lahar')
       CALL eval_rheology_lahar_explicit( &
            h, T, u, v, mod_vel, mod_hor_vel, &
            mod_vel2, rhom, inv_rhom, red_grav, &
            grav_coeff, alphas, curvature_term, B_second_xx, &
            B_second_xy, B_second_yy, use_pore_pressure, exc_pore_pres, moving, &
            use_stochastic, stochastic_value, rate )
    CASE('rheology_empirical')
       CALL eval_rheology_empirical_explicit( &
            h, T, u, v, mod_vel, mod_hor_vel, &
            mod_vel2, rhom, inv_rhom, red_grav, &
            grav_coeff, alphas, curvature_term, B_second_xx, &
            B_second_xy, B_second_yy, use_pore_pressure, exc_pore_pres, moving, &
            use_stochastic, stochastic_value, rate )
    CASE('rheology_quadratic')
       CALL eval_rheology_quadratic_explicit( &
            h, T, u, v, mod_vel, mod_hor_vel, &
            mod_vel2, rhom, inv_rhom, red_grav, &
            grav_coeff, alphas, curvature_term, B_second_xx, &
            B_second_xy, B_second_yy, use_pore_pressure, exc_pore_pres, moving, &
            use_stochastic, stochastic_value, rate )
    CASE('rheology_newmodel')
       CALL eval_rheology_newmodel_explicit( &
            h, T, u, v, mod_vel, mod_hor_vel, &
            mod_vel2, rhom, inv_rhom, red_grav, &
            grav_coeff, alphas, curvature_term, B_second_xx, &
            B_second_xy, B_second_yy, use_pore_pressure, exc_pore_pres, moving, &
            use_stochastic, stochastic_value, rate )
    CASE('rheology_froude')
       CALL eval_rheology_froude_explicit( &
            h, T, u, v, mod_vel, mod_hor_vel, &
            mod_vel2, rhom, inv_rhom, red_grav, &
            grav_coeff, alphas, curvature_term, B_second_xx, &
            B_second_xy, B_second_yy, use_pore_pressure, exc_pore_pres, moving, &
            use_stochastic, stochastic_value, rate )
    CASE('rheology_lucas')
       CALL eval_rheology_lucas_explicit( &
            h, T, u, v, mod_vel, mod_hor_vel, &
            mod_vel2, rhom, inv_rhom, red_grav, &
            grav_coeff, alphas, curvature_term, B_second_xx, &
            B_second_xy, B_second_yy, use_pore_pressure, exc_pore_pres, moving, &
            use_stochastic, stochastic_value, rate )
    CASE('rheology_mui')
       CALL eval_rheology_mui_explicit( &
            h, T, u, v, mod_vel, mod_hor_vel, &
            mod_vel2, rhom, inv_rhom, red_grav, &
            grav_coeff, alphas, curvature_term, B_second_xx, &
            B_second_xy, B_second_yy, use_pore_pressure, exc_pore_pres, moving, &
            use_stochastic, stochastic_value, rate )
    CASE('rheology_mui_voellmy')
       CALL eval_rheology_mui_voellmy_explicit( &
            h, T, u, v, mod_vel, mod_hor_vel, &
            mod_vel2, rhom, inv_rhom, red_grav, &
            grav_coeff, alphas, curvature_term, B_second_xx, &
            B_second_xy, B_second_yy, use_pore_pressure, exc_pore_pres, moving, &
            use_stochastic, stochastic_value, rate )
    CASE DEFAULT
       CALL fatal_error('Invalid rheology name in explicit dispatch')
    END SELECT
  END SUBROUTINE eval_rheology_explicit

  SUBROUTINE eval_rheology_semi_implicit(model_name, h,          &
       u, v, mod_vel, mod_hor_vel, rhom,    &
       red_grav, grav_coeff, alphas, curvature_term, &
       B_second_xx, B_second_xy, B_second_yy, use_pore_pressure,          &
       exc_pore_pres, use_stochastic, stochastic_value, rate)

    CHARACTER(*), INTENT(IN) :: model_name
    REAL(wp), INTENT(IN) :: h, u, v, mod_vel
    REAL(wp), INTENT(IN) :: mod_hor_vel, rhom, red_grav
    REAL(wp), INTENT(IN) :: grav_coeff, alphas(:)
    LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure, use_stochastic
    REAL(wp), INTENT(IN) :: B_second_xx, B_second_xy, B_second_yy
    REAL(wp), INTENT(IN) :: exc_pore_pres, stochastic_value
    REAL(wp), INTENT(OUT) :: rate(2)

    REAL(wp) :: curvature_term_value

    rate = 0.0_wp
    SELECT CASE(TRIM(model_name))
    CASE('rheology_voellmy')
       CALL eval_rheology_voellmy_semi_implicit(h, u, &
            v, mod_vel, rhom, red_grav,             &
            grav_coeff, curvature_term, B_second_xx, B_second_xy, &
            B_second_yy, use_pore_pressure, exc_pore_pres,           &
            rate(1), rate(2))
    CASE('rheology_lava')
       CALL eval_rheology_lava_semi_implicit(u, v,   &
            mod_vel, rate(1), rate(2))
    CASE('rheology_lahar')
       CALL eval_rheology_lahar_semi_implicit(h, u,  &
            v, mod_hor_vel, rhom,                  &
            SUM(alphas), mod_vel .GT. 0.0_wp, rate(1), rate(2))
    CASE('rheology_froude')
       CALL eval_rheology_froude_semi_implicit(h, u, v,   &
            mod_vel, mod_hor_vel, rhom, red_grav,      &
            grav_coeff, 0.0_wp, use_stochastic, stochastic_value,  &
            rate(1), rate(2))
    CASE('rheology_lucas')
       CALL eval_rheology_lucas_semi_implicit(h, u, v,    &
            mod_vel, mod_hor_vel, rhom, red_grav,      &
            grav_coeff, 0.0_wp, use_stochastic, stochastic_value,  &
            rate(1), rate(2))
    CASE('rheology_mui')
       curvature_term_value = 0.0_wp
       IF (curvature_term) curvature_term_value = curvature_acceleration(   &
            B_second_xx, B_second_xy, B_second_yy, u, v)
       CALL eval_rheology_mui_semi_implicit(h, u, v,     &
            mod_vel, rhom, red_grav, grav_coeff,  &
            alphas, curvature_term_value, use_pore_pressure,      &
            exc_pore_pres, rate(1), rate(2))
    CASE('rheology_mui_voellmy')
       curvature_term_value = 0.0_wp
       IF (curvature_term) curvature_term_value = curvature_acceleration(   &
            B_second_xx, B_second_xy, B_second_yy, u, v)
       CALL eval_rheology_mui_voellmy_semi_implicit(h, u,          &
            v, mod_vel, rhom, red_grav,            &
            grav_coeff, alphas, curvature_term_value,    &
            use_pore_pressure, exc_pore_pres, rate(1), rate(2))
    CASE('none', 'rheology_plastic', 'rheology_empirical',       &
         'rheology_quadratic', 'rheology_newmodel')
       CONTINUE
    CASE DEFAULT
       CALL fatal_error('Invalid rheology name in semi-implicit momentum-rate dispatch')
    END SELECT

  END SUBROUTINE eval_rheology_semi_implicit

  SUBROUTINE eval_rheology_implicit(model_name, h, T,   &
       rhom, inv_rhom, u, v, mod_vel,      &
       mod_hor_vel, mod_vel2, red_grav, alphas,    &
       moving, use_stochastic, stochastic_value, rate)

    CHARACTER(*), INTENT(IN) :: model_name
    COMPLEX(wp), INTENT(IN) :: h, T, rhom
    COMPLEX(wp), INTENT(IN) :: inv_rhom, u, v
    COMPLEX(wp), INTENT(IN) :: mod_vel, mod_hor_vel, mod_vel2
    COMPLEX(wp), INTENT(IN) :: red_grav, alphas(:)
    LOGICAL, INTENT(IN) :: moving, use_stochastic
    REAL(wp), INTENT(IN) :: stochastic_value
    COMPLEX(wp), INTENT(OUT) :: rate(2)

    rate = CMPLX(0.0_wp, 0.0_wp, wp)
    SELECT CASE(TRIM(model_name))
    CASE('rheology_voellmy')
       CALL eval_rheology_voellmy_implicit(u, v,      &
            rhom, red_grav, mod_vel, mod_vel2, moving,  &
            use_stochastic, stochastic_value, rate(1), rate(2))
    CASE('rheology_mui_voellmy')
       CALL eval_rheology_mui_voellmy_implicit(u, v,  &
            rhom, red_grav, mod_vel, mod_vel2, moving,  &
            use_stochastic, stochastic_value, rate(1), rate(2))
    CASE('rheology_plastic')
       CALL eval_rheology_plastic_implicit(u, v,      &
            mod_vel, rhom, moving, rate(1), rate(2))
    CASE('rheology_lava')
       CALL eval_rheology_lava_implicit(h, T,        &
            rhom, u, v, moving,                &
            rate(1), rate(2))
    CASE('rheology_lahar')
       CALL eval_rheology_lahar_implicit(h, u,        &
            v, mod_vel, mod_hor_vel, T, rhom, &
            inv_rhom, alphas, moving, rate(1), rate(2))
    CASE('rheology_empirical')
       CALL eval_rheology_empirical_implicit(h, u,    &
            v, mod_hor_vel, mod_vel, rhom, moving,   &
            rate(1), rate(2))
    CASE('rheology_newmodel')
       CALL eval_rheology_newmodel_implicit(u,v, &
            rhom,rate(1),rate(2))
    CASE('rheology_quadratic')
       CALL eval_rheology_quadratic_implicit(u, v,   &
            mod_hor_vel, mod_vel2, rhom, moving,       &
            rate(1), rate(2))
    CASE('none', 'rheology_froude', 'rheology_lucas', 'rheology_mui')
       CONTINUE
    CASE DEFAULT
       CALL fatal_error('Invalid rheology name in implicit momentum-rate dispatch')
    END SELECT

  END SUBROUTINE eval_rheology_implicit

  SUBROUTINE solve_rheology_analytic(model_name, dt,    &
       h, T, rhom, speed_ratio, velocity_before, &
       red_grav, grav_coeff, alphas,                &
       curvature_acceleration, use_pore_pressure, exc_pore_pres,      &
       use_stochastic, stochastic_value, velocity_after, model_supported)

    CHARACTER(*), INTENT(IN) :: model_name
    REAL(wp), INTENT(IN) :: dt, h, T
    REAL(wp), INTENT(IN) :: rhom, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2), red_grav
    REAL(wp), INTENT(IN) :: grav_coeff, alphas(:)
    REAL(wp), INTENT(IN) :: curvature_acceleration, exc_pore_pres
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
       CALL solve_rheology_voellmy_analytic(dt, h,    &
            speed_ratio, velocity_before, rhom, red_grav,  &
            grav_coeff, curvature_acceleration, use_pore_pressure, &
            exc_pore_pres, use_stochastic, stochastic_value,         &
            velocity_after)
    CASE ('rheology_plastic')
       CALL solve_rheology_plastic_analytic(dt, h,    &
            speed_ratio, velocity_before, velocity_after)
    CASE ('rheology_lava')
       CALL solve_rheology_lava_analytic(dt, h,      &
            T, rhom, speed_ratio, velocity_before,     &
            velocity_after)
    CASE ('rheology_lahar')
       CALL solve_rheology_lahar_analytic(dt, h,     &
            T, SUM(alphas), rhom, speed_ratio,&
            velocity_before, velocity_after)
    CASE ('rheology_empirical')
       CALL solve_rheology_empirical_analytic(dt, h, &
            speed_ratio, velocity_before, velocity_after)
    CASE ('rheology_quadratic')
       CALL solve_rheology_quadratic_analytic(dt, h, &
            speed_ratio, velocity_before, velocity_after)
    CASE ('rheology_froude')
       CALL solve_rheology_froude_analytic(dt, h,    &
            speed_ratio, velocity_before, red_grav,                  &
            grav_coeff, curvature_acceleration, use_stochastic,   &
            stochastic_value, velocity_after)
    CASE ('rheology_lucas')
       CALL solve_rheology_lucas_analytic(dt, h,     &
            speed_ratio, velocity_before, red_grav,                  &
            grav_coeff, curvature_acceleration, use_stochastic,   &
            stochastic_value, velocity_after)
    CASE ('rheology_mui')
       CALL solve_rheology_mui_analytic(dt, h,       &
            speed_ratio, velocity_before, rhom, red_grav, &
            grav_coeff, alphas, curvature_acceleration,   &
            use_pore_pressure, exc_pore_pres, velocity_after)
    CASE ('rheology_mui_voellmy')
       CALL solve_rheology_mui_voellmy_analytic(dt,           &
            h, speed_ratio, velocity_before, rhom,      &
            red_grav, grav_coeff, alphas,           &
            curvature_acceleration, use_pore_pressure, exc_pore_pres, &
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
  !> \param[in] u Cartesian depth-averaged x velocity [m s^-1].
  !> \param[in] v Cartesian depth-averaged y velocity [m s^-1].
  !> \return Signed bed-curvature acceleration [m s^-2].
  !******************************************************************************
  PURE FUNCTION curvature_acceleration(Bsecond_xx, Bsecond_xy, Bsecond_yy,  &
       u, v) RESULT(acceleration)

    REAL(wp), INTENT(IN) :: Bsecond_xx, Bsecond_xy, Bsecond_yy
    REAL(wp), INTENT(IN) :: u, v
    REAL(wp) :: acceleration

    acceleration = Bsecond_xx*u**2 +                              &
         2.0_wp*Bsecond_xy*u*v + Bsecond_yy*v**2

  END FUNCTION curvature_acceleration

  !******************************************************************************
  !> \brief Evaluate the selected implicit rheology momentum rate.
  !> \details All state values are explicit arguments so the dispatch has no
  !> per-cell workspace allocation or hidden cell-state conversion.
  !> \param[in] model_name Selected rheology identifier.
  !> \param[in] h Flow depth [m].
  !> \param[in] T Mixture T [K].
  !> \param[in] rhom Mixture density [kg m^-3].
  !> \param[in] inv_rhom Reciprocal mixture density [m^3 kg^-1].
  !> \param[in] u x velocity [m s^-1].
  !> \param[in] v y velocity [m s^-1].
  !> \param[in] mod_vel Total flow-mod_vel magnitude [m s^-1].
  !> \param[in] mod_hor_vel Horizontal flow-mod_vel magnitude [m s^-1].
  !> \param[in] mod_vel2 Squared flow-mod_vel magnitude [m^2 s^-2].
  !> \param[in] red_grav Reduced gravitational acceleration [m s^-2].
  !> \param[in] alphas Solid volume fractions, one per class.
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
  !> \param[in] h Flow depth [m].
  !> \param[in] u x velocity [m s^-1].
  !> \param[in] v y velocity [m s^-1].
  !> \param[in] mod_vel Total flow-mod_vel magnitude [m s^-1].
  !> \param[in] mod_hor_vel Horizontal flow-mod_vel magnitude [m s^-1].
  !> \param[in] rhom Mixture density [kg m^-3].
  !> \param[in] red_grav Reduced gravitational acceleration [m s^-2].
  !> \param[in] grav_coeff Effective gravity multiplier [dimensionless].
  !> \param[in] alphas Solid volume fractions, one per class.
  !> \param[in] curvature_term Enable curvature acceleration for applicable laws.
  !> \param[in] B_second_xx Bed curvature in xx [m^-1].
  !> \param[in] B_second_xy Bed curvature in xy [m^-1].
  !> \param[in] B_second_yy Bed curvature in yy [m^-1].
  !> \param[in] use_pore_pressure Enable the pore-pressure correction.
  !> \param[in] exc_pore_pres Excess pore pressure [Pa].
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
  !> \param[in] dt Exact-flow interval [s].
  !> \param[in] h Flow depth [m].
  !> \param[in] T Flow T [K].
  !> \param[in] rhom Mixture density [kg m^-3].
  !> \param[in] speed_ratio Total-to-horizontal mod_vel ratio.
  !> \param[in] velocity_before Horizontal velocity before integration [m s^-1].
  !> \param[in] red_grav Reduced gravity [m s^-2].
  !> \param[in] grav_coeff Gravity coefficient used by rheology laws.
  !> \param[in] alphas Per-class solid volume fractions.
  !> \param[in] curvature_acceleration Signed bed-curvature acceleration [m s^-2].
  !> \param[in] use_pore_pressure Whether excess pore pressure is active.
  !> \param[in] exc_pore_pres Excess pore pressure.
  !> \param[in] use_stochastic Whether stochastic forcing is active.
  !> \param[in] stochastic_value Local stochastic rheology value.
  !> \param[out] velocity_after Horizontal velocity after integration [m s^-1].
  !> \param[out] model_supported True when model_name has an analytic velocity law.
  !******************************************************************************




END MODULE rheology_2d
