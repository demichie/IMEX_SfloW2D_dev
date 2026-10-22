!********************************************************************************
!> \brief Temperature-dependent viscous and yield-stress lava rheology.
!> \details The implicit viscous drag uses a temperature-dependent kinematic
!> viscosity and a depth-averaged shear factor:
!> \f[
!> \gamma=\frac{3\nu_{\rm ref}}{h_*}
!>   \exp[-b(T-T_{\rm ref})],\qquad
!> \boldsymbol{S}_{v}=-\rho_m\gamma\boldsymbol{u}_h,
!> \f]
!> where \f$h_*=\max(h,10^{-10})\f$ in the depth branch selected by the
!> implementation. The semi-implicit yield-stress source is
!> \f$\boldsymbol{S}_{y}=-\tau_0\boldsymbol{u}_h/|\boldsymbol{u}|\f$.
!> \ref eval_lava_implicit_source evaluates the complex-valued viscous source;
!> \ref eval_lava_semi_impl_source evaluates the real-valued yield source.
!> \ref nu_ref, \ref visc_par and \ref tau0 are the model parameters. Thermal
!> cooling is owned by the thermal-source module, not by this rheology.
!> \ref solve_lava_analytic_velocity integrates the frozen affine velocity
!> source exactly with temperature and other non-velocity inputs explicit.
!********************************************************************************

MODULE rheology_lava_2d

  USE parameters_2d, ONLY : wp
  USE rheology_parameters_2d, ONLY : T_ref
  USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs,  &
       fatal_analytic_input_error, solve_exact_velocity_flow, EXACT_FLOW_AFFINE

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nu_ref, visc_par, tau0
  PUBLIC :: eval_lava_implicit_source, eval_lava_semi_impl_source
  PUBLIC :: solve_lava_analytic_velocity

  REAL(wp) :: nu_ref
  REAL(wp) :: visc_par
  REAL(wp) :: tau0

CONTAINS

  !******************************************************************************
  !> \brief Integrate the lava-rheology velocity source exactly.
  !> \details Holds depth, temperature, density, and speed ratio fixed. The
  !> affine speed-decay law is integrated exactly, with stopping at zero.
  !> \param integration_step Exact-flow interval, normally dt * a_ii.
  !> \param flow_depth Explicit flow depth.
  !> \param temperature Explicit temperature.
  !> \param mixture_density Explicit mixture density.
  !> \param speed_ratio Explicit total-to-horizontal speed ratio.
  !> \param velocity_before Input horizontal velocity.
  !> \param velocity_after Returned horizontal velocity.
  !******************************************************************************
  SUBROUTINE solve_lava_analytic_velocity(integration_step, flow_depth,     &
       temperature, mixture_density, speed_ratio, velocity_before,           &
       velocity_after)

    REAL(wp), INTENT(IN) :: integration_step, flow_depth, temperature
    REAL(wp), INTENT(IN) :: mixture_density, speed_ratio
    REAL(wp), INTENT(IN) :: velocity_before(2)
    REAL(wp), INTENT(OUT) :: velocity_after(2)

    REAL(wp), PARAMETER :: h_threshold = 1.0E-10_wp
    REAL(wp) :: speed_before, speed_after, effective_depth, gamma
    REAL(wp) :: linear_drag, yield_drag

    velocity_after = velocity_before
    speed_before = SQRT(SUM(velocity_before**2))
    CALL validate_analytic_velocity_inputs(integration_step, flow_depth,    &
         speed_ratio, speed_before, 'lava')
    IF (integration_step .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
    IF (mixture_density .LE. 0.0_wp) THEN
       CALL fatal_analytic_input_error('Nonpositive mixture density in lava analytic solver')
    END IF

    effective_depth = MAX(flow_depth, h_threshold)
    gamma = 3.0_wp * nu_ref / effective_depth *                             &
         EXP(-visc_par * (temperature - T_ref))
    linear_drag = gamma / flow_depth
    yield_drag = tau0 / (mixture_density * flow_depth * speed_ratio)
    CALL solve_exact_velocity_flow(EXACT_FLOW_AFFINE,                       &
         [yield_drag, linear_drag], integration_step, speed_before,         &
         speed_after)
    velocity_after = (speed_after / speed_before) * velocity_before

  END SUBROUTINE solve_lava_analytic_velocity

  !******************************************************************************
  !> \brief Evaluate the implicit temperature-dependent viscous drag.
  !> \details Computes the depth-regularized drag coefficient and returns the
  !> x/y momentum sources for a moving state. Complex inputs preserve
  !> temperature, depth and velocity sensitivities in the implicit Jacobian.
  !> \param flow_depth Flow depth.
  !> \param temperature Lava temperature.
  !> \param mixture_density Mixture density.
  !> \param velocity_x x component of horizontal velocity.
  !> \param velocity_y y component of horizontal velocity.
  !> \param speed Speed magnitude (part of the common source interface).
  !> \param moving Whether the state is moving.
  !> \param source_x Returned x momentum source.
  !> \param source_y Returned y momentum source.
  !******************************************************************************
  SUBROUTINE eval_lava_implicit_source(flow_depth, temperature, mixture_density, &
       velocity_x, velocity_y, speed, moving, source_x, source_y)

    COMPLEX(wp), INTENT(IN) :: flow_depth, temperature, mixture_density
    COMPLEX(wp), INTENT(IN) :: velocity_x, velocity_y, speed
    LOGICAL, INTENT(IN) :: moving
    COMPLEX(wp), INTENT(OUT) :: source_x, source_y

    COMPLEX(wp) :: gamma
    REAL(wp), PARAMETER :: h_threshold = 1.0E-10_wp

    source_x = CMPLX(0.0_wp, 0.0_wp, wp)
    source_y = CMPLX(0.0_wp, 0.0_wp, wp)
    IF (REAL(flow_depth, wp) .GT. h_threshold) THEN
       gamma = 3.0_wp * nu_ref / flow_depth * EXP(-visc_par * (temperature - T_ref))
    ELSE
       gamma = 3.0_wp * nu_ref / h_threshold * EXP(-visc_par * (temperature - T_ref))
    END IF

    IF (moving) THEN
       source_x = -mixture_density * gamma * velocity_x
       source_y = -mixture_density * gamma * velocity_y
    END IF

  END SUBROUTINE eval_lava_implicit_source

  !******************************************************************************
  !> \brief Evaluate the semi-implicit lava yield-stress source.
  !> \details Returns \f$-\tau_0(u_x,u_y)/|\boldsymbol{u}|\f$ for positive
  !> speed and zero otherwise.
  !> \param velocity_x x component of horizontal velocity.
  !> \param velocity_y y component of horizontal velocity.
  !> \param speed Speed magnitude.
  !> \param source_x Returned x momentum source.
  !> \param source_y Returned y momentum source.
  !******************************************************************************
  SUBROUTINE eval_lava_semi_impl_source(velocity_x, velocity_y, speed, source_x, source_y)

    REAL(wp), INTENT(IN) :: velocity_x, velocity_y, speed
    REAL(wp), INTENT(OUT) :: source_x, source_y

    source_x = 0.0_wp
    source_y = 0.0_wp
    IF (speed .LE. 0.0_wp) RETURN

    source_x = -tau0 * velocity_x / speed
    source_y = -tau0 * velocity_y / speed

  END SUBROUTINE eval_lava_semi_impl_source

END MODULE rheology_lava_2d
