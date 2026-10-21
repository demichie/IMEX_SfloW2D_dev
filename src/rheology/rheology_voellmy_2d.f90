!********************************************************************************
!> \brief Voellmy-Salm basal-friction model.
!> \details The model combines a Coulomb basal contribution and a
!> velocity-squared turbulent contribution. The implicit turbulent source is
!> \f[
!> \boldsymbol{S}_t=-\frac{\rho_m g'}{\xi}|\boldsymbol{u}|^2
!>                  \frac{\boldsymbol{u}_h}{|\boldsymbol{u}|}.
!> \f]
!> The semi-implicit Coulomb stress is
!> \f[
!> \tau_b=\max\!\left(0,\,
!> \mu\rho_m h(g'+a_c)\sqrt{G}
!> -\mathbf{1}_{p}\mu\sqrt{G}\,p_e\right),\qquad
!> \boldsymbol{S}_b=-\tau_b\boldsymbol{u}_h/|\boldsymbol{u}|,
!> \f]
!> where \f$a_c\f$ is the optional curvature acceleration and the pore
!> pressure subtraction is enabled only when requested. The module provides
!> \ref eval_voellmy_implicit_source for the turbulent term and
!> \ref eval_voellmy_semi_impl_source for the Coulomb term; \ref mu and
!> \ref xi are their respective friction parameters. When enabled,
!> \f$a_c=c_{xx}u_x^2+2c_{xy}u_xu_y+c_{yy}u_y^2\f$.
!********************************************************************************

MODULE rheology_voellmy_2d

  USE parameters_2d, ONLY : wp
  USE rheology_parameters_2d, ONLY : xi

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mu, eval_voellmy_implicit_source, eval_voellmy_semi_impl_source

  REAL(wp) :: mu

CONTAINS

  !******************************************************************************
  !> \brief Evaluate the implicit turbulent Voellmy source.
  !> \details Computes the velocity-squared stress using \f$\xi\f$ and returns
  !> its x/y components opposite to the velocity direction. Optional
  !> stochastic variation of \f$\xi\f$ is bounded to [1, 10^5].
  !> \param velocity_x x component of horizontal velocity.
  !> \param velocity_y y component of horizontal velocity.
  !> \param mixture_density Mixture density \f$\rho_m\f$.
  !> \param reduced_gravity Reduced gravity \f$g'\f$.
  !> \param speed Speed magnitude.
  !> \param speed_squared Squared speed used in the turbulent stress.
  !> \param moving Whether the state is moving.
  !> \param use_stochastic_xi Enable the stochastic \f$\xi\f$ variation.
  !> \param stochastic_value Stochastic increment to \f$\xi\f$.
  !> \param source_x Returned x momentum source.
  !> \param source_y Returned y momentum source.
  !******************************************************************************
  SUBROUTINE eval_voellmy_implicit_source(velocity_x, velocity_y, mixture_density, &
       reduced_gravity, speed, speed_squared, moving, use_stochastic_xi,          &
       stochastic_value, source_x, source_y)

    COMPLEX(wp), INTENT(IN) :: velocity_x, velocity_y
    COMPLEX(wp), INTENT(IN) :: mixture_density, reduced_gravity
    COMPLEX(wp), INTENT(IN) :: speed, speed_squared
    LOGICAL, INTENT(IN) :: moving, use_stochastic_xi
    REAL(wp), INTENT(IN) :: stochastic_value
    COMPLEX(wp), INTENT(OUT) :: source_x, source_y

    REAL(wp) :: xi_effective
    COMPLEX(wp) :: turbulent_stress

    source_x = CMPLX(0.0_wp, 0.0_wp, wp)
    source_y = CMPLX(0.0_wp, 0.0_wp, wp)
    IF (.NOT. moving) RETURN

    IF (use_stochastic_xi) THEN
       xi_effective = xi + stochastic_value
       IF (xi_effective .LT. 1.0_wp) xi_effective = 1.0_wp
       IF (xi_effective .GT. 1.0E5_wp) xi_effective = 1.0E5_wp
    ELSE
       xi_effective = xi
    END IF

    turbulent_stress = mixture_density * reduced_gravity / xi_effective * speed_squared
    source_x = -turbulent_stress * (velocity_x / speed)
    source_y = -turbulent_stress * (velocity_y / speed)

  END SUBROUTINE eval_voellmy_implicit_source

  !******************************************************************************
  !> \brief Evaluate the semi-implicit Coulomb Voellmy source.
  !> \details Builds the basal stress from flow depth, reduced gravity and
  !> optional curvature acceleration, clips negative stress to zero, and
  !> optionally reduces it by excess pore pressure. The returned x/y sources
  !> oppose the horizontal flow direction.
  !> \param flow_depth Flow depth.
  !> \param velocity_x x component of horizontal velocity.
  !> \param velocity_y y component of horizontal velocity.
  !> \param speed Speed magnitude.
  !> \param mixture_density Mixture density.
  !> \param reduced_gravity Reduced gravity.
  !> \param gravity_coefficient Gravity multiplier.
  !> \param curvature_term Enable curvature acceleration.
  !> \param curvature_xx xx curvature coefficient.
  !> \param curvature_xy xy curvature coefficient.
  !> \param curvature_yy yy curvature coefficient.
  !> \param use_pore_pressure Enable pore-pressure stress reduction.
  !> \param excess_pore_pressure Excess pore pressure.
  !> \param source_x Returned x momentum source.
  !> \param source_y Returned y momentum source.
  !******************************************************************************
  SUBROUTINE eval_voellmy_semi_impl_source(flow_depth, velocity_x, velocity_y, &
       speed, mixture_density, reduced_gravity, gravity_coefficient,            &
       curvature_term, curvature_xx, curvature_xy, curvature_yy,                &
       use_pore_pressure, excess_pore_pressure, source_x, source_y)

    REAL(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y, speed
    REAL(wp), INTENT(IN) :: mixture_density, reduced_gravity, gravity_coefficient
    LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure
    REAL(wp), INTENT(IN) :: curvature_xx, curvature_xy, curvature_yy
    REAL(wp), INTENT(IN) :: excess_pore_pressure
    REAL(wp), INTENT(OUT) :: source_x, source_y

    REAL(wp) :: curvature_acceleration, friction_stress

    source_x = 0.0_wp
    source_y = 0.0_wp
    IF (speed .LE. 0.0_wp) RETURN

    curvature_acceleration = 0.0_wp
    IF (curvature_term) THEN
       curvature_acceleration = curvature_xx * velocity_x**2                &
            + 2.0_wp * curvature_xy * velocity_x * velocity_y                &
            + curvature_yy * velocity_y**2
    END IF

    friction_stress = mu * mixture_density *                                  &
         (reduced_gravity + curvature_acceleration) * flow_depth              &
         * SQRT(gravity_coefficient)
    friction_stress = MAX(0.0_wp, friction_stress)
    IF (use_pore_pressure) friction_stress = MAX(0.0_wp, friction_stress      &
         - mu * SQRT(gravity_coefficient) * excess_pore_pressure)

    source_x = -friction_stress * velocity_x / speed
    source_y = -friction_stress * velocity_y / speed

  END SUBROUTINE eval_voellmy_semi_impl_source

END MODULE rheology_voellmy_2d
