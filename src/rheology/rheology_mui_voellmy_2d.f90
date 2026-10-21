!********************************************************************************
!> \brief Combined regularized \f$\mu(I)\f$ and Voellmy model (model 12).
!> \details The semi-implicit basal term uses the regularized granular law
!> \f[
!> I=\frac{d\dot\gamma}{\sqrt{\sigma'_n/\rho_p}},\quad
!> \mu(I)=\frac{\mu_s I_0+\mu_2 I+\mu_{I,\infty}I^2}{I_0+I},\quad
!> \tau_b=\mu(I)\sigma'_n,
!> \f]
!> with \f$\dot\gamma=(5/2)|\boldsymbol{u}|/(h_{\rm eff}\sqrt{G})\f$ and
!> \f$\sigma'_n=\max(\epsilon,\sqrt{G}\max(0,\rho_m(g'+a_c)h_{\rm eff}
!> -\mathbf{1}_p p_e))\f$. Here \f$h_{\rm eff}\f$ is depth regularized by
!> machine epsilon, and \f$\mathbf{1}_p\f$ enables pore-pressure reduction.
!> The basal source is
!> \f$\boldsymbol{S}_b=-\tau_b\boldsymbol{u}_h/|\boldsymbol{u}|\f$. The
!> implicit turbulent term is
!> \f$\boldsymbol{S}_t=-\rho_mg'|\boldsymbol{u}|^2
!> \boldsymbol{u}_h/(\xi|\boldsymbol{u}|)\f$. The module exposes
!> \ref eval_mui_voellmy_source for the regularized basal part and
!> \ref eval_mui_voellmy_implicit_source for the turbulent part. The granular
!> source uses Sauter particle diameter and average particle density when
!> multiple solid classes are present.
!********************************************************************************

MODULE rheology_mui_voellmy_2d

  USE parameters_2d, ONLY : wp, n_solid
  USE constitutive_parameters_2d, ONLY : diam_s, rho_s
  USE rheology_parameters_2d, ONLY : I_0, mu_2, mu_s, muI_inf, xi
  USE state_conversion_2d, ONLY : sauter_diameter, average_density_solids

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_mui_voellmy_source, eval_mui_voellmy_implicit_source

CONTAINS

  !******************************************************************************
  !> \brief Evaluate the implicit Voellmy turbulent source.
  !> \details Computes the velocity-squared stress with \f$\xi\f$ and returns
  !> its x/y components. Optional stochastic variation of \f$\xi\f$ is bounded
  !> to [1, 10^5]; a stationary state returns zero.
  !> \param velocity_x x component of horizontal velocity.
  !> \param velocity_y y component of horizontal velocity.
  !> \param mixture_density Mixture density.
  !> \param reduced_gravity Reduced gravity.
  !> \param speed Speed magnitude.
  !> \param speed_squared Squared speed.
  !> \param moving Whether the state is moving.
  !> \param use_stochastic_xi Enable stochastic \f$\xi\f$ variation.
  !> \param stochastic_value Stochastic increment to \f$\xi\f$.
  !> \param source_x Returned x momentum source.
  !> \param source_y Returned y momentum source.
  !******************************************************************************
  SUBROUTINE eval_mui_voellmy_implicit_source(velocity_x, velocity_y,        &
       mixture_density, reduced_gravity, speed, speed_squared, moving,       &
       use_stochastic_xi, stochastic_value, source_x, source_y)

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

  END SUBROUTINE eval_mui_voellmy_implicit_source

  !******************************************************************************
  !> \brief Evaluate model 12's regularized \f$\mu(I)\f$ basal source.
  !> \details Computes the inertial number and effective normal stress,
  !> including optional pore-pressure reduction, then returns the basal
  !> friction source in x and y. Stationary states return zero.
  !> \param flow_depth Flow depth.
  !> \param velocity_x x component of horizontal velocity.
  !> \param velocity_y y component of horizontal velocity.
  !> \param speed Speed magnitude.
  !> \param mixture_density Mixture density.
  !> \param reduced_gravity Reduced gravity.
  !> \param gravity_coefficient Gravity multiplier.
  !> \param solid_fraction Solid volume fractions by class.
  !> \param curvature_acceleration Curvature contribution to vertical stress.
  !> \param use_pore_pressure Enable pore-pressure reduction.
  !> \param excess_pore_pressure Excess pore pressure.
  !> \param source_x Returned x momentum source.
  !> \param source_y Returned y momentum source.
  !******************************************************************************

  SUBROUTINE eval_mui_voellmy_source( flow_depth, velocity_x, velocity_y,    &
       speed, mixture_density, reduced_gravity, gravity_coefficient,         &
       solid_fraction, curvature_acceleration, use_pore_pressure,             &
       excess_pore_pressure, source_x, source_y )

    REAL(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y, speed
    REAL(wp), INTENT(IN) :: mixture_density, reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient
    REAL(wp), INTENT(IN) :: solid_fraction(n_solid)
    REAL(wp), INTENT(IN) :: curvature_acceleration
    LOGICAL, INTENT(IN) :: use_pore_pressure
    REAL(wp), INTENT(IN) :: excess_pore_pressure
    REAL(wp), INTENT(OUT) :: source_x, source_y

    REAL(wp) :: particle_diameter
    REAL(wp) :: particle_density
    REAL(wp) :: effective_depth
    REAL(wp) :: effective_vertical_stress
    REAL(wp) :: effective_normal_stress
    REAL(wp) :: shear_rate
    REAL(wp) :: inertial_number
    REAL(wp) :: friction_coefficient
    REAL(wp) :: friction_stress

    source_x = 0.0_wp
    source_y = 0.0_wp
    IF ( speed .LE. 0.0_wp ) RETURN

    IF ( n_solid .GT. 1 ) THEN
       particle_diameter = sauter_diameter(solid_fraction)
       particle_density = average_density_solids(solid_fraction)
    ELSE
       particle_diameter = diam_s(1)
       particle_density = rho_s(1)
    END IF

    effective_depth = flow_depth
    IF ( effective_depth .LT. EPSILON(1.0_wp) ) THEN
       effective_depth = EPSILON(1.0_wp)
    END IF

    effective_vertical_stress = mixture_density *                          &
         ( reduced_gravity + curvature_acceleration ) * effective_depth
    IF ( use_pore_pressure ) THEN
       effective_vertical_stress = effective_vertical_stress -              &
            excess_pore_pressure
    END IF
    effective_vertical_stress = MAX(effective_vertical_stress, 0.0_wp)

    effective_normal_stress = MAX(0.0_wp, effective_vertical_stress *       &
         SQRT(gravity_coefficient))
    IF ( effective_normal_stress .LT. EPSILON(1.0_wp) ) THEN
       effective_normal_stress = EPSILON(1.0_wp)
    END IF

    shear_rate = 5.0_wp/2.0_wp * speed /                                    &
         (effective_depth * SQRT(gravity_coefficient))
    inertial_number = particle_diameter * shear_rate /                      &
         SQRT(effective_normal_stress / particle_density)
    friction_coefficient = (mu_s * I_0 + mu_2 * inertial_number +            &
         muI_inf * inertial_number**2) / (I_0 + inertial_number)

    friction_stress = friction_coefficient * effective_normal_stress
    source_x = source_x - friction_stress * (velocity_x / speed)
    source_y = source_y - friction_stress * (velocity_y / speed)

  END SUBROUTINE eval_mui_voellmy_source

END MODULE rheology_mui_voellmy_2d
