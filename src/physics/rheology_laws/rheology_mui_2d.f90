!********************************************************************************
!> \brief Regularized \f$\mu(I)\f$ granular basal-friction model (model 11).
!> \details The inertial number is formed from particle diameter \f$d\f$,
!> shear rate \f$\dot\gamma=(5/2)|\boldsymbol{u}|/(h_{\rm eff}\sqrt{G})\f$
!> and effective normal stress
!> \f$\sigma'_n=\max(\epsilon,\sqrt{G}\max(0,\rho_m(g'+a_c)h_{\rm eff}
!> -\mathbf{1}_p p_e))\f$, where \f$h_{\rm eff}\f$ is depth regularized by
!> machine epsilon and \f$\mathbf{1}_p\f$ enables pore-pressure reduction:
!> \f[
!> I=\frac{d\dot\gamma}{\sqrt{\sigma'_n/\rho_p}},\qquad
!> \mu(I)=\frac{\mu_s I_0+\mu_2 I+\mu_{I,\infty}I^2}{I_0+I}.
!> \f]
!> The basal stress is
!> \f$\tau_b=\mu(I)\sigma'_n\f$ and the returned source is
!> \f$\boldsymbol{S}=-\tau_b\boldsymbol{u}_h/|\boldsymbol{u}|\f$.
!> \ref eval_rheology_mui_semi_implicit evaluates the real-valued x/y momentum rates,
!> choosing a Sauter diameter and average solid density for mixtures with
!> multiple solid classes. \ref solve_rheology_mui_analytic solves its
!> rational mod_vel-dependent frozen-coefficient source in time.
!********************************************************************************
MODULE rheology_mui_2d
    USE diagnostics_2d, ONLY : fatal_error
    USE parameters_2d, ONLY : wp, n_solid
    USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
    USE constitutive_parameters_2d, ONLY : diam_s, rho_s
    USE rheology_parameters_2d, ONLY : I_0, mu_2, mu_s, muI_inf
    USE state_conversion_2d, ONLY : sauter_diameter, average_density_solids
    USE rheology_analytic_utils_2d, ONLY : validate_analytic_velocity_inputs,                      &
            fatal_analytic_input_error, solve_exact_velocity_flow,                                 &
            EXACT_FLOW_MUI
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_rheology_mui_explicit
    PUBLIC :: eval_rheology_mui_semi_implicit
    PUBLIC :: eval_rheology_mui_implicit
    PUBLIC :: solve_rheology_mui_analytic
    ! Uniform source treatment entry points. Existing kernels remain unchanged.
CONTAINS
    !> \param[in] state Shared real or complex physical state.
    !> \brief Evaluate the explicit mui contribution to the horizontal momentum source.
    !> \details Compute the momentum source under the regularized mu(I) friction law, including
    !> effective stress and granular inertial-number dependence, using the treatment-specific
    !> velocity and effective-stress inputs. Outputs are conservative x/y momentum source rates.
    !> \param[in] mod_vel Three-dimensional speed [m s^-1].
    !> \param[in] grav_coeff Local dimensionless gravitational correction.
    !> \param[in] curvature_term Enable curvature correction.
    !> \param[in] B_second_xx Topographic second derivative in x [m^-1].
    !> \param[in] B_second_xy Mixed topographic derivative [m^-1].
    !> \param[in] B_second_yy Topographic second derivative in y [m^-1].
    !> \param[in] use_pore_pressure Enable pore-pressure correction.
    !> \param[out] rate Two-component conservative momentum source.
    SUBROUTINE eval_rheology_mui_explicit(state, &
            mod_vel, &
            grav_coeff, &
            curvature_term, &
            B_second_xx, &
            B_second_xy, &
            B_second_yy, &
            use_pore_pressure, &
            rate)
        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: mod_vel, grav_coeff
        REAL(wp), INTENT(IN) :: B_second_xx, B_second_xy, B_second_yy
        LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure
        REAL(wp), INTENT(OUT) :: rate(2)
        REAL(wp) :: sx, sy, curvature_accel
        COMPLEX(wp) :: ix, iy
        TYPE(complex_physical_state_type) :: cs
        curvature_accel=0.0_wp
        IF (curvature_term) curvature_accel=B_second_xx*state%u**2 +                                     &
                2.0_wp*B_second_xy*state%u*state%v+B_second_yy*state%v**2
        CALL eval_rheology_mui_semi_implicit(state, mod_vel, grav_coeff, curvature_accel, &
             use_pore_pressure, sx, sy)
        CALL eval_rheology_mui_implicit(cs, ix, iy)
        rate=[sx,sy]+REAL([ix,iy],wp)
    END SUBROUTINE eval_rheology_mui_explicit

    !> \param[in] state Shared real or complex physical state.
    !> \brief Evaluate the semi-implicit mui contribution to the horizontal momentum source.
    !> \details Compute the momentum source under the regularized mu(I) friction law, including
    !> effective stress and granular inertial-number dependence, using the treatment-specific
    !> velocity and effective-stress inputs. Outputs are conservative x/y momentum source rates.
    !> \param[in] mod_vel Three-dimensional speed [m s^-1].
    !> \param[in] grav_coeff Local dimensionless gravitational correction.
    !> \param[in] curvature_acceleration Normal acceleration due to topographic curvature [m s^-2].
    !> \param[in] use_pore_pressure Enable pore-pressure correction.
    !> \param[out] rate_x x momentum source.
    !> \param[out] rate_y y momentum source.
    SUBROUTINE eval_rheology_mui_semi_implicit(state, mod_vel, grav_coeff, curvature_acceleration, &
            use_pore_pressure, rate_x, rate_y)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: mod_vel, grav_coeff, curvature_acceleration
        LOGICAL, INTENT(IN) :: use_pore_pressure
        REAL(wp), INTENT(OUT) :: rate_x, rate_y
        REAL(wp) :: particle_diameter
        REAL(wp) :: particle_density
        REAL(wp) :: effective_depth
        REAL(wp) :: effective_vertical_stress
        REAL(wp) :: effective_normal_stress
        REAL(wp) :: shear_rate
        REAL(wp) :: inertial_number
        REAL(wp) :: friction_coefficient
        REAL(wp) :: friction_stress

        ASSOCIATE(h => state%h, &
                  u => state%u, &
                  v => state%v, &
                  rhom => state%rho_m, &
                  red_grav => state%red_grav, &
                  solid_fraction => state%alphas, &
                  exc_pore_pres => state%exc_pore_pres)
        rate_x = 0.0_wp
        rate_y = 0.0_wp
        IF ( mod_vel .LE. 0.0_wp ) RETURN
        IF ( n_solid .GT. 1 ) THEN
            particle_diameter = sauter_diameter(solid_fraction)
            particle_density = average_density_solids(solid_fraction)

        ELSE
            particle_diameter = diam_s(1)
            particle_density = rho_s(1)

        END IF

        effective_depth = h
        IF ( effective_depth .LT. EPSILON(1.0_wp) ) THEN
            effective_depth = EPSILON(1.0_wp)

        END IF

        effective_vertical_stress = rhom *                                                         &
                ( red_grav + curvature_acceleration ) * effective_depth
        IF ( use_pore_pressure ) THEN
            effective_vertical_stress = effective_vertical_stress -                                &
                    exc_pore_pres

        END IF

        effective_vertical_stress = MAX(effective_vertical_stress, 0.0_wp)
        effective_normal_stress = MAX(0.0_wp, effective_vertical_stress *                          &
                SQRT(grav_coeff))
        IF ( effective_normal_stress .LT. EPSILON(1.0_wp) ) THEN
            effective_normal_stress = EPSILON(1.0_wp)

        END IF

        shear_rate = 5.0_wp/2.0_wp * mod_vel /                                                     &
                (effective_depth * SQRT(grav_coeff))
        inertial_number = particle_diameter * shear_rate /                                         &
                SQRT(effective_normal_stress / particle_density)
        friction_coefficient = (mu_s * I_0 + mu_2 * inertial_number +                              &
                muI_inf * inertial_number**2) / (I_0 + inertial_number)
        friction_stress = friction_coefficient * effective_normal_stress
        rate_x = rate_x - friction_stress * (u / mod_vel)
        rate_y = rate_y - friction_stress * (v / mod_vel)
    
        END ASSOCIATE
    END SUBROUTINE eval_rheology_mui_semi_implicit

    !> \param[in] state Shared real or complex physical state.
    !> \brief Evaluate the implicit mui contribution to the horizontal momentum source.
    !> \details Return the contribution assigned to the implicit treatment. This implementation
    !> supplies a zero source when the model has no term in this treatment.
    !> \param[out] rate_x x momentum source.
    !> \param[out] rate_y y momentum source.
    SUBROUTINE eval_rheology_mui_implicit(state, rate_x, rate_y)

        TYPE(complex_physical_state_type), INTENT(IN) :: state
        COMPLEX(wp), INTENT(OUT) :: rate_x, rate_y

        rate_x = CMPLX(0.0_wp,0.0_wp,wp)
        rate_y = CMPLX(0.0_wp,0.0_wp,wp)
    END SUBROUTINE eval_rheology_mui_implicit

    !> \param[in] state Shared real physical state before integration.
    !> \brief Integrate the mui velocity source in time.
    !> \details Advance the scalar speed under the regularized mu(I) friction law, including
    !> effective stress and granular inertial-number dependence with model coefficients fixed over
    !> dt; the horizontal direction is retained in the returned velocity vector.
    !> \param[in] dt Source integration interval [s].
    !> \param[in] speed_ratio Ratio of total to horizontal speed.
    !> \param[in] grav_coeff Local dimensionless gravitational correction.
    !> \param[in] curvature_acceleration Normal acceleration due to topographic curvature [m s^-2].
    !> \param[in] use_pore_pressure Enable pore-pressure correction.
    !> \param[out] velocity_after Horizontal velocity vector after integration [m s^-1].
    SUBROUTINE solve_rheology_mui_analytic(dt, state, speed_ratio, grav_coeff, curvature_acceleration, &
            use_pore_pressure, velocity_after)

        TYPE(physical_state_type), INTENT(IN) :: state
        REAL(wp), INTENT(IN) :: dt, speed_ratio, grav_coeff, curvature_acceleration
        LOGICAL, INTENT(IN) :: use_pore_pressure
        REAL(wp), INTENT(OUT) :: velocity_after(2)
        REAL(wp) :: velocity_before(2)
        REAL(wp) :: speed_before, speed_after, particle_diameter, particle_density
        REAL(wp) :: effective_depth, vertical_stress, normal_stress
        REAL(wp) :: inertial_slope, basal_factor, sqrt_gravity

        ASSOCIATE(h => state%h, &
                  rhom => state%rho_m, &
                  red_grav => state%red_grav, &
                  solid_fraction => state%alphas, &
                  exc_pore_pres => state%exc_pore_pres)
        velocity_before = [state%u,state%v]
        velocity_after = velocity_before
        speed_before = SQRT(SUM(velocity_before**2))
        CALL validate_analytic_velocity_inputs(dt, h,                                              &
                speed_ratio, speed_before, 'mu(I)')
        IF (dt .EQ. 0.0_wp .OR. speed_before .EQ. 0.0_wp) RETURN
        IF (rhom .LE. 0.0_wp .OR. grav_coeff .LE. 0.0_wp) THEN
            CALL fatal_analytic_input_error(                                                       &
                'Nonpositive density or gravity multiplier in mu(I) analytic solver')

        END IF

        IF (n_solid .GT. 1) THEN
            particle_diameter = sauter_diameter(solid_fraction)
            particle_density = average_density_solids(solid_fraction)

        ELSE
            particle_diameter = diam_s(1)
            particle_density = rho_s(1)

        END IF

        IF (particle_diameter .LE. 0.0_wp .OR. particle_density .LE. 0.0_wp) THEN
            CALL fatal_analytic_input_error(                                                       &
                'Nonpositive particle property in mu(I) analytic solver')

        END IF

        IF (I_0 .LT. 0.0_wp) THEN
            CALL fatal_analytic_input_error('Negative I_0 in mu(I) analytic solver')

        END IF

        effective_depth = MAX(h, EPSILON(1.0_wp))
        vertical_stress = rhom *                                                                   &
                (red_grav + curvature_acceleration) * effective_depth
        IF (use_pore_pressure) vertical_stress = vertical_stress -                                 &
                exc_pore_pres
        vertical_stress = MAX(vertical_stress, 0.0_wp)
        normal_stress = MAX(EPSILON(1.0_wp),                                                       &
                vertical_stress * SQRT(grav_coeff))
        sqrt_gravity = SQRT(grav_coeff)
        inertial_slope = 2.5_wp * particle_diameter * speed_ratio /                                &
                (effective_depth * sqrt_gravity *                                                  &
                SQRT(normal_stress/particle_density))
        basal_factor = normal_stress / (rhom*h*speed_ratio)
        CALL solve_exact_velocity_flow(EXACT_FLOW_MUI,                                             &
                [basal_factor*mu_s*I_0, basal_factor*mu_2*inertial_slope,                          &
                basal_factor*muI_inf*inertial_slope**2, I_0, inertial_slope],                      &
                dt, speed_before, speed_after)
        velocity_after = (speed_after / speed_before) * velocity_before
    
        END ASSOCIATE
    END SUBROUTINE solve_rheology_mui_analytic

    ! Unsupported physical treatment, rejected by the source catalog.
    ! Valid zero component of the configured Newton treatment.
END MODULE rheology_mui_2d