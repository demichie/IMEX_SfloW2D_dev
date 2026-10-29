!> @brief Evaluate erosion process contributions to mass, energy and bed evolution.
!> @details Explicit treatment is supported; other treatment interfaces reject
!> unsupported configurations. Shared resource budgets belong to the coordinator.
MODULE mass_exchange_erosion_2d
    USE parameters_2d, ONLY : idx_h, idx_T
    USE parameters_2d, ONLY : wp, n_solid, n_vars, idx_u, idx_v, slope_correction_flag, &
              erodible_deposit_flag, mass_exchange_rates_type, reset_mass_exchange_rates
    USE diagnostics_2d, ONLY : fatal_error
    USE state_conversion_2d, ONLY : mixt_var
    USE constitutive_parameters_2d, ONLY : erosion_coeff, erodible_porosity, &
              erodible_fract, alphastot_min, rho_s, rho_c_sub, sp_heat_s, T_erodible
    IMPLICIT NONE
    PRIVATE

    PUBLIC :: eval_mass_exchange_erosion_explicit
    PUBLIC :: eval_mass_exchange_erosion_semi_implicit
    PUBLIC :: eval_mass_exchange_erosion_implicit
    PUBLIC :: solve_mass_exchange_erosion_analytic
    PUBLIC :: eval_erosion_rates
    ! Four-entry physical treatment contract; only EXPLICIT is supported.

CONTAINS

    !> @brief Evaluate explicit erosion from erodible substrate.
    !> @param[in] qp Local primitive flow state; component order follows model configuration.
    !> @param[in] Bx Bed elevation gradient in the x direction [dimensionless].
    !> @param[in] By Bed elevation gradient in the y direction [dimensionless].
    !> @param[in] alpha Total dispersed-solid volume fraction [dimensionless].
    !> @param[in] erodible Available erodible solid inventory by species [m].
    !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_erosion_explicit(qp,Bx,By,alpha,erodible,dt,r)
        REAL(wp), INTENT(IN) :: qp(n_vars+2), Bx, By, alpha
        REAL(wp), INTENT(IN) :: erodible(n_solid), dt
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
        REAL(wp) :: h, u, v, T_flow
        REAL(wp) :: normal_speed, vertical_velocity, substrate_temperature
        REAL(wp) :: Ri,rhom,rhoc,redgrav,cp_sub,cpm

        CALL reset_mass_exchange_rates(r)
        IF (alpha<=alphastot_min) RETURN
        h=qp(idx_h)
        u=qp(idx_u)
        v=qp(idx_v)
        T_flow=qp(idx_T)
        vertical_velocity=0.0_wp
        IF (slope_correction_flag) vertical_velocity=u*Bx+v*By
        normal_speed=SQRT(u*u+v*v+vertical_velocity*vertical_velocity)

        ! Erodible inventory limits species withdrawal before constructing mass and bed terms.
        CALL eval_erosion_volume_kernel(h,normal_speed,alpha,erodible,dt, &
                  r%solid_erosion_rate,r%carrier_erosion_rate)
        r%solid_mass_rate=rho_s*r%solid_erosion_rate
        r%mixture_mass_rate=SUM(r%solid_mass_rate)+rho_c_sub*r%carrier_erosion_rate
        CALL mixt_var(qp,Ri,rhom,rhoc,redgrav,cp_sub,cpm)
        substrate_temperature=T_erodible
        IF (erodible_deposit_flag) substrate_temperature=T_flow
        r%energy_rate=substrate_temperature*(SUM(sp_heat_s*r%solid_mass_rate)+ &
                  rho_c_sub*cp_sub*r%carrier_erosion_rate)
        r%dB_dt=-SUM(r%solid_erosion_rate)/(1.0_wp-erodible_porosity)
    END SUBROUTINE

    !> @brief Reject unsupported semi-implicit treatment for erosion.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] qp Local primitive flow state; component order follows model configuration.
    !> @param[in] Bx Bed elevation gradient in the x direction [dimensionless].
    !> @param[in] By Bed elevation gradient in the y direction [dimensionless].
    !> @param[in] alpha Total dispersed-solid volume fraction [dimensionless].
    !> @param[in] erodible Available erodible solid inventory by species [m].
    !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_erosion_semi_implicit(qp,Bx,By,alpha,erodible,dt,r)
        REAL(wp), INTENT(IN) :: qp(n_vars+2), Bx, By, alpha
        REAL(wp), INTENT(IN) :: erodible(n_solid), dt
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('erosion: unsupported physical treatment')
    END SUBROUTINE eval_mass_exchange_erosion_semi_implicit

    !> @brief Reject unsupported implicit treatment for erosion from erodible substrate.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] qp Local primitive flow state; component order follows model configuration.
    !> @param[in] Bx Bed elevation gradient in the x direction [dimensionless].
    !> @param[in] By Bed elevation gradient in the y direction [dimensionless].
    !> @param[in] alpha Total dispersed-solid volume fraction [dimensionless].
    !> @param[in] erodible Available erodible solid inventory by species [m].
    !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_erosion_implicit(qp,Bx,By,alpha,erodible,dt,r)
        REAL(wp), INTENT(IN) :: qp(n_vars+2), Bx, By, alpha
        REAL(wp), INTENT(IN) :: erodible(n_solid), dt
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('erosion: unsupported physical treatment')
    END SUBROUTINE eval_mass_exchange_erosion_implicit

    !> @brief Reject unsupported analytic treatment for erosion from erodible substrate.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] qp Local primitive flow state; component order follows model configuration.
    !> @param[in] Bx Bed elevation gradient in the x direction [dimensionless].
    !> @param[in] By Bed elevation gradient in the y direction [dimensionless].
    !> @param[in] alpha Total dispersed-solid volume fraction [dimensionless].
    !> @param[in] erodible Available erodible solid inventory by species [m].
    !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE solve_mass_exchange_erosion_analytic(qp,Bx,By,alpha,erodible,dt,r)
        REAL(wp), INTENT(IN) :: qp(n_vars+2), Bx, By, alpha
        REAL(wp), INTENT(IN) :: erodible(n_solid), dt
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('erosion: unsupported physical treatment')
    END SUBROUTINE solve_mass_exchange_erosion_analytic

    !> @brief Compute erosion volumes under the finite erodible-bed inventory constraint.
    !> @param[in] r_h Mobile flow thickness [m].
    !> @param[in] mod_vel Slope-corrected flow speed [m/s].
    !> @param[in] alphas_tot Total solid volume fraction.
    !> @param[in] erodible Available erodible solid inventory by species [m].
    !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
    !> @param[out] erosion_rate Solid erosion volume fluxes by species [m/s].
    !> @param[out] carrier_erosion_rate Companion carrier erosion volume flux [m/s].
    SUBROUTINE eval_erosion_volume_kernel(r_h, mod_vel, alphas_tot, erodible, dt, &
              erosion_rate, carrier_erosion_rate)
        REAL(wp), INTENT(IN) :: r_h, mod_vel, alphas_tot, erodible(n_solid), dt
        REAL(wp), INTENT(OUT) :: erosion_rate(n_solid), carrier_erosion_rate

        IF (erosion_coeff .GT. 0.0_wp) THEN
              CALL eval_erosion_rates(r_h, mod_vel, alphas_tot, erosion_coeff,       &
                        erodible_porosity, erodible_fract(1:n_solid),                    &
                        erodible(1:n_solid), dt, erosion_rate(1:n_solid),                &
                        carrier_erosion_rate)

        ELSE

              erosion_rate = 0.0_wp
              carrier_erosion_rate = 0.0_wp
        END IF

    END SUBROUTINE eval_erosion_volume_kernel

    !> @brief Apply the empirical erosion flux and cap solid removal by available bed material.
    !> @param[in] flow_depth Mobile flow thickness [m].
    !> @param[in] speed Flow speed used by the erosion closure [m/s].
    !> @param[in] total_solid_fraction Total solid volume fraction.
    !> @param[in] erosion_coefficient Empirical coefficient in the erosion closure.
    !> @param[in] erodible_porosity Porosity of the erodible substrate [dimensionless].
    !> @param[in] erodible_fraction Species composition of the erodible bed [dimensionless].
    !> @param[in] erodible_inventory Available substrate inventory by species [m].
    !> @param[in] time_step Physical time-step duration [s].
    !> @param[out] solid_erosion_rate Capped solid erosion volume fluxes by species [m/s].
    !> @param[out] carrier_erosion_rate Companion carrier erosion volume flux [m/s].
    PURE SUBROUTINE eval_erosion_rates(flow_depth, speed, total_solid_fraction, &
              erosion_coefficient, erodible_porosity, erodible_fraction,            &
              erodible_inventory, time_step, solid_erosion_rate,                    &
              carrier_erosion_rate)

        REAL(wp), INTENT(IN) :: flow_depth, speed, total_solid_fraction
        REAL(wp), INTENT(IN) :: erosion_coefficient, erodible_porosity
        REAL(wp), INTENT(IN) :: erodible_fraction(:), erodible_inventory(:)
        REAL(wp), INTENT(IN) :: time_step
        REAL(wp), INTENT(OUT) :: solid_erosion_rate(:)
        REAL(wp), INTENT(OUT) :: carrier_erosion_rate

        REAL(wp) :: total_erosion, total_solid_erosion

        IF (erosion_coefficient .GT. 0.0_wp) THEN
              total_erosion = erosion_coefficient * speed * flow_depth *            &
                        (1.0_wp - total_solid_fraction)
              total_solid_erosion = total_erosion * (1.0_wp - erodible_porosity)
              solid_erosion_rate = erodible_fraction * total_solid_erosion

        ELSE

              solid_erosion_rate = 0.0_wp
        END IF

        ! Restrict each species independently: the bed cannot supply more than it contains.
        solid_erosion_rate = MAX(0.0_wp, MIN(solid_erosion_rate,                 &
                  erodible_inventory / time_step))
        total_solid_erosion = SUM(solid_erosion_rate)
        total_erosion = total_solid_erosion / (1.0_wp - erodible_porosity)
        carrier_erosion_rate = total_erosion * erodible_porosity

    END SUBROUTINE eval_erosion_rates

END MODULE mass_exchange_erosion_2d
