!> @brief Evaluate deposition process contributions to mass, energy and bed evolution.
!> @details Explicit treatment is supported; other treatment interfaces reject
!> unsupported configurations. Shared resource budgets belong to the coordinator.
MODULE mass_exchange_deposition_2d
    USE parameters_2d, ONLY : wp, n_solid, n_add_gas, gas_flag, sutherland_flag,                                       &
            mass_exchange_rates_type, reset_mass_exchange_rates
    USE constitutive_parameters_2d, ONLY : alphastot_min, kin_visc_c, muRef_Suth, Tref_Suth, S_mu, pres, &
        sp_gas_const_a, settling_flag, grav, collective_settling_flag, A_drag, B_drag
    USE solid_transport_parameters_2d, ONLY : rho_s, sp_heat_s, diam_s
    USE substrate_parameters_2d, ONLY : coeff_porosity

    USE diagnostics_2d, ONLY : fatal_error
    USE state_types, ONLY : physical_state_real, physical_state_complex
    USE rheology_lahar_2d, ONLY : eval_lahar_fluid_viscosity
    USE rheology_parameters_2d, ONLY : rheology_config_type

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: eval_mass_exchange_deposition_explicit
    PUBLIC :: eval_mass_exchange_deposition_semi_implicit
    PUBLIC :: eval_mass_exchange_deposition_implicit
    PUBLIC :: solve_mass_exchange_deposition_analytic
    PUBLIC :: settling_velocity
    PUBLIC :: eval_settling_fluid_properties
    PUBLIC :: eval_solid_deposition_rate
    PUBLIC :: eval_lahar_deposition_rate
    PUBLIC :: eval_hindered_settling_deposition_rate
        ! Four-entry physical treatment contract; only EXPLICIT is supported.

CONTAINS

        !> @brief Evaluate explicit sediment settling and deposition.
        !> @param[in] state Local real physical state.
        !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
        !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_deposition_explicit(config, state,dt,r)
            TYPE(rheology_config_type), INTENT(IN) :: config

            TYPE(physical_state_real), INTENT(IN) :: state
            REAL(wp), INTENT(IN) :: dt
            TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
            REAL(wp) :: alpha
            REAL(wp) :: gas_phase_mass_fraction

            CALL reset_mass_exchange_rates(r)
            alpha=SUM(state%alphas)
            IF (alpha<=alphastot_min) RETURN
            gas_phase_mass_fraction=0.0_wp

            IF (gas_flag) THEN
                gas_phase_mass_fraction=1.0_wp-SUM(state%xs)
                gas_phase_mass_fraction=MAX(0.0_wp,gas_phase_mass_fraction)

            END IF

            ! Evaluate settling before producing the conservative species, energy and bed terms.
            CALL eval_deposition_volume_kernel(config,state%h, state%T, state%alphas, alpha,                           &
                state%rho_m, state%rho_c, dt,                                                                          &
                r%solid_deposition_rate)
            r%carrier_deposition_rate=coeff_porosity*SUM(r%solid_deposition_rate)
            r%solid_mass_rate=-rho_s*r%solid_deposition_rate
            r%mixture_mass_rate=SUM(r%solid_mass_rate)-state%rho_c*r%carrier_deposition_rate
            r%energy_rate=state%T*(SUM(sp_heat_s*r%solid_mass_rate)-                                                   &
                        state%rho_c*state%cp_c*r%carrier_deposition_rate)
            IF (gas_flag .AND. n_add_gas>0 .AND. gas_phase_mass_fraction>0.0_wp)                                       &
                        r%add_gas_mass_rate= -state%rho_c*r%carrier_deposition_rate*                                   &
                        state%xg/gas_phase_mass_fraction
            r%dB_dt=SUM(r%solid_deposition_rate)+r%carrier_deposition_rate
    END SUBROUTINE eval_mass_exchange_deposition_explicit

        !> @brief Reject unsupported semi-implicit treatment for deposition.
        !> @details This entry point raises fatal_error; the rates are reset first.
        !> @param[in] state Local real physical state.
        !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
        !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_deposition_semi_implicit(state,dt,r)

            TYPE(physical_state_real), INTENT(IN) :: state
            REAL(wp), INTENT(IN) :: dt
            TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

            CALL reset_mass_exchange_rates(r)
            CALL fatal_error('deposition: unsupported physical treatment')
    END SUBROUTINE eval_mass_exchange_deposition_semi_implicit

        !> @brief Reject unsupported implicit treatment for sediment settling and deposition.
        !> @details This entry point raises fatal_error; the rates are reset first.
        !> @param[in] state Local real physical state.
        !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
        !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_deposition_implicit(state,dt,r)

            TYPE(physical_state_complex), INTENT(IN) :: state
            REAL(wp), INTENT(IN) :: dt
            TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

            CALL reset_mass_exchange_rates(r)
            CALL fatal_error('deposition: unsupported physical treatment')
    END SUBROUTINE eval_mass_exchange_deposition_implicit

        !> @brief Reject unsupported analytic treatment for sediment settling and deposition.
        !> @details This entry point raises fatal_error; the rates are reset first.
        !> @param[in] state Local real physical state.
        !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
        !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE solve_mass_exchange_deposition_analytic(state,dt,r)

            TYPE(physical_state_real), INTENT(IN) :: state
            REAL(wp), INTENT(IN) :: dt
            TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r

            CALL reset_mass_exchange_rates(r)
            CALL fatal_error('deposition: unsupported physical treatment')
    END SUBROUTINE solve_mass_exchange_deposition_analytic

        !> @brief Compute per-species settling and choose carrier-fluid properties.
        !> @param[in] r_h Mobile flow thickness [m].
        !> @param[in] r_T Local flow temperature [K].
        !> @param[in] r_alphas Solid volume fractions by species.
        !> @param[in] alphas_tot Total solid volume fraction.
        !> @param[in] r_rho_m Mixture density [kg/m3].
        !> @param[in] r_rho_c Carrier density [kg/m3].
        !> @param[in] dt Physical step duration [s]; must be positive where used as a divisor.
        !> @param[out] deposition_rate Solid deposition thickness flux for each species [m/s].
    SUBROUTINE eval_deposition_volume_kernel(config, r_h, r_T, r_alphas, alphas_tot, r_rho_m,                          &
        r_rho_c, dt, deposition_rate)
            TYPE(rheology_config_type), INTENT(IN) :: config

            REAL(wp), INTENT(IN) :: r_h, r_T, r_alphas(n_solid), alphas_tot
            REAL(wp), INTENT(IN) :: r_rho_m, r_rho_c, dt
            REAL(wp), INTENT(OUT) :: deposition_rate(n_solid)
            REAL(wp) :: alpha_max, hind_exp, kin_visc_c_local
            REAL(wp) :: dyn_visc_c, rho_c, rhoc, inv_kin_visc, settling_vel
            INTEGER :: i_solid

            alpha_max = 0.6_wp
            hind_exp = 4.65_wp
            kin_visc_c_local = kin_visc_c
            deposition_rate = 0.0_wp

            ! For the gas carrier, use temperature-dependent viscosity and ideal-gas density.
            IF ( gas_flag .AND. sutherland_flag ) THEN

                    dyn_visc_c = muRef_Suth * ( r_T / Tref_Suth )**1.5_wp *                                            &
                                ( Tref_Suth + S_mu ) / ( r_T + S_mu )

                    rho_c = pres / ( sp_gas_const_a * r_T )
                    kin_visc_c_local = dyn_visc_c / rho_c

            END IF

            CALL eval_settling_fluid_properties(config, r_T, alphas_tot, r_rho_m, r_rho_c, kin_visc_c_local,           &
                rhoc, inv_kin_visc)

            DO i_solid=1,n_solid

                    IF ( ( r_alphas(i_solid) .GT. 0.0_wp ) .AND. ( settling_flag ) ) THEN

                            settling_vel = settling_velocity( diam_s(i_solid) , rho_s(i_solid) ,                       &
                                        rhoc , inv_kin_visc )

                            deposition_rate(i_solid) = eval_solid_deposition_rate(                                     &
                                        config%rheology_model, r_alphas(i_solid),                                      &
                                            settling_vel, r_h, dt,                                                     &
                                        alphas_tot, alpha_max, hind_exp)

                            IF ( deposition_rate(i_solid) .LT. 0.0_wp ) THEN
                                    CALL fatal_error('negative solid deposition rate')

                            END IF

                    END IF

            END DO

    END SUBROUTINE eval_deposition_volume_kernel

        !> @brief Evaluate the terminal settling speed with continuous high-Re drag.
        !> @details Schiller-Naumann (with the optional collective correction)
        !>          is retained up to the physical particle Reynolds number 1000.
        !>          Beyond that point C_D is held at C_D(1000), continuously.
        !>          This is a constitutive choice; it replaces the former
        !>          discontinuous shortcut C_D=1 at the *initial* Reynolds.
        !> @param[in] diam Solid particle diameter [m].
        !> @param[in] rhos Solid particle density [kg/m3].
        !> @param[in] rhoc Carrier density [kg/m3].
        !> @param[in] inv_kin_visc Inverse carrier kinematic viscosity [s/m2].
        REAL(wp) FUNCTION settling_velocity(diam, rhos, rhoc, inv_kin_visc)

            REAL(wp), INTENT(IN) :: diam, rhos, rhoc, inv_kin_visc
            REAL(wp), PARAMETER :: transition_re=1000.0_wp
            REAL(wp), PARAMETER :: solver_tolerance=1.0E-11_wp
            REAL(wp) :: drag_transition, initial_re, initial_re_squared
            REAL(wp) :: const_part, reynolds, lower_re, upper_re
            REAL(wp) :: residual, derivative, candidate, reynolds_squared
            INTEGER :: iteration

            IF (rhos <= rhoc) THEN
                settling_velocity=0.0_wp
                RETURN
            END IF

            const_part=SQRT((4.0_wp/3.0_wp)*(rhos/rhoc-1.0_wp)*diam*grav)
            initial_re=diam*const_part*inv_kin_visc
            IF (initial_re<=0.0_wp) THEN
                settling_velocity=0.0_wp
                RETURN
            END IF

            ! C_D(Re) is Schiller-Naumann, optionally including collective
            ! drag, through the physical Reynolds number Re=1000.
            ! Above that point C_D stays constant at its transition value.
            drag_transition=24.0_wp*(1.0_wp+0.15_wp*transition_re**0.687_wp)/transition_re
            IF (collective_settling_flag) drag_transition=drag_transition+                                             &
                A_drag*(LOG(transition_re)-B_drag)
            IF (drag_transition<=0.0_wp) CALL fatal_error('settling: nonpositive drag at Re=1000')

            ! Re^2*C_D(Re)=Re0^2: compare the Reynolds of the solution,
            ! not the original C_D=1 estimate. The upper branch is explicit.
            IF (initial_re>=transition_re*SQRT(drag_transition)) THEN
                settling_velocity=const_part/SQRT(drag_transition)
                RETURN
            END IF

            ! Solve F(Re)=Re^2*C_D(Re)-Re0^2=0 on [0,1000].
            ! In this interval F(0)<0 and F(1000)>0. Schiller-Naumann
            ! has a unique root; the collective term is evaluated using its
            ! original logarithmic correction without approximating drag.
            initial_re_squared=initial_re*initial_re
            lower_re=0.0_wp
            upper_re=transition_re

            ! The Stokes and inertial asymptotes give a cheap starting point.
            ! Safeguarding below guarantees the root stays bracketed even if
            ! this estimate is inaccurate for collective settling.
            reynolds=MIN(transition_re,initial_re_squared/24.0_wp,                                                     &
                (initial_re_squared/3.6_wp)**(1.0_wp/1.687_wp))
            reynolds=MAX(reynolds,TINY(reynolds))

            ! Safeguarded Newton: use the exact derivative of Re^2*C_D(Re),
            ! reject nonphysical or out-of-bracket Newton steps, and bisect
            ! when necessary. The rate and all updates remain strictly local.
            DO iteration=1,64
                reynolds_squared=reynolds*reynolds
                residual=24.0_wp*reynolds*(1.0_wp+0.15_wp*reynolds**0.687_wp)                                          &
                    -initial_re_squared
                derivative=24.0_wp+6.0732_wp*reynolds**0.687_wp
                IF (collective_settling_flag) THEN
                    residual=residual+A_drag*reynolds_squared*(LOG(reynolds)-B_drag)
                    derivative=derivative+A_drag*reynolds*(2.0_wp*(LOG(reynolds)-B_drag)+1.0_wp)
                END IF
                IF (ABS(residual)<=solver_tolerance*initial_re_squared) EXIT
                IF (residual>0.0_wp) THEN
                    upper_re=reynolds
                ELSE
                    lower_re=reynolds
                END IF
                candidate=0.5_wp*(lower_re+upper_re)
                IF (derivative>0.0_wp) THEN
                    ! Newton is accepted only when it preserves the bracket.
                    ! This also avoids division by a nonpositive derivative.
                    IF (reynolds-residual/derivative>lower_re .AND.                                                    &
                        reynolds-residual/derivative<upper_re) candidate=reynolds-residual/derivative
                END IF
                reynolds=candidate
            END DO
            IF (iteration>64) CALL fatal_error('settling: safeguarded Newton did not converge')

            settling_velocity=reynolds/(diam*inv_kin_visc)
        END FUNCTION settling_velocity

        !> @brief Select carrier density and inverse kinematic viscosity for settling.
        !> @param[in] model_name Canonical rheology closure identifier.
        !> @param[in] temperature Material temperature [K].
        !> @param[in] total_solid_fraction Total solid volume fraction.
        !> @param[in] mixture_density Mixture mass density [kg/m3].
        !> @param[in] carrier_density_fallback Carrier density from mixture properties [kg/m3].
        !> @param[in] kinematic_viscosity_fallback Default carrier kinematic viscosity [m2/s].
        !> @param[out] carrier_density Carrier density selected for the settling calculation
        !> [kg/m3].
        !> @param[out] inverse_kinematic_viscosity Inverse selected kinematic viscosity [s/m2].
    SUBROUTINE eval_settling_fluid_properties(config, temperature, total_solid_fraction, mixture_density,              &
            carrier_density_fallback, kinematic_viscosity_fallback, carrier_density, inverse_kinematic_viscosity)
            TYPE(rheology_config_type), INTENT(IN) :: config

            REAL(wp), INTENT(IN) :: temperature, total_solid_fraction, mixture_density
            REAL(wp), INTENT(IN) :: carrier_density_fallback
            REAL(wp), INTENT(IN) :: kinematic_viscosity_fallback
            REAL(wp), INTENT(OUT) :: carrier_density, inverse_kinematic_viscosity
            REAL(wp) :: fluid_viscosity

            carrier_density = carrier_density_fallback

            IF (config%rheology_model .EQ. 'lahar') THEN
                    fluid_viscosity = eval_lahar_fluid_viscosity(config,                                               &
                        temperature,                                                                                   &
                                total_solid_fraction)
                    carrier_density = mixture_density
                    inverse_kinematic_viscosity = mixture_density / fluid_viscosity

            ELSE

                    inverse_kinematic_viscosity = 1.0_wp / kinematic_viscosity_fallback

            END IF

    END SUBROUTINE eval_settling_fluid_properties

        !> @brief Dispatch the selected solid settling/deposition closure.
        !> @param[in] model_name Canonical rheology closure identifier.
        !> @param[in] solid_fraction Volume fraction of the settling solid species.
        !> @param[in] settling_velocity Species terminal settling speed [m/s].
        !> @param[in] flow_depth Mobile flow thickness [m].
        !> @param[in] time_step Physical time-step duration [s].
        !> @param[in] total_solid_fraction Total solid volume fraction.
        !> @param[in] maximum_solid_fraction Packing concentration used by the hindrance closure.
        !> @param[in] hindered_settling_exponent Exponent controlling settling suppression at high
        !> concentration.
        !> @return Physical rate from the indicated constitutive closure.
    PURE FUNCTION eval_solid_deposition_rate(model_name, solid_fraction,                                               &
                    settling_velocity, flow_depth, time_step, total_solid_fraction,                                    &
                    maximum_solid_fraction, hindered_settling_exponent) RESULT(rate)

            CHARACTER(*), INTENT(IN) :: model_name
            REAL(wp), INTENT(IN) :: solid_fraction, settling_velocity
            REAL(wp), INTENT(IN) :: flow_depth, time_step, total_solid_fraction
            REAL(wp), INTENT(IN) :: maximum_solid_fraction, hindered_settling_exponent
            REAL(wp) :: rate

            IF (model_name .EQ. 'lahar') THEN
                    rate = eval_lahar_deposition_rate(solid_fraction, settling_velocity,                               &
                                flow_depth, time_step)

            ELSE

                    rate = eval_hindered_settling_deposition_rate(solid_fraction,                                      &
                                settling_velocity, flow_depth, time_step, total_solid_fraction,                        &
                                maximum_solid_fraction, hindered_settling_exponent)

            END IF

    END FUNCTION eval_solid_deposition_rate

        !> @brief Limit lahar settling flux by the mobile solid inventory per time step.
        !> @param[in] solid_fraction Volume fraction of the settling solid species.
        !> @param[in] settling_velocity Species terminal settling speed [m/s].
        !> @param[in] flow_depth Mobile flow thickness [m].
        !> @param[in] time_step Physical time-step duration [s].
        !> @return Physical rate from the indicated constitutive closure.
    PURE FUNCTION eval_lahar_deposition_rate(solid_fraction, settling_velocity,                                        &
                    flow_depth, time_step) RESULT(deposition_rate)

            REAL(wp), INTENT(IN) :: solid_fraction, settling_velocity
            REAL(wp), INTENT(IN) :: flow_depth, time_step
            REAL(wp) :: deposition_rate

            deposition_rate = MIN(solid_fraction * settling_velocity,                                                  &
                        flow_depth * solid_fraction / time_step)

    END FUNCTION eval_lahar_deposition_rate

        !> @brief Apply concentration-dependent settling hindrance and inventory limiting.
        !> @param[in] solid_fraction Volume fraction of the settling solid species.
        !> @param[in] settling_velocity Species terminal settling speed [m/s].
        !> @param[in] flow_depth Mobile flow thickness [m].
        !> @param[in] time_step Physical time-step duration [s].
        !> @param[in] total_solid_fraction Total solid volume fraction.
        !> @param[in] maximum_solid_fraction Packing concentration used by the hindrance closure.
        !> @param[in] hindered_settling_exponent Exponent controlling settling suppression at high
        !> concentration.
        !> @return Physical rate from the indicated constitutive closure.
    PURE FUNCTION eval_hindered_settling_deposition_rate(solid_fraction,                                               &
                    settling_velocity, flow_depth, time_step, total_solid_fraction,                                    &
                    maximum_solid_fraction, hindered_settling_exponent) RESULT(rate)

            REAL(wp), INTENT(IN) :: solid_fraction, settling_velocity
            REAL(wp), INTENT(IN) :: flow_depth, time_step, total_solid_fraction
            REAL(wp), INTENT(IN) :: maximum_solid_fraction, hindered_settling_exponent
            REAL(wp) :: rate

            rate = solid_fraction * settling_velocity
            rate = rate * (1.0_wp - MIN(1.0_wp, total_solid_fraction /                                                 &
                        maximum_solid_fraction))**hindered_settling_exponent
            rate = MIN(rate, flow_depth * solid_fraction / time_step)

    END FUNCTION eval_hindered_settling_deposition_rate

END MODULE mass_exchange_deposition_2d
