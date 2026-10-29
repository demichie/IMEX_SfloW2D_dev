!> \brief Compute solid settling and deposition rates.
!> \details The explicit treatment preserves the established local fluid
!> properties and settling closures; unsupported methods are rejected.
MODULE mass_exchange_deposition_2d
  USE parameters_2d
  USE constitutive_parameters_2d

  USE diagnostics_2d, ONLY : fatal_error
  USE rheology_lahar_2d, ONLY : eval_lahar_fluid_viscosity
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



  SUBROUTINE eval_mass_exchange_deposition_explicit(h,T,alphas,alpha,rho_m,rho_c,cp_c,dt,r)
    REAL(wp), INTENT(IN) :: h,T,alphas(n_solid),alpha,rho_m,rho_c,cp_c,dt
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL eval_deposition_volume_kernel(h,T,alphas,alpha,rho_m,rho_c, &
         dt,r%solid_deposition_rate)
    r%continuous_phase_deposition_rate=coeff_porosity*SUM(r%solid_deposition_rate)
    r%solid_mass_rate=-rho_s*r%solid_deposition_rate
    r%mixture_mass_rate=SUM(r%solid_mass_rate)-rho_c*r%continuous_phase_deposition_rate
    r%energy_rate=T*(SUM(sp_heat_s*r%solid_mass_rate)- &
         rho_c*cp_c*r%continuous_phase_deposition_rate)
    r%dB_dt=SUM(r%solid_deposition_rate)+r%continuous_phase_deposition_rate
  END SUBROUTINE

  SUBROUTINE eval_mass_exchange_deposition_semi_implicit(h,T,alphas,alpha,rho_m,rho_c,cp_c,dt,r)
    REAL(wp), INTENT(IN) :: h,T,alphas(n_solid),alpha,rho_m,rho_c,cp_c,dt
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('deposition: unsupported physical treatment')
  END SUBROUTINE eval_mass_exchange_deposition_semi_implicit

  SUBROUTINE eval_mass_exchange_deposition_implicit(h,T,alphas,alpha,rho_m,rho_c,cp_c,dt,r)
    REAL(wp), INTENT(IN) :: h,T,alphas(n_solid),alpha,rho_m,rho_c,cp_c,dt
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('deposition: unsupported physical treatment')
  END SUBROUTINE eval_mass_exchange_deposition_implicit

  SUBROUTINE solve_mass_exchange_deposition_analytic(h,T,alphas,alpha,rho_m,rho_c,cp_c,dt,r)
    REAL(wp), INTENT(IN) :: h,T,alphas(n_solid),alpha,rho_m,rho_c,cp_c,dt
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('deposition: unsupported physical treatment')
  END SUBROUTINE solve_mass_exchange_deposition_analytic

  SUBROUTINE eval_deposition_volume_kernel(r_h, r_T, r_alphas, alphas_tot, &
       r_rho_m, r_rho_c, dt, deposition_term)
    REAL(wp), INTENT(IN) :: r_h, r_T, r_alphas(n_solid), alphas_tot
    REAL(wp), INTENT(IN) :: r_rho_m, r_rho_c, dt
    REAL(wp), INTENT(OUT) :: deposition_term(n_solid)
    REAL(wp) :: alpha_max, hind_exp, kin_visc_c_local
    REAL(wp) :: dyn_visc_c, rho_c, rhoc, inv_kin_visc, settling_vel
    INTEGER :: i_solid

    alpha_max = 0.6_wp
    hind_exp = 4.65_wp
    kin_visc_c_local = kin_visc_c
    deposition_term = 0.0_wp

    IF ( gas_flag .AND. sutherland_flag ) THEN

       dyn_visc_c = muRef_Suth * ( r_T / Tref_Suth )**1.5_wp *               &
            ( Tref_Suth + S_mu ) / ( r_T + S_mu )

       rho_c = pres / ( sp_gas_const_a * r_T )
       kin_visc_c_local = dyn_visc_c / rho_c

    END IF

    CALL eval_settling_fluid_properties(rheology_model, r_T, alphas_tot,  &
         r_rho_m, r_rho_c, kin_visc_c_local, rhoc, inv_kin_visc)

    DO i_solid=1,n_solid

       IF ( ( r_alphas(i_solid) .GT. 0.0_wp ) .AND. ( settling_flag ) ) THEN

          settling_vel = settling_velocity( diam_s(i_solid) , rho_s(i_solid) ,  &
               rhoc , inv_kin_visc )

          deposition_term(i_solid) = eval_solid_deposition_rate(              &
               rheology_model, r_alphas(i_solid), settling_vel, r_h, dt,  &
               alphas_tot, alpha_max, hind_exp)

          IF ( deposition_term(i_solid) .LT. 0.0_wp ) THEN
             CALL fatal_error('negative solid deposition rate')

          END IF

       END IF

    END DO
  END SUBROUTINE eval_deposition_volume_kernel

  REAL(wp) FUNCTION settling_velocity(diam, rhos, rhoc, inv_kin_visc)

    REAL(wp), INTENT(IN) :: diam
    REAL(wp), INTENT(IN) :: rhos
    REAL(wp), INTENT(IN) :: rhoc
    REAL(wp), INTENT(IN) :: inv_kin_visc

    REAL(wp) :: Rey
    REAL(wp) :: inv_sqrt_C_D
    INTEGER :: i
    REAL(wp) :: const_part
    REAL(wp) :: set_vel_old

    inv_sqrt_C_D = 1.0_wp

    IF (rhos <= rhoc) THEN
       settling_velocity = 0.0_wp
       RETURN
    END IF

    const_part = SQRT((4.0_wp / 3.0_wp) * (rhos / rhoc - 1.0_wp) * diam * grav)

    settling_velocity = const_part * inv_sqrt_C_D
    Rey = diam * settling_velocity * inv_kin_visc

    IF (Rey <= 1000.0_wp) THEN
       C_D_loop: DO i = 1, 20
          set_vel_old = settling_velocity

          IF (collective_settling_flag) THEN
             inv_sqrt_C_D = SQRT(1.0_wp / (24.0_wp * (1.0_wp +               &
                  0.15_wp * Rey**(0.687_wp)) / Rey +                          &
                  A_drag * (LOG(Rey) - B_drag)))
          ELSE
             inv_sqrt_C_D = SQRT(Rey / (24.0_wp * (1.0_wp +                 &
                  0.15_wp * Rey**(0.687_wp))))
          END IF

          settling_velocity = const_part * inv_sqrt_C_D

          IF (ABS(set_vel_old - settling_velocity) / set_vel_old             &
               < 1.0E-6_wp) THEN
             RETURN
          END IF

          Rey = diam * settling_velocity * inv_kin_visc
       END DO C_D_loop
    END IF

    RETURN

  END FUNCTION settling_velocity

  SUBROUTINE eval_settling_fluid_properties(model_name, temperature,           &
       total_solid_fraction, mixture_density, carrier_density_fallback,     &
       kinematic_viscosity_fallback, carrier_density,                       &
       inverse_kinematic_viscosity)

    CHARACTER(*), INTENT(IN) :: model_name
    REAL(wp), INTENT(IN) :: temperature, total_solid_fraction, mixture_density
    REAL(wp), INTENT(IN) :: carrier_density_fallback
    REAL(wp), INTENT(IN) :: kinematic_viscosity_fallback
    REAL(wp), INTENT(OUT) :: carrier_density, inverse_kinematic_viscosity
    REAL(wp) :: fluid_viscosity

    carrier_density = carrier_density_fallback

    IF (model_name .EQ. 'lahar') THEN
       fluid_viscosity = eval_lahar_fluid_viscosity(temperature,            &
            total_solid_fraction)
       carrier_density = mixture_density
       inverse_kinematic_viscosity = mixture_density / fluid_viscosity
    ELSE
       inverse_kinematic_viscosity = 1.0_wp / kinematic_viscosity_fallback
    END IF

  END SUBROUTINE eval_settling_fluid_properties

  PURE FUNCTION eval_solid_deposition_rate(model_name, solid_fraction,         &
       settling_velocity, flow_depth, time_step, total_solid_fraction,       &
       maximum_solid_fraction, hindered_settling_exponent) RESULT(rate)

    CHARACTER(*), INTENT(IN) :: model_name
    REAL(wp), INTENT(IN) :: solid_fraction, settling_velocity
    REAL(wp), INTENT(IN) :: flow_depth, time_step, total_solid_fraction
    REAL(wp), INTENT(IN) :: maximum_solid_fraction, hindered_settling_exponent
    REAL(wp) :: rate

    IF (model_name .EQ. 'lahar') THEN
       rate = eval_lahar_deposition_rate(solid_fraction, settling_velocity,   &
            flow_depth, time_step)
    ELSE
       rate = eval_hindered_settling_deposition_rate(solid_fraction,         &
            settling_velocity, flow_depth, time_step, total_solid_fraction,  &
            maximum_solid_fraction, hindered_settling_exponent)
    END IF

  END FUNCTION eval_solid_deposition_rate

  PURE FUNCTION eval_lahar_deposition_rate(solid_fraction, settling_velocity, &
       flow_depth, time_step) RESULT(deposition_rate)

    REAL(wp), INTENT(IN) :: solid_fraction, settling_velocity
    REAL(wp), INTENT(IN) :: flow_depth, time_step
    REAL(wp) :: deposition_rate

    deposition_rate = MIN(solid_fraction * settling_velocity,                 &
         flow_depth * solid_fraction / time_step)

  END FUNCTION eval_lahar_deposition_rate

  PURE FUNCTION eval_hindered_settling_deposition_rate(solid_fraction,       &
       settling_velocity, flow_depth, time_step, total_solid_fraction,       &
       maximum_solid_fraction, hindered_settling_exponent) RESULT(rate)

    REAL(wp), INTENT(IN) :: solid_fraction, settling_velocity
    REAL(wp), INTENT(IN) :: flow_depth, time_step, total_solid_fraction
    REAL(wp), INTENT(IN) :: maximum_solid_fraction, hindered_settling_exponent
    REAL(wp) :: rate

    rate = solid_fraction * settling_velocity
    rate = rate * (1.0_wp - MIN(1.0_wp, total_solid_fraction /               &
         maximum_solid_fraction))**hindered_settling_exponent
    rate = MIN(rate, flow_depth * solid_fraction / time_step)

  END FUNCTION eval_hindered_settling_deposition_rate
END MODULE mass_exchange_deposition_2d
