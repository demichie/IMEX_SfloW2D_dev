!> Compute solid settling and deposition rates, preserving legacy fluid properties.
MODULE mass_exchange_deposition_2d
  USE parameters_2d
  USE constitutive_parameters_2d
  USE mass_exchange_laws_2d, ONLY : eval_settling_fluid_properties, &
       settling_velocity, eval_solid_deposition_rate
  USE diagnostics_2d, ONLY : fatal_error
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_mass_exchange_deposition
CONTAINS
  SUBROUTINE eval_mass_exchange_deposition(r_h, r_T, r_alphas, alphas_tot, &
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

    CALL eval_settling_fluid_properties(rheology_model_id, r_T, alphas_tot,  &
         r_rho_m, r_rho_c, kin_visc_c_local, rhoc, inv_kin_visc)

    DO i_solid=1,n_solid

       IF ( ( r_alphas(i_solid) .GT. 0.0_wp ) .AND. ( settling_flag ) ) THEN

          settling_vel = settling_velocity( diam_s(i_solid) , rho_s(i_solid) ,  &
               rhoc , inv_kin_visc )

          deposition_term(i_solid) = eval_solid_deposition_rate(              &
               rheology_model_id, r_alphas(i_solid), settling_vel, r_h, dt,  &
               alphas_tot, alpha_max, hind_exp)

          IF ( deposition_term(i_solid) .LT. 0.0_wp ) THEN
             CALL fatal_error('negative solid deposition rate')

          END IF

       END IF

    END DO
  END SUBROUTINE eval_mass_exchange_deposition
END MODULE mass_exchange_deposition_2d
