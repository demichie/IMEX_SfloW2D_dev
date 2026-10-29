!> Erosion process: limit solid and carrier entrainment by available bed inventory.
MODULE mass_exchange_erosion_2d
  USE parameters_2d, ONLY : wp, n_solid, n_vars, slope_correction_flag, &
       erodible_deposit_flag, mass_exchange_rates_type, reset_mass_exchange_rates
  USE diagnostics_2d, ONLY : fatal_error
  USE state_conversion_2d, ONLY : mixt_var
  USE constitutive_parameters_2d, ONLY : erosion_coeff, erodible_porosity, erodible_fract

  USE constitutive_parameters_2d, ONLY: rho_s, rho_c_sub, sp_heat_s, T_erodible
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_mass_exchange_erosion_explicit
  PUBLIC :: eval_mass_exchange_erosion_semi_implicit
  PUBLIC :: eval_mass_exchange_erosion_implicit
  PUBLIC :: solve_mass_exchange_erosion_analytic
  PUBLIC :: eval_erosion_rates
  ! Four-entry physical treatment contract; only EXPLICIT is supported.
CONTAINS



  SUBROUTINE eval_mass_exchange_erosion_explicit(h,u,v,Bx,By,alpha,erodible,dt,T_flow,qp,r)
    REAL(wp), INTENT(IN) :: h,u,v,Bx,By,alpha,erodible(n_solid),dt,T_flow
    REAL(wp), INTENT(IN) :: qp(n_vars+2)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    REAL(wp) :: normal_speed,vertical_velocity,substrate_temperature
    REAL(wp) :: Ri,rhom,rhoc,redgrav,cp_sub,cpm
    CALL reset_mass_exchange_rates(r)
    vertical_velocity=0.0_wp
    IF (slope_correction_flag) vertical_velocity=u*Bx+v*By
    normal_speed=SQRT(u*u+v*v+vertical_velocity*vertical_velocity)
    CALL eval_erosion_volume_kernel(h,normal_speed,alpha,erodible,dt, &
         r%solid_erosion_rate,r%continuous_phase_erosion_rate)
    r%solid_mass_rate=rho_s*r%solid_erosion_rate
    r%mixture_mass_rate=SUM(r%solid_mass_rate)+rho_c_sub*r%continuous_phase_erosion_rate
    CALL mixt_var(qp,Ri,rhom,rhoc,redgrav,cp_sub,cpm)
    substrate_temperature=T_erodible
    IF (erodible_deposit_flag) substrate_temperature=T_flow
    r%energy_rate=substrate_temperature*(SUM(sp_heat_s*r%solid_mass_rate)+ &
         rho_c_sub*cp_sub*r%continuous_phase_erosion_rate)
    r%dB_dt=-SUM(r%solid_erosion_rate)/(1.0_wp-erodible_porosity)
  END SUBROUTINE

  SUBROUTINE eval_mass_exchange_erosion_semi_implicit(h,u,v,Bx,By,alpha,erodible,dt,T_flow,qp,r)
    REAL(wp), INTENT(IN) :: h,u,v,Bx,By,alpha,erodible(n_solid),dt,T_flow
    REAL(wp), INTENT(IN) :: qp(n_vars+2)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('erosion: unsupported physical treatment')
  END SUBROUTINE eval_mass_exchange_erosion_semi_implicit

  SUBROUTINE eval_mass_exchange_erosion_implicit(h,u,v,Bx,By,alpha,erodible,dt,T_flow,qp,r)
    REAL(wp), INTENT(IN) :: h,u,v,Bx,By,alpha,erodible(n_solid),dt,T_flow
    REAL(wp), INTENT(IN) :: qp(n_vars+2)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('erosion: unsupported physical treatment')
  END SUBROUTINE eval_mass_exchange_erosion_implicit

  SUBROUTINE solve_mass_exchange_erosion_analytic(h,u,v,Bx,By,alpha,erodible,dt,T_flow,qp,r)
    REAL(wp), INTENT(IN) :: h,u,v,Bx,By,alpha,erodible(n_solid),dt,T_flow
    REAL(wp), INTENT(IN) :: qp(n_vars+2)
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
    CALL reset_mass_exchange_rates(r)
    CALL fatal_error('erosion: unsupported physical treatment')
  END SUBROUTINE solve_mass_exchange_erosion_analytic

  SUBROUTINE eval_erosion_volume_kernel(r_h, mod_vel, alphas_tot, erodible, dt, &
       erosion_rate, continuous_phase_erosion_rate)
    REAL(wp), INTENT(IN) :: r_h, mod_vel, alphas_tot, erodible(n_solid), dt
    REAL(wp), INTENT(OUT) :: erosion_rate(n_solid), continuous_phase_erosion_rate
    IF (erosion_coeff .GT. 0.0_wp) THEN
       CALL eval_erosion_rates(r_h, mod_vel, alphas_tot, erosion_coeff,       &
            erodible_porosity, erodible_fract(1:n_solid),                    &
            erodible(1:n_solid), dt, erosion_rate(1:n_solid),                &
            continuous_phase_erosion_rate)
    ELSE
       erosion_rate = 0.0_wp
       continuous_phase_erosion_rate = 0.0_wp
    END IF
  END SUBROUTINE eval_erosion_volume_kernel

  PURE SUBROUTINE eval_erosion_rates(flow_depth, speed, total_solid_fraction, &
       erosion_coefficient, erodible_porosity, erodible_fraction,            &
       erodible_inventory, time_step, solid_erosion_rate,                    &
       continuous_phase_erosion_rate)

    REAL(wp), INTENT(IN) :: flow_depth, speed, total_solid_fraction
    REAL(wp), INTENT(IN) :: erosion_coefficient, erodible_porosity
    REAL(wp), INTENT(IN) :: erodible_fraction(:), erodible_inventory(:)
    REAL(wp), INTENT(IN) :: time_step
    REAL(wp), INTENT(OUT) :: solid_erosion_rate(:)
    REAL(wp), INTENT(OUT) :: continuous_phase_erosion_rate

    REAL(wp) :: total_erosion, total_solid_erosion

    IF (erosion_coefficient .GT. 0.0_wp) THEN
       total_erosion = erosion_coefficient * speed * flow_depth *            &
            (1.0_wp - total_solid_fraction)
       total_solid_erosion = total_erosion * (1.0_wp - erodible_porosity)
       solid_erosion_rate = erodible_fraction * total_solid_erosion
    ELSE
       solid_erosion_rate = 0.0_wp
    END IF

    solid_erosion_rate = MAX(0.0_wp, MIN(solid_erosion_rate,                 &
         erodible_inventory / time_step))
    total_solid_erosion = SUM(solid_erosion_rate)
    total_erosion = total_solid_erosion / (1.0_wp - erodible_porosity)
    continuous_phase_erosion_rate = total_erosion * erodible_porosity

  END SUBROUTINE eval_erosion_rates
END MODULE mass_exchange_erosion_2d
