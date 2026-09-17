!********************************************************************************
!> \brief State conversion and directly derived physical quantities
!********************************************************************************
MODULE state_conversion_2d

  USE constitutive_parameters_2d

  USE parameters_2d, ONLY : wp
  USE parameters_2d, ONLY : n_vars, n_solid, n_add_gas, n_quad
  USE parameters_2d, ONLY : energy_flag, liquid_flag, gas_flag, alpha_flag,     &
       stoch_transport_flag, pore_pressure_flag, sutherland_flag

  USE parameters_2d, ONLY : idx_alfas_first, idx_alfas_last, idx_addGas_first,  &
       idx_addGas_last, idx_stoch, idx_pore, idx_u, idx_v

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: r_phys_var, c_phys_var
  PUBLIC :: qc_to_qp, qp_to_qc, qp_to_qp2
  PUBLIC :: mixt_var, eval_sp_heat
  PUBLIC :: avg_profiles_mix
  PUBLIC :: u_log_profile, u_log_profile_scalar, u_log_profile_array
  PUBLIC :: alphas_exp_profile, alphas_exp_profile_scalar
  PUBLIC :: alphas_exp_profile_array
  PUBLIC :: dynamic_pressure, sauter_diameter, average_density_solids
  PUBLIC :: settling_velocity

  INTERFACE u_log_profile
     MODULE PROCEDURE u_log_profile_scalar
     MODULE PROCEDURE u_log_profile_array
  END INTERFACE u_log_profile

  INTERFACE alphas_exp_profile
     MODULE PROCEDURE alphas_exp_profile_scalar
     MODULE PROCEDURE alphas_exp_profile_array
  END INTERFACE alphas_exp_profile

CONTAINS

  FUNCTION u_log_profile_scalar(b,z)

    REAL(wp) :: u_log_profile_scalar
    REAL(wp), INTENT(IN) :: b
    REAL(wp), INTENT(IN) :: z

    u_log_profile_scalar = LOG( b*z + 1.0_wp )

  END FUNCTION u_log_profile_scalar

  FUNCTION u_log_profile_array(b,z)

    REAL(wp) :: u_log_profile_array(n_quad)
    REAL(wp), INTENT(IN) :: b
    REAL(wp), INTENT(IN) :: z(n_quad)

    u_log_profile_array = LOG( b*z + 1.0_wp )

  END FUNCTION u_log_profile_array

  FUNCTION alphas_exp_profile_scalar(a,z)

    REAL(wp) :: alphas_exp_profile_scalar
    REAL(wp), INTENT(IN) :: a
    REAL(wp), INTENT(IN) :: z

    alphas_exp_profile_scalar = EXP( a*z )

  END FUNCTION alphas_exp_profile_scalar

  FUNCTION alphas_exp_profile_array(a,z)

    REAL(wp) :: alphas_exp_profile_array(n_quad)
    REAL(wp), INTENT(IN) :: a
    REAL(wp), INTENT(IN) :: z(n_quad)

    alphas_exp_profile_array = EXP( a*z )

  END FUNCTION alphas_exp_profile_array

  FUNCTION dynamic_pressure(rho_c, alphas, normalizing_coeff_alpha , a,  &
       u, normalizing_coeff_u , b, z)

    !> dynamic pressure
    REAL(wp) :: dynamic_pressure


    !> density of the carrier phase
    REAL(wp), INTENT(IN) :: rho_c
    !> depth-averaged volumetric fractions of solid phases
    REAL(wp), INTENT(IN) :: alphas(n_solid)
    !> normalizing coefficients of the solid profile functions
    REAL(wp), INTENT(IN) :: normalizing_coeff_alpha(n_solid)
    !> shape parameters of the solid profile function
    REAL(wp), INTENT(IN) :: a(n_solid)
    ! depth-averaged velocity
    REAL(wp), INTENT(IN) :: u
    !> normalizing coefficients of the velocity profile function
    REAL(wp), INTENT(IN) :: normalizing_coeff_u
    !> shape parameter of the velocity profile function
    REAL(wp), INTENT(IN) :: b
    !> elevation at which the dynamic pressure is computed
    REAL(wp), INTENT(IN) :: z

    !> volumetric fractions of solid at z
    REAL(wp) :: alphas_z(n_solid)
    !> mixture density at z
    REAL(wp) :: rhom_z
    !> mixture velocity at z
    REAL(wp) :: u_z

    !> volumetric fraction of carrier phase at z
    REAL(wp) :: alphac_z

    !> loop counter
    INTEGER :: i_solid

    DO i_solid = 1,n_solid

       alphas_z(i_solid) = alphas(i_solid) *                                    &
            normalizing_coeff_alpha(i_solid) * alphas_exp_profile(a(i_solid),z)

    END DO

    alphac_z = 1.0_wp - SUM(alphas_z)

    rhom_z = alphac_z * rho_c + SUM( rho_s * alphas_z )

    u_z = ( u * normalizing_coeff_u * u_log_profile(b,z) )

    dynamic_pressure = 0.5_wp * rhom_z * u_z**2

  END FUNCTION dynamic_pressure

  !> Function that calculates the Sauter diameter
  FUNCTION sauter_diameter(alpha_solids)
      !> Sauter diameter
      REAL(wp) :: sauter_diameter

      !> depth-averaged volumetric fractions of solid phases
      REAL(wp), INTENT(IN) :: alpha_solids(n_solid)

      IF ( SUM( alpha_solids ) .LT. EPSILON(1.0_wp) ) THEN
         sauter_diameter = SUM( diam_s ) / n_solid
      ELSE
         sauter_diameter = SUM( alpha_solids ) / SUM( alpha_solids /  diam_s )
      END IF

  END FUNCTION sauter_diameter

   !> Function that calculates the average density of solid phases
  FUNCTION average_density_solids(alpha_solids)
      !> Average density
      REAL(wp) :: average_density_solids

      !> depth-averaged volumetric fractions of solid phases
      REAL(wp), INTENT(IN) :: alpha_solids(n_solid)

      IF ( SUM( alpha_solids ) .LT. EPSILON(1.0_wp) ) THEN
         average_density_solids = SUM( rho_s ) / n_solid
      ELSE
         average_density_solids = SUM( alpha_solids * rho_s ) / SUM( alpha_solids )
      END IF

  END FUNCTION average_density_solids

  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the conservative local variables qj
  !> the local physical variables  (\f$h,u,v,\alpha_s,\rho_m,T,\alpha_l \f$).
  !> \param[in]    r_qj        real conservative variables
  !> \param[out]   r_h         real-value flow thickness
  !> \param[out]   r_u         real-value flow x-velocity
  !> \param[out]   r_v         real-value flow y-velocity
  !> \param[out]   r_alphas    real-value solid volume fractions
  !> \param[out]   r_rho_m     real-value flow density
  !> \param[out]   r_T         real-value flow temperature
  !> \param[out]   r_alphal    real-value liquid volume fraction
  !> \param[out]   r_red_grav  real-value reduced gravity
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !> \date 2019/12/13
  !******************************************************************************

  SUBROUTINE r_phys_var(r_qj , r_h , r_u , r_v , r_alphas , r_rho_m , r_T ,     &
       r_alphal , r_alphag , r_red_grav , p_dyn , r_Zs , r_exc_pore_pres)

    USE geometry_2d, ONLY : lambertw , lambertw0 , lambertwm1
    USE geometry_2d, ONLY : z_quad , w_quad

    USE parameters_2d, ONLY : eps_sing , eps_sing4 , vertical_profiles_flag
    IMPLICIT none

    REAL(wp), INTENT(IN) :: r_qj(n_vars)       !< real-value conservative var
    REAL(wp), INTENT(OUT) :: r_h               !< real-value flow thickness
    REAL(wp), INTENT(OUT) :: r_u               !< real-value x-velocity
    REAL(wp), INTENT(OUT) :: r_v               !< real-value y-velocity
    REAL(wp), INTENT(OUT) :: r_alphas(n_solid) !< real-value solid volume fracts
    REAL(wp), INTENT(OUT) :: r_rho_m           !< real-value mixture density
    REAL(wp), INTENT(OUT) :: r_T               !< real-value temperature
    REAL(wp), INTENT(OUT) :: r_alphal          !< real-value liquid volume fract
    REAL(wp), INTENT(OUT) :: r_alphag(n_add_gas) !< real-value gas volume fracts
    REAL(wp), INTENT(OUT) :: r_red_grav        !< real-value reduced gravity
    REAL(wp), INTENT(OUT) :: p_dyn
    REAL(wp), INTENT(OUT) :: r_Zs              !< real-value stochastic variable
    REAL(wp), INTENT(OUT) :: r_exc_pore_pres   !< real-value pore pressure

    REAL(wp) :: r_inv_rhom
    REAL(wp) :: r_xs(n_solid)     !< real-value solid mass fractions
    REAL(wp) :: r_xg(n_add_gas)     !< real-value additional gas mass fractions
    REAL(wp) :: r_xs_tot

    REAL(wp) :: r_Ri            !< real-value Richardson number
    REAL(wp) :: r_xl            !< real-value liquid mass fraction
    REAL(wp) :: r_xc            !< real-value carrier phase mass fraction
    REAL(wp) :: r_alphac        !< real-value carrier phase volume fraction
    REAL(wp) :: r_sp_heat_c     !< real-value specific heat of carrier phase
    REAL(wp) :: r_sp_heat_mix   !< real-value specific heat of mixture
    REAL(wp) :: r_sp_gas_const_c!< real-value gas constant of carrier phase
    REAL(wp) :: r_rho_c         !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_inv_rho_c
    REAL(wp) :: r_inv_rho_g(n_add_gas)    !< add. gas density reciprocal

    REAL(wp) :: inv_qj1

    REAL(wp) :: rhos_alfas(n_solid)

    REAL(wp) :: rhos_alfas_tot_u
    REAL(wp) :: rhos_alfas_tot_v

    REAL(wp) :: settling_vel(n_solid)

    REAL(wp) :: inv_kin_visc

    !REAL(wp) :: a_crit_rel
    REAL(wp) :: h_rel

    REAL(wp) :: a,b,c,d

    REAL(wp) :: h0_rel
    REAL(wp) :: h0_rel_1
    REAL(wp) :: h0_rel_2

    REAL(wp) :: normalizing_coeff_u

    REAL(wp) :: h0
    REAL(wp) :: u_rel0

    REAL(wp) :: uRho_avg
    REAL(wp) :: uRho_avg_new

    REAL(wp) :: u_avg_guess
    REAL(wp) :: u_avg_new

    REAL(wp) :: rhom_avg

    REAL(wp) :: x0,x1,x2

    INTEGER :: i_aitken
    REAL(wp) :: abs_tol , rel_tol
    REAL(wp) :: denominator
    REAL(wp) :: aitkenX
    REAL(wp) :: lambda

    INTEGER :: i_solid

    REAL(wp) :: z(n_quad)
    REAL(wp) :: w(n_quad)

    REAL(wp) :: u_log_avg

    REAL(wp) :: r_sp_heat_c_by_xc

    REAL(wp) :: dyn_visc_c
    REAL(wp) :: kin_visc_c_local

    ! Optional transported quantities must have a defined value even when
    ! their corresponding model is disabled.
    r_alphal = 0.0_wp
    r_alphag = 0.0_wp
    r_Zs = 0.0_wp
    r_exc_pore_pres = 0.0_wp
    r_xl = 0.0_wp

    ! compute solid mass fractions
    IF ( r_qj(1) .GT. EPSILON(1.0_wp) ) THEN

       inv_qj1 = 1.0_wp / r_qj(1)

       r_xs(1:n_solid) = r_qj(idx_alfas_first:idx_alfas_last) * inv_qj1

       IF ( SUM( r_qj(idx_alfas_first:idx_alfas_last) ) .EQ. r_qj(1) ) THEN

          r_xs(1:n_solid) = r_xs(1:n_solid) / SUM( r_xs(1:n_solid) )

       END IF

       IF ( n_add_gas .GT. 0 ) r_xg(1:n_add_gas) =                              &
            r_qj(idx_addGas_first:idx_addGas_last) * inv_qj1

       IF ( stoch_transport_flag ) r_Zs = r_qj(idx_stoch) * inv_qj1

       IF ( pore_pressure_flag ) r_exc_pore_pres = r_qj(idx_pore) * inv_qj1

    ELSE

       r_h = 0.0_wp
       r_u = 0.0_wp
       r_v = 0.0_wp
       r_alphas = 0.0_wp
       r_rho_m = rho_a_amb
       r_T = T_ambient
       r_alphal = 0.0_wp
       r_alphag = 0.0_wp
       r_red_grav = 0.0_wp
       r_rho_c = rho_a_amb
       p_dyn = 0.0_wp
       r_Zs = 0.0_wp
       r_exc_pore_pres = 0.0_wp

       RETURN

    END IF

    r_xs_tot = SUM(r_xs)

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! compute liquid mass fraction
       r_xl = r_qj(n_vars) * inv_qj1

       ! compute carrier phase (gas) mass fraction
       r_xc =  1.0_wp - r_xs_tot - r_xl

       ! compute specific heat of gas phase (weighted average of specific heat
       ! of gas components, with weights given by mass fractions)
       r_sp_heat_c_by_xc = ( ( r_xc - SUM( r_xg(1:n_add_gas) ) ) * sp_heat_a +  &
            DOT_PRODUCT( r_xg(1:n_add_gas) , sp_heat_g(1:n_add_gas) ) )

       ! specific heat of the mixutre: mass average of sp. heat pf phases
       r_sp_heat_mix = DOT_PRODUCT( r_xs(1:n_solid) , sp_heat_s(1:n_solid) )    &
            + r_xl * sp_heat_l + r_xc * r_sp_heat_c_by_xc

    ELSE

       ! compute carrier phase (gas or liquid) mass fraction
       r_xc = 1.0_wp - r_xs_tot

       IF ( gas_flag ) THEN

          r_sp_heat_c_by_xc = ( ( r_xc - SUM( r_xg(1:n_add_gas) ) ) * sp_heat_a &
               + DOT_PRODUCT( r_xg(1:n_add_gas) , sp_heat_g(1:n_add_gas) ) )

       ELSE

          r_sp_heat_c_by_xc = sp_heat_l * r_xc

       END IF

       ! specific heaf of the mixutre: mass average of sp. heat pf phases
       r_sp_heat_mix = DOT_PRODUCT( r_xs(1:n_solid) , sp_heat_s(1:n_solid) )    &
            + r_sp_heat_c_by_xc

    END IF

    ! compute temperature from energy
    IF ( r_qj(1) .GT. eps_sing ) THEN

       IF ( energy_flag ) THEN

          r_T = ( r_qj(4) - 0.5_wp * ( r_qj(2)**2 + r_qj(3)**2 ) * inv_qj1 ) /  &
               ( r_qj(1) * r_sp_heat_mix )

       ELSE

          r_T = r_qj(4) / ( r_qj(1) * r_sp_heat_mix )

       END IF

       IF ( r_T .LE. 0.0_wp ) r_T = T_ambient

    ELSE

       r_T = T_ambient

    END IF

    IF ( gas_flag ) THEN

       ! carrier phase is gas
       IF ( r_xc .GT. EPSILON(1.0_wp) ) THEN

          r_sp_gas_const_c = ( ( r_xc - SUM( r_xg(1:n_add_gas) ) ) * sp_gas_const_a&
               + DOT_PRODUCT( r_xg(1:n_add_gas) , sp_gas_const_g(1:n_add_gas) ) )  &
               / r_xc

       ELSE

          r_sp_gas_const_c = sp_gas_const_a

       END IF

       r_rho_c =  pres / ( r_sp_gas_const_c * r_T )
       r_inv_rho_c = r_sp_gas_const_c * r_T * inv_pres

       r_inv_rho_g(1:n_add_gas) = sp_gas_const_g(1:n_add_gas) * r_T * inv_pres

    ELSE

       r_rho_c = rho_l
       r_inv_rho_c = inv_rho_l
       sp_heat_c = sp_heat_l

    END IF

    ! the liquid contribution (if present) is added below
    r_inv_rhom = DOT_PRODUCT( r_xs(1:n_solid) , inv_rho_s(1:n_solid) )          &
         + r_xc * r_inv_rho_c

    IF ( gas_flag .AND. liquid_flag ) THEN

       r_inv_rhom = r_inv_rhom + r_xl * inv_rho_l

    END IF

    ! mixture density
    r_rho_m = 1.0_wp / r_inv_rhom

    IF ( gas_flag .AND. liquid_flag ) THEN

       r_alphal = r_xl * r_rho_m * inv_rho_l

    END IF

    ! convert from mass fraction to volume fraction
    r_alphas(1:n_solid) = r_xs(1:n_solid) * r_rho_m * inv_rho_s(1:n_solid)

    ! convert from mass fraction to volume fraction
    r_alphag(1:n_add_gas) = r_xg(1:n_add_gas) * r_rho_m                         &
         * r_inv_rho_g(1:n_add_gas)

    ! convert from mass fraction to volume fraction
    r_alphac = r_xc * r_rho_m * r_inv_rho_c

    r_h = r_qj(1) * r_inv_rhom

    ! reduced gravity
    r_red_grav = ( r_rho_m - rho_a_amb ) * r_inv_rhom * grav

    kin_visc_c_local = kin_visc_c

    IF ( vertical_profiles_flag ) THEN

       rhos_alfas_tot_u = r_xs_tot * r_qj(2) / r_h
       rhos_alfas_tot_v = r_xs_tot * r_qj(3) / r_h

       rhos_alfas(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid)

       IF ( gas_flag .AND. sutherland_flag ) THEN

          dyn_visc_c = muRef_Suth * ( r_T / Tref_Suth )**1.5_wp *               &
               ( Tref_Suth + S_mu ) / ( r_T + S_mu )

          kin_visc_c_local = dyn_visc_c * r_inv_rho_c

       END IF

       ! Viscosity read from input file [m2 s-1]
       inv_kin_visc = 1.0_wp / kin_visc_c_local

       DO i_solid=1,n_solid

          settling_vel(i_solid) = settling_velocity( diam_s(i_solid) ,          &
               rho_s(i_solid) , r_rho_c , inv_kin_visc )

       END DO


       ! The profile parameters depend on h/k_s, not on the absolute value of h
       h_rel = r_h / k_s

       IF ( h_rel .GT. H_crit_rel ) THEN

          ! we search for h0_rel such that the average integral between 0 and
          ! h_rel is equal to 1
          ! For h_rel > H_crit_rel this integral is the sum of two pieces:
          ! integral between 0 and h0_rel of the log profile
          ! integral between h0_rel and h_rel of the costant profile

          a = h_rel * vonK / SQRT(friction_factor)
          b = 1.0_wp / 30.0_wp + h_rel
          c = 30.0_wp

          ! solve b*log(c*z+1)-z=a for z
          d =  a / b - 1.0_wp / ( b*c )

          h0_rel_1 = -b*lambertw0( -EXP(d)/(b*c) ) - 1.0_wp / c
          h0_rel_2 = -b*lambertwm1( -EXP(d)/(b*c) ) - 1.0_wp / c
          h0_rel = MIN( h0_rel_1 , h0_rel_2)

       ELSE

          ! when h_rel <= H_crit_rel we have only the log profile and we have to
          ! rescale it in order to have the integral between o and h_rel equal to
          ! 1
          h0_rel = h_rel

       END IF

       h0 = h0_rel * k_s

       b = 30.0_wp / k_s

       ! Quadrature points and weights for the interval [0;h0]
       z = 0.5_wp * h0 * ( z_quad + 1.0_wp )
       w = 0.5_wp * h0 * w_quad

       u_log_avg = ( SUM( w * u_log_profile(b,z) ) + u_log_profile(b,h0)*(r_h-h0) )  &
            / r_h

       !u_log_avg = ( SUM( w * LOG( b*z + 1.0_wp ) ) + LOG( b*h0 + 1.0_wp )*(r_h-h0) )  &
       !     / r_h

       normalizing_coeff_u = 1.0_wp / u_log_avg

       ! velocity at h0
       u_rel0 = normalizing_coeff_u * u_log_profile(b,h0)
       ! u_rel0 = normalizing_coeff_u * LOG( b*h0 + 1.0_wp )

       uRho_avg = SQRT( r_qj(2)**2 + r_qj(3)**2 ) / r_h

       u_avg_guess = uRho_avg / r_rho_m
       x0 = u_avg_guess

       ! loop to compute the average velocity from average rho*alpha and average
       ! uRho ( = 1/h*int( u*rhog*alphag + sum[u*rhos(i)*alphas(i)] ) )

       rel_tol = 1.e-8
       abs_tol = 1.e-8

       aitken_loop:DO i_aitken=1,10

          x0 = u_avg_guess

          CALL avg_profiles_mix( r_h , settling_vel , rhos_alfas(1:n_solid) ,   &
               u_avg_guess , h0 , b , u_rel0 , r_rho_c , rhom_avg ,             &
               uRho_avg_new , p_dyn )

          u_avg_new = u_avg_guess * uRho_avg / ( uRho_avg_new)

          x1 = u_avg_new

          CALL avg_profiles_mix( r_h , settling_vel , rhos_alfas(1:n_solid) ,   &
               u_avg_new , h0 , b , u_rel0 , r_rho_c , rhom_avg ,               &
               uRho_avg_new , p_dyn )

          u_avg_new = u_avg_new * uRho_avg / ( uRho_avg_new)

          x2 = u_avg_new

          IF (x1 .NE.  x0) THEN

             lambda = ABS((x2 - x1)/(x1 - x0))

          END IF

          denominator = (x2 - x1) - (x1 - x0)

          IF ( ABS(denominator) .LT. 0.1*abs_tol ) EXIT aitken_loop

          aitkenX = x2 - ( (x2 - x1)**2 ) / denominator

          u_avg_new = aitkenX

          IF ( ( ABS(u_avg_guess-u_avg_new)/u_avg_guess < rel_tol ) .OR.        &
               ( ABS(u_avg_guess-u_avg_new) < abs_tol ) ) THEN

             EXIT aitken_loop

          END IF

          u_avg_guess = u_avg_new


       END DO aitken_loop

       r_u = u_avg_new * r_qj(2) / ( SQRT( r_qj(2)**2 + r_qj(3)**2 ) )
       r_v = u_avg_new * r_qj(3) / ( SQRT( r_qj(2)**2 + r_qj(3)**2 ) )

    ELSE

       ! velocity components
       IF ( r_qj(1) .GT. eps_sing ) THEN

          r_u = r_qj(2) * inv_qj1
          r_v = r_qj(3) * inv_qj1

       ELSE

          r_u = SQRT(2.0_wp) * r_qj(1) * r_qj(2) / SQRT( r_qj(1)**4 + eps_sing4 )
          r_v = SQRT(2.0_wp) * r_qj(1) * r_qj(3) / SQRT( r_qj(1)**4 + eps_sing4 )

       END IF

    END IF

    p_dyn = 0.5 * r_rho_m * ( r_u**2 + r_v**2 )

    ! Richardson number
    IF ( ( r_u**2 + r_v**2 ) .GT. 0.0_wp ) THEN

       r_Ri = r_red_grav * r_h / ( r_u**2 + r_v**2 )

    ELSE

       r_Ri = 0.0_wp

    END IF

    RETURN

  END SUBROUTINE r_phys_var


  SUBROUTINE avg_profiles_mix( h , settling_vel , rho_alphas_avg , u_guess ,    &
       h0 , b , u_rel0 , rho_c , rhom_avg , uRho_avg_new , p_dyn )

    USE geometry_2d, ONLY : z_quad , w_quad

    USE geometry_2d, ONLY : calcei , gaulegf

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: h
    REAL(wp), INTENT(IN) :: settling_vel(n_solid)
    REAL(wp), INTENT(IN) :: rho_alphas_avg(n_solid)

    REAL(wp), INTENT(IN) :: u_guess
    REAL(wp), INTENT(IN) :: h0
    REAL(wp), INTENT(IN) :: b
    REAL(wp), INTENT(IN) :: u_rel0
    REAL(wp), INTENT(IN) :: rho_c
    REAL(wp), INTENT(OUT) :: rhom_avg
    REAL(wp), INTENT(OUT) :: uRho_avg_new
    REAL(wp), INTENT(OUT) :: p_dyn

    !> Shear velocity computed from u_guess
    REAL(wp) :: shear_vel

    !>  Rouse numbers for the particle classes
    REAL(wp) :: Rouse_no(n_solid)

    !> array for depth-averaged value of rho*u(z)*C(z)
    REAL(wp) :: rho_u_alphas(n_solid)


    REAL(wp) :: rho_alphas(n_solid)
    REAL(wp) :: alphas(n_solid)

    INTEGER :: i_solid

    REAL(wp) :: normalizing_coeff_u

    REAL(wp) :: a(n_solid)

    REAL(wp) :: alphas_exp_avg
    REAL(wp) :: u_log_avg

    REAL(wp) :: normalizing_coeff_alpha(n_solid)
    REAL(wp) :: alphas_rel0
    REAL(wp) :: y

    REAL(wp) :: int_def1, int_def2

    REAL(wp) :: epsilon_s
    REAL(wp) :: a_coeff

    INTEGER ( kind = 4 ) :: i
    ! REAL(wp) :: x,ei

    REAL(wp) :: z(n_quad)
    REAL(wp) :: w(n_quad)
    REAL(wp) :: int_quad

    REAL(wp) :: rhom_z
    REAL(wp) :: u_z

    REAL(wp) :: z_test

    ! Shear velocity computed from u_guess
    shear_vel = u_guess * SQRT(friction_factor)

    ! Rouse numbers for the particle classes
    DO i_solid=1,n_solid

       IF ( shear_vel .GT. 0.0_wp ) THEN

          Rouse_no(i_solid) = settling_vel(i_solid) / ( vonK * shear_vel )

       ELSE

          Rouse_no(i_solid) = 0.0_wp

       END IF

    END DO

    ! Quadrature points and weights for the interval [0;h0]
    z = 0.5_wp * h0 * ( z_quad + 1.0_wp )
    w = 0.5_wp * h0 * w_quad

    epsilon_s = Sc * shear_vel * vonK * ( ( h0/6.0_wp ) + ( k_s / 60.0_wp ) )
    a_coeff = - vonK * shear_vel / epsilon_s

    u_log_avg = ( SUM( w * u_log_profile(b,z) ) + u_log_profile(b,h0)*(h-h0) )  &
         / h

    normalizing_coeff_u = 1.0_wp / u_log_avg

    DO i_solid=1,n_solid

       a(i_solid) = a_coeff * Rouse_no(i_solid)

       ! depth-average value of exp(a*x)
       alphas_exp_avg = ( SUM( w * alphas_exp_profile(a(i_solid),z) ) +         &
            alphas_exp_profile(a(i_solid),h0)*(h-h0) ) / h

       normalizing_coeff_alpha(i_solid) = 1.0_wp / alphas_exp_avg

       int_quad = SUM( w * ( alphas_exp_profile(a(i_solid),z) *                 &
            u_log_profile(b,z) ) )

       ! integral of alfa_rel_i*u between in the boundary layer
       int_def1 = normalizing_coeff_u * normalizing_coeff_alpha(i_solid) *      &
            int_quad

       ! relative concentration alphas_rel at depth h0 (from the bottom)
       ! alphas_rel is defined as alphas(z)/alphas_avg
       alphas_rel0 = normalizing_coeff_alpha(i_solid) *                         &
            alphas_exp_profile(a(i_solid),h0)

       ! integral of alfa_rel_i*u in the free-stream layer
       int_def2 =  ( h - h0 ) * u_rel0 * alphas_rel0

       ! we add the contribution of the integral of the constant region, we
       ! average by dividing by h and we multiply by the density of solid and
       ! average concentration and by u_guess.
       rho_u_alphas(i_solid) = rho_alphas_avg(i_solid) * u_guess *              &
            ( int_def1 + int_def2 ) / h

    END DO

    ! we add the contribution of the gas phase to the depth-averaged mixture
    ! density. This value should be equal to that used to compute the input
    ! values rhoalphas_avg.
    rhom_avg = rho_c + SUM( ( rho_s - rho_c ) / rho_s * rho_alphas_avg )

    ! we add the contribution of the gas phase to the mixture depth-averaged
    ! momentum
    uRho_avg_new = ( u_guess*rho_c + SUM((rho_s-rho_c) / rho_s * rho_u_alphas) )


    IF ( z_dyn .GT. 0.0_wp ) THEN

       z_test = MIN(z_dyn,h0)

       alphas = rho_alphas_avg / rho_s

       p_dyn = dynamic_pressure(rho_c, alphas, normalizing_coeff_alpha , a, &
            u_guess, normalizing_coeff_u , b, z_test)

    ELSE

       p_dyn = 0.5_wp * rhom_avg * uRho_avg_new

    END IF

  END SUBROUTINE avg_profiles_mix

  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the conservative local variables qj
  !> the local physical variables  (\f$h,u,v,T,\rho_m,red grav,\alpha_s \f$).
  !> \param[in]    c_qj      complex conservative variables
  !> \param[out]   h         complex-value flow thickness
  !> \param[out]   u         complex-value flow x-velocity
  !> \param[out]   v         complex-value flow y-velocity
  !> \param[out]   T         complex-value flow temperature
  !> \param[out]   rho_m     complex-value flow density
  !> \param[out]   alphas    complex-value solid volume fractions
  !> \param[out]   inv_rhom  complex-value mixture density reciprocal
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !> \date 2019/12/13
  !******************************************************************************

  SUBROUTINE c_phys_var( c_qj , h , u , v , T , rho_m , alphas , alphag ,       &
       inv_rhom , Zs , exc_pore_pres )

    USE COMPLEXIFY
    USE parameters_2d, ONLY : eps_sing , eps_sing4
    IMPLICIT none

    COMPLEX(wp), INTENT(IN) :: c_qj(n_vars)
    COMPLEX(wp), INTENT(OUT) :: h               !< height [m]
    COMPLEX(wp), INTENT(OUT) :: u               !< velocity (x direction) [m s-1]
    COMPLEX(wp), INTENT(OUT) :: v               !< velocity (y direction) [m s-1]
    COMPLEX(wp), INTENT(OUT) :: T               !< temperature [K]
    COMPLEX(wp), INTENT(OUT) :: rho_m           !< mixture density [kg m-3]
    COMPLEX(wp), INTENT(OUT) :: alphas(n_solid) !< sediment volume fractions
    COMPLEX(wp), INTENT(OUT) :: alphag(n_add_gas) !< additional-gas volume fractions
    COMPLEX(wp), INTENT(OUT) :: inv_rhom        !< 1/mixture density [kg-1 m3]
    COMPLEX(wp), INTENT(OUT) :: Zs              !< stochastic variable
    COMPLEX(wp), INTENT(OUT) :: exc_pore_pres   !< excess pore pressure

    COMPLEX(wp) :: xs(n_solid)             !< sediment mass fractions
    COMPLEX(wp) :: xg(n_add_gas)           !< additional gas comp. mass fractions
    COMPLEX(wp) :: xs_tot                  !< sum of solid mass fraction
    COMPLEX(wp) :: xl                      !< liquid mass fraction
    COMPLEX(wp) :: xc                      !< carrier phase mass fraction
    COMPLEX(wp) :: sp_heat_c               !< Specific heat of carrier phase
    COMPLEX(wp) :: sp_heat_mix             !< Specific heat of mixture
    COMPLEX(wp) :: sp_gas_const_c          !< Gas constant of carrier phase
    COMPLEX(wp) :: inv_cqj1                !< reciprocal of 1st cons. variable
    COMPLEX(wp) :: inv_rho_c               !< carrier phase density reciprocal
    COMPLEX(wp) :: inv_rho_g(n_add_gas)    !< add. gas density reciprocal

    Zs = CMPLX(0.0_wp,0.0_wp,wp)
    exc_pore_pres = CMPLX(0.0_wp,0.0_wp,wp)

    ! compute solid mass fractions
    IF ( REAL(c_qj(1)) .GT.  EPSILON(1.0_wp) ) THEN

       inv_cqj1 = 1.0_wp / c_qj(1)
       xs(1:n_solid) = c_qj(idx_alfas_first:idx_alfas_last) * inv_cqj1

       xg(1:n_add_gas) = c_qj(idx_addGas_first:idx_addGas_last) * inv_cqj1

       IF ( stoch_transport_flag ) Zs = c_qj(idx_stoch) * inv_cqj1

       IF ( pore_pressure_flag) exc_pore_pres = c_qj(idx_pore) * inv_cqj1

    ELSE

       h = CMPLX(0.0_wp,0.0_wp,wp)
       u = CMPLX(0.0_wp,0.0_wp,wp)
       v = CMPLX(0.0_wp,0.0_wp,wp)
       T = CMPLX(T_ambient,0.0_wp,wp)
       rho_m = CMPLX(rho_a_amb,0.0_wp,wp)
       alphas = CMPLX(0.0_wp,0.0_wp,wp)
       alphag = CMPLX(0.0_wp,0.0_wp,wp)
       inv_rhom = 1.0_wp / rho_m
       Zs = 0.0_wp
       exc_pore_pres = 0.0_wp

       RETURN

    END IF

    xs_tot = SUM(xs)

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! compute liquid mass fraction
       xl = c_qj(n_vars) * inv_cqj1

       ! compute carrier phase (gas) mass fraction
       xc = 1.0_wp - xs_tot - xl

       sp_heat_c = ( ( xc - SUM( xg(1:n_add_gas) ) ) * sp_heat_a +              &
            DOT_PRODUCT( xg(1:n_add_gas) , sp_heat_g(1:n_add_gas) ) ) / xc

       ! specific heat of the mixutre: mass average of sp. heat pf phases
       sp_heat_mix = DOT_PRODUCT( xs(1:n_solid) , sp_heat_s(1:n_solid) )        &
            + xl * sp_heat_l + xc * sp_heat_c

    ELSE

       ! compute carrier phase (gas or liquid) mass fraction
       xc = 1.0_wp - xs_tot

       IF ( gas_flag ) THEN

          sp_heat_c = ( ( xc - SUM( xg(1:n_add_gas) ) ) * sp_heat_a +           &
               DOT_PRODUCT( xg(1:n_add_gas) , sp_heat_g(1:n_add_gas) ) ) / xc

       ELSE

          sp_heat_c = CMPLX(sp_heat_l,0.0_wp,wp)

       END IF

       ! specific heaf of the mixutre: mass average of sp. heat pf phases
       sp_heat_mix = DOT_PRODUCT( xs(1:n_solid) , sp_heat_s(1:n_solid) )        &
            + xc * sp_heat_c

    END IF

    ! compute temperature from energy
    IF ( REAL(c_qj(1)) .GT. eps_sing ) THEN

       IF ( energy_flag ) THEN

          T = ( c_qj(4) - 0.5_wp * ( c_qj(2)**2 + c_qj(3)**2 ) * inv_cqj1 ) /   &
               ( c_qj(1) * sp_heat_mix )

       ELSE

          T = c_qj(4) / ( c_qj(1) * sp_heat_mix )

       END IF

       IF ( REAL(T) .LE. 0.0_wp ) T = CMPLX(T_ambient,0.0_wp,wp)

    ELSE

       T = CMPLX(T_ambient,0.0_wp,wp)

    END IF

    IF ( gas_flag ) THEN

       ! carrier phase is gas
       sp_gas_const_c = ( ( xc - SUM( xg(1:n_add_gas) ) ) * sp_gas_const_a      &
            + DOT_PRODUCT( xg(1:n_add_gas) , sp_gas_const_g(1:n_add_gas) ) )    &
            / xc

       inv_rho_c = sp_gas_const_c * T * inv_pres

       inv_rho_g(1:n_add_gas) = sp_gas_const_g(1:n_add_gas) * T * inv_pres

    ELSE

       inv_rho_c = CMPLX(inv_rho_l,0.0_wp,wp)

    END IF

    inv_rhom = DOT_PRODUCT( xs(1:n_solid) , c_inv_rho_s(1:n_solid) )            &
         + xc * inv_rho_c

    IF ( gas_flag .AND. liquid_flag ) inv_rhom = inv_rhom + xl * inv_rho_l

    rho_m = 1.0_wp / inv_rhom

    ! convert from mass fraction to volume fraction
    alphas(1:n_solid) = rho_m * xs(1:n_solid) * c_inv_rho_s(1:n_solid)

    ! convert from mass fraction to volume fraction
    alphag(1:n_add_gas) = rho_m * xg(1:n_add_gas) * inv_rho_g(1:n_add_gas)

    h = c_qj(1) * inv_rhom

    ! velocity components
    IF ( REAL( c_qj(1) ) .GT. eps_sing ) THEN

       u = c_qj(2) * inv_cqj1
       v = c_qj(3) * inv_cqj1

    ELSE

       u = SQRT(2.0_wp) * c_qj(1) * c_qj(2) / SQRT( c_qj(1)**4 + eps_sing4 )
       v = SQRT(2.0_wp) * c_qj(1) * c_qj(3) / SQRT( c_qj(1)**4 + eps_sing4 )

    END IF

    RETURN

  END SUBROUTINE c_phys_var


  !******************************************************************************
  !> \brief Mixture variables
  !
  !> This subroutine evaluates from the physical real-value local variables qpj,
  !> some mixture variable.
  !> \param[in]    qpj          real-valued physical variables
  !> \param[out]   r_Ri         real-valued Richardson number
  !> \param[out]   r_rho_m      real-valued mixture density
  !> \param[out]   r_rho_c      real-valued carrier phase density
  !> \param[out]   r_red_grav   real-valued reduced gravity
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !> \date 10/10/2019
  !******************************************************************************

  SUBROUTINE mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,sp_heat_flag,         &
       r_sp_heat_c,r_sp_heat_mix)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2) !< real-value physical variables
    REAL(wp), INTENT(OUT) :: r_Ri         !< real-value Richardson number
    REAL(wp), INTENT(OUT) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(wp), INTENT(OUT) :: r_rho_c !< real-value carrier phase density [kg/m3]
    REAL(wp), INTENT(OUT) :: r_red_grav   !< real-value reduced gravity
    LOGICAL, INTENT(IN) :: sp_heat_flag
    REAL(wp), INTENT(OUT) :: r_sp_heat_c
    REAL(wp), INTENT(OUT) :: r_sp_heat_mix


    REAL(wp) :: r_u                       !< real-value x-velocity
    REAL(wp) :: r_v                       !< real-value y-velocity
    REAL(wp) :: r_h                       !< real-value flow thickness
    REAL(wp) :: r_alphas(n_solid)         !< real-value solid volume fractions
    REAL(wp) :: r_alphag(n_add_gas)       !< real-value add.gas volume fractions
    REAL(wp) :: r_rho_a                   !< real-value atm.gas density
    REAL(wp) :: r_rho_g(n_add_gas)        !< real-value add.gas densities
    REAL(wp) :: r_T                       !< real-value temperature [K]
    REAL(wp) :: r_alphal                  !< real-value liquid volume fraction
    REAL(wp) :: r_alphac

    REAL(wp) :: alphas_tot                !< total solid fraction

    REAL(wp) :: r_inv_rhom

    r_h = qpj(1)

    IF ( qpj(1) .LE. EPSILON(1.0_wp) ) THEN

       r_red_grav = 0.0_wp
       r_rho_m = rho_a_amb
       r_rho_c = rho_a_amb

       IF ( liquid_flag ) THEN

          r_sp_heat_c = sp_heat_l
          r_sp_heat_mix = sp_heat_l

       ELSE

          r_sp_heat_c = sp_heat_a
          r_sp_heat_mix = sp_heat_a


       END IF

       r_Ri = 0.0_wp

       RETURN

    END IF

    r_u = qpj(idx_u)
    r_v = qpj(idx_v)
    r_T = qpj(4)

    IF ( alpha_flag ) THEN

       r_alphas(1:n_solid) = qpj(idx_alfas_first:idx_alfas_last)
       r_alphag(1:n_add_gas) = qpj(idx_addGas_first:idx_addGas_last)

    ELSE

       r_alphas(1:n_solid) = qpj(idx_alfas_first:idx_alfas_last) / qpj(1)
       IF ( n_add_gas .GT. 0 ) r_alphag(1:n_add_gas) =                       &
            qpj(idx_addGas_first:idx_addGas_last) / qpj(1)

    END IF

    alphas_tot = SUM(r_alphas)

    r_alphal = 0.0_wp

    IF ( gas_flag .AND. liquid_flag ) THEN

       IF ( alpha_flag ) THEN

          r_alphal = qpj(n_vars)

       ELSE

          r_alphal = qpj(n_vars) / qpj(1)

       END IF

    END IF

    ! carrier phase volume fraction
    r_alphac = 1.0_wp - alphas_tot - r_alphal

    IF ( gas_flag ) THEN

       ! continuous phase is gas
       r_rho_a =  pres / ( sp_gas_const_a * r_T )
       r_rho_g(1:n_add_gas) = pres / ( sp_gas_const_g(1:n_add_gas) * r_T )

       r_rho_c = ( ( 1.0_wp - r_alphal - alphas_tot - SUM(r_alphag) ) * r_rho_a &
            + DOT_PRODUCT( r_alphag(1:n_add_gas) , r_rho_g(1:n_add_gas) ) )     &
            / ( 1.0_wp - r_alphal - alphas_tot )

    ELSE

       ! continuous phase is liquid
       r_rho_c = rho_l

    END IF

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! density of mixture of carrier (gas), liquid and solids
       r_rho_m = ( 1.0_wp - alphas_tot - r_alphal ) * r_rho_c                   &
            + DOT_PRODUCT( r_alphas , rho_s ) + r_alphal * rho_l

    ELSE

       ! density of mixture of carrier phase and solids
       r_rho_m = ( 1.0_wp - alphas_tot ) * r_rho_c + DOT_PRODUCT( r_alphas ,    &
            rho_s )

    END IF

    r_inv_rhom = 1.0_wp / r_rho_m

    ! reduced gravity
    r_red_grav = ( r_rho_m - rho_a_amb ) / r_rho_m * grav

    ! Richardson number
    IF ( ( r_u**2 + r_v**2 ) .GT. 0.0_wp ) THEN

       r_Ri = MIN(1.E15_wp,r_red_grav * r_h / ( r_u**2 + r_v**2 ))

    ELSE

       r_Ri = 0.0_wp

    END IF

    IF ( sp_heat_flag ) THEN

       CALL eval_sp_heat( r_alphal , r_alphas , r_alphag, r_rho_g , r_inv_rhom ,&
            r_sp_heat_c , r_sp_heat_mix )

    END IF

    RETURN

  END SUBROUTINE mixt_var

  !******************************************************************************
  !> \brief Specific heat
  !
  !> This subroutine evaluates the specific heat of the carrier phase and of the
  !> mixture.
  !> \param[in]    r_alphal        real-value liquid volume fraction
  !> \param[in]    r_alphas        real-value solid volume fraction
  !> \param[in]    r_alphag        real-value add.gas volume fraction
  !> \param[in]    rho_g           real-value gas density
  !> \param[in]    r_inv_rhom      real-value mixture density reciprocal
  !> \param[out]   r_sp_heat_c     real-valued carrier phase specific heat
  !> \param[out]   r_sp_heat_mix   real-valued mixture specific heat
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !> \date 2021/07/09
  !******************************************************************************

  SUBROUTINE eval_sp_heat( r_alphal , r_alphas , r_alphag , rho_g ,r_inv_rhom , &
       r_sp_heat_c , r_sp_heat_mix )

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: r_alphal          !< real-value liquid volume fraction
    REAL(wp), INTENT(IN) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp), INTENT(IN) :: r_alphag(n_add_gas)!< real-value add.gas volume fractions
    ! REAL(wp), INTENT(IN) :: r_alphac    !< real-value carrier phase volume fraction

    REAL(wp), INTENT(IN) :: rho_g(n_add_gas)!< real-value add.gas densities [kg/m3]
    ! REAL(wp), INTENT(IN) :: r_rho_c   !< real-value carrier phase density [kg/m3]
    REAL(wp), INTENT(IN) :: r_inv_rhom        !< real-value mixture density [kg/m3]

    REAL(wp), INTENT(OUT) :: r_sp_heat_c
    REAL(wp), INTENT(OUT) :: r_sp_heat_mix

    REAL(wp) :: r_xl              !< real-value liquid mass fraction
    REAL(wp) :: r_xc              !< real-value carrier phase mass fraction

    REAL(wp) :: r_xs(n_solid)     !< real-value solid mass fractions
    REAL(wp) :: r_xg(n_add_gas)   !< real-value add.gas mass fractions

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! liquid mass fraction
       r_xl = r_alphal * rho_l * r_inv_rhom

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid) * r_inv_rhom

       ! additional gas mass fractions
       r_xg(1:n_add_gas) = r_alphag(1:n_add_gas) * rho_g(1:n_add_gas) * r_inv_rhom

       ! carrier (gas) mass fraction
       r_xc = 1.0_wp - ( r_xl + SUM(r_xs(1:n_solid) ) )

       ! specific heat of gas (mass. avg. of sp.heat of gas components)

       IF ( r_xc .GT. EPSILON(1.0_wp) ) THEN

          r_sp_heat_c = ( ( r_xc - SUM( r_xg(1:n_add_gas) ) ) * sp_heat_a +        &
               DOT_PRODUCT( r_xg(1:n_add_gas) , sp_heat_g(1:n_add_gas) ) ) / r_xc

       ELSE

          r_sp_heat_c = sp_heat_a

       END IF

       ! mass averaged mixture specific heat
       r_sp_heat_mix =  DOT_PRODUCT( r_xs , sp_heat_s ) + r_xl * sp_heat_l      &
            + r_xc * r_sp_heat_c

    ELSE

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid) * r_inv_rhom

       ! additional gas mass fractions
       r_xg(1:n_add_gas) = r_alphag(1:n_add_gas) * rho_g(1:n_add_gas)           &
            * r_inv_rhom

       ! carrier (gas or liquid) mass fraction
       r_xc = 1.0_wp - SUM( r_xs(1:n_solid) )

       r_sp_heat_c = 0.0_wp

       IF ( gas_flag ) THEN

          IF ( r_xc .GT. EPSILON(1.0_wp) ) THEN

             r_sp_heat_c = ( ( r_xc - SUM( r_xg(1:n_add_gas) ) ) * sp_heat_a +     &
                  DOT_PRODUCT( r_xg(1:n_add_gas) , sp_heat_g(1:n_add_gas) ) ) / r_xc

          ELSE

             r_sp_heat_c = sp_heat_a

          END IF

       ELSE

          r_sp_heat_c = sp_heat_l

       END IF

       ! mass averaged mixture specific heat
       r_sp_heat_mix =  DOT_PRODUCT( r_xs , sp_heat_s ) + r_xc * r_sp_heat_c

    END IF

    RETURN

  END SUBROUTINE eval_sp_heat

  !******************************************************************************
  !> \brief Conservative to physical variables
  !
  !> This subroutine evaluates from the conservative variables qc the
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h \f$
  !> - qp(2) = \f$ hu \f$
  !> - qp(3) = \f$ hv \f$
  !> - qp(4) = \f$ T \f$
  !> - qp(idx_alfas_first:idx_alfas_last) = \f$ alphas(1:n_solid) \f$
  !> - qp(idx_addGas_first:idx_addGas_last) = \f$ alphas(1:n_add_gas) \f$
  !> - qp(n_vars) = \f$ alphal \f$
  !> - qp(idx_u) = \f$ u \f$
  !> - qp(idx_v) = \f$ v \f$
  !> .
  !> The physical variables are those used for the linear reconstruction at the
  !> cell interfaces. Limiters are applied to the reconstructed slopes.
  !> \param[in]     qc     local conservative variables
  !> \param[out]    qp     local physical variables
  !> \param[out]    p_dyn  local dynamic pressure
  !
  !> \date 2019/11/11
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE qc_to_qp(qc,qp,p_dyn)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qc(n_vars)
    REAL(wp), INTENT(OUT) :: qp(n_vars+2)
    REAL(wp), INTENT(OUT) :: p_dyn

    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(wp) :: r_T               !< real-value temperature [K]
    REAL(wp) :: r_alphal          !< real-value liquid volume fraction
    REAL(wp) :: r_alphag(n_add_gas) !< real-value add. gas volume fractions
    REAL(wp) :: r_red_grav
    REAL(wp) :: r_Zs             !< real-value stochastic variable
    REAL(wp) :: r_exc_pore_pres  !< real-value pore pressure

    CALL r_phys_var( qc , r_h , r_u , r_v , r_alphas , r_rho_m , r_T ,          &
         r_alphal , r_alphag , r_red_grav , p_dyn , r_Zs , r_exc_pore_pres )

    qp(1) = r_h

    qp(2) = r_h*r_u
    qp(3) = r_h*r_v

    qp(4) = r_T

    IF ( alpha_flag ) THEN

       qp(idx_alfas_first:idx_alfas_last) = r_alphas(1:n_solid)
       qp(idx_addGas_first:idx_addGas_last) = r_alphag(1:n_add_gas)
       IF ( gas_flag .AND. liquid_flag ) qp(n_vars) = r_alphal

    ELSE

       qp(idx_alfas_first:idx_alfas_last) = r_alphas(1:n_solid) * r_h
       qp(idx_addGas_first:idx_addGas_last) = r_alphag(1:n_add_gas) * r_h
       IF ( gas_flag .AND. liquid_flag ) qp(n_vars) = r_alphal * r_h

    END IF

    IF ( stoch_transport_flag) qp(idx_stoch) = r_Zs

    IF ( pore_pressure_flag) qp(idx_pore) = r_exc_pore_pres

    qp(idx_u) = r_u
    qp(idx_v) = r_v

    RETURN

  END SUBROUTINE qc_to_qp

  !******************************************************************************
  !> \brief Physical to conservative variables
  !
  !> This subroutine evaluates the conservative real_value variables qc from the
  !> array of real_valued physical variables qp:\n
  !> - qp(1) = \f$ h \f$
  !> - qp(2) = \f$ h*u \f$
  !> - qp(3) = \f$ h*v \f$
  !> - qp(4) = \f$ T \f$
  !> - qp(idx_alfas_first:idx_alfas_last) = \f$ alphas(1:n_s) \f$
  !> - qp(idx_addGas_first:idx_addGas_last) = \f$ alphas(1:n_g) \f$
  !> - qp(n_vars) = \f$ alphal \f$
  !> - qp(idx_u) = \f$ u \f$
  !> - qp(idx_v) = \f$ v \f$
  !> .
  !> \param[in]    qp      physical variables
  !> \param[in]    B       local topography
  !> \param[out]   qc      conservative variables
  !
  !> \date 2019/11/18
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE qp_to_qc(qp,qc)

    USE geometry_2d, ONLY : z_quad, w_quad

    USE parameters_2d, ONLY : vertical_profiles_flag
    USE geometry_2d, ONLY : gaulegf
    USE geometry_2d, ONLY : lambertw,lambertw0,lambertwm1

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qp(n_vars+2)
    REAL(wp), INTENT(OUT) :: qc(n_vars)

    REAL(wp) :: r_sp_heat_mix
    REAL(wp) :: sum_sl

    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity
    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_hu              !< real-value volumetric x-flow
    REAL(wp) :: r_hv              !< real-value volumetric y-flow
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_alphag(n_add_gas)!< real-value add.gas volume fractions
    REAL(wp) :: r_xl              !< real-value liquid mass fraction
    REAL(wp) :: r_xc              !< real-value carrier phase mass fraction
    REAL(wp) :: r_T               !< real-value temperature [K]
    REAL(wp) :: r_alphal          !< real-value liquid volume fraction
    REAL(wp) :: r_alphac          !< real-value carrier phase volume fraction
    REAL(wp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(wp) :: r_rho_c           !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_rho_a           !< real-value atm.gas density [kg/m3]
    REAL(wp) :: r_rho_g(n_add_gas)!< real-value add.gas densities [kg/m3]
    REAL(wp) :: r_xs(n_solid)     !< real-value solid mass fractions
    REAL(wp) :: r_xg(n_add_gas)   !< real-value add.gas mass fractions

    REAL(wp) :: r_Zs              !< real-value stochastic variable
    REAL(wp) :: r_exc_pore_pres   !< real-value pore pressure

    REAL(wp) :: r_alphas_rhos(n_solid)
    REAL(wp) :: r_alphag_rhog(n_add_gas)
    REAL(wp) :: alphas_tot

    REAL(wp) :: r_sp_heat_c

    REAL(wp) :: r_inv_rhom

    REAL(wp) :: rho_u_alphas(n_solid)
    REAL(wp) :: uRho_avg

    REAL(wp) :: u_rel0
    REAL(wp) :: normalizing_coeff_u
    REAL(wp) :: shear_vel

    REAL(wp) :: a , b , c , d

    REAL(wp) :: z(n_quad) , w(n_quad)

    REAL(wp) :: inv_kin_visc
    REAL(wp) :: settling_vel

    REAL(wp) :: Rouse_no(n_solid)

    INTEGER :: i_solid

    REAL(wp) :: r_w
    REAL(wp) :: mod_vel
    REAL(wp) :: log_term_h0
    REAL(wp) :: alphas_exp_avg
    REAL(wp) :: int_quad
    REAL(wp) :: int_def1 , int_def2
    REAL(wp) :: h_rel , h0_rel , h0
    REAL(wp) :: h0_rel_1
    REAL(wp) :: h0_rel_2

    REAL(wp) :: epsilon_s
    REAL(wp) :: a_coeff
    REAL(wp) :: alphas_rel0
    real(wp) :: normalizing_coeff_alpha

    REAL(wp) :: u_log_avg

    REAL(wp) :: dyn_visc_c
    REAL(wp) :: kin_visc_c_local

    r_xl = 0.0_wp
    r_Zs = 0.0_wp
    r_exc_pore_pres = 0.0_wp

    r_h = qp(1)

    IF ( r_h .GT. EPSILON(1.0_wp) ) THEN

       r_hu = qp(2)
       r_hv = qp(3)

       r_u = qp(idx_u)
       r_v = qp(idx_v)

    ELSE

       qc(1:n_vars) = 0.0_wp
       RETURN

    END IF

    r_T  = qp(4)

    IF ( alpha_flag ) THEN

       r_alphas(1:n_solid) = qp(idx_alfas_first:idx_alfas_last)
       r_alphag(1:n_add_gas) = qp(idx_addGas_first:idx_addGas_last)

    ELSE

       r_alphas(1:n_solid) = qp(idx_alfas_first:idx_alfas_last) / qp(1)
       r_alphag(1:n_add_gas) = qp(idx_addGas_first:idx_addGas_last) / qp(1)

    END IF

    alphas_tot = SUM(r_alphas)

    r_alphas_rhos(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid)

    r_alphal = 0.0_wp

    IF ( gas_flag .AND. liquid_flag ) THEN

       IF ( alpha_flag ) THEN

          r_alphal = qp(n_vars)

       ELSE

          r_alphal = qp(n_vars) / qp(1)

       END IF

    END IF

    IF ( gas_flag ) THEN

       ! continuous phase is air
       r_rho_a =  pres / ( sp_gas_const_a * r_T )
       r_rho_g(1:n_add_gas) = pres / ( sp_gas_const_g(1:n_add_gas) * r_T )
       r_alphag_rhog(1:n_add_gas) = r_alphag(1:n_add_gas) * r_rho_g(1:n_add_gas)

       r_rho_c = ( ( 1.0_wp - r_alphal - alphas_tot - SUM(r_alphag) ) * r_rho_a &
            + SUM( r_alphag_rhog(1:n_add_gas) ) ) /                             &
            ( 1.0_wp - r_alphal - alphas_tot )

    ELSE

       ! carrier phase is liquid
       r_rho_c = rho_l

    END IF

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! check and correction on dispersed phases volume fractions
       IF ( ( alphas_tot + r_alphal ) .GT. 1.0_wp ) THEN

          sum_sl = alphas_tot + r_alphal
          r_alphas(1:n_solid) = r_alphas(1:n_solid) / sum_sl
          r_alphal = r_alphal / sum_sl

       ELSEIF ( ( alphas_tot + r_alphal ) .LT. 0.0_wp ) THEN

          r_alphas(1:n_solid) = 0.0_wp
          r_alphal = 0.0_wp

       END IF

       ! carrier phase volume fraction
       r_alphac = 1.0_wp - alphas_tot - r_alphal

       ! volume averaged mixture density: carrier (gas) + solids + liquid
       r_rho_m = r_alphac * r_rho_c + SUM( r_alphas_rhos(1:n_solid) )           &
            + r_alphal * rho_l

       r_inv_rhom = 1.0_wp / r_rho_m

       ! liquid mass fraction
       r_xl = r_alphal * rho_l * r_inv_rhom

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas_rhos(1:n_solid) * r_inv_rhom

       ! additional gas mass fractions
       r_xg(1:n_add_gas) = r_alphag_rhog(1:n_add_gas) * r_inv_rhom

       ! carrier (gas) mass fraction
       r_xc = r_alphac * r_rho_c * r_inv_rhom

       ! specific heat of gas (mass. avg. of sp.heat of gas components)
       IF ( r_xc .GT. EPSILON(1.0_wp) ) THEN

          r_sp_heat_c = ( ( r_xc - SUM( r_xg(1:n_add_gas) ) ) * sp_heat_a +        &
               DOT_PRODUCT( r_xg(1:n_add_gas) , sp_heat_g(1:n_add_gas) ) ) / r_xc

       ELSE

          r_sp_heat_c = sp_heat_a

       END IF

       ! mass averaged mixture specific heat
       r_sp_heat_mix =  DOT_PRODUCT( r_xs , sp_heat_s ) + r_xl * sp_heat_l      &
            + r_xc * r_sp_heat_c

    ELSE

       ! mixture of carrier phase ( gas or liquid ) and solid

       ! check and corrections on dispersed phases
       IF ( alphas_tot .GT. 1.0_wp ) THEN

          r_alphas(1:n_solid) = r_alphas(1:n_solid) / alphas_tot

       ELSEIF ( alphas_tot .LT. 0.0_wp ) THEN

          r_alphas(1:n_solid) = 0.0_wp

       END IF

       ! carrier (gas or liquid) volume fraction
       r_alphac = 1.0_wp - alphas_tot

       ! volume averaged mixture density: carrier (gas or liquid) + solids
       r_rho_m = r_alphac * r_rho_c + SUM( r_alphas_rhos(1:n_solid) )

       r_inv_rhom = 1.0_wp / r_rho_m

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas_rhos(1:n_solid) * r_inv_rhom

       ! additional gas mass fractions
       r_xg(1:n_add_gas) = r_alphag_rhog(1:n_add_gas) * r_inv_rhom

       ! carrier (gas or liquid) mass fraction
       r_xc = r_alphac * r_rho_c * r_inv_rhom

       IF ( gas_flag ) THEN

          IF ( r_xc .GT. EPSILON(1.0_wp) ) THEN

             r_sp_heat_c = ( ( r_xc - SUM( r_xg(1:n_add_gas) ) ) * sp_heat_a +     &
                  DOT_PRODUCT( r_xg(1:n_add_gas) , sp_heat_g(1:n_add_gas) ) ) / r_xc

          ELSE

             r_sp_heat_c = sp_heat_a

          END IF

       ELSE

          r_sp_heat_c = sp_heat_l

       END IF

       ! mass averaged mixture specific heat
       r_sp_heat_mix =  DOT_PRODUCT( r_xs , sp_heat_s ) + r_xc * r_sp_heat_c

    END IF

    IF ( stoch_transport_flag) r_Zs = qp(idx_stoch)

    IF ( pore_pressure_flag ) r_exc_pore_pres = qp(idx_pore)

    qc(1) = r_rho_m * r_h

    IF ( vertical_profiles_flag ) THEN

       r_w = 0.0_wp

       mod_vel = SQRT( r_u**2 + r_v**2 + r_w**2 )

       shear_vel = SQRT( friction_factor ) * mod_vel

       kin_visc_c_local = kin_visc_c

       IF ( gas_flag .AND. sutherland_flag ) THEN

          dyn_visc_c = muRef_Suth * ( r_T / Tref_Suth )**1.5_wp *               &
               ( Tref_Suth + S_mu ) / ( r_T + S_mu )

          kin_visc_c_local = dyn_visc_c / r_rho_c

       END IF

       ! Viscosity read from input file [m2 s-1]
       inv_kin_visc = 1.0_wp / kin_visc_c_local

       DO i_solid=1,n_solid

          settling_vel = settling_velocity( diam_s(i_solid) , rho_s(i_solid) ,  &
               r_rho_c , inv_kin_visc )

          IF ( shear_vel .GT. 0.0_wp ) THEN

             Rouse_no(i_solid) = settling_vel / ( vonK * shear_vel )

          ELSE

             Rouse_no(i_solid) = 0.0_wp

          END IF

       END DO

       ! The profile parameters depend on h/k_s, not on the absolute value of h.
       h_rel = r_h / k_s

       IF ( h_rel .GT. H_crit_rel ) THEN

          ! we search for h0_rel such that the average integral between 0 and
          ! h_rel is equal to 1
          ! For h_rel > H_crit_rel this integral is the sum of two pieces:
          ! integral between 0 and h0_rel of the log profile
          ! integral between h0_rel and h_rel of the costant profile

          a = h_rel * vonK / SQRT(friction_factor)
          b = 1.0_wp / 30.0_wp + h_rel
          c = 30.0_wp

          ! solve b*log(c*z+1)-z=a for z
          d = a/b - 1.0_wp / (b*c)

          h0_rel_1 = -b*lambertw0( -EXP(d)/(b*c) ) - 1.0_wp / c
          h0_rel_2 = -b*lambertwm1( -EXP(d)/(b*c) ) - 1.0_wp / c
          h0_rel = MIN( h0_rel_1 , h0_rel_2)

       ELSE

          ! when h_rel <= H_crit_rel we have only the log profile and we have to
          ! rescale it in order to have the integral between o and h_rel equal to
          ! 1
          h0_rel = h_rel

       END IF

       h0 = h0_rel * k_s

       b = 30.0_wp / k_s

       log_term_h0 = u_log_profile(b,h0)**2

       epsilon_s = Sc * shear_vel * vonK * (( h0 / 6.0_wp ) + ( k_s / 60.0_wp ))
       a_coeff = - vonK * shear_vel / epsilon_s

       z = 0.5_wp * h0 * ( z_quad + 1.0_wp )
       w = 0.5_wp * h0 * w_quad

       u_log_avg = ( SUM( w * u_log_profile(b,z) ) + u_log_profile(b,h0) *      &
            ( r_h -h0 ) ) / r_h

       normalizing_coeff_u = 1.0_wp / u_log_avg

       ! relative velocty at h0
       u_rel0 = normalizing_coeff_u * u_log_profile(b,h0)

       DO i_solid = 1,n_solid

          a = a_coeff * Rouse_no(i_solid)

          alphas_exp_avg = ( SUM( w * alphas_exp_profile(a,z) ) +               &
               alphas_exp_profile(a,h0)*(r_h-h0) ) / r_h

          normalizing_coeff_alpha = 1.0_wp / alphas_exp_avg

          int_quad = SUM( w * ( alphas_exp_profile(a,z) * u_log_profile(b,z) ) )

          ! integral of alfa_i*u between in the boundary layer
          int_def1 = ( mod_vel * normalizing_coeff_u ) *                        &
               ( r_alphas(i_solid) * normalizing_coeff_alpha ) * int_quad

          ! relative concentration alphas_rel at h0
          alphas_rel0 = normalizing_coeff_alpha * alphas_exp_profile(a,h0)

          ! integral of alfa_i*u in the free-stream layer
          int_def2 =  ( r_h - h0 ) * ( mod_vel * u_rel0 ) *                     &
               ( alphas_rel0 * r_alphas(i_solid) )

          ! we add the contribution of the integral of the constant region, we
          ! average by dividing by h and we multiply by the density of solid and
          ! average concentration and by u_guess.
          rho_u_alphas(i_solid) = rho_s(i_solid) * ( int_def1 + int_def2 ) / r_h

       END DO

       ! we add the contribution of the gas phase to the mixture depth-averaged
       ! momentum
       uRho_avg = ( mod_vel * r_rho_c + SUM((rho_s - r_rho_c) / rho_s *         &
            rho_u_alphas) )

       qc(2) = r_h * uRho_avg * r_u / mod_vel
       qc(3) = r_h * uRho_avg * r_v / mod_vel

    ELSE

       qc(2) = r_rho_m * r_hu
       qc(3) = r_rho_m * r_hv

    END IF

    IF ( energy_flag ) THEN

       IF ( r_h .GT. 0.0_wp ) THEN

          ! total energy (internal and kinetic)
          qc(4) = r_h * r_rho_m * ( r_sp_heat_mix * r_T                         &
               + 0.5_wp * ( r_u**2 + r_v**2 ) )

       ELSE

          qc(4) = 0.0_wp

       END IF

    ELSE

       ! internal energy
       qc(4) = r_h * r_rho_m * r_sp_heat_mix * r_T

    END IF

    qc(idx_alfas_first:idx_alfas_last) = r_xs * qc(1)
    qc(idx_addGas_first:idx_addGas_last) = r_xg * qc(1)

    IF ( stoch_transport_flag ) qc(idx_stoch) = r_Zs * qc(1)

    IF ( pore_pressure_flag ) qc(idx_pore) = r_exc_pore_pres * qc(1)

    IF ( gas_flag .AND. liquid_flag ) qc(n_vars) = r_xl * qc(1)

    RETURN

  END SUBROUTINE qp_to_qc

  !******************************************************************************
  !> \brief Additional Physical variables
  !
  !> This subroutine evaluates from the physical local variables qpj, the two
  !> additional local variables qp2j = (h+B,u,v).
  !> \param[in]    qpj    real-valued physical variables
  !> \param[in]    Bj     real-valued local topography
  !> \param[out]   qp2j   real-valued physical variables
  !> @author
  !> Mattia de' Michieli Vitturi
  !> \date 10/10/2019
  !******************************************************************************

  SUBROUTINE qp_to_qp2(qpj,Bj,qp2j)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: Bj
    REAL(wp), INTENT(OUT) :: qp2j(3)

    qp2j(1) = qpj(1) + Bj

    IF ( qpj(1) .LE. 0.0_wp ) THEN

       qp2j(2) = 0.0_wp
       qp2j(3) = 0.0_wp

    ELSE

       qp2j(2) = qpj(2)/qpj(1)
       qp2j(3) = qpj(3)/qpj(1)

    END IF

    RETURN

  END SUBROUTINE qp_to_qp2

  !------------------------------------------------------------------------------
  !> Settling velocity function
  !
  !> This subroutine compute the settling velocity of the particles, as a
  !> function of diameter, density of particles and carrier phase and viscosity.
  !> \date 2019/11/11
  !> \param[in]    diam          particle diameter
  !> \param[in]    rhos          particle density
  !> \param[in]    rhoc          carrier phase density
  !> \param[in]    inv_kin_visc  reciprocal of kinetic viscosity
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !------------------------------------------------------------------------------

  REAL(wp) FUNCTION settling_velocity(diam,rhos,rhoc,inv_kin_visc)

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: diam          !< particle diameter [m]
    REAL(wp), INTENT(IN) :: rhos          !< particle density [kg/m3]
    REAL(wp), INTENT(IN) :: rhoc          !< carrier phase density [kg/m3]
    REAL(wp), INTENT(IN) :: inv_kin_visc  !< carrier phase viscosity reciprocal

    REAL(wp) :: Rey           !< Reynolds number
    REAL(wp) :: inv_sqrt_C_D  !< Reciprocal of sqrt of Drag coefficient

    ! loop variables
    INTEGER :: i              !< loop counter for iterative procedure
    REAL(wp) :: const_part    !< term not changing in iterative procedure
    REAL(wp) :: inv_sqrt_C_D_old  !< previous iteration sqrt of drag coefficient
    REAL(wp) :: set_vel_old   !< previous iteration settling velocity

    ! INTEGER :: dig          !< order of magnitude of settling velocity

    inv_sqrt_C_D = 1.0_wp

    IF ( rhos .LE. rhoc ) THEN

       settling_velocity = 0.0_wp
       RETURN

    END IF

    const_part =  SQRT( (4.0_wp/3.0_wp) * ( rhos / rhoc - 1.0_wp ) * diam * grav )

    settling_velocity = const_part * inv_sqrt_C_D

    Rey = diam * settling_velocity * inv_kin_visc

    IF ( Rey .LE. 1000.0_wp ) THEN

       C_D_loop:DO i=1,20

          set_vel_old = settling_velocity
          inv_sqrt_C_D_old = inv_sqrt_C_D

          IF ( collective_settling_flag ) THEN

             ! Optional collective-settling drag law from the development branch
             inv_sqrt_C_D = SQRT( 1.0_wp / ( 24.0_wp * ( 1.0_wp +              &
                  0.15_wp*Rey**(0.687_wp) ) / Rey +                             &
                  A_drag * ( LOG(Rey) - B_drag ) ) )

          ELSE

             ! Original Schiller-Naumann drag law
             inv_sqrt_C_D = SQRT( Rey / ( 24.0_wp * ( 1.0_wp +                 &
                  0.15_wp*Rey**(0.687_wp) ) ) )

          END IF

          settling_velocity = const_part * inv_sqrt_C_D

          IF ( ABS( set_vel_old - settling_velocity ) / set_vel_old             &
               .LT. 1.0E-6_wp ) THEN

!!$             ! round to first three significative digits
!!$             dig = FLOOR(LOG10(set_vel_old))
!!$             settling_velocity = 10.0_wp**(dig-3)                            &
!!$                  * FLOOR( 10.0_wp**(-dig+3)*set_vel_old )

             RETURN

          END IF

          Rey = diam * settling_velocity * inv_kin_visc

       END DO C_D_loop

    END IF

    RETURN

  END FUNCTION settling_velocity

END MODULE state_conversion_2d
