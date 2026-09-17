!********************************************************************************
!> \brief State conversion and directly derived physical quantities
!********************************************************************************
MODULE state_conversion_2d

  USE constitutive_parameters_2d

  USE parameters_2d, ONLY : wp
  USE parameters_2d, ONLY : n_vars, n_solid, n_add_gas
  USE parameters_2d, ONLY : energy_flag, liquid_flag, gas_flag, alpha_flag,     &
       stoch_transport_flag, pore_pressure_flag, sutherland_flag

  USE parameters_2d, ONLY : idx_alfas_first, idx_alfas_last, idx_addGas_first,  &
       idx_addGas_last, idx_stoch, idx_pore, idx_u, idx_v

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: r_phys_var, c_phys_var
  PUBLIC :: qc_to_qp, qp_to_qc, qp_to_qp2
  PUBLIC :: mixt_var, eval_sp_heat
  PUBLIC :: sauter_diameter, average_density_solids
  PUBLIC :: settling_velocity

CONTAINS

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

  SUBROUTINE r_phys_var(qj, h, u, v, alphas, rho_m, T, alphal, alphag,          &
       red_grav, p_dyn, Zs, exc_pore_pres)

    USE parameters_2d, ONLY : eps_sing, eps_sing4

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: qj(n_vars)
    REAL(wp), INTENT(OUT) :: h, u, v, T, rho_m
    REAL(wp), INTENT(OUT) :: alphas(n_solid), alphag(n_add_gas)
    REAL(wp), INTENT(OUT) :: alphal, red_grav, p_dyn, Zs, exc_pore_pres

    REAL(wp) :: inv_rhom, rho_c, inv_rho_c, xs_tot
    REAL(wp) :: xs(n_solid), xg(n_add_gas)
    REAL(wp) :: xl, xc, inv_qj1
    REAL(wp) :: carrier_sp_heat_mass, sp_heat_mix
    REAL(wp) :: sp_gas_const_c, inv_rho_g(n_add_gas)
    LOGICAL :: is_wet

    ! The shared fragment returns immediately for a dry state, so initialize
    ! the REAL-only derived outputs before entering it.
    red_grav = 0.0_wp
    p_dyn = 0.0_wp

    INCLUDE 'state_conversion_phys_var.inc'

    ! Preserve the historical update of the carrier specific heat.
    IF ( .NOT. gas_flag ) sp_heat_c = sp_heat_l

    red_grav = ( rho_m - rho_a_amb ) * inv_rhom * grav

    p_dyn = 0.5_wp * rho_m * ( u**2 + v**2 )

  END SUBROUTINE r_phys_var


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

  SUBROUTINE c_phys_var(qj, h, u, v, T, rho_m, alphas, alphag, inv_rhom, Zs,    &
       exc_pore_pres)

    USE parameters_2d, ONLY : eps_sing, eps_sing4

    IMPLICIT NONE

    COMPLEX(wp), INTENT(IN) :: qj(n_vars)
    COMPLEX(wp), INTENT(OUT) :: h, u, v, T, rho_m
    COMPLEX(wp), INTENT(OUT) :: alphas(n_solid), alphag(n_add_gas)
    COMPLEX(wp), INTENT(OUT) :: inv_rhom, Zs, exc_pore_pres

    COMPLEX(wp) :: alphal, rho_c, inv_rho_c, xs_tot
    COMPLEX(wp) :: xs(n_solid), xg(n_add_gas)
    COMPLEX(wp) :: xl, xc, inv_qj1
    COMPLEX(wp) :: carrier_sp_heat_mass, sp_heat_mix
    COMPLEX(wp) :: sp_gas_const_c, inv_rho_g(n_add_gas)
    LOGICAL :: is_wet

    INCLUDE 'state_conversion_phys_var.inc'

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

    qc(2) = r_rho_m * r_hu
    qc(3) = r_rho_m * r_hv

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
