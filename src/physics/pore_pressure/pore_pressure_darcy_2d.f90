!> Darcy pore-pressure diffusion and optional interstitial carrier-gas loss.
!> The four model interfaces return conservative rates; the analytic interface
!> integrates the coupled frozen-state mass and pressure system.
!> Gas loss is also still used by the legacy mass-exchange transaction while
!> shared carrier-inventory and packing limits are being migrated.

MODULE pore_pressure_darcy_2d

  USE parameters_2d, ONLY : wp, n_solid, gas_flag, gas_loss_flag, pore_pressure_flag, &
       n_eqns, n_add_gas, idx_totMassEqn, idx_uEqn, idx_vEqn, idx_engyEqn, &
       idx_poreEqn, idx_addGasEqn_first, idx_solidEqn_first, &
       idx_solidEqn_last, sutherland_flag
  USE COMPLEXIFY
  USE diagnostics_2d, ONLY : fatal_error
  USE constitutive_parameters_2d, ONLY : alpha_trans, N_inh, f_inhibit_mode,    &
       pascal_coeff, hydraulic_permeability, kin_visc_c, sp_heat_a,             &
       sp_gas_const_a, pres, diam_s, sphericity_s, maximum_solid_packing,       &
       dynamic_permeability_flag, muRef_Suth, Tref_Suth, S_mu
  USE geometry_2d, ONLY : pi_g
  USE state_conversion_2d, ONLY : sauter_diameter, physical_state_type, &
       complex_physical_state_type, complex_state_from_real, mixt_var_state, &
       eval_mixture_heat_capacity
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_pore_pressure_darcy_explicit
  PUBLIC :: eval_pore_pressure_darcy_semi_implicit
  PUBLIC :: eval_pore_pressure_darcy_implicit
  PUBLIC :: solve_pore_pressure_darcy_analytic

CONTAINS

  !> The model returns a single conservative source vector, including
  !> diffusion and (when enabled) gas withdrawal. Coordination of shared
  !> carrier limits belongs to the source executor, not this model.
  !> Assemble the complete real conservative rate without repeating physics.
  SUBROUTINE eval_pore_pressure_darcy_explicit(state,rate,include_gas_loss)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(OUT) :: rate(n_eqns)
    LOGICAL, INTENT(IN), OPTIONAL :: include_gas_loss
    TYPE(complex_physical_state_type) :: complex_state
    REAL(wp) :: semi_rate(n_eqns)
    COMPLEX(wp) :: implicit_rate(n_eqns)

    rate=0.0_wp
    IF (.NOT.pore_pressure_flag) RETURN
    CALL eval_pore_pressure_darcy_semi_implicit(state,semi_rate,include_gas_loss)
    CALL complex_state_from_real(state,complex_state)
    CALL eval_pore_pressure_darcy_implicit(complex_state,implicit_rate,include_gas_loss)
    rate=semi_rate+REAL(implicit_rate,wp)
  END SUBROUTINE eval_pore_pressure_darcy_explicit

  !> No separately semi-implicit Darcy contribution.
  SUBROUTINE eval_pore_pressure_darcy_semi_implicit(state,rate,include_gas_loss)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(OUT) :: rate(n_eqns)
    LOGICAL, INTENT(IN), OPTIONAL :: include_gas_loss
    rate=0.0_wp
  END SUBROUTINE eval_pore_pressure_darcy_semi_implicit

  !> Complete complex Darcy rate, including mass-pressure coupling.
  SUBROUTINE eval_pore_pressure_darcy_implicit(state,rate,include_gas_loss)
    TYPE(complex_physical_state_type), INTENT(IN) :: state
    COMPLEX(wp), INTENT(OUT) :: rate(n_eqns)
    LOGICAL, INTENT(IN), OPTIONAL :: include_gas_loss
    COMPLEX(wp) :: mass_rate,cp_c,cp_m,xc,gas_fraction
    COMPLEX(wp) :: diffusion_coefficient,inhibition
    LOGICAL :: active
    LOGICAL :: do_gas_loss

    rate=CMPLX(0.0_wp,0.0_wp,wp)
    IF (.NOT.pore_pressure_flag) RETURN
    IF (REAL(state%h,wp)>0.0_wp .AND. &
        REAL(state%exc_pore_pres,wp)>=0.0_wp) THEN
      CALL eval_pore_pressure_diffusion_factors(state%h,state%T, &
           state%alphas,diffusion_coefficient,inhibition,active)
      IF (active) rate(idx_poreEqn)= &
           -state%rho_m*(pi_g/2.0_wp)**2*diffusion_coefficient* &
           state%exc_pore_pres*inhibition / &
           MERGE(state%h,CMPLX(1.0e-10_wp,0.0_wp,wp), &
           REAL(state%h,wp)>1.0e-10_wp)
    END IF
    do_gas_loss=gas_flag .AND. gas_loss_flag
    IF (PRESENT(include_gas_loss)) do_gas_loss=do_gas_loss .AND. include_gas_loss
    IF (.NOT.do_gas_loss) RETURN
    IF (REAL(SUM(state%alphas),wp)>=maximum_solid_packing) RETURN
    IF (REAL(state%exc_pore_pres,wp)<=0.0_wp) RETURN
    ! With unit carrier density, Darcy velocity equals the mass-loss flux.
    CALL pore_gas_loss_complex(state%h,state%T,state%alphas, &
         SUM(state%alphas),CMPLX(1.0_wp,0.0_wp,wp), &
         state%exc_pore_pres,mass_rate)
    mass_rate=-mass_rate
    CALL eval_mixture_heat_capacity(state%xs,state%xg,xc,cp_c,cp_m)
    rate(idx_totMassEqn)=mass_rate
    rate(idx_uEqn)=mass_rate*state%u
    rate(idx_vEqn)=mass_rate*state%v
    rate(idx_engyEqn)=mass_rate*state%T*cp_c
    IF (n_add_gas>0) THEN
      gas_fraction=1.0_wp-SUM(state%xs)
      IF (REAL(gas_fraction,wp)>0.0_wp) &
           rate(idx_addGasEqn_first:idx_addGasEqn_first+n_add_gas-1)= &
           mass_rate*state%xg/gas_fraction
    END IF
    rate(idx_poreEqn)=rate(idx_poreEqn)+mass_rate*state%exc_pore_pres
  END SUBROUTINE eval_pore_pressure_darcy_implicit

  !> Integrate frozen-coefficient Darcy directly and return the full
  !> conservative mean rate for the selected substep.
  SUBROUTINE solve_pore_pressure_darcy_analytic( &
       state,q_before,interval,rate,admissible,include_gas_loss)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: interval,q_before(n_eqns)
    REAL(wp), INTENT(OUT) :: rate(n_eqns)
    LOGICAL, INTENT(OUT) :: admissible
    REAL(wp) :: q_after(n_eqns)
    LOGICAL, INTENT(IN), OPTIONAL :: include_gas_loss
    LOGICAL :: do_gas_loss
    COMPLEX(wp) :: diffusion, inhibition
    REAL(wp) :: lambda, A, loss, mass_after, pressure_after, pressure_q
    REAL(wp) :: rho_c, cp_c, Ri, rhom, redgrav, cpm
    REAL(wp) :: gas_inventory, gas_fraction
    REAL(wp) :: exponential_factor, pressure_integral, coefficient
    LOGICAL :: active
    COMPLEX(wp) :: solid_fractions(n_solid)
    REAL(wp), PARAMETER :: h_min=1.0e-10_wp
    REAL(wp), PARAMETER :: gas_h_min=1.0e-5_wp

    q_after=q_before
    rate=0.0_wp
    admissible=.FALSE.
    IF (.NOT.pore_pressure_flag) THEN
      admissible=.TRUE.
      RETURN
    END IF
    IF (interval<0.0_wp .OR. .NOT.ieee_is_finite(interval)) RETURN
    IF (q_before(idx_totMassEqn)<=0.0_wp) RETURN
    IF (state%h<=0.0_wp .OR. state%T<=0.0_wp) THEN
      admissible=.TRUE.
      RETURN
    END IF
    solid_fractions=CMPLX(state%alphas,0.0_wp,wp)
    CALL eval_pore_pressure_diffusion_factors(CMPLX(state%h,0.0_wp,wp), &
         CMPLX(state%T,0.0_wp,wp),solid_fractions,diffusion,inhibition,active)
    lambda=0.0_wp
    IF (active) lambda=REAL((pi_g/2.0_wp)**2*diffusion*inhibition,wp) / &
         (state%h*MAX(h_min,state%h))
    A=0.0_wp
    do_gas_loss=gas_flag .AND. gas_loss_flag
    IF (PRESENT(include_gas_loss)) do_gas_loss=do_gas_loss .AND. include_gas_loss
    IF (do_gas_loss .AND. &
         SUM(state%alphas)<maximum_solid_packing) THEN
      CALL mixt_var_state(state,Ri,rhom,rho_c,redgrav,cp_c,cpm)
      IF (state%exc_pore_pres>0.0_wp) THEN
        ! This velocity is linear in p_e with the frozen closures.
        ! Evaluate at p_e=1 Pa to obtain the mass-loss coefficient A.
        BLOCK
          COMPLEX(wp) :: gas_loss_coefficient
          CALL pore_gas_loss_complex(CMPLX(state%h,0.0_wp,wp), &
               CMPLX(state%T,0.0_wp,wp),CMPLX(state%alphas,0.0_wp,wp), &
               CMPLX(SUM(state%alphas),0.0_wp,wp), &
               CMPLX(rho_c,0.0_wp,wp),CMPLX(1.0_wp,0.0_wp,wp), &
               gas_loss_coefficient)
          A=rho_c*REAL(gas_loss_coefficient,wp)
        END BLOCK
      END IF
    END IF
    mass_after=q_before(idx_totMassEqn)
    pressure_after=state%exc_pore_pres
    pressure_q=mass_after*pressure_after
    loss=0.0_wp
    IF (.NOT.ieee_is_finite(lambda) .OR. .NOT.ieee_is_finite(A) .OR. &
        .NOT.ieee_is_finite(state%exc_pore_pres)) RETURN
    IF (lambda<0.0_wp .OR. A<0.0_wp) RETURN
    admissible=.TRUE.
    IF (state%exc_pore_pres>0.0_wp .AND. interval>0.0_wp) THEN
      IF (.NOT.ieee_is_finite(q_before(idx_totMassEqn)) .OR. &
          .NOT.ieee_is_finite(state%exc_pore_pres) .OR. &
          .NOT.ieee_is_finite(interval) .OR. &
          .NOT.ieee_is_finite(lambda) .OR. &
          .NOT.ieee_is_finite(A)) RETURN
      IF (q_before(idx_totMassEqn)<=0.0_wp .OR. interval<0.0_wp .OR. &
          lambda<0.0_wp .OR. A<0.0_wp) RETURN
  
      ! The existing v2/dev diffusion closure does not dissipate
      ! negative excess pore pressure; preserve that guard.
  
      exponential_factor=EXP(-lambda*interval)
      pressure_after=state%exc_pore_pres*exponential_factor
  
      coefficient=0.0_wp
      IF (do_gas_loss) coefficient=A
      IF (lambda>0.0_wp) THEN
        ! Stable evaluation of (1-exp(-lambda*dt))/lambda,
        ! including small lambda*dt without cancellation.
        IF (ABS(lambda*interval)<1.0e-5_wp) THEN
          pressure_integral=interval*(1.0_wp-0.5_wp*lambda*interval+ &
               (lambda*interval)**2/6.0_wp- &
               (lambda*interval)**3/24.0_wp)
        ELSE
          pressure_integral=(1.0_wp-exponential_factor)/lambda
        END IF
      ELSE
        pressure_integral=interval
      END IF
  
      loss=coefficient*state%exc_pore_pres*pressure_integral
      mass_after=q_before(idx_totMassEqn)-loss
      IF (.NOT.ieee_is_finite(mass_after)) RETURN
      pressure_q=mass_after*pressure_after
  
    END IF
    gas_inventory=q_before(idx_totMassEqn)
    gas_inventory=gas_inventory-SUM(q_before( &
         idx_solidEqn_first:idx_solidEqn_last))
    IF (gas_inventory<0.0_wp) THEN
      admissible=.FALSE.
      RETURN
    END IF
    ! Independent Darcy inventory limiter: do not withdraw more carrier
    ! mass than the frozen state contains. Preserve a positive total mass.
    IF (loss>0.0_wp) THEN
      loss=MIN(loss,MAX(0.0_wp,gas_inventory), &
           MAX(0.0_wp,q_before(idx_totMassEqn)- &
           EPSILON(1.0_wp)*MAX(1.0_wp,q_before(idx_totMassEqn))))
      mass_after=q_before(idx_totMassEqn)-loss
      pressure_q=mass_after*pressure_after
    END IF
    IF (n_add_gas>0 .AND. loss>0.0_wp) THEN
      IF (ANY(q_before(idx_addGasEqn_first: &
           idx_addGasEqn_first+n_add_gas-1)<0.0_wp)) THEN
        admissible=.FALSE.
        RETURN
      END IF
    END IF
    q_after(idx_totMassEqn)=mass_after
    q_after(idx_poreEqn)=pressure_q
    q_after(idx_uEqn)=q_before(idx_uEqn)-state%u*loss
    q_after(idx_vEqn)=q_before(idx_vEqn)-state%v*loss
    IF (loss>0.0_wp) THEN
      ! No gas is removed from solid-class conservative masses.
      q_after(idx_engyEqn)=q_before(idx_engyEqn)- &
           state%T*cp_c*loss
      IF (n_add_gas>0) THEN
        gas_fraction=loss/gas_inventory
        q_after(idx_addGasEqn_first: &
             idx_addGasEqn_first+n_add_gas-1)= &
             q_before(idx_addGasEqn_first: &
             idx_addGasEqn_first+n_add_gas-1)*(1.0_wp-gas_fraction)
      END IF
    END IF
    IF (interval>0.0_wp) rate=(q_after-q_before)/interval
  END SUBROUTINE solve_pore_pressure_darcy_analytic

  ! Internal Darcy closures and other auxiliary source calculations.







  !******************************************************************************
  !> \brief Evaluate the pore-gas loss rate from Darcy outflow.
  !> \details Computes an upward pore-gas velocity from local permeability,
  !>          kinematic viscosity, carrier density, flow depth and excess
  !>          pressure, then multiplies it by the configured solid-fraction
  !>          inhibition factor. The velocity is
  !>          \f$v_D=\frac{\pi k p_e}{2\nu_g\rho_c\max(10^{-5}\,\mathrm{m},h)}\f$,
  !>          and \f$\texttt{gas\_loss}=Iv_D\f$. The dynamic-permeability and Sutherland
  !>          options use the same mixture closures as pore-pressure diffusion.
  !> \param[in] flow_depth Flow depth h [m].
  !> \param[in] temperature Gas temperature T [K].
  !> \param[in] solid_fractions Per-class solid volume fractions.
  !> \param[in] total_solid_fraction Sum of solid volume fractions.
  !> \param[in] carrier_density Carrier-phase density [kg m^-3].
  !> \param[in] excess_pressure Excess pore pressure p_e [Pa].
  !> \param[out] gas_loss Inhibited upward carrier-gas loss velocity [m s^-1].
  !******************************************************************************


  SUBROUTINE pore_gas_loss_complex(flow_depth,temperature,solid_fractions, &
       total_solid_fraction,carrier_density,excess_pressure,gas_loss)
    COMPLEX(wp), INTENT(IN) :: flow_depth,temperature
    COMPLEX(wp), INTENT(IN) :: solid_fractions(n_solid)
    COMPLEX(wp), INTENT(IN) :: total_solid_fraction,carrier_density,excess_pressure
    COMPLEX(wp), INTENT(OUT) :: gas_loss
    COMPLEX(wp) :: diffusion,inhibition,mobility,depth
    LOGICAL :: active

    gas_loss=CMPLX(0.0_wp,0.0_wp,wp)
    CALL eval_pore_pressure_diffusion_factors(flow_depth,temperature, &
         solid_fractions,diffusion,inhibition,active,mobility)
    IF (.NOT.active) RETURN
    IF (REAL(carrier_density,wp)<=0.0_wp) RETURN
    depth=MERGE(flow_depth,CMPLX(1.0e-5_wp,0.0_wp,wp), &
         REAL(flow_depth,wp)>1.0e-5_wp)
    gas_loss=(pi_g/2.0_wp)*mobility*excess_pressure*inhibition/ &
         (carrier_density*depth)
  END SUBROUTINE pore_gas_loss_complex

  !******************************************************************************
  !> \brief Compute pressure-diffusion and inhibition factors for a local mixture.
  !> \details Computes gas compressibility and density from the ideal-gas
  !>          closure, pore volume from phi=1-sum(alpha_i), and the diffusion
  !>          coefficient
  !>          \f$D=k/(\phi\nu_g\rho_gC_g)\f$. Permeability may be fixed
  !>          or evaluated from the Kozeny-Carman mixture closure; the Sutherland
  !>          viscosity law is used when configured. The inhibition output
  !>          implements the selected OFF, STEP, STATIC or DYNAMIC closure.
  !> \param[in] flow_depth Complex-valued flow depth h [m].
  !> \param[in] temperature Complex-valued gas/mixture temperature T [K].
  !> \param[in] solid_fractions Complex-valued per-class solid volume fractions.
  !> \param[out] diffusion_coefficient Pressure-diffusion coefficient D [m^2 s^-1].
  !> \param[out] inhibition Dimensionless solid-fraction inhibition factor I.
  !> \param[out] active True when the state has carrier/pore volume for diffusion;
  !>                    otherwise both returned factors are their neutral defaults.
  !> \note Real parts select piecewise branches while smooth closure algebra retains
  !>       complex perturbations for complex-step differentiation.
  !******************************************************************************
  SUBROUTINE eval_pore_pressure_diffusion_factors(flow_depth, temperature,    &
       solid_fractions, diffusion_coefficient, inhibition, active, mobility)

    COMPLEX(wp), INTENT(IN) :: flow_depth, temperature
    COMPLEX(wp), INTENT(IN) :: solid_fractions(n_solid)
    COMPLEX(wp), INTENT(OUT) :: diffusion_coefficient, inhibition
    LOGICAL, INTENT(OUT) :: active
    COMPLEX(wp), INTENT(OUT), OPTIONAL :: mobility

    COMPLEX(wp) :: gas_compressibility, gas_density, porosity
    COMPLEX(wp) :: dynamic_viscosity, local_permeability
    COMPLEX(wp) :: local_kinematic_viscosity, sauter_diam
    COMPLEX(wp) :: mean_sphericity, x, smoothstep
    COMPLEX(wp) :: normalized_fraction, term, width, dynamic_transition
    REAL(wp) :: gamma_gas
    INTEGER :: n

    diffusion_coefficient = CMPLX(0.0_wp, 0.0_wp, wp)
    inhibition = CMPLX(1.0_wp, 0.0_wp, wp)
    active = .FALSE.
    IF (PRESENT(mobility)) mobility=CMPLX(0.0_wp,0.0_wp,wp)

    IF (REAL(1.0_wp-SUM(solid_fractions),wp) .LE. EPSILON(1.0_wp)) RETURN
    IF (dynamic_permeability_flag) THEN
       IF (REAL(SUM(solid_fractions),wp) .LE. EPSILON(1.0_wp)) RETURN
    END IF
    active = .TRUE.

    local_permeability = CMPLX(hydraulic_permeability, 0.0_wp, wp)
    local_kinematic_viscosity = CMPLX(kin_visc_c, 0.0_wp, wp)
    gamma_gas = sp_heat_a / (sp_heat_a - sp_gas_const_a)
    gas_compressibility = 1.0_wp / (gamma_gas * pres)
    gas_density = pres / (sp_gas_const_a * temperature)
    porosity = 1.0_wp - SUM(solid_fractions)

    IF (dynamic_permeability_flag) THEN
       IF (n_solid .GT. 1) THEN
          sauter_diam = SUM(solid_fractions) / SUM(solid_fractions / diam_s)
          mean_sphericity = SUM(solid_fractions * diam_s**2 * sphericity_s) / &
               SUM(solid_fractions * diam_s**2)
       ELSE
          sauter_diam = diam_s(1)
          mean_sphericity = sphericity_s(1)
       END IF
       local_permeability = (porosity**3.0_wp *                             &
            (sauter_diam * mean_sphericity)**2.0_wp) /                      &
            (150.0_wp * (1.0_wp - porosity)**2.0_wp)
    END IF

    IF (gas_flag .AND. sutherland_flag) THEN
       dynamic_viscosity = muRef_Suth * (temperature / Tref_Suth)**1.5_wp *  &
            (Tref_Suth + S_mu) / (temperature + S_mu)
       local_kinematic_viscosity = dynamic_viscosity / gas_density
    END IF

    IF (PRESENT(mobility)) mobility=local_permeability/local_kinematic_viscosity
    diffusion_coefficient = local_permeability /                            &
         (porosity * local_kinematic_viscosity * gas_density *               &
         gas_compressibility)

    SELECT CASE (TRIM(ADJUSTL(f_inhibit_mode)))
    CASE ('OFF')
       inhibition = 1.0_wp

    CASE ('STEP')
       IF (REAL(SUM(solid_fractions),wp) .LT. maximum_solid_packing) THEN
          inhibition = 1.0_wp
       ELSE
          inhibition = 0.0_wp
       END IF

    CASE ('STATIC')
       IF (N_inh .LT. 0) THEN
          WRITE(*,*) 'ERROR: N_inh must be >= 0. N_inh=', N_inh
          STOP
       END IF
       x = SUM(solid_fractions)
       width = maximum_solid_packing - alpha_trans
       normalized_fraction = (x - alpha_trans) / width
       IF (REAL(normalized_fraction,wp) .LT. 0.0_wp) normalized_fraction = 0.0_wp
       IF (REAL(normalized_fraction,wp) .GT. 1.0_wp) normalized_fraction = 1.0_wp
       smoothstep = CMPLX(0.0_wp, 0.0_wp, wp)
       DO n = 0, N_inh
          term = pascal_coeff(n) * normalized_fraction**(N_inh + n + 1)
          smoothstep = smoothstep + term
       END DO
       inhibition = 1.0_wp - smoothstep

    CASE ('DYNAMIC')
       dynamic_transition = 10.0_wp**(-0.28_wp) * flow_depth**0.12_wp
       IF (REAL(dynamic_transition,wp) .GE. maximum_solid_packing) THEN
          dynamic_transition = 0.99_wp * maximum_solid_packing
       END IF
       IF (N_inh .LT. 0) THEN
          WRITE(*,*) 'ERROR: N_inh must be >= 0. N_inh=', N_inh
          STOP
       END IF
       x = SUM(solid_fractions)
       width = maximum_solid_packing - dynamic_transition
       normalized_fraction = (x - dynamic_transition) / width
       IF (REAL(normalized_fraction,wp) .LT. 0.0_wp) normalized_fraction = 0.0_wp
       IF (REAL(normalized_fraction,wp) .GT. 1.0_wp) normalized_fraction = 1.0_wp
       smoothstep = CMPLX(0.0_wp, 0.0_wp, wp)
       DO n = 0, N_inh
          term = pascal_coeff(n) * normalized_fraction**(N_inh + n + 1)
          smoothstep = smoothstep + term
       END DO
       inhibition = 1.0_wp - smoothstep
    END SELECT

  END SUBROUTINE eval_pore_pressure_diffusion_factors



END MODULE pore_pressure_darcy_2d
