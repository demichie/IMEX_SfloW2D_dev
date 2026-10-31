!********************************************************************************
!> \brief Excess pore-pressure diffusion and gas-loss source laws.
!> \details Evaluates Darcy-type excess-pressure diffusion, pore-gas loss,
!>          pressure sources induced by inlet injection and mass exchange, and
!>          the constitutive diffusion/inhibition factors shared by those
!>          calculations. The diffusion routine accepts complex states for
!>          Jacobian differentiation. Separate real-valued routines provide
!>          both the exact frozen-state ODE flow and the closed-form algebraic
!>          DIRK stage solve.
!********************************************************************************

MODULE pore_pressure_darcy_2d

  USE parameters_2d, ONLY : wp, n_solid, gas_flag, gas_loss_flag, pore_pressure_flag, &
       n_eqns, n_add_gas, idx_totMassEqn, idx_uEqn, idx_vEqn, idx_engyEqn, &
       idx_poreEqn, idx_addGasEqn_first, sutherland_flag
  USE COMPLEXIFY
  USE diagnostics_2d, ONLY : fatal_error
  USE constitutive_parameters_2d, ONLY : alpha_trans, N_inh, f_inhibit_mode,    &
       pascal_coeff, hydraulic_permeability, kin_visc_c, sp_heat_a,             &
       sp_gas_const_a, pres, diam_s, sphericity_s, maximum_solid_packing,       &
       dynamic_permeability_flag, muRef_Suth, Tref_Suth, S_mu
  USE geometry_2d, ONLY : pi_g
  USE state_conversion_2d, ONLY : sauter_diameter, physical_state_type, &
       complex_physical_state_type, complex_state_from_real, mixt_var_state
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: apply_frozen_darcy_withdrawal
  PUBLIC :: eval_pore_pressure_darcy_explicit
  PUBLIC :: eval_pore_pressure_darcy_semi_implicit
  PUBLIC :: eval_pore_pressure_darcy_implicit
  PUBLIC :: solve_pore_pressure_darcy_analytic
  PUBLIC :: solve_frozen_darcy_mass_pressure
  PUBLIC :: eval_pore_pressure_dissipation_explicit
  PUBLIC :: eval_pore_pressure_dissipation_semi_implicit
  PUBLIC :: eval_pore_pressure_dissipation_implicit
  PUBLIC :: solve_pore_pressure_dissipation_analytic
  PUBLIC :: advance_pore_pressure_diffusion
  PUBLIC :: eval_pore_pressure_gas_loss
  PUBLIC :: eval_pore_pressure_inlet_source
  PUBLIC :: eval_pore_pressure_mass_exchange_source

  ! Uniform treatment surface for this source family.

CONTAINS

  SUBROUTINE eval_pore_pressure_dissipation_explicit(flow_depth, temperature, &
       mixture_density, solid_fractions, excess_pressure, rate)
    REAL(wp), INTENT(IN) :: flow_depth, temperature, mixture_density
    REAL(wp), INTENT(IN) :: solid_fractions(n_solid), excess_pressure
    REAL(wp), INTENT(OUT) :: rate
    REAL(wp) :: si_rate
    COMPLEX(wp) :: im_rate

    CALL eval_pore_pressure_dissipation_semi_implicit(si_rate)
    CALL eval_pore_pressure_dissipation_implicit(CMPLX(flow_depth,0.0_wp,wp), &
         CMPLX(temperature,0.0_wp,wp), CMPLX(mixture_density,0.0_wp,wp), &
         CMPLX(solid_fractions,0.0_wp,wp), CMPLX(excess_pressure,0.0_wp,wp), im_rate)
    rate = si_rate + REAL(im_rate,wp)
  END SUBROUTINE eval_pore_pressure_dissipation_explicit

  SUBROUTINE eval_pore_pressure_dissipation_semi_implicit(rate)
    REAL(wp), INTENT(OUT) :: rate
    rate = 0.0_wp
  END SUBROUTINE eval_pore_pressure_dissipation_semi_implicit

  SUBROUTINE eval_pore_pressure_dissipation_implicit(flow_depth, temperature,      &
       mixture_density, solid_fractions, excess_pressure, source)

    COMPLEX(wp), INTENT(IN) :: flow_depth
    COMPLEX(wp), INTENT(IN) :: temperature
    COMPLEX(wp), INTENT(IN) :: mixture_density
    COMPLEX(wp), INTENT(IN) :: solid_fractions(n_solid)
    COMPLEX(wp), INTENT(IN) :: excess_pressure
    COMPLEX(wp), INTENT(OUT) :: source

    COMPLEX(wp) :: diffusion_coefficient
    COMPLEX(wp) :: inhibition
    LOGICAL :: active
    REAL(wp), PARAMETER :: h_threshold = 1.0E-10_wp

    source = CMPLX(0.0_wp, 0.0_wp, wp)

    IF (.NOT. (REAL(flow_depth,wp) .GT. 0.0_wp)) RETURN
    IF (.NOT. (REAL(excess_pressure,wp) .GE. 0.0_wp)) RETURN

    CALL eval_pore_pressure_diffusion_factors(flow_depth, temperature,         &
         solid_fractions, diffusion_coefficient, inhibition, active)
    IF (.NOT. active) RETURN

    source = -mixture_density * (pi_g / 2.0_wp)**2 * diffusion_coefficient   &
         / MAX(h_threshold, flow_depth) * excess_pressure * inhibition

  END SUBROUTINE eval_pore_pressure_dissipation_implicit

  SUBROUTINE solve_pore_pressure_dissipation_analytic(excess_pressure, integration_step,   &
       flow_depth, temperature, solid_fractions, pressure_after)

    REAL(wp), INTENT(IN) :: excess_pressure, integration_step, flow_depth
    REAL(wp), INTENT(IN) :: temperature, solid_fractions(n_solid)
    REAL(wp), INTENT(OUT) :: pressure_after

    COMPLEX(wp) :: diffusion_coefficient, inhibition
    COMPLEX(wp) :: complex_solid_fractions(n_solid)
    REAL(wp) :: decay_rate, decay_exponent
    LOGICAL :: active
    REAL(wp), PARAMETER :: h_threshold = 1.0E-10_wp

    IF (.NOT. ieee_is_finite(excess_pressure) .OR.                           &
         .NOT. ieee_is_finite(integration_step) .OR.                         &
         .NOT. ieee_is_finite(flow_depth) .OR.                               &
         .NOT. ieee_is_finite(temperature) .OR.                              &
         .NOT. ALL(ieee_is_finite(solid_fractions))) THEN
       CALL fatal_error('Nonfinite input to pore-pressure diffusion integrator')
    END IF
    IF (integration_step .LT. 0.0_wp) THEN
       CALL fatal_error('Negative integration step in pore-pressure diffusion integrator')
    END IF

    pressure_after = excess_pressure
    IF (integration_step .EQ. 0.0_wp .OR. excess_pressure .LE. 0.0_wp .OR.   &
         flow_depth .LE. 0.0_wp) RETURN
    IF (REAL(1.0_wp-SUM(solid_fractions),wp) .LE. EPSILON(1.0_wp)) RETURN
    IF (dynamic_permeability_flag) THEN
       IF (REAL(SUM(solid_fractions),wp) .LE. EPSILON(1.0_wp)) RETURN
    END IF
    IF (temperature .LE. 0.0_wp) THEN
       CALL fatal_error('Nonpositive temperature in pore-pressure diffusion integrator')
    END IF

    complex_solid_fractions = CMPLX(solid_fractions, 0.0_wp, wp)
    CALL eval_pore_pressure_diffusion_factors(                               &
         CMPLX(flow_depth, 0.0_wp, wp), CMPLX(temperature, 0.0_wp, wp),      &
         complex_solid_fractions, diffusion_coefficient, inhibition, active)
    IF (.NOT. active) RETURN

    decay_rate = REAL((pi_g / 2.0_wp)**2 * diffusion_coefficient * inhibition,&
         wp) / MAX(h_threshold, flow_depth) / flow_depth
    IF (.NOT. ieee_is_finite(decay_rate) .OR. decay_rate .LT. 0.0_wp) THEN
       CALL fatal_error('Invalid decay rate in pore-pressure diffusion integrator')
    END IF
    IF (decay_rate .EQ. 0.0_wp) RETURN

    IF (integration_step .GT. LOG(HUGE(1.0_wp)) / decay_rate) THEN
       pressure_after = 0.0_wp
       RETURN
    END IF
    decay_exponent = integration_step * decay_rate
    pressure_after = excess_pressure * EXP(-decay_exponent)

  END SUBROUTINE solve_pore_pressure_dissipation_analytic


  !******************************************************************************
  !> \brief Evaluate the excess pore-pressure diffusion source.
  !> \details For wet cells with nonnegative excess pressure and active pore
  !>          volume, returns
  !>          \f$S_p=-\rho_m(\pi/2)^2D p_e I/\max(h_{\mathrm{threshold}},h)\f$.
  !>          D is the pressure-diffusion coefficient and I is the configured
  !>          solid-fraction inhibition factor. Complex arguments preserve
  !>          derivatives for the implicit Jacobian.
  !> \param[in] flow_depth Flow depth h [m].
  !> \param[in] temperature Gas/mixture temperature T [K].
  !> \param[in] mixture_density Mixture density rho_m [kg m^-3].
  !> \param[in] solid_fractions Per-class solid volume fractions.
  !> \param[in] excess_pressure Nonnegative excess pore pressure p_e [Pa].
  !> \param[out] source Conservative pore-pressure diffusion source S_p
  !>                   [kg^2 m^-3 s^-3].
  !> \note Returns zero for dry cells, negative excess pressure, absent pore volume, or
  !>       inactive diffusion closure.
  !******************************************************************************

  !******************************************************************************
  !> \brief Solve exact-flow excess pore pressure over a stage interval.
  !> \details With mixture density, depth, temperature and solid fractions
  !>          frozen, the transported variable
  !>          \f$Q_p=\rho_m h p_e\f$ satisfies
  !>          \f$\frac{dQ_p}{dt}=-\lambda Q_p\f$, where
  !>          \f$\lambda=(\pi/2)^2DI/[h\max(h_{\mathrm{threshold}},h)]\f$.
  !>          Therefore p_e is multiplied by
  !>          \f$\exp(-\lambda\,\texttt{integration\_step})\f$.
  !>          Dry states, nonpositive excess pressure, absent carrier phases and
  !>          inactive diffusion closures are returned unchanged.
  !> \param[in] excess_pressure Initial excess pore pressure p_e [Pa].
  !> \param[in] integration_step Effective stage interval, normally dt*a_ii [s].
  !> \param[in] flow_depth Frozen flow depth h [m].
  !> \param[in] temperature Frozen temperature T [K].
  !> \param[in] solid_fractions Frozen per-class solid volume fractions.
  !> \param[out] pressure_after Exact-flow excess pore pressure after the interval [Pa].
  !> \note Requires finite inputs and a nonnegative integration step. A negative step,
  !>       nonfinite input, or nonpositive temperature in an active state is fatal.
  !******************************************************************************

  !******************************************************************************
  !> \brief Backward-compatible name for the exact-flow pore-pressure solver.
  !> \param[in] excess_pressure Initial excess pore pressure [Pa].
  !> \param[in] integration_step Integration interval [s].
  !> \param[in] flow_depth Flow depth [m].
  !> \param[in] temperature Mixture temperature [K].
  !> \param[in] solid_fractions Solid volume fractions, one per class.
  !> \param[out] pressure_after Exact-flow excess pressure after the interval [Pa].
  !******************************************************************************
  SUBROUTINE advance_pore_pressure_diffusion(excess_pressure, integration_step,&
       flow_depth, temperature, solid_fractions, pressure_after)

    REAL(wp), INTENT(IN) :: excess_pressure, integration_step, flow_depth
    REAL(wp), INTENT(IN) :: temperature, solid_fractions(n_solid)
    REAL(wp), INTENT(OUT) :: pressure_after

    CALL solve_pore_pressure_dissipation_analytic(excess_pressure, integration_step,     &
         flow_depth, temperature, solid_fractions, pressure_after)

  END SUBROUTINE advance_pore_pressure_diffusion

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
  SUBROUTINE eval_pore_pressure_gas_loss(flow_depth, temperature,              &
       solid_fractions, total_solid_fraction, carrier_density, excess_pressure,&
       gas_loss)

    REAL(wp), INTENT(IN) :: flow_depth
    REAL(wp), INTENT(IN) :: temperature
    REAL(wp), INTENT(IN) :: solid_fractions(n_solid)
    REAL(wp), INTENT(IN) :: total_solid_fraction
    REAL(wp), INTENT(IN) :: carrier_density
    REAL(wp), INTENT(IN) :: excess_pressure
    REAL(wp), INTENT(OUT) :: gas_loss

    REAL(wp) :: inhibition
    REAL(wp) :: velocity_loss
    REAL(wp) :: sauter_diam
    REAL(wp) :: mean_sphericity
    REAL(wp) :: x
    REAL(wp) :: smoothstep
    REAL(wp) :: normalized_fraction
    REAL(wp) :: term
    REAL(wp) :: width
    REAL(wp) :: dynamic_transition
    REAL(wp) :: dynamic_viscosity
    REAL(wp) :: gas_density
    REAL(wp) :: local_permeability
    REAL(wp) :: local_kinematic_viscosity
    INTEGER :: n

    gas_loss = 0.0_wp
    inhibition = 1.0_wp
    local_permeability = hydraulic_permeability
    local_kinematic_viscosity = kin_visc_c

    SELECT CASE (TRIM(ADJUSTL(f_inhibit_mode)))
    CASE ('OFF')
       inhibition = 1.0_wp

    CASE ('STEP')
       IF (total_solid_fraction .LT. maximum_solid_packing) THEN
          inhibition = 1.0_wp
       ELSE
          inhibition = 0.0_wp
       END IF

    CASE ('STATIC')
       IF (N_inh .LT. 0) THEN
          WRITE(*,*) 'ERROR: N_inh must be >= 0. N_inh=', N_inh
          STOP
       END IF
       x = total_solid_fraction
       width = maximum_solid_packing - alpha_trans
       normalized_fraction = (x - alpha_trans) / width
       IF (normalized_fraction .LT. 0.0_wp) normalized_fraction = 0.0_wp
       IF (normalized_fraction .GT. 1.0_wp) normalized_fraction = 1.0_wp
       smoothstep = 0.0_wp
       DO n = 0, N_inh
          term = pascal_coeff(n) * normalized_fraction**(N_inh + n + 1.0_wp)
          smoothstep = smoothstep + term
       END DO
       inhibition = 1.0_wp - smoothstep

    CASE ('DYNAMIC')
       dynamic_transition = 10.0_wp**(-0.28_wp) *                           &
            MAX(flow_depth, EPSILON(1.0_wp))**0.12_wp
       IF (dynamic_transition .GE. maximum_solid_packing) THEN
          dynamic_transition = 0.99_wp * maximum_solid_packing
       END IF
       IF (N_inh .LT. 0) THEN
          WRITE(*,*) 'ERROR: N_inh must be >= 0. N_inh=', N_inh
          STOP
       END IF
       x = total_solid_fraction
       width = maximum_solid_packing - dynamic_transition
       normalized_fraction = (x - dynamic_transition) / width
       IF (normalized_fraction .LT. 0.0_wp) normalized_fraction = 0.0_wp
       IF (normalized_fraction .GT. 1.0_wp) normalized_fraction = 1.0_wp
       smoothstep = 0.0_wp
       DO n = 0, N_inh
          term = pascal_coeff(n) * normalized_fraction**(N_inh + n + 1.0_wp)
          smoothstep = smoothstep + term
       END DO
       inhibition = 1.0_wp - smoothstep
    END SELECT

    IF (dynamic_permeability_flag .AND.                                     &
         (total_solid_fraction .GT. EPSILON(1.0_wp))) THEN
       IF (n_solid .GT. 1) THEN
          sauter_diam = sauter_diameter(solid_fractions)
          mean_sphericity = DOT_PRODUCT(solid_fractions / diam_s, sphericity_s) &
               / SUM(solid_fractions / diam_s)
       ELSE
          sauter_diam = diam_s(1)
          mean_sphericity = sphericity_s(1)
       END IF
       local_permeability = ((1.0_wp-total_solid_fraction)**3.0_wp *         &
            (sauter_diam * mean_sphericity)**2.0_wp) /                       &
            (150.0_wp * total_solid_fraction**2.0_wp)
    END IF

    IF (gas_flag .AND. sutherland_flag) THEN
       dynamic_viscosity = muRef_Suth * (temperature / Tref_Suth)**1.5_wp *  &
            (Tref_Suth + S_mu) / (temperature + S_mu)
       gas_density = pres / (sp_gas_const_a * temperature)
       local_kinematic_viscosity = dynamic_viscosity / gas_density
    END IF

    velocity_loss = local_permeability / (local_kinematic_viscosity *         &
         carrier_density) / MAX(1.0e-5_wp,flow_depth) * 0.5_wp * pi_g *       &
         excess_pressure
    gas_loss = velocity_loss * inhibition

  END SUBROUTINE eval_pore_pressure_gas_loss

  !******************************************************************************
  !> \brief Compute the excess pore-pressure source associated with inlet injection.
  !> \details Returns the time-weighted source
  !>          \f$S_p=c_t\dot{H}\rho_m(q_1f_pg'+p_e)\f$,
  !>          where \f$q_1=\rho_mh\f$,
  !>          combining the hydrostatic pressure carried by injected mixture
  !>          with the injected excess-pressure fraction.
  !> \param[in] time_coefficient Dimensionless temporal activation or step-averaged weight.
  !> \param[in] q1 Local mass-per-area/carrier state quantity used by the inlet closure.
  !> \param[in] pressure_fraction Configured fraction f_p of the injected hydrostatic pressure.
  !> \param[in] injection_rate Equivalent injected layer-thickness rate Hdot [m s^-1].
  !> \param[in] mixture_density Density of the injected mixture rho_m [kg m^-3].
  !> \param[in] reduced_gravity Buoyancy-reduced gravitational acceleration g' [m s^-2].
  !> \param[in] excess_pressure Local excess pore pressure p_e [Pa].
  !> \return Conservative excess-pore-pressure source [kg^2 m^-3 s^-3].
  !******************************************************************************
  PURE FUNCTION eval_pore_pressure_inlet_source(time_coefficient, q1,          &
       pressure_fraction, injection_rate, mixture_density, reduced_gravity,    &
       excess_pressure) RESULT(source)

    REAL(wp), INTENT(IN) :: time_coefficient
    REAL(wp), INTENT(IN) :: q1
    REAL(wp), INTENT(IN) :: pressure_fraction
    REAL(wp), INTENT(IN) :: injection_rate
    REAL(wp), INTENT(IN) :: mixture_density
    REAL(wp), INTENT(IN) :: reduced_gravity
    REAL(wp), INTENT(IN) :: excess_pressure
    REAL(wp) :: source

    source = time_coefficient * (q1 * pressure_fraction * injection_rate *    &
         mixture_density * reduced_gravity + excess_pressure * injection_rate &
         * mixture_density)

  END FUNCTION eval_pore_pressure_inlet_source

  !******************************************************************************
  !> \brief Compute the excess-pressure source induced by a mass exchange.
  !> \details Applies the local proportional closure S_p = S_m*p_e.
  !> \param[in] mass_source Conservative mass source S_m for the exchanged carrier phase.
  !> \param[in] excess_pressure Local excess pore pressure p_e [Pa].
  !> \return Excess-pore-pressure equation source [kg^2 m^-3 s^-3].
  !******************************************************************************
  PURE FUNCTION eval_pore_pressure_mass_exchange_source(mass_source,           &
       excess_pressure) RESULT(source)

    REAL(wp), INTENT(IN) :: mass_source
    REAL(wp), INTENT(IN) :: excess_pressure
    REAL(wp) :: source

    source = mass_source * excess_pressure

  END FUNCTION eval_pore_pressure_mass_exchange_source

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
       solid_fractions, diffusion_coefficient, inhibition, active)

    COMPLEX(wp), INTENT(IN) :: flow_depth, temperature
    COMPLEX(wp), INTENT(IN) :: solid_fractions(n_solid)
    COMPLEX(wp), INTENT(OUT) :: diffusion_coefficient, inhibition
    LOGICAL, INTENT(OUT) :: active

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

    diffusion_coefficient = local_permeability /                            &
         (porosity * local_kinematic_viscosity * gas_density *               &
         gas_compressibility)

    SELECT CASE (TRIM(ADJUSTL(f_inhibit_mode)))
    CASE ('OFF')
       inhibition = 1.0_wp

    CASE ('STEP')
       IF (SUM(solid_fractions) .LT. maximum_solid_packing) THEN
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


  ! Valid zero contribution for the semi-implicit part of Newton.

  SUBROUTINE solve_frozen_darcy_mass_pressure(mass_before, pressure_before, &
       interval, dissipation_rate, loss_coefficient, gas_loss_enabled, &
       mass_after, pressure_after, conservative_pressure_after, &
       gas_mass_lost, admissible)
    REAL(wp), INTENT(IN) :: mass_before, pressure_before, interval
    REAL(wp), INTENT(IN) :: dissipation_rate, loss_coefficient
    LOGICAL, INTENT(IN) :: gas_loss_enabled
    REAL(wp), INTENT(OUT) :: mass_after, pressure_after
    REAL(wp), INTENT(OUT) :: conservative_pressure_after, gas_mass_lost
    LOGICAL, INTENT(OUT) :: admissible
    REAL(wp) :: exponential_factor, pressure_integral, coefficient

    mass_after=mass_before
    pressure_after=pressure_before
    conservative_pressure_after=mass_before*pressure_before
    gas_mass_lost=0.0_wp
    admissible=.FALSE.

    IF (.NOT.ieee_is_finite(mass_before) .OR. &
        .NOT.ieee_is_finite(pressure_before) .OR. &
        .NOT.ieee_is_finite(interval) .OR. &
        .NOT.ieee_is_finite(dissipation_rate) .OR. &
        .NOT.ieee_is_finite(loss_coefficient)) RETURN
    IF (mass_before<=0.0_wp .OR. interval<0.0_wp .OR. &
        dissipation_rate<0.0_wp .OR. loss_coefficient<0.0_wp) RETURN

    ! The existing v2/dev diffusion closure does not dissipate
    ! negative excess pore pressure; preserve that guard.
    IF (pressure_before<=0.0_wp .OR. interval==0.0_wp) THEN
      admissible=.TRUE.
      RETURN
    END IF

    exponential_factor=EXP(-dissipation_rate*interval)
    pressure_after=pressure_before*exponential_factor

    coefficient=0.0_wp
    IF (gas_loss_enabled) coefficient=loss_coefficient
    IF (dissipation_rate>0.0_wp) THEN
      ! Stable evaluation of (1-exp(-lambda*dt))/lambda,
      ! including small lambda*dt without cancellation.
      IF (ABS(dissipation_rate*interval)<1.0e-5_wp) THEN
        pressure_integral=interval*(1.0_wp-0.5_wp*dissipation_rate*interval+ &
             (dissipation_rate*interval)**2/6.0_wp- &
             (dissipation_rate*interval)**3/24.0_wp)
      ELSE
        pressure_integral=(1.0_wp-exponential_factor)/dissipation_rate
      END IF
    ELSE
      pressure_integral=interval
    END IF

    gas_mass_lost=coefficient*pressure_before*pressure_integral
    mass_after=mass_before-gas_mass_lost
    admissible=(mass_after>0.0_wp .AND. ieee_is_finite(mass_after))
    IF (.NOT.admissible) RETURN
    conservative_pressure_after=mass_after*pressure_after
  END SUBROUTINE solve_frozen_darcy_mass_pressure


  !> Conservative withdrawal using properties frozen at the start of the
  !> analytic Darcy substep. The lost material is carrier gas alone.
  !> gas_mass_fractions are normalized WITHIN the gas phase, not the mixture.
  !> The caller must supply the energy per unit mass of removed gas in the
  !> convention of the conservative energy equation (including kinetic energy
  !> if that equation stores it).
  SUBROUTINE apply_frozen_darcy_withdrawal(mass_lost, u, v, gas_energy, &
       gas_mass_fractions, momentum_x, momentum_y, energy, gas_masses, admissible)
    REAL(wp), INTENT(IN) :: mass_lost, u, v, gas_energy
    REAL(wp), INTENT(IN) :: gas_mass_fractions(:)
    REAL(wp), INTENT(INOUT) :: momentum_x, momentum_y, energy
    REAL(wp), INTENT(INOUT) :: gas_masses(:)
    LOGICAL, INTENT(OUT) :: admissible
    REAL(wp) :: withdrawals(SIZE(gas_masses)), tolerance

    admissible=.FALSE.
    IF (SIZE(gas_mass_fractions)/=SIZE(gas_masses)) RETURN
    IF (.NOT.ieee_is_finite(mass_lost) .OR. mass_lost<0.0_wp) RETURN
    IF (.NOT.ieee_is_finite(u) .OR. .NOT.ieee_is_finite(v) .OR. &
        .NOT.ieee_is_finite(gas_energy)) RETURN
    IF (.NOT.ALL(ieee_is_finite(gas_masses)) .OR. &
        .NOT.ALL(ieee_is_finite(gas_mass_fractions))) RETURN
    IF (ANY(gas_mass_fractions<0.0_wp)) RETURN
    IF (mass_lost>0.0_wp) THEN
      tolerance=1.0e-10_wp
      IF (ABS(SUM(gas_mass_fractions)-1.0_wp)>tolerance) RETURN
    END IF
    withdrawals=gas_mass_fractions*mass_lost
    IF (ANY(withdrawals>gas_masses+ &
         1.0e-12_wp*MAX(1.0_wp,ABS(gas_masses)))) RETURN
    ! Perform all feasibility checks before touching any conservative variable.
    momentum_x=momentum_x-u*mass_lost
    momentum_y=momentum_y-v*mass_lost
    energy=energy-gas_energy*mass_lost
    gas_masses=MAX(0.0_wp,gas_masses-withdrawals)
    admissible=.TRUE.
  END SUBROUTINE apply_frozen_darcy_withdrawal


  !> The model returns a single conservative source vector, including
  !> diffusion and (when enabled) gas withdrawal. Coordination of shared
  !> carrier limits belongs to the source executor, not this model.
  SUBROUTINE eval_pore_pressure_darcy_explicit(state,rate)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(OUT) :: rate(n_eqns)
    COMPLEX(wp) :: pressure_rate
    REAL(wp) :: loss_volume, rho_c, cp_c, Ri, rhom, redgrav, cpm
    REAL(wp) :: mass_rate, gas_mass_fraction
    TYPE(complex_physical_state_type) :: cs
    rate=0.0_wp
    IF (.NOT.pore_pressure_flag) RETURN
    CALL complex_state_from_real(state,cs)
    CALL eval_pore_pressure_dissipation_implicit(cs%h,cs%T,cs%rho_m, &
         cs%alphas,cs%exc_pore_pres,pressure_rate)
    rate(idx_poreEqn)=REAL(pressure_rate,wp)
    IF (.NOT.(gas_flag .AND. gas_loss_flag)) RETURN
    IF (SUM(state%alphas)>=maximum_solid_packing) RETURN
    CALL mixt_var_state(state,Ri,rhom,rho_c,redgrav,cp_c,cpm)
    loss_volume=0.0_wp
    CALL eval_pore_pressure_gas_loss(state%h,state%T,state%alphas, &
         SUM(state%alphas),rho_c,state%exc_pore_pres,loss_volume)
    mass_rate=-rho_c*MAX(0.0_wp,loss_volume)
    rate(idx_totMassEqn)=mass_rate
    rate(idx_uEqn)=mass_rate*state%u
    rate(idx_vEqn)=mass_rate*state%v
    rate(idx_engyEqn)=mass_rate*state%T*cp_c
    IF (n_add_gas>0) THEN
      gas_mass_fraction=MAX(0.0_wp,1.0_wp-SUM(state%xs))
      IF (gas_mass_fraction>0.0_wp) &
        rate(idx_addGasEqn_first:idx_addGasEqn_first+n_add_gas-1)= &
          mass_rate*state%xg/gas_mass_fraction
    END IF
    rate(idx_poreEqn)=rate(idx_poreEqn)+mass_rate*state%exc_pore_pres
  END SUBROUTINE eval_pore_pressure_darcy_explicit

  SUBROUTINE eval_pore_pressure_darcy_semi_implicit(state,rate)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(OUT) :: rate(n_eqns)
    rate=0.0_wp
  END SUBROUTINE eval_pore_pressure_darcy_semi_implicit

  SUBROUTINE eval_pore_pressure_darcy_implicit(state,rate)
    TYPE(complex_physical_state_type), INTENT(IN) :: state
    COMPLEX(wp), INTENT(OUT) :: rate(n_eqns)
    rate=CMPLX(0.0_wp,0.0_wp,wp)
    IF (.NOT.pore_pressure_flag) RETURN
    CALL eval_pore_pressure_dissipation_implicit(state%h,state%T, &
         state%rho_m,state%alphas,state%exc_pore_pres,rate(idx_poreEqn))
  END SUBROUTINE eval_pore_pressure_darcy_implicit

  !> The full coupled analytic map is a separate existing kernel. Until
  !> the source registry supports multi-equation pressure groups, this
  !> state-level adapter preserves the original pressure-only analytic step.
  SUBROUTINE solve_pore_pressure_darcy_analytic(state,interval,pressure_after)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: interval
    REAL(wp), INTENT(OUT) :: pressure_after
    pressure_after=state%exc_pore_pres
    IF (.NOT.pore_pressure_flag) RETURN
    CALL solve_pore_pressure_dissipation_analytic(state%exc_pore_pres, &
         interval,state%h,state%T,state%alphas,pressure_after)
  END SUBROUTINE solve_pore_pressure_darcy_analytic

END MODULE pore_pressure_darcy_2d
