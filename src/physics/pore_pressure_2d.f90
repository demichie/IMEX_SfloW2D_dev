!********************************************************************************
!> \brief Excess pore-pressure source terms and their local closures.
!********************************************************************************

MODULE pore_pressure_2d

  USE parameters_2d, ONLY : wp, n_solid, gas_flag, sutherland_flag
  USE COMPLEXIFY
  USE diagnostics_2d, ONLY : fatal_error
  USE constitutive_parameters_2d, ONLY : alpha_trans, N_inh, f_inhibit_mode,    &
       pascal_coeff, hydraulic_permeability, kin_visc_c, sp_heat_a,             &
       sp_gas_const_a, pres, diam_s, sphericity_s, maximum_solid_packing,       &
       dynamic_permeability_flag, muRef_Suth, Tref_Suth, S_mu
  USE geometry_2d, ONLY : pi_g
  USE state_conversion_2d, ONLY : sauter_diameter
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_pore_pressure_diffusion_source
  PUBLIC :: advance_pore_pressure_diffusion
  PUBLIC :: eval_pore_pressure_gas_loss
  PUBLIC :: eval_pore_pressure_inlet_source
  PUBLIC :: eval_pore_pressure_mass_exchange_source

CONTAINS

  SUBROUTINE eval_pore_pressure_diffusion_source(flow_depth, temperature,      &
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

  END SUBROUTINE eval_pore_pressure_diffusion_source

  !******************************************************************************
  !> \brief Advance excess pore pressure analytically over one DIRK integration step.
  !> \details With mixture density, depth, temperature and solid fractions
  !>          frozen, the transported variable Q_p = rho_m*h*p_e satisfies
  !>          dQ_p/dt = -lambda*Q_p,
  !>          where lambda = (pi/2)^2*D*I/(h*max(h_threshold,h)). Therefore the
  !>          excess pressure is multiplied by exp(-lambda*integration_step).
  !>          Existing dry, negative-pressure and no-carrier states are unchanged.
  !> \param integration_step Effective DIRK integration step (dt * a_ii).
  !******************************************************************************
  SUBROUTINE advance_pore_pressure_diffusion(excess_pressure, integration_step,&
       flow_depth, temperature, solid_fractions, pressure_after)

    REAL(wp), INTENT(IN) :: excess_pressure, integration_step, flow_depth
    REAL(wp), INTENT(IN) :: temperature, solid_fractions(n_solid)
    REAL(wp), INTENT(OUT) :: pressure_after

    COMPLEX(wp) :: diffusion_coefficient, inhibition
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

    CALL eval_pore_pressure_diffusion_factors(                               &
         CMPLX(flow_depth, 0.0_wp, wp), CMPLX(temperature, 0.0_wp, wp),      &
         CMPLX(solid_fractions, 0.0_wp, wp), diffusion_coefficient,          &
         inhibition, active)
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

  END SUBROUTINE advance_pore_pressure_diffusion

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

  PURE FUNCTION eval_pore_pressure_mass_exchange_source(mass_source,           &
       excess_pressure) RESULT(source)

    REAL(wp), INTENT(IN) :: mass_source
    REAL(wp), INTENT(IN) :: excess_pressure
    REAL(wp) :: source

    source = mass_source * excess_pressure

  END FUNCTION eval_pore_pressure_mass_exchange_source

END MODULE pore_pressure_2d
