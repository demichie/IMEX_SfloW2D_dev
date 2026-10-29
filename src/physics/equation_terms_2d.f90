!********************************************************************************
!> \brief Equation partitioning, inertial transport and hydrostatic flux operators.
!> \details Configurable source assembly and injection window coefficients now
!> reside in source_terms_2d and mass_injection_coefficients_2d, respectively.
!********************************************************************************

MODULE equation_terms_2d

  USE diagnostics_2d, ONLY : fatal_error
  USE constitutive_parameters_2d
  USE state_conversion_2d

  USE equation_metadata_2d, ONLY : equation_partition_type
  USE parameters_2d, ONLY : wp, &
       n_eqns, &
       n_vars, &
       thermal_cooling_flag, &
       liquid_flag, &
       gas_flag, &
       pore_pressure_flag, &
       idx_pore, &
       idx_u, &
       idx_v, &
       idx_solidEqn_first, &
       idx_solidEqn_last, &
       idx_addGasEqn_first, &
       idx_addGasEqn_last

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_local_speeds_direction
  PUBLIC :: eval_inertial_flux, eval_hydrostatic_coefficient
  PUBLIC :: limit_component_mass_flux

CONTAINS

  !******************************************************************************
  !> \brief Evaluate conservative transport without the hydrostatic pressure term.
  !>
  !> Every retained conservative component is transported by the normal face
  !> velocity.  Hydrostatic pressure is represented exclusively by the path
  !> contribution and must not appear in this flux.
  !> The flux is \f$\mathbf{F}_n=u_n\mathbf{q}_c\f$ for h>0 and is zero for a
  !> dry state.
  !>
  !> \param[in] qcj Local conservative state vector associated with qpj.
  !> \param[in] qpj Local physical state of length n_vars+2 with component mass fractions.
  !> \param[in] dir Cartesian flux direction: 1 for x, 2 for y.
  !> \param[out] flux Conservative flux vector, one entry per balance equation.
  !******************************************************************************

  SUBROUTINE eval_inertial_flux(qcj,qpj,dir,flux)

    REAL(wp), INTENT(IN) :: qcj(n_vars)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    INTEGER, INTENT(IN) :: dir
    REAL(wp), INTENT(OUT) :: flux(n_eqns)

    REAL(wp) :: normal_velocity

    flux = 0.0_wp
    IF ( qpj(1) .LE. EPSILON(1.0_wp) ) RETURN

    SELECT CASE ( dir )
    CASE ( 1 )
       normal_velocity = qpj(idx_u)
    CASE ( 2 )
       normal_velocity = qpj(idx_v)
    CASE DEFAULT
       ERROR STOP 'eval_inertial_flux: invalid direction'
    END SELECT

    flux = normal_velocity * qcj(1:n_eqns)

  END SUBROUTINE eval_inertial_flux

  !******************************************************************************
  !> \brief Derive reduced gravity and the pressure coefficient Gamma from a physical state.
  !>
  !> \f$\Gamma=\rho_mg'\f$ is recomputed from the final reconstructed primitive state;
  !> no derived thermodynamic quantity is reconstructed independently.
  !>
  !> \param[in] qpj Local physical state of length n_vars+2 with component mass fractions.
  !> \param[out] reduced_gravity Buoyancy-reduced gravitational acceleration g prime [m s^-2].
  !> \param[out] gamma Hydrostatic pressure coefficient \f$\Gamma=\rho_mg'\f$
  !>                   [kg m^-2 s^-2].
  !******************************************************************************

  SUBROUTINE eval_hydrostatic_coefficient(qpj,reduced_gravity,gamma)

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(OUT) :: reduced_gravity, gamma

    REAL(wp) :: Richardson, rho_m, rho_c
    REAL(wp) :: sp_heat_c, sp_heat_mix

    CALL mixt_var(qpj,Richardson,rho_m,rho_c,reduced_gravity,sp_heat_c,       &
         sp_heat_mix)
    gamma = rho_m * reduced_gravity

  END SUBROUTINE eval_hydrostatic_coefficient


  !> \brief Shared characteristic-speed calculation for a Cartesian direction.
  !> \details Uses the shallow-flow eigenvalue bounds
  !>          \f$\lambda_{\min,\max}=u_n\mp\sqrt{Gg'h}\f$. If the radicand is
  !>          negative, both bounds are set to u_n; a zero radicand also gives
  !>          equal bounds.
  !> \param[in] qpj Local physical state of length n_vars+2: depth, temperature, component
  !>                fractions and appended x/y velocities.
  !> \param[in] grav_coeff Large-slope gravity factor G at the state location.
  !> \param[in] normal_velocity_index Index of the normal velocity in qpj (idx_u for x,
  !>                                  idx_v for y).
  !> \param[out] vel_min Lower characteristic speed [m s^-1].
  !> \param[out] vel_max Upper characteristic speed [m s^-1].
  !> \note normal_velocity_index must identify one of the two Cartesian velocity components.
  SUBROUTINE eval_local_speeds_direction(qpj,grav_coeff,normal_velocity_index, &
       vel_min,vel_max)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2), grav_coeff
    INTEGER, INTENT(IN) :: normal_velocity_index
    REAL(wp), INTENT(OUT) :: vel_min, vel_max
    REAL(wp) :: r_Ri, r_rho_m, r_rho_c, r_red_grav
    REAL(wp) :: r_sp_heat_c, r_sp_heat_mix, r_h, normal_velocity, celerity

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,r_sp_heat_c,          &
         r_sp_heat_mix)
    r_h = qpj(1)
    normal_velocity = qpj(normal_velocity_index)

    IF ( r_red_grav * r_h .LT. 0.0_wp ) THEN
       vel_min = normal_velocity
       vel_max = normal_velocity
    ELSE
       celerity = SQRT( grav_coeff * r_red_grav * r_h )
       vel_min = normal_velocity - celerity
       vel_max = normal_velocity + celerity
    END IF
  END SUBROUTINE eval_local_speeds_direction

  !> \brief Bound component-mass transport by the signed total-mixture mass flux.
  !>
  !> Bound the signed sum of all independently transported component-mass
  !> fluxes by the signed total-mixture mass flux. A single common factor
  !> preserves the relative component proportions.
  !> If the signed component sum exceeds the total flux in its direction, each
  !> component flux is multiplied by |F_total/F_components|; all other flux
  !> components are left unchanged.
  !>
  !> \param[in,out] flux Conservative flux vector [equation]; its solid, additional-gas and,
  !>                  when both gas and liquid phases are enabled, liquid components may be
  !>                  scaled in place.

  SUBROUTINE limit_component_mass_flux(flux)

    REAL(wp), INTENT(INOUT) :: flux(n_eqns)
    REAL(wp) :: component_flux, scale

    component_flux = SUM(flux(idx_solidEqn_first:idx_solidEqn_last))          &
         + SUM(flux(idx_addGasEqn_first:idx_addGasEqn_last))

    IF ((flux(1) .GT. 0.0_wp .AND. component_flux .GT. flux(1)) .OR.          &
        (flux(1) .LT. 0.0_wp .AND. component_flux .LT. flux(1))) THEN
       scale = flux(1) / component_flux
       flux(idx_solidEqn_first:idx_solidEqn_last) =                           &
            scale * flux(idx_solidEqn_first:idx_solidEqn_last)
       flux(idx_addGasEqn_first:idx_addGasEqn_last) =                         &
            scale * flux(idx_addGasEqn_first:idx_addGasEqn_last)
    END IF

  END SUBROUTINE limit_component_mass_flux

END MODULE equation_terms_2d
