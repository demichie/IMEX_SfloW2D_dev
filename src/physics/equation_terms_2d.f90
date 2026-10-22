!********************************************************************************
!> \brief Physical fluxes and source terms for the depth-averaged equations.
!> \details Provides the equation-partition metadata, inertial transport flux,
!>          characteristic wave-speed bounds, hydrostatic coefficients, bed
!>          curvature sources, and explicit, implicit, semi-implicit and
!>          mass-exchange source vectors. It also evaluates inlet states and
!>          exact step-averaged activation coefficients for periodic, interval
!>          and fissure sources. Rheology, pore-pressure and thermal laws are
!>          delegated to their owning physics modules. Spatial path assembly
!>          and time integration remain in their dedicated modules.
!>          The module retains private per-step source-window coefficients
!>          between set_source_time_coefficients and eval_expl_terms.
!********************************************************************************

MODULE equation_terms_2d

  USE diagnostics_2d, ONLY : fatal_error

  USE constitutive_parameters_2d
  USE state_conversion_2d
  USE rheology_froude_2d, ONLY : eval_froude_source
  USE rheology_lucas_2d, ONLY : eval_lucas_source
  USE rheology_mui_2d, ONLY : eval_mui_source
  USE rheology_mui_voellmy_2d, ONLY : eval_mui_voellmy_source,               &
       eval_mui_voellmy_implicit_source
  USE pore_pressure_2d, ONLY : eval_pore_pressure_diffusion_source,            &
       eval_pore_pressure_gas_loss, eval_pore_pressure_inlet_source,            &
       eval_pore_pressure_mass_exchange_source
  USE thermal_sources_2d, ONLY : eval_inlet_thermal_energy_source,              &
       eval_lava_cooling_source, eval_mass_exchange_thermal_energy,              &
       eval_vaporization_sources
  USE rheology_voellmy_2d, ONLY : eval_voellmy_implicit_source,              &
       eval_voellmy_semi_impl_source
  USE rheology_plastic_2d, ONLY : eval_plastic_implicit_source
  USE rheology_lava_2d, ONLY : eval_lava_implicit_source,                    &
       eval_lava_semi_impl_source
  USE rheology_lahar_2d, ONLY : eval_lahar_implicit_source,                  &
       eval_lahar_semi_impl_source, eval_lahar_fluid_viscosity,              &
       eval_lahar_deposition_rate
  USE rheology_empirical_2d, ONLY : eval_empirical_implicit_source
  USE rheology_quadratic_2d, ONLY : eval_quadratic_implicit_source

  USE equation_metadata_2d, ONLY : equation_partition_type,                   &
       IM_GROUP_METHOD_ANALYTIC, IM_GROUP_METHOD_NEWTON

  USE parameters_2d, ONLY : wp, dry_thickness_tolerance
  USE parameters_2d, ONLY : n_eqns, n_vars, n_solid, n_add_gas,                &
       n_stoch_vars, n_pore_vars
  USE parameters_2d, ONLY : rheology_flag, rheology_model_id,                   &
       rheology_id_voellmy, rheology_id_plastic, rheology_id_lava, rheology_id_lahar,       &
       rheology_id_empirical, rheology_id_quadratic, rheology_id_froude, rheology_id_lucas, &
       rheology_id_mui, rheology_id_mui_voellmy,                                     &
       liquid_flag, gas_flag, slope_correction_flag,                            &
       curvature_term_flag, stochastic_flag,                  &
       stoch_transport_flag, pore_pressure_flag, sutherland_flag

  USE parameters_2d, ONLY : idx_h, idx_hu, idx_hv, idx_T, idx_solid_first,      &
       idx_solid_last, idx_add_gas_first, idx_add_gas_last, idx_stoch, idx_pore,  &
       idx_u, idx_v

  USE parameters_2d, ONLY : idx_totMassEqn, idx_uEqn, idx_vEqn, idx_engyEqn,    &
       idx_solidEqn_first, idx_solidEqn_last, idx_addGasEqn_first,              &
       idx_addGasEqn_last, idx_stochEqn, idx_poreEqn

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_problem_param
  PUBLIC :: eval_local_speeds_direction
  PUBLIC :: eval_inertial_flux, eval_hydrostatic_coefficient
  PUBLIC :: limit_component_mass_flux
  PUBLIC :: curvature_acceleration, eval_curvature_momentum_source
  PUBLIC :: eval_expl_terms
  PUBLIC :: eval_implicit_terms, eval_nh_semi_impl_terms
  PUBLIC :: eval_mass_exchange_terms, eval_source_bdry
  PUBLIC :: set_source_time_coefficients

  ! Time coefficients of the sources averaged over the current step, set by
  ! set_source_time_coefficients once the step length is known. A negative
  ! value means they have not been set, and eval_expl_terms uses the value
  ! at its time argument.
  REAL(wp) :: bottom_coeff_step = -1.0_wp
  REAL(wp) :: bottom_vel_step = -1.0_wp
  REAL(wp) :: lateral_coeff_step = -1.0_wp
  REAL(wp), ALLOCATABLE :: fissure_coeff_step(:)

CONTAINS

  !******************************************************************************
  !> \brief Contract the filtered bed Hessian with the Cartesian velocity.
  !>
  !> The Cartesian velocity components are contracted with the cell-centred
  !> Hessian obtained from the retained quadratic least-squares bed fit.
  !>
  !> \param[in] Bsecond_xx Filtered bed Hessian xx entry [m^-1].
  !> \param[in] Bsecond_xy Filtered bed Hessian xy entry [m^-1].
  !> \param[in] Bsecond_yy Filtered bed Hessian yy entry [m^-1].
  !> \param[in] velocity_x Cartesian depth-averaged x velocity [m s^-1].
  !> \param[in] velocity_y Cartesian depth-averaged y velocity [m s^-1].
  !> \return Signed bed-curvature acceleration
  !>         \f$a_{\mathrm{curv}}=B_{xx}u^2+2B_{xy}uv+B_{yy}v^2\f$ [m s^-2].
  !******************************************************************************

  PURE FUNCTION curvature_acceleration( Bsecond_xx, Bsecond_xy, Bsecond_yy,   &
       velocity_x, velocity_y ) RESULT(acceleration)

    REAL(wp), INTENT(IN) :: Bsecond_xx, Bsecond_xy, Bsecond_yy
    REAL(wp), INTENT(IN) :: velocity_x, velocity_y
    REAL(wp) :: acceleration

    acceleration = Bsecond_xx * velocity_x**2                               &
         + 2.0_wp * Bsecond_xy * velocity_x * velocity_y                    &
         + Bsecond_yy * velocity_y**2

  END FUNCTION curvature_acceleration

  !******************************************************************************
  !> \brief Evaluate the bed-curvature contribution to both momentum equations.
  !>
  !> Implements \f$\mathbf{S}_{\mathrm{curv}}=-G_cM
  !> a_{\mathrm{curv}}\nabla B\f$, independently of the hydrostatic path.
  !>
  !> \param[in] Bprime_x Filtered bed derivative dB/dx [dimensionless].
  !> \param[in] Bprime_y Filtered bed derivative dB/dy [dimensionless].
  !> \param[in] Bsecond_xx Filtered bed Hessian xx entry [m^-1].
  !> \param[in] Bsecond_xy Filtered bed Hessian xy entry [m^-1].
  !> \param[in] Bsecond_yy Filtered bed Hessian yy entry [m^-1].
  !> \param[in] grav_coefficient Large-slope gravity factor G, not the gravitational acceleration.
  !> \param[in] mixture_mass Mixture mass per unit horizontal area, rho_m*h [kg m^-2].
  !> \param[in] velocity_x Cartesian depth-averaged x velocity [m s^-1].
  !> \param[in] velocity_y Cartesian depth-averaged y velocity [m s^-1].
  !> \param[out] source_x Curvature contribution to the x-momentum source.
  !> \param[out] source_y Curvature contribution to the y-momentum source.
  !******************************************************************************

  PURE SUBROUTINE eval_curvature_momentum_source( Bprime_x, Bprime_y,         &
       Bsecond_xx, Bsecond_xy, Bsecond_yy, grav_coefficient,                 &
       mixture_mass, velocity_x, velocity_y, source_x, source_y )

    REAL(wp), INTENT(IN) :: Bprime_x, Bprime_y
    REAL(wp), INTENT(IN) :: Bsecond_xx, Bsecond_xy, Bsecond_yy
    REAL(wp), INTENT(IN) :: grav_coefficient
    REAL(wp), INTENT(IN) :: mixture_mass, velocity_x, velocity_y
    REAL(wp), INTENT(OUT) :: source_x, source_y
    REAL(wp) :: acceleration

    acceleration = curvature_acceleration( Bsecond_xx, Bsecond_xy,          &
         Bsecond_yy, velocity_x, velocity_y )
    source_x = -grav_coefficient * mixture_mass * acceleration * Bprime_x
    source_y = -grav_coefficient * mixture_mass * acceleration * Bprime_y

  END SUBROUTINE eval_curvature_momentum_source

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
    IF (gas_flag .AND. liquid_flag) component_flux = component_flux           &
         + flux(n_vars)

    IF ((flux(1) .GT. 0.0_wp .AND. component_flux .GT. flux(1)) .OR.          &
        (flux(1) .LT. 0.0_wp .AND. component_flux .LT. flux(1))) THEN
       scale = flux(1) / component_flux
       flux(idx_solidEqn_first:idx_solidEqn_last) =                           &
            scale * flux(idx_solidEqn_first:idx_solidEqn_last)
       flux(idx_addGasEqn_first:idx_addGasEqn_last) =                         &
            scale * flux(idx_addGasEqn_first:idx_addGasEqn_last)
       IF (gas_flag .AND. liquid_flag) flux(n_vars) = scale * flux(n_vars)
    END IF

  END SUBROUTINE limit_component_mass_flux

  !******************************************************************************
  !> \brief Classify the configured equations as explicit or implicit.
  !
  !> Builds the compact equation partition consumed by the local implicit
  !> solver. The two momentum equations are implicit; lava temperature is
  !> implicit for the lava rheology, pore pressure is implicit when enabled,
  !> and transported solid-component equations remain explicit. The implicit
  !> groups are ordered as ANALYTIC lava energy, ANALYTIC pore pressure, and
  !> NEWTON momentum, omitting inactive groups and keeping the numbering dense.
  !> \date 07/09/2012
  !>
  !> \param[in,out] equation_partition Receives the implicit mask, compact
  !>                                   explicit/implicit maps and group methods.
  !>
  !> \note Reads configured rheology and pore-pressure flags. Group methods are
  !>       metadata only here; this routine does not dispatch a solve by group.
  !>       Reports the number of implicit equations.
  !******************************************************************************

  SUBROUTINE init_problem_param( equation_partition )

    USE parameters_2d, ONLY : pore_pressure_flag, rheology_model_id,        &
         rheology_id_lava, imex_velocity_group, imex_energy_group,           &
         imex_pore_pressure_group, imex_group_method
    USE equation_metadata_2d, ONLY : resolve_imex_group_method
    IMPLICIT NONE

    TYPE(equation_partition_type), INTENT(INOUT) :: equation_partition

    LOGICAL :: implicit_mask(n_eqns)
    INTEGER :: IM_group(n_eqns)
    INTEGER :: IM_group_method(3), configured_group_method(3)
    INTEGER :: n_IM_groups, configured_group, compact_group
    LOGICAL :: group_is_active(3)

    implicit_mask(1:n_eqns) = .FALSE.
    implicit_mask(2) = .TRUE.
    implicit_mask(3) = .TRUE.

    ! Temperature
    IF ( rheology_model_id .EQ. rheology_id_lava ) THEN

       implicit_mask(4) = .TRUE.

    END IF

    ! Transported solid-component masses are explicit equations.
    implicit_mask(idx_solidEqn_first:idx_solidEqn_last) = .FALSE.

    IF ( pore_pressure_flag ) THEN

       implicit_mask(idx_pore) = .TRUE.

    END IF

    IM_group = 0
    IM_group_method = 0
    configured_group_method = 0
    n_IM_groups = 0

    IF (MIN(imex_velocity_group, imex_energy_group,                        &
         imex_pore_pressure_group) .LT. 1 .OR.                            &
         MAX(imex_velocity_group, imex_energy_group,                       &
         imex_pore_pressure_group) .GT. 3) THEN
       CALL fatal_error('IMEX group IDs must be integers from 1 to 3')
    END IF

    DO configured_group = 1, 3
       CALL resolve_imex_group_method(imex_group_method(configured_group), &
            configured_group_method(configured_group))
       IF (configured_group_method(configured_group) .LT. 0) THEN
          CALL fatal_error('IMEX_GROUP_METHOD entries must be ANALYTIC or NEWTON')
       END IF
    END DO

    ! Group IDs express the requested sweep order; repeated IDs combine terms.
    IF ( implicit_mask(4) ) THEN
       IM_group(4) = imex_energy_group
    END IF

    IF ( pore_pressure_flag ) THEN
       IM_group(idx_pore) = imex_pore_pressure_group
    END IF

    IF ( implicit_mask(2) .OR. implicit_mask(3) ) THEN
       IM_group(2:3) = imex_velocity_group
    END IF

    group_is_active = .FALSE.
    DO configured_group = 1, 3
       group_is_active(configured_group) = ANY(IM_group .EQ. configured_group)
    END DO

    DO configured_group = 1, 3
       IF (.NOT. group_is_active(configured_group)) CYCLE
       n_IM_groups = n_IM_groups + 1
       IM_group_method(n_IM_groups) = configured_group_method(configured_group)
       DO compact_group = 1, n_eqns
          IF (IM_group(compact_group) .EQ. configured_group)                &
               IM_group(compact_group) = n_IM_groups
       END DO
    END DO

    CALL equation_partition%initialize( implicit_mask, IM_group,             &
         IM_group_method(1:n_IM_groups) )

    WRITE(*,*) 'Implicit equations =',equation_partition%n_implicit

  RETURN

  END SUBROUTINE init_problem_param


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
  !> \brief Evaluate cell-local explicit sources, including configured inlet and thermal terms.
  !> \details Returns the right-hand-side source vector for explicit cell-local
  !>          physics. It includes the velocity-dependent bed-curvature term
  !>          and any active bottom-radial, fissural or lateral-radial inlet
  !>          sources, including their mass, momentum, thermal, component and
  !>          pore-pressure contributions. Hydrostatic pressure and ordinary
  !>          bed-slope forcing are represented by the spatial path operator.
  !>          Curvature is suppressed below dry_thickness_tolerance without
  !>          discarding positive conservative mass or disabling injection.
  !>          If step-averaged source coefficients have been set, those
  !>          coefficients replace pointwise window evaluation.
  !> \date 2019/12/13
  !> @author Mattia de' Michieli Vitturi
  !> \param[in] Bprimej_x Filtered cell bed derivative dB/dx [dimensionless].
  !> \param[in] Bprimej_y Filtered cell bed derivative dB/dy [dimensionless].
  !> \param[in] Bsecondj_xx Filtered cell bed derivative d^2B/dx^2 [m^-1].
  !> \param[in] Bsecondj_xy Filtered mixed bed derivative d^2B/(dx dy) [m^-1].
  !> \param[in] Bsecondj_yy Filtered cell bed derivative d^2B/dy^2 [m^-1].
  !> \param[in] grav_coeff Large-slope gravity factor G at the supplied state location.
  !> \param[in] qpj Local physical state of length n_vars+2 with component mass fractions.
  !> \param[out] expl_term Explicit cell-source vector, with the right-hand-side sign convention.
  !> \param[in] time Current simulation time [s].
  !> \param[in] cell_fract_jk Fraction of the current cell covered by the bottom radial source.
  !> \param[in] cell_fissure_fract Current cell coverage fractions for all configured fissures.
  !> \param[in] lat_arc_perim_jk Active emission-arc length within the current cell [m].
  !> \param[in] lat_n_x_jk Mean x component of the emitting arc normal.
  !> \param[in] lat_n_y_jk Mean y component of the emitting arc normal.
  !> \param[in] cell_area_jk Horizontal area of the current cell [m^2].
  !******************************************************************************

  SUBROUTINE eval_expl_terms( Bprimej_x, Bprimej_y, Bsecondj_xx , Bsecondj_xy , &
       Bsecondj_yy, grav_coeff, qpj, expl_term, time, cell_fract_jk,             &
       cell_fissure_fract,                                                       &
       lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, cell_area_jk )

    USE parameters_2d, ONLY : vel_source , T_source , xs_source , xg_source,    &
         xl_source , time_param , bottom_radial_source_flag,                    &
         pore_pressure_flag , pore_pres_fract ,                                &
         n_intervals , t_intervals , vel_intervals ,                            &
         radial_source_flag , h_source, bottom_fissural_source_flag,            &
         n_fissures, linear_vel_fissures, T_fissures, time_param_fissures,      &
         xs_fissures, xg_fissures, xl_fissures

    USE geometry_2d, ONLY : pi_g , bottom_source_area_ratio


    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y
    REAL(wp), INTENT(IN) :: Bsecondj_xx
    REAL(wp), INTENT(IN) :: Bsecondj_xy
    REAL(wp), INTENT(IN) :: Bsecondj_yy
    REAL(wp), INTENT(IN) :: grav_coeff

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(OUT) :: expl_term(n_eqns)

    REAL(wp), INTENT(IN) :: time
    REAL(wp), INTENT(IN) :: cell_fract_jk
      REAL(wp), INTENT(IN) :: cell_fissure_fract(:)

    ! Lateral radial-source per-cell geometry. Pass 0 to disable the lateral
    ! injection.
    REAL(wp), INTENT(IN) :: lat_arc_perim_jk
    REAL(wp), INTENT(IN) :: lat_n_x_jk
    REAL(wp), INTENT(IN) :: lat_n_y_jk
    REAL(wp), INTENT(IN) :: cell_area_jk

    REAL(wp) :: r_h          !< real-value flow thickness
    REAL(wp) :: r_u          !< real-value x-velocity
    REAL(wp) :: r_v          !< real-value y-velocity
    REAL(wp) :: r_Ri         !< real-value Richardson number
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity

    REAL(wp) :: r_sp_heat_mix !< real_value mixture specific heat
    REAL(wp) :: r_sp_heat_c

    REAL(wp) :: t_rem
    REAL(wp) :: t_coeff
    REAL(wp) :: h_dot

    REAL(wp) :: qp_source(n_vars+2)

    REAL(wp) :: exc_pore_pres

    REAL(wp) :: q1

    REAL(wp) :: vel_local
    INTEGER :: i_int
      INTEGER :: i_fissure


    expl_term(1:n_eqns) = 0.0_wp
    q1 = 0.0_wp

      IF ( ( qpj(1) .LE. EPSILON(1.0_wp) ) .AND. ( cell_fract_jk .EQ. 0.0_wp )    &
             .AND. ( lat_arc_perim_jk .EQ. 0.0_wp ) .AND.                          &
             ALL(cell_fissure_fract .EQ. 0.0_wp) ) RETURN

    ! Gravity terms - only meaningful when the cell has mass.
    IF ( qpj(1) .GT. EPSILON(1.0_wp) ) THEN

       r_h = qpj(1)
       r_u = qpj(idx_u)
       r_v = qpj(idx_v)

       CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,                      &
            r_sp_heat_c,r_sp_heat_mix)

       q1 = r_h * r_rho_m

       IF ( curvature_term_flag .AND. r_h .GT. dry_thickness_tolerance ) THEN

          ! Hydrostatic pressure and the ordinary bed-slope source are already
          ! represented by the HP path. Curvature remains a separate,
          ! velocity-dependent cell source using the fitted Hessian and G_c.
          ! Below the shared dry threshold the auxiliary velocity is not a
          ! resolved physical velocity: do not square it in this source.
          ! This guard neither discards mass nor disables external injection.
          CALL eval_curvature_momentum_source( Bprimej_x, Bprimej_y,         &
               Bsecondj_xx, Bsecondj_xy, Bsecondj_yy, grav_coeff, q1,        &
               r_u, r_v, expl_term(2), expl_term(3) )

       END IF

       ! The hydrostatic path/source has no thermal-energy component.
       expl_term(4) = 0.0_wp

    END IF

    ! ----------- ADDITIONAL EXPLICIT TERMS FOR BOTTOM RADIAL SOURCE ------------

    IF ( bottom_radial_source_flag .AND. ( cell_fract_jk .GT. 0.0_wp ) ) THEN

    IF ( n_intervals .GT. 0 ) THEN

       ! --- Variable velocity intervals mode ---
       ! Find which interval the current time falls in (reverse search)
       vel_local = 0.0_wp

       DO i_int = n_intervals, 1, -1

          IF ( time .GE. t_intervals(i_int) ) THEN

             vel_local = vel_intervals(i_int)
             EXIT

          END IF

       END DO

       IF ( bottom_vel_step .GE. 0.0_wp ) vel_local = bottom_vel_step

       ! vel_local = 0 means source is off
       IF ( vel_local .GT. 0.0_wp ) THEN

          t_coeff = 1.0_wp

       ELSE

          t_coeff = 0.0_wp

       END IF

       h_dot = cell_fract_jk * vel_local * bottom_source_area_ratio

    ELSE

       ! --- Original pulse mode (backward compatible) ---
       t_rem = MOD( time + time_param(4) , time_param(1) )

       IF ( time_param(3) .EQ. 0.0_wp ) THEN

          IF ( t_rem .LE. time_param(2) ) THEN

             t_coeff = 1.0_wp

          ELSE

             t_coeff = 0.0_wp

          END IF

       ELSE

          IF ( t_rem .LE. time_param(3) ) THEN

             t_coeff = ( t_rem / time_param(3) )

          ELSEIF ( t_rem .LE. time_param(2) - time_param(3) ) THEN

             t_coeff = 1.0_wp

          ELSEIF ( t_rem .LE. time_param(2) ) THEN

             t_coeff = 1.0_wp - ( t_rem - time_param(2) + time_param(3) ) /     &
                  time_param(3)

          ELSE

             t_coeff = 0.0_wp

          END IF

       END IF

       IF ( bottom_coeff_step .GE. 0.0_wp ) t_coeff = bottom_coeff_step

       h_dot = cell_fract_jk * vel_source * bottom_source_area_ratio

    END IF

    qp_source = 0.0_wp
    qp_source(1) = 1.0_wp
    qp_source(2) = 0.0_wp
    qp_source(3) = 0.0_wp
    qp_source(4) = t_source

    qp_source(idx_solid_first:idx_solid_last) = xs_source(1:n_solid)
    qp_source(idx_add_gas_first:idx_add_gas_last) = xg_source(1:n_add_gas)
    IF ( gas_flag .AND. liquid_flag ) qp_source(n_vars) = xl_source

    ! Source term transport stochastic equation
    IF ( stoch_transport_flag) qp_source(idx_stoch) = 0.0_wp

    IF ( pore_pressure_flag ) qp_source(idx_poreEqn) = 0.0_wp

    qp_source(idx_u) = 0.0_wp
    qp_source(idx_v) = 0.0_wp


    CALL mixt_var(qp_source,r_Ri,r_rho_m,r_rho_c,r_red_grav,                   &
         r_sp_heat_c,r_sp_heat_mix)

    expl_term(1) = expl_term(1) + t_coeff * h_dot * r_rho_m
    expl_term(2) = expl_term(2) + 0.0_wp
    expl_term(3) = expl_term(3) + 0.0_wp

    expl_term(4) = expl_term(4) + eval_inlet_thermal_energy_source(t_coeff,    &
         h_dot, r_rho_m, r_sp_heat_mix, t_source)

    ! source terms for the solid equations
    expl_term(idx_solidEqn_first:idx_solidEqn_last) =                          &
         expl_term(idx_solidEqn_first:idx_solidEqn_last) + t_coeff             &
         * h_dot * r_rho_m * xs_source(1:n_solid)

    ! source terms for the additional gas equations
    expl_term(idx_addGasEqn_first:idx_addGasEqn_last) =                         &
         expl_term(idx_addGasEqn_first:idx_addGasEqn_last) + t_coeff            &
         * h_dot * r_rho_m * xg_source(1:n_add_gas)

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! source term for the liquid phase
       expl_term(n_vars) = expl_term(n_vars) + t_coeff * h_dot * r_rho_m       &
            * xl_source

    END IF

    IF ( pore_pressure_flag ) THEN

       exc_pore_pres = qpj(idx_pore)

       ! we multiply the pore pressure inlet rate by the rate for q1
       ! the units of this source term are: kg^2 m^-3 s^-3
       expl_term(idx_poreEqn) = expl_term(idx_poreEqn) +                       &
            eval_pore_pressure_inlet_source(t_coeff, q1, pore_pres_fract,      &
            h_dot, r_rho_m, r_red_grav, exc_pore_pres)

    END IF

    END IF   ! bottom_radial_source_flag

    ! Each fissure uses its own clipped cell coverage, pulse and temperature,
    ! and the source composition converted at that temperature.
    IF ( bottom_fissural_source_flag ) THEN

       DO i_fissure = 1, n_fissures

          IF ( cell_fissure_fract(i_fissure) .LE. 0.0_wp ) CYCLE

          t_rem = MOD(time + time_param_fissures(4,i_fissure),                &
               time_param_fissures(1,i_fissure))
          IF (time_param_fissures(3,i_fissure) .EQ. 0.0_wp) THEN
             IF (t_rem .LE. time_param_fissures(2,i_fissure)) THEN
                t_coeff = 1.0_wp
             ELSE
                t_coeff = 0.0_wp
             END IF
          ELSEIF (t_rem .LT. time_param_fissures(3,i_fissure)) THEN
             t_coeff = t_rem / time_param_fissures(3,i_fissure)
          ELSEIF (t_rem .LE. time_param_fissures(2,i_fissure) -              &
               time_param_fissures(3,i_fissure)) THEN
             t_coeff = 1.0_wp
          ELSEIF (t_rem .LE. time_param_fissures(2,i_fissure)) THEN
             t_coeff = 1.0_wp - (t_rem - time_param_fissures(2,i_fissure) +  &
                  time_param_fissures(3,i_fissure)) /                        &
                  time_param_fissures(3,i_fissure)
          ELSE
             t_coeff = 0.0_wp
          END IF

          IF ( ALLOCATED( fissure_coeff_step ) )                              &
               t_coeff = fissure_coeff_step(i_fissure)

          IF (t_coeff .LE. 0.0_wp) CYCLE

          h_dot = cell_fissure_fract(i_fissure) * linear_vel_fissures(i_fissure)
          qp_source = 0.0_wp
          qp_source(1) = 1.0_wp
          qp_source(4) = T_fissures(i_fissure)
          qp_source(idx_solid_first:idx_solid_last) =                          &
               xs_fissures(1:n_solid,i_fissure)
          qp_source(idx_add_gas_first:idx_add_gas_last) =                      &
               xg_fissures(1:n_add_gas,i_fissure)
          IF (gas_flag .AND. liquid_flag) qp_source(n_vars) =                  &
               xl_fissures(i_fissure)
          IF (stoch_transport_flag) qp_source(idx_stoch) = 0.0_wp
          IF (pore_pressure_flag) qp_source(idx_poreEqn) = 0.0_wp
          qp_source(idx_u) = 0.0_wp
          qp_source(idx_v) = 0.0_wp

          CALL mixt_var(qp_source, r_Ri, r_rho_m, r_rho_c, r_red_grav,       &
               r_sp_heat_c, r_sp_heat_mix)

          expl_term(1) = expl_term(1) + t_coeff * h_dot * r_rho_m
          expl_term(4) = expl_term(4) + eval_inlet_thermal_energy_source(      &
               t_coeff, h_dot, r_rho_m, r_sp_heat_mix, T_fissures(i_fissure))
          expl_term(idx_solidEqn_first:idx_solidEqn_last) =                   &
               expl_term(idx_solidEqn_first:idx_solidEqn_last) +              &
               t_coeff * h_dot * r_rho_m * xs_fissures(1:n_solid,i_fissure)
          expl_term(idx_addGasEqn_first:idx_addGasEqn_last) =                  &
               expl_term(idx_addGasEqn_first:idx_addGasEqn_last) +             &
               t_coeff * h_dot * r_rho_m * xg_fissures(1:n_add_gas,i_fissure)

          IF (gas_flag .AND. liquid_flag) THEN
             expl_term(n_vars) = expl_term(n_vars) + t_coeff * h_dot         &
                  * r_rho_m * xl_fissures(i_fissure)
          END IF

          IF (pore_pressure_flag) THEN
             exc_pore_pres = qpj(idx_pore)
             expl_term(idx_poreEqn) = expl_term(idx_poreEqn) +                &
                  eval_pore_pressure_inlet_source(t_coeff, q1,                 &
                  pore_pres_fract, h_dot, r_rho_m, r_red_grav, exc_pore_pres)
          END IF

       END DO

    END IF

    ! ----------- LATERAL RADIAL SOURCE (volume formulation) ------------------
    !
    ! Replaces the prior Dirichlet-BC implementation. For each in-arc ring
    ! cell, inject mass + momentum at a rate that integrates to the user-
    ! specified MFR. cell_arc_perim is rescaled in init_source so the
    ! per-cell rate rho_source * h_source * vel_source * cell_arc_perim_jk
    ! / cell_area_jk sums to MFR over all in-arc cells.

    IF ( radial_source_flag .AND. ( lat_arc_perim_jk .GT. 0.0_wp ) .AND.       &
         ( h_source .GT. 0.0_wp ) ) THEN

       ! Time gate: keep the current-main semantics used by eval_source_bdry.
       ! TIME_PARAM(4) is an absolute cut-off for the lateral radial source.
       IF ( time .GE. time_param(4) ) THEN

          t_coeff = 0.0_wp

       ELSE

          t_rem = MOD( time , time_param(1) )
          t_coeff = 0.0_wp

          IF ( time_param(3) .EQ. 0.0_wp ) THEN
             IF ( t_rem .LE. time_param(2) ) t_coeff = 1.0_wp
          ELSE
             IF ( t_rem .LT. time_param(3) ) THEN
                t_coeff = 0.5_wp * ( 1.0_wp - COS( pi_g * t_rem /              &
                     time_param(3) ) )
             ELSEIF ( t_rem .LE. ( time_param(2) - time_param(3) ) ) THEN
                t_coeff = 1.0_wp
             ELSEIF ( t_rem .LE. time_param(2) ) THEN
                t_coeff = 0.5_wp * ( 1.0_wp + COS( pi_g * ( ( t_rem -          &
                     time_param(2) ) / time_param(3) + 1.0_wp ) ) )
             END IF
          END IF

       END IF

       IF ( lateral_coeff_step .GE. 0.0_wp ) t_coeff = lateral_coeff_step

       IF ( t_coeff .GT. 0.0_wp ) THEN

          qp_source(:) = 0.0_wp
          qp_source(1) = 1.0_wp
          qp_source(4) = T_source

          qp_source(idx_solid_first:idx_solid_last) = xs_source(1:n_solid)
          qp_source(idx_add_gas_first:idx_add_gas_last) = xg_source(1:n_add_gas)
          IF ( gas_flag .AND. liquid_flag ) qp_source(n_vars) = xl_source

          IF ( stoch_transport_flag ) qp_source(idx_stoch) = 0.0_wp
          IF ( pore_pressure_flag )   qp_source(idx_poreEqn) = 0.0_wp

          qp_source(idx_u) = 0.0_wp
          qp_source(idx_v) = 0.0_wp

          CALL mixt_var( qp_source, r_Ri, r_rho_m, r_rho_c, r_red_grav,        &
               r_sp_heat_c, r_sp_heat_mix )

          ! Per-cell injection rate (m/s) - integrates to MFR over all cells.
          h_dot = h_source * vel_source * lat_arc_perim_jk / cell_area_jk

          ! Mass equation.
          expl_term(1) = expl_term(1) + t_coeff * h_dot * r_rho_m

          ! Momentum equations - emit along the cell's outward unit normal.
          expl_term(2) = expl_term(2) + t_coeff * h_dot * r_rho_m * vel_source &
               * lat_n_x_jk
          expl_term(3) = expl_term(3) + t_coeff * h_dot * r_rho_m * vel_source &
               * lat_n_y_jk

       ! Thermal-energy equation.
          expl_term(4) = expl_term(4) + eval_inlet_thermal_energy_source(     &
               t_coeff, h_dot, r_rho_m, r_sp_heat_mix, T_source)

          ! Solid transport.
          expl_term(idx_solidEqn_first:idx_solidEqn_last) =                   &
               expl_term(idx_solidEqn_first:idx_solidEqn_last) +              &
               t_coeff * h_dot * r_rho_m * xs_source(1:n_solid)

          ! Additional gases.

          expl_term(idx_addGasEqn_first:idx_addGasEqn_last) =                  &
               expl_term(idx_addGasEqn_first:idx_addGasEqn_last) +             &
               t_coeff * h_dot * r_rho_m * xg_source(1:n_add_gas)

          IF ( gas_flag .AND. liquid_flag ) THEN
             expl_term(n_vars) = expl_term(n_vars) +                           &
                  t_coeff * h_dot * r_rho_m * xl_source
          END IF

          IF ( pore_pressure_flag ) THEN
             exc_pore_pres = qpj(idx_pore)
             expl_term(idx_poreEqn) = expl_term(idx_poreEqn) +                &
                  eval_pore_pressure_inlet_source(t_coeff, q1,                 &
                  pore_pres_fract, h_dot, r_rho_m, r_red_grav, exc_pore_pres)
          END IF

       END IF

    END IF

    RETURN

  END SUBROUTINE eval_expl_terms

  !******************************************************************************
  !> \brief Cache exact step-mean coefficients for configured inlet sources.
  !> \details Computes the mean source activation over [time,time+dt] rather
  !>          than sampling a ramp at a single time. Bottom radial and fissure
  !>          windows use their configured pulse phases; the lateral radial
  !>          window also honors its absolute cutoff and cosine ramps. Variable
  !>          bottom-source velocity intervals are averaged separately.
  !>          The resulting private coefficients are consumed by
  !>          eval_expl_terms during the current step.
  !>
  !> \param[in] time Time at the start of the step [s].
  !> \param[in] dt Length of the step [s].
  !> \note A nonpositive dt leaves the cached coefficients unchanged. No value is returned.
  !******************************************************************************

  SUBROUTINE set_source_time_coefficients( time , dt )

    USE parameters_2d, ONLY : time_param , bottom_radial_source_flag ,          &
         radial_source_flag , n_intervals , bottom_fissural_source_flag ,      &
         n_fissures , time_param_fissures

    REAL(wp), INTENT(IN) :: time
    REAL(wp), INTENT(IN) :: dt

    INTEGER :: i_fissure

    IF ( dt .LE. 0.0_wp ) RETURN

    IF ( bottom_radial_source_flag ) THEN

       IF ( n_intervals .GT. 0 ) THEN

          bottom_vel_step = interval_mean_velocity( time , dt )

       ELSE

          bottom_coeff_step = window_mean( time + time_param(4) , dt ,         &
               time_param(1) , time_param(2) , time_param(3) , .FALSE. ,      &
               HUGE(1.0_wp) )

       END IF

    END IF

    IF ( bottom_fissural_source_flag ) THEN

       IF ( .NOT. ALLOCATED( fissure_coeff_step ) )                            &
            ALLOCATE( fissure_coeff_step(n_fissures) )

       DO i_fissure = 1, n_fissures

          fissure_coeff_step(i_fissure) = window_mean( time +                 &
               time_param_fissures(4,i_fissure) , dt ,                        &
               time_param_fissures(1,i_fissure) ,                             &
               time_param_fissures(2,i_fissure) ,                             &
               time_param_fissures(3,i_fissure) , .FALSE. , HUGE(1.0_wp) )

       END DO

    END IF

    ! The lateral radial source has cosine ramps and stops at TIME_PARAM(4)
    IF ( radial_source_flag ) lateral_coeff_step = window_mean( time , dt ,     &
         time_param(1) , time_param(2) , time_param(3) , .TRUE. ,              &
         time_param(4) )

  END SUBROUTINE set_source_time_coefficients

  !******************************************************************************
  !> \brief Mean of a periodic source window over [v0, v0+dt].
  !>
  !> In each period the coefficient rises over the ramp (linearly, or as a half
  !> cosine), is 1 until act-ramp, falls to 0 at act and stays 0 to the end of
  !> the period; it is 0 from cutoff on. With ramp = 0 the window is on for
  !> 0 <= u <= act, as in the point evaluation.
  !>
  !> \param[in] v0 Start of the step, phase included [s].
  !> \param[in] dt Length of the step [s].
  !> \param[in] per Period [s].
  !> \param[in] act Active time in each period [s].
  !> \param[in] ramp Ramp time [s].
  !> \param[in] cosine Cosine ramps when true, linear ramps otherwise.
  !> \param[in] cutoff Time from which the window is 0 [s].
  !> \return Mean coefficient over the step.
  !> \note Assumes per > 0 and dt > 0; the result is clipped to [0,1].
  !******************************************************************************

  REAL(wp) FUNCTION window_mean( v0 , dt , per , act , ramp , cosine , cutoff )

    REAL(wp), INTENT(IN) :: v0 , dt , per , act , ramp , cutoff
    LOGICAL, INTENT(IN) :: cosine

    REAL(wp) :: v1 , n0 , n1 , u0 , u1 , total

    v1 = MIN( v0 + dt , cutoff )

    IF ( v1 .LE. v0 ) THEN

       window_mean = 0.0_wp
       RETURN

    END IF

    ! A window with no ramp that is active over the whole period is always on
    IF ( ( ramp .LE. 0.0_wp ) .AND. ( act .GE. per ) .AND.                    &
         ( v0 + dt .LE. cutoff ) ) THEN

       window_mean = 1.0_wp
       RETURN

    END IF

    n0 = AINT( v0 / per )
    n1 = AINT( v1 / per )
    u0 = v0 - n0 * per
    u1 = v1 - n1 * per

    IF ( n1 .EQ. n0 ) THEN

       ! A step that stays on the plateau gets exactly 1
       IF ( ( u0 .GE. ramp ) .AND. ( u1 .LE. act - ramp ) .AND.                &
            ( v0 + dt .LE. cutoff ) ) THEN

          window_mean = 1.0_wp
          RETURN

       END IF

       total = window_integral( u0 , u1 , act , ramp , cosine )

    ELSE

       total = window_integral( u0 , per , act , ramp , cosine )               &
            + ( n1 - n0 - 1.0_wp ) * ( act - ramp )                            &
            + window_integral( 0.0_wp , u1 , act , ramp , cosine )

    END IF

    window_mean = MIN( MAX( total / dt , 0.0_wp ) , 1.0_wp )

  END FUNCTION window_mean

  !******************************************************************************
  !> \brief Integral of one period of a source window between a and b.
  !>
  !> Each segment is integrated in closed form, with the cosine ramps written as
  !> sin/cos products so that short intervals lose no precision.
  !>
  !> \param[in] a Lower limit, 0 <= a <= b <= period [s].
  !> \param[in] b Upper limit [s].
  !> \param[in] act Active time [s].
  !> \param[in] ramp Ramp time [s].
  !> \param[in] cosine Cosine ramps when true, linear ramps otherwise.
  !> \return Integral of the coefficient [s].
  !> \note The limits are expected to satisfy 0 <= a <= b. A nonpositive ramp selects
  !>       a rectangular window.
  !******************************************************************************

  REAL(wp) FUNCTION window_integral( a , b , act , ramp , cosine )

    USE geometry_2d, ONLY : pi_g

    REAL(wp), INTENT(IN) :: a , b , act , ramp
    LOGICAL, INTENT(IN) :: cosine

    REAL(wp) :: x0 , x1

    window_integral = 0.0_wp

    IF ( ramp .LE. 0.0_wp ) THEN

       x1 = MIN( b , act )
       IF ( x1 .GT. a ) window_integral = x1 - a
       RETURN

    END IF

    ! ramp up
    x0 = MAX( a , 0.0_wp )
    x1 = MIN( b , ramp )

    IF ( x1 .GT. x0 ) THEN

       IF ( cosine ) THEN
          window_integral = window_integral + 0.5_wp * ( x1 - x0 )             &
               - ramp / pi_g * COS( 0.5_wp * pi_g * ( x0 + x1 ) / ramp )      &
               * SIN( 0.5_wp * pi_g * ( x1 - x0 ) / ramp )
       ELSE
          window_integral = window_integral                                    &
               + 0.5_wp * ( x1 - x0 ) * ( x1 + x0 ) / ramp
       END IF

    END IF

    ! plateau
    x0 = MAX( a , ramp )
    x1 = MIN( b , act - ramp )
    IF ( x1 .GT. x0 ) window_integral = window_integral + ( x1 - x0 )

    ! ramp down
    x0 = MAX( a , act - ramp )
    x1 = MIN( b , act )

    IF ( x1 .GT. x0 ) THEN

       IF ( cosine ) THEN
          window_integral = window_integral + 0.5_wp * ( x1 - x0 )             &
               + ramp / pi_g * COS( 0.5_wp * pi_g * ( x0 + x1 - 2.0_wp *       &
               ( act - ramp ) ) / ramp ) * SIN( 0.5_wp * pi_g * ( x1 - x0 )    &
               / ramp )
       ELSE
          window_integral = window_integral                                    &
               + 0.5_wp * ( x1 - x0 ) * ( 2.0_wp * act - x0 - x1 ) / ramp
       END IF

    END IF

  END FUNCTION window_integral

  !******************************************************************************
  !> \brief Mean velocity of the bottom radial source intervals over a step.
  !>
  !> Interval i holds VEL_INTERVALS(i) from T_INTERVALS(i) to T_INTERVALS(i+1),
  !> the last one to the end of the run; before T_INTERVALS(1) the source is
  !> off. A step inside one interval gets its velocity exactly.
  !>
  !> \param[in] t0 Start of the step [s].
  !> \param[in] dt Length of the step [s].
  !> \return Mean velocity over the step [m s^-1].
  !> \note Uses the configured interval start times and velocities; a positive dt is required.
  !******************************************************************************

  REAL(wp) FUNCTION interval_mean_velocity( t0 , dt )

    USE parameters_2d, ONLY : n_intervals , t_intervals , vel_intervals

    REAL(wp), INTENT(IN) :: t0 , dt

    REAL(wp) :: t1 , x0 , x1 , total
    INTEGER :: i_int

    t1 = t0 + dt
    total = 0.0_wp

    DO i_int = 1, n_intervals

       x0 = MAX( t0 , t_intervals(i_int) )

       IF ( i_int .LT. n_intervals ) THEN
          x1 = MIN( t1 , t_intervals(i_int+1) )
       ELSE
          x1 = t1
       END IF

       IF ( x1 .GT. x0 ) THEN

          IF ( ( x0 .EQ. t0 ) .AND. ( x1 .EQ. t1 ) ) THEN
             interval_mean_velocity = vel_intervals(i_int)
             RETURN
          END IF

          total = total + vel_intervals(i_int) * ( x1 - x0 )

       END IF

    END DO

    interval_mean_velocity = total / dt

  END FUNCTION interval_mean_velocity


  !******************************************************************************
  !> \brief Evaluate stiff cell-local sources in the real or complex state path.
  !> \details Dispatches the configured implicit rheology law and lava cooling
  !>          source, and evaluates pore-pressure diffusion when enabled. The
  !>          selected result is an equation-sized right-hand-side vector.
  !>          Supply either the real state/output pair or the complex
  !>          state/output pair. The complex path preserves smooth complex-step
  !>          derivatives for the Newton Jacobian; piecewise branches are
  !>          selected using real state values.
  !> \date 01/06/2012
  !> @author Mattia de' Michieli Vitturi
  !> \param[in] Bprimej_x Filtered cell bed derivative dB/dx [dimensionless].
  !> \param[in] Bprimej_y Filtered cell bed derivative dB/dy [dimensionless].
  !> \param[in] Zij Local effective stochastic fluctuation supplied to the rheology.
  !> \param[in] c_qj Optional complex conservative cell state [variable], used for complex-step
  !>                 source derivatives.
  !> \param[out] c_nh_term_impl Optional complex implicit-source vector [equation] corresponding
  !>                            to c_qj.
  !> \param[in] r_qj Optional real conservative cell state [variable].
  !> \param[out] r_nh_term_impl Optional real implicit-source vector [equation] corresponding
  !>                            to r_qj.
  !>
  !> \note Provide exactly one complete input/output pair: c_qj with c_nh_term_impl, or r_qj
  !>       with r_nh_term_impl. Zij is the effective stochastic fluctuation passed to rheology
  !>       closures. Pore diffusion is evaluated independently of the rheology flag; dry or
  !>       negative-pressure states, absent pore volume, and dynamic permeability without
  !>       solids contribute no pressure relaxation.
  !******************************************************************************

  SUBROUTINE eval_implicit_terms( Bprimej_x, Bprimej_y, Zij, c_qj,              &
       c_nh_term_impl, r_qj , r_nh_term_impl )

    USE COMPLEXIFY

    USE parameters_2d, ONLY : pore_pressure_flag
    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y
    REAL(wp), INTENT(IN):: Zij
    COMPLEX(wp), INTENT(IN), OPTIONAL :: c_qj(n_vars)
    COMPLEX(wp), INTENT(OUT), OPTIONAL :: c_nh_term_impl(n_eqns)
    REAL(wp), INTENT(IN), OPTIONAL :: r_qj(n_vars)
    REAL(wp), INTENT(OUT), OPTIONAL :: r_nh_term_impl(n_eqns)

    COMPLEX(wp) :: h                       !< height [m]
    COMPLEX(wp) :: u                       !< velocity (x direction) [m/s]
    COMPLEX(wp) :: v                       !< velocity (y direction) [m/s]
    COMPLEX(wp) :: w                       !< velocity (z direction) [m/s]
    COMPLEX(wp) :: T                       !< temperature [K]
    COMPLEX(wp) :: rho_m                   !< mixture density [kg/m3]
    COMPLEX(wp) :: alphas(n_solid)         !< sediment volume fractions
    COMPLEX(wp) :: alphag(n_add_gas)       !< add. gas volume fractions
    COMPLEX(wp) :: inv_rho_m               !< 1/mixture density [kg-1 m3]

    COMPLEX(wp) :: qj(n_vars)
    COMPLEX(wp) :: nh_term(n_eqns)
    COMPLEX(wp) :: source_term(n_eqns)

    COMPLEX(wp) :: mod_vel
    COMPLEX(wp) :: mod_vel_hor
    COMPLEX(wp) :: mod_vel2
    INTEGER :: i
    LOGICAL :: moving

    COMPLEX(wp) :: Zs

    COMPLEX(wp) :: exc_pore_pres

    COMPLEX(wp) :: red_grav

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       qj = c_qj

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       DO i = 1,n_vars

          qj(i) = CMPLX( r_qj(i),0.0_wp,wp )

       END DO

    ELSE

      WRITE(*,*) 'Constitutive, eval_implicit_terms: problem with arguments'
       STOP

    END IF

    ! initialize the source terms
    source_term(1:n_eqns) = CMPLX(0.0_wp,0.0_wp,wp)
    IF (rheology_flag .OR. pore_pressure_flag) THEN
       CALL c_phys_var(qj,h,u,v,T,rho_m,alphas,alphag,inv_rho_m,Zs,             &
            exc_pore_pres)
    END IF

    IF (rheology_flag) THEN

       red_grav = ( rho_m - rho_a_amb ) / rho_m * grav

       IF ( slope_correction_flag ) THEN

          w = u * Bprimej_x + v * Bprimej_y

       ELSE

          w = CMPLX( 0.0_wp , 0.0_wp , wp )

       END IF

       mod_vel2 = u**2 + v**2 + w**2
       mod_vel = SQRT( mod_vel2 )
       mod_vel_hor = SQRT( u**2 + v**2 )
       ! A purely imaginary velocity perturbation about rest still carries
       ! the derivative of linear viscous laws. Branching on REAL(V) alone
       ! would erase it. This does not assign a derivative to Coulomb's cusp.
       moving = (REAL(mod_vel,wp) .NE. 0.0_wp) .OR. (AIMAG(mod_vel) .NE. 0.0_wp)

       IF (rheology_model_id .EQ. rheology_id_voellmy) THEN
          CALL eval_voellmy_implicit_source(u, v, rho_m, red_grav, mod_vel,   &
               mod_vel2, moving, stochastic_flag, Zij, source_term(2),       &
               source_term(3))
       ELSEIF (rheology_model_id .EQ. rheology_id_mui_voellmy) THEN
          CALL eval_mui_voellmy_implicit_source(u, v, rho_m, red_grav,        &
               mod_vel, mod_vel2, moving, stochastic_flag, Zij,              &
               source_term(2), source_term(3))
       ELSEIF (rheology_model_id .EQ. rheology_id_plastic) THEN
          CALL eval_plastic_implicit_source(u, v, mod_vel, rho_m, moving,     &
               source_term(2), source_term(3))
       ELSEIF (rheology_model_id .EQ. rheology_id_lava) THEN
          CALL eval_lava_implicit_source(h, T, rho_m, u, v, mod_vel, moving,  &
               source_term(2), source_term(3))
          source_term(4) = eval_lava_cooling_source(T, radiative_term_coeff,   &
               convective_term_coeff, T_env)
       ELSEIF (rheology_model_id .EQ. rheology_id_lahar) THEN
          CALL eval_lahar_implicit_source(h, u, v, mod_vel, mod_vel_hor, T,   &
               rho_m, inv_rho_m, alphas, moving, source_term(2),             &
               source_term(3))
       ELSEIF (rheology_model_id .EQ. rheology_id_empirical) THEN
          CALL eval_empirical_implicit_source(h, u, v, mod_vel_hor, mod_vel,  &
               rho_m, moving, source_term(2), source_term(3))
       ELSEIF (rheology_model_id .EQ. rheology_id_quadratic) THEN
          CALL eval_quadratic_implicit_source(u, v, mod_vel, mod_vel_hor,     &
               mod_vel2, rho_m, moving, source_term(2), source_term(3))
       ENDIF

    ENDIF

    IF ( pore_pressure_flag ) CALL eval_pore_pressure_diffusion_source(h, T,   &
         rho_m, alphas, exc_pore_pres, source_term(idx_poreEqn))

    nh_term = source_term

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       c_nh_term_impl = nh_term

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       r_nh_term_impl = REAL( nh_term )

    END IF

    RETURN

  END SUBROUTINE eval_implicit_terms

  !******************************************************************************
  !> \brief Evaluate semi-implicit rheology and geometric source contributions.
  !> \details Dispatches the semi-implicit part of the configured rheology law
  !>          for the supplied physical state and returns it in the momentum
  !>          equations of an otherwise zero equation-sized source vector.
  !>          The selected law may use bed slope, bed curvature, reduced
  !>          gravity, pore pressure, solid fractions and the effective
  !>          stochastic variable. A dry cell or exact rest state contributes
  !>          zero, avoiding undefined velocity norms or divisions.
  !> \date 20/01/2018
  !> @author Mattia de' Michieli Vitturi
  !> \param[in] Bprimej_x Filtered cell bed derivative dB/dx [dimensionless].
  !> \param[in] Bprimej_y Filtered cell bed derivative dB/dy [dimensionless].
  !> \param[in] Bsecondj_xx Filtered cell bed derivative d^2B/dx^2 [m^-1].
  !> \param[in] Bsecondj_xy Filtered mixed bed derivative d^2B/(dx dy) [m^-1].
  !> \param[in] Bsecondj_yy Filtered cell bed derivative d^2B/dy^2 [m^-1].
  !> \param[in] grav_coeff Large-slope gravity factor G at the supplied state location.
  !> \param[in] qpj Local physical state of length n_vars+2 with component mass fractions.
  !> \param[out] nh_semi_impl_term Cell-source vector handled by the semi-implicit part of the time
  !>                               scheme.
  !> \param[in] Zj Local effective stochastic fluctuation supplied to the rheology.
  !> \note The Froude-dependent law (model 9) requires positive reduced gravity.
  !>       Neutral or buoyant mixtures have no compressive buoyancy load in this law and receive
  !>       no basal friction; their Froude denominator must not be evaluated.
  !******************************************************************************

  SUBROUTINE eval_nh_semi_impl_terms( Bprimej_x , Bprimej_y , Bsecondj_xx ,     &
       Bsecondj_xy , Bsecondj_yy , grav_coeff , qpj , nh_semi_impl_term ,       &
       Zj )

    USE parameters_2D, ONLY: pore_pressure_flag

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y
    REAL(wp), INTENT(IN) :: Bsecondj_xx
    REAL(wp), INTENT(IN) :: Bsecondj_xy
    REAL(wp), INTENT(IN) :: Bsecondj_yy
    REAL(wp), INTENT(IN) :: grav_coeff

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: Zj ! value stochastic process

    REAL(wp), INTENT(OUT) :: nh_semi_impl_term(n_eqns)

    REAL(wp) :: source_term(n_eqns)

    REAL(wp) :: mod_vel
    REAL(wp) :: mod_hor_vel

    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity
    REAL(wp) :: r_w               !< real_value z-velocity
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(wp) :: r_alphal          !< real-value liquid volume fraction
    REAL(wp) :: r_alphag(n_add_gas) !< real-value add. gas volume fractions
    REAL(wp) :: r_red_grav
    REAL(wp) :: r_rho_c
    REAL(wp) :: r_Ri
    ! mu(I) rheology variables
    REAL(wp) :: centr_force_term
    REAL(wp) :: rheology_source_x
    REAL(wp) :: rheology_source_y

    REAL(wp) :: r_sp_heat_c
    REAL(wp) :: r_sp_heat_mix

    REAL(wp) :: exc_pore_pres       !< excess pore pressure

    nh_semi_impl_term(1:n_eqns) = 0.0_wp

    ! A dry cell carries no friction source.
    IF ( qpj(1) .LE. EPSILON(1.0_wp) ) RETURN

    ! initialize and evaluate the forces terms
    source_term(1:n_eqns) = 0.0_wp
    centr_force_term = 0.0_wp

    rheology_if:IF (rheology_flag) THEN

       r_h = qpj(1)
       r_u = qpj(idx_u)
       r_v = qpj(idx_v)

       ! Every SI law below is zero at exact rest. Return before norms and
       ! divisions so optimized builds with arithmetic traps cannot speculate
       ! a reciprocal of zero from a later moving-state branch.
       IF ((r_u .EQ. 0.0_wp) .AND. (r_v .EQ. 0.0_wp)) RETURN

       CALL primitive_to_volume_fractions(qpj, r_alphas, r_alphag, r_alphal)

       CALL mixt_var(qpj, r_Ri, r_rho_m, r_rho_c, r_red_grav,                  &
            r_sp_heat_c, r_sp_heat_mix)


       IF ( slope_correction_flag ) THEN

          r_w = r_u * Bprimej_x + r_v * Bprimej_y

       ELSE

          r_w = 0.0_wp

       END IF

       mod_vel = SQRT( r_u**2 + r_v**2 + r_w**2 )

       mod_hor_vel = SQRT( r_u**2 + r_v**2 )

       rheology_model_if:IF (rheology_model_id .EQ. rheology_id_voellmy) THEN
          exc_pore_pres = 0.0_wp
          IF (pore_pressure_flag) exc_pore_pres = qpj(idx_pore)
          CALL eval_voellmy_semi_impl_source(r_h, r_u, r_v, mod_vel, r_rho_m, &
               r_red_grav, grav_coeff, curvature_term_flag, Bsecondj_xx,      &
               Bsecondj_xy, Bsecondj_yy, pore_pressure_flag, exc_pore_pres,   &
               source_term(2), source_term(3))

       ELSEIF (rheology_model_id .EQ. rheology_id_lava) THEN
          CALL eval_lava_semi_impl_source(r_u, r_v, mod_vel, source_term(2),   &
               source_term(3))

       ELSEIF (rheology_model_id .EQ. rheology_id_lahar) THEN
          CALL eval_lahar_semi_impl_source(r_h, r_u, r_v, mod_hor_vel,        &
               r_rho_m, SUM(r_alphas), mod_vel .GT. 0.0_wp,                  &
               source_term(2), source_term(3))

       ELSEIF ( rheology_model_id .EQ. rheology_id_froude ) THEN
          CALL eval_froude_source(r_h, r_u, r_v, mod_vel, mod_hor_vel,       &
               r_rho_m, r_red_grav, grav_coeff, centr_force_term,            &
               stochastic_flag, Zj, rheology_source_x, rheology_source_y)
          source_term(2) = rheology_source_x
          source_term(3) = rheology_source_y

      ELSEIF ( rheology_model_id .EQ. rheology_id_lucas ) THEN
          CALL eval_lucas_source(r_h, r_u, r_v, mod_vel, mod_hor_vel,        &
               r_rho_m, r_red_grav, grav_coeff, centr_force_term,            &
               stochastic_flag, Zj, rheology_source_x, rheology_source_y)
          source_term(2) = rheology_source_x
          source_term(3) = rheology_source_y

      ELSEIF ( rheology_model_id .EQ. rheology_id_mui ) THEN
         IF ( mod_vel .GT. 0.0_wp ) THEN ! v>0 (avoid div by 0)

            IF ( curvature_term_flag ) THEN
               centr_force_term = curvature_acceleration( Bsecondj_xx,       &
                    Bsecondj_xy, Bsecondj_yy, r_u, r_v )
            ELSE
               centr_force_term = 0.0_wp
            END IF

            exc_pore_pres = 0.0_wp
            IF ( pore_pressure_flag ) exc_pore_pres = qpj(idx_pore)
            CALL eval_mui_source(r_h, r_u, r_v, mod_vel, r_rho_m,          &
                 r_red_grav, grav_coeff, r_alphas, centr_force_term,       &
                 pore_pressure_flag, exc_pore_pres, rheology_source_x,     &
                 rheology_source_y)
            source_term(2) = rheology_source_x
            source_term(3) = rheology_source_y

         END IF


       ELSEIF ( rheology_model_id .EQ. rheology_id_mui_voellmy ) THEN

         IF ( mod_vel .GT. 0.0_wp ) THEN

            IF ( curvature_term_flag ) THEN
               centr_force_term = curvature_acceleration( Bsecondj_xx,       &
                    Bsecondj_xy, Bsecondj_yy, r_u, r_v )
            ELSE
                centr_force_term = 0.0_wp
            END IF

            exc_pore_pres = 0.0_wp
            IF ( pore_pressure_flag ) exc_pore_pres = qpj(idx_pore)
            CALL eval_mui_voellmy_source(r_h, r_u, r_v, mod_vel, r_rho_m,   &
                 r_red_grav, grav_coeff, r_alphas, centr_force_term,       &
                 pore_pressure_flag, exc_pore_pres, rheology_source_x,     &
                 rheology_source_y)
            source_term(2) = rheology_source_x
            source_term(3) = rheology_source_y

         END IF

       ENDIF rheology_model_if

    ENDIF rheology_if

    nh_semi_impl_term = source_term

    RETURN

  END SUBROUTINE eval_nh_semi_impl_terms

  !******************************************************************************
  !> \brief Compute limited erosion, deposition, entrainment and associated equation sources.
  !> \details Evaluates per-class solid erosion and settling/deposition, carrier
  !>          erosion and loss, optional air entrainment, vaporization,
  !>          thermal and pore-pressure sources, and the corresponding
  !>          conservative equation source. Erosion is limited by the available
  !>          erodible inventory over dt; deposition and carrier loss are
  !>          limited by the mass available in the flow and by the maximum
  !>          packing constraint. The bed-elevation rate is
  !>          \f$\dot{B}=(\sum_iD_i-\sum_iE_i)/(1-\phi_e)\f$.
  !>          If the flow is dry or its total solid fraction is at/below the
  !>          configured cutoff, every output is zero.
  !> \date 20/01/2018
  !> @author Mattia de' Michieli Vitturi
  !> \param[in] qpj Local physical state of length n_vars+2 with component mass fractions.
  !> \param[in] B_zone Connected-water-region indicator used by the exchange laws.
  !> \param[in] B_prime_x Filtered bed derivative dB/dx [dimensionless].
  !> \param[in] B_prime_y Filtered bed derivative dB/dy [dimensionless].
  !> \param[in] erodible Available erodible solid inventory expressed as thickness per class [m].
  !> \param[in] dt Time increment [s].
  !> \param[out] erosion_term Per-solid erosion rates after inventory limiting [m s^-1].
  !> \param[out] deposition_term Per-solid deposition rates [m s^-1].
  !> \param[out] continuous_phase_erosion_term Carrier-phase volume rate entrained with erosion
  !>                                           from the substrate [m s^-1].
  !> \param[out] continuous_phase_loss_term Carrier-phase volume loss rate, including configured
  !>                                        loss and any active pore-gas outflow [m s^-1].
  !> \param[out] eqns_term Conservative source vector [equation] for the cell-local mass-exchange
  !>                        update; includes mass, momentum, energy and component sources.
  !> \param[out] topo_term Proposed cell bed-elevation rate, positive for bed growth [m s^-1].
  !> \note B_zone selects the connected-water zone used by vaporization. Substrate temperature
  !>       is selected locally so parallel calls never mutate shared material configuration.
  !******************************************************************************

  SUBROUTINE eval_mass_exchange_terms( qpj , B_zone , B_prime_x , B_prime_y ,   &
       erodible , dt , erosion_term , deposition_term ,                         &
       continuous_phase_erosion_term , continuous_phase_loss_term , eqns_term , &
       topo_term  )

    USE parameters_2d, ONLY : erodible_deposit_flag , liquid_vaporization_flag ,&
         pore_pressure_flag , gas_loss_flag

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(IN) :: B_prime_x
    REAL(wp), INTENT(IN) :: B_prime_y
    REAL(wp), INTENT(IN) :: erodible(n_solid)
    REAL(wp), INTENT(IN) :: dt

    REAL(wp), INTENT(OUT) :: erosion_term(n_solid)
    REAL(wp), INTENT(OUT) :: deposition_term(n_solid)
    REAL(wp), INTENT(OUT) :: continuous_phase_erosion_term
    REAL(wp), INTENT(OUT) :: continuous_phase_loss_term
    REAL(wp), INTENT(OUT) :: eqns_term(n_eqns)
    REAL(wp), INTENT(OUT) :: topo_term


    REAL(wp) :: mod_vel
    REAL(wp) :: mod_vel2

    REAL(wp) :: hind_exp
    REAL(wp) :: alpha_max

    INTEGER :: i_solid

    REAL(wp) :: r_h          !< real-value flow thickness
    REAL(wp) :: r_u          !< real-value x-velocity
    REAL(wp) :: r_v          !< real-value y-velocity
    REAL(wp) :: r_W          !< real-value z-velocity
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_alphag(n_add_gas) !< real-value add.gas volume fractions
    REAL(wp) :: r_alphal          !< real-value liquid volume fraction
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_T          !< real-value mixture temperature [K]
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(wp) :: r_Zs         !< real-value stochastic variable
    REAL(wp) :: r_exc_pore_pres   !< real-value excess p pres

    REAL(wp) :: tot_erosion  !< total erosion rate [m/s]
    REAL(wp) :: tot_solid_erosion !< total solid erosion rate [m/s]

    REAL(wp) :: alphas_tot   !< total solid volume fraction

    REAL(wp) :: fluid_visc
    REAL(wp) :: inv_kin_visc
    REAL(wp) :: rhoc
    REAL(wp) :: dyn_visc_c
    REAL(wp) :: rho_c
    REAL(wp) :: kin_visc_c_local
    REAL(wp) :: r_Ri
    REAL(wp) :: r_red_grav

    !> Hindered settling velocity (units: m s-1 )
    REAL(wp) :: settling_vel

    REAL(wp) :: r_sp_heat_c
    REAL(wp) :: r_sp_heat_mix

    REAL(wp) :: entr_coeff
    REAL(wp) :: air_entr

    REAL(wp) :: dep_tot
    REAL(wp) :: ers_tot
    REAL(wp) :: rho_dep_tot
    REAL(wp) :: rho_ers_tot

    REAL(wp) :: T_erodible_local !< cell-local substrate temperature [K]
    REAL(wp) :: mass_vap_rate
    REAL(wp) :: vapor_energy_source

    REAL(wp) :: pore_pressure_term

    erosion_term(1:n_solid) = 0.0_wp
    deposition_term(1:n_solid) = 0.0_wp
    continuous_phase_erosion_term = 0.0_wp
    continuous_phase_loss_term = 0.0_wp
    eqns_term(1:n_eqns) = 0.0_wp
    topo_term = 0.0_wp
    r_Zs = 0.0_wp
    r_exc_pore_pres = 0.0_wp
    pore_pressure_term = 0.0_wp
    kin_visc_c_local = kin_visc_c

    IF ( qpj(1) .LE. epsilon(1.0_wp) ) THEN

       RETURN

    END IF

    ! parameters for Michaels and Bolger (1962) sedimentation correction
    alpha_max = 0.6_wp
    hind_exp = 4.65_wp

    IF ( qpj(1) .LE. EPSILON(1.0_wp) ) RETURN

    r_h = qpj(1)
    r_u = qpj(idx_u)
    r_v = qpj(idx_v)

    CALL primitive_to_volume_fractions(qpj, r_alphas, r_alphag, r_alphal)

    alphas_tot = SUM(r_alphas)

    ! Below the configured cutoff no conservative or bed source is applied.
    ! Return before proposing any erosion too: otherwise the caller would
    ! consume substrate inventory while eqns_term and topo_term stay zero.
    IF ( alphas_tot .LE. alphastot_min ) RETURN

    IF ( stoch_transport_flag ) r_Zs = qpj(idx_stoch)

    IF ( pore_pressure_flag ) r_exc_pore_pres = qpj(idx_pore)

    IF ( slope_correction_flag ) THEN

       r_w = r_u * B_prime_x + r_v * B_prime_y

    ELSE

       r_w = 0.0_wp

    END IF

    mod_vel2 = r_u**2 + r_v**2 + r_w**2
    mod_vel = SQRT( mod_vel2 )

    IF ( erosion_coeff .GT. 0.0_wp ) THEN

       ! empirical formulation (see Fagents & Baloga 2006, Eq. 5)
       ! here we use the solid volume fraction instead of relative density
       ! This term has units: m s-1
       tot_erosion = erosion_coeff * mod_vel * r_h * ( 1.0_wp-alphas_tot )

       tot_solid_erosion = tot_erosion * ( 1.0_wp - erodible_porosity )

       erosion_term(1:n_solid) = erodible_fract(1:n_solid)  * tot_solid_erosion

    ELSE

       tot_solid_erosion = 0.0_wp

       erosion_term(1:n_solid) = 0.0_wp

    END IF

    ! Limit the deposition during a single time step
    erosion_term(1:n_solid) = MAX(0.0_wp,MIN( erosion_term(1:n_solid),          &
         erodible(1:n_solid) / dt ) )

    tot_solid_erosion = SUM( erosion_term(1:n_solid) )
    tot_erosion = tot_solid_erosion / ( 1.0_wp - erodible_porosity )

    continuous_phase_erosion_term = tot_erosion * erodible_porosity

    r_T = qpj(4)

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,r_sp_heat_c,             &
         r_sp_heat_mix)

    IF ( rheology_model_id .EQ. rheology_id_lahar ) THEN
       fluid_visc = eval_lahar_fluid_viscosity(r_T, alphas_tot)
       ! Kinematic viscosity [m2 s-1]
       inv_kin_visc = r_rho_m / fluid_visc
       ! Continuous phase density used for the settling velocity
       rhoc = r_rho_m

    ELSE

       IF ( gas_flag .AND. sutherland_flag ) THEN

          dyn_visc_c = muRef_Suth * ( r_T / Tref_Suth )**1.5_wp *               &
               ( Tref_Suth + S_mu ) / ( r_T + S_mu )

          rho_c = pres / ( sp_gas_const_a * r_T )
          kin_visc_c_local = dyn_visc_c / rho_c

       END IF

       ! Viscosity read from input file [m2 s-1]
       inv_kin_visc = 1.0_wp / kin_visc_c_local
       ! Continuous phase density used for the settling velocity
       rhoc = r_rho_c

    END IF

    DO i_solid=1,n_solid

       IF ( ( r_alphas(i_solid) .GT. 0.0_wp ) .AND. ( settling_flag ) ) THEN

          settling_vel = settling_velocity( diam_s(i_solid) , rho_s(i_solid) ,  &
               rhoc , inv_kin_visc )

          IF (rheology_model_id .EQ. rheology_id_lahar) THEN
             deposition_term(i_solid) = eval_lahar_deposition_rate(           &
                  r_alphas(i_solid), settling_vel, r_h, dt)
          ELSE
             deposition_term(i_solid) = r_alphas(i_solid) * settling_vel

             ! Michaels and Bolger (1962) sedimentation correction accounting
             ! for hindered settling due to the presence of particles
             deposition_term(i_solid) = deposition_term(i_solid) *              &
                  ( 1.0_wp - MIN( 1.0_wp , alphas_tot / alpha_max ) )**hind_exp

             ! Limit deposition to the available amount in this time step.
             deposition_term(i_solid) = MIN(deposition_term(i_solid),           &
                  r_h * r_alphas(i_solid) / dt)
          END IF

          IF ( deposition_term(i_solid) .LT. 0.0_wp ) THEN
             CALL fatal_error('negative solid deposition rate')

          END IF

       END IF

    END DO

    IF ( liquid_flag ) THEN

       ! set the rate of loss of continuous phase
       continuous_phase_loss_term = loss_rate

    ELSE

       continuous_phase_loss_term = 0.0_wp

    END IF

    ! add the loss associated with solid deposition
    continuous_phase_loss_term =  continuous_phase_loss_term +                  &
         coeff_porosity * SUM( deposition_term(1:n_solid) )

    IF ( pore_pressure_flag .AND. gas_loss_flag ) THEN

       IF ( alphas_tot .LT. maximum_solid_packing ) THEN

         CALL eval_pore_pressure_gas_loss(r_h, r_T, r_alphas, alphas_tot,      &
              r_rho_c, r_exc_pore_pres, pore_pressure_term)

         continuous_phase_loss_term = continuous_phase_loss_term +             &
               pore_pressure_term

       END IF

    END IF

    ! limit the loss accountaing for maximum solid packing
    continuous_phase_loss_term = MIN( continuous_phase_loss_term ,              &
         r_h * MAX( 0.0_wp , maximum_solid_packing - SUM(r_alphas) ) / dt )

    ! Limit to the available carrier after reserving the pore volume of the
    ! remaining solids. An already exhausted reserve is zero, not negative:
    ! a nonnegative requested loss must not become a spurious carrier gain.
    continuous_phase_loss_term = MIN( continuous_phase_loss_term ,              &
         MAX(0.0_wp, ( r_h * ( 1.0_wp - SUM(r_alphas) ) - coeff_porosity *      &
         ( r_h * SUM(r_alphas) - dt * SUM(deposition_term(1:n_solid)) ) ) / dt) )


    ! Re-erodible deposits use the current flow temperature by model policy,
    ! but the configured substrate temperature is shared read-only state.
    T_erodible_local = T_erodible
    IF ( erodible_deposit_flag ) T_erodible_local = r_T

    IF ( entrainment_flag .AND.  ( r_h .GT. 0.0_wp ) .AND.                      &
         ( r_Ri .GT. 0.0_wp ) ) THEN

       entr_coeff = 0.075_wp / SQRT( 1.0_wp + 718.0_wp * r_Ri**2.4_wp )

       air_entr = entr_coeff * SQRT(mod_vel2)

    ELSE

       air_entr = 0.0_wp

    END IF

    eqns_term(1:n_eqns) = 0.0_wp

    ! solid total volume deposition rate [m s-1]
    dep_tot = SUM( deposition_term )
    ! solid total volume deposition rate [m s-1]
    ers_tot = SUM( erosion_term )

    ! solid total mass deposition rate [kg m-2 s-1]
    rho_dep_tot = DOT_PRODUCT( rho_s , deposition_term )
    ! solid total mass erosion rate [kg m-2 s-1]
    rho_ers_tot = DOT_PRODUCT( rho_s , erosion_term )

    ! total mass equation source term [kg m-2 s-1]:
    ! deposition, erosion and entrainment are considered
    eqns_term(1) = rho_a_amb * air_entr + rho_ers_tot - rho_dep_tot             &
         + rho_c_sub * continuous_phase_erosion_term                            &
         - r_rho_c * continuous_phase_loss_term

    ! x-momenutm equation source term [kg m-1 s-2]:
    ! only deposition contribute to change in momentum, erosion does not carry
    ! any momentum inside the flow
    eqns_term(2) = - r_u * ( rho_dep_tot + r_rho_c * continuous_phase_loss_term )

    ! y-momentum equation source term [kg m-1 s-2]:
    ! only deposition contribute to change in momentum, erosion does not carry
    ! any momentum inside the flow
    eqns_term(3) = - r_v * ( rho_dep_tot + r_rho_c * continuous_phase_loss_term )

    ! Thermal-energy equation source term [kg s-3]:
    ! deposition, erosion and entrainment are considered
    eqns_term(4) = eval_mass_exchange_thermal_energy(r_T, rho_s, sp_heat_s,    &
         deposition_term, erosion_term, r_rho_c, r_sp_heat_c,                  &
         continuous_phase_loss_term, T_erodible_local, rho_c_sub,              &
         continuous_phase_erosion_term, T_ambient, sp_heat_a, rho_a_amb,       &
         air_entr)

    ! solid phase mass equation source term [kg m-2 s-1]:
    ! due to solid erosion and deposition
    eqns_term(idx_solidEqn_first:idx_solidEqn_last) = rho_s(1:n_solid) *        &
         ( erosion_term(1:n_solid) - deposition_term(1:n_solid) )

    IF ( gas_flag .AND. liquid_vaporization_flag .AND. ( B_zone .NE. 0 ) ) THEN

       CALL eval_vaporization_sources(r_T, rho_s, sp_heat_s, deposition_term,  &
            gamma_steam, sp_heat_g(1), mass_vap_rate, vapor_energy_source)

       eqns_term(1) = eqns_term(1) + mass_vap_rate
       eqns_term(4) = eqns_term(4) + vapor_energy_source
       eqns_term(5+n_solid) = eqns_term(5+n_solid) + mass_vap_rate

    END IF

    IF ( stoch_transport_flag) THEN

       ! Evaluate the equation term related to the noise transport equation
       ! (if transport flag is false nothing is done)
       eqns_term(idx_stochEqn) = eqns_term(1) * r_Zs

    END IF

    IF ( pore_pressure_flag ) THEN

       ! Equation for q1*exc_pore_press
       eqns_term(idx_poreEqn) = eval_pore_pressure_mass_exchange_source(       &
            eqns_term(1), r_exc_pore_pres)

    END IF

    ! erodible layer thickness source terms [m s-1]:
    ! due to erosion and deposition of solid+continuous phase
    topo_term = ( dep_tot - ers_tot ) / ( 1.0_wp - erodible_porosity )

    RETURN

  END SUBROUTINE eval_mass_exchange_terms


  !******************************************************************************
  !> \brief Construct the physical inlet state for an active radial-source pulse.
  !> \details Evaluates the configured periodic pulse at time and returns the
  !>          inlet state in the canonical q/qp layout: depth, depth-weighted
  !>          x/y velocity, temperature, component fractions and appended
  !>          x/y velocities. During an inactive pulse, depth, momentum and
  !>          composition remain zero while temperature is set to T_source.
  !>          For an active pulse with coefficient c_t, depth and velocity are
  !>          scaled as h=c_t*h_source and u,v=sqrt(c_t)*vel_source*(vect_x,
  !>          vect_y), so depth-weighted momentum scales as c_t^(3/2).
  !> \date 2019/12/01
  !> @author Mattia de' Michieli Vitturi
  !> \param[in] time Current simulation time [s].
  !> \param[in] vect_x X component of the inlet direction vector.
  !> \param[in] vect_y Y component of the inlet direction vector.
  !> \param[out] source_bdry Inlet state of length n_vars+2 in the canonical physical-state
  !>                         layout; composition is populated only when the pulse is active.
  !> \note Returns the temperature-only inactive state when time reaches the configured cutoff.
  !******************************************************************************

  SUBROUTINE eval_source_bdry( time, vect_x , vect_y , source_bdry )

    USE parameters_2d, ONLY : h_source , vel_source , T_source , xs_source,    &
         xg_source, xl_source, time_param

    USE geometry_2d, ONLY : pi_g

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: time
    REAL(wp), INTENT(IN) :: vect_x
    REAL(wp), INTENT(IN) :: vect_y
    REAL(wp), INTENT(OUT) :: source_bdry(n_vars+2)

    REAL(wp) :: t_rem
    REAL(wp) :: t_coeff

    source_bdry = 0.0_wp
    source_bdry(4) = T_source

    IF ( time .GE. time_param(4) ) THEN

       RETURN

    END IF

    t_rem = MOD( time , time_param(1) )

    t_coeff = 0.0_wp

    IF ( time_param(3) .EQ. 0.0_wp ) THEN

       IF ( t_rem .LE. time_param(2) ) t_coeff = 1.0_wp

    ELSE

       IF ( t_rem .LT. time_param(3) ) THEN

          t_coeff = 0.5_wp * ( 1.0_wp - COS( pi_g * t_rem / time_param(3) ) )

       ELSEIF ( t_rem .LE. ( time_param(2) - time_param(3) ) ) THEN

          t_coeff = 1.0_wp

       ELSEIF ( t_rem .LE. time_param(2) ) THEN

          t_coeff = 0.5_wp * ( 1.0_wp + COS( pi_g * ( ( t_rem - time_param(2) ) &
               / time_param(3) + 1.0_wp ) ) )

       END IF

    END IF

    ! Preserve the common dry-state convention while the periodic source is
    ! inactive. In particular, no composition is attached to zero thickness.
    IF (t_coeff .LE. EPSILON(1.0_wp)) RETURN

    ! The exponents of t_coeff are such that Ri does not depend on t_coeff
    source_bdry(1) = t_coeff * h_source
    source_bdry(2) = t_coeff**1.5_wp * h_source * vel_source * vect_x
    source_bdry(3) = t_coeff**1.5_wp * h_source * vel_source * vect_y
    source_bdry(4) = T_source

    source_bdry(idx_solid_first:idx_solid_last) = xs_source(1:n_solid)
    source_bdry(idx_add_gas_first:idx_add_gas_last) = xg_source(1:n_add_gas)
    IF ( gas_flag .AND. liquid_flag ) source_bdry(n_vars) = xl_source

    source_bdry(idx_u) = t_coeff**0.5_wp * vel_source * vect_x
    source_bdry(idx_v) = t_coeff**0.5_wp * vel_source * vect_y

    RETURN

  END SUBROUTINE eval_source_bdry

END MODULE equation_terms_2d
