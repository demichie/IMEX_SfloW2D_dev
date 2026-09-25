!********************************************************************************
!> \brief Model-specific equation terms
!********************************************************************************
MODULE equation_terms_2d

  USE diagnostics_2d, ONLY : fatal_error

  USE constitutive_parameters_2d
  USE state_conversion_2d

  USE equation_metadata_2d, ONLY : equation_partition_type

  USE parameters_2d, ONLY : wp, sp, tolh
  USE parameters_2d, ONLY : n_eqns, n_vars, n_solid, n_add_gas,                &
       n_stoch_vars, n_pore_vars
  USE parameters_2d, ONLY : rheology_flag, rheology_model,                      &
       liquid_flag, gas_flag, slope_correction_flag,                            &
       curvature_term_flag, stochastic_flag, mean_field_flag,                  &
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
  PUBLIC :: eval_local_speeds_x, eval_local_speeds_y
  PUBLIC :: eval_fluxes
  PUBLIC :: eval_inertial_flux, eval_hydrostatic_coefficient
  PUBLIC :: limit_component_mass_flux
  PUBLIC :: eval_expl_terms, integrate_friction_term
  PUBLIC :: eval_implicit_terms, eval_nh_semi_impl_terms
  PUBLIC :: eval_mass_exchange_terms, eval_source_bdry

CONTAINS

  !> Bound the signed sum of all independently transported component-mass
  !> fluxes by the signed total-mixture mass flux. A single common factor
  !> preserves the relative component proportions.
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
  !> \brief Initialization of relaxation flags
  !
  !> This subroutine set the number and the flags of the non-hyperbolic
  !> terms.
  !> \date 07/09/2012
  !******************************************************************************

  SUBROUTINE init_problem_param( equation_partition )

    USE parameters_2d, ONLY : pore_pressure_flag
    IMPLICIT NONE

    TYPE(equation_partition_type), INTENT(INOUT) :: equation_partition

    LOGICAL :: implicit_mask(n_eqns)

    implicit_mask(1:n_eqns) = .FALSE.
    implicit_mask(2) = .TRUE.
    implicit_mask(3) = .TRUE.

    ! Temperature
    IF ( rheology_model .EQ. 3 ) THEN

       implicit_mask(4) = .TRUE.

    END IF

    ! Transported solid-component masses are explicit equations.
    implicit_mask(idx_solidEqn_first:idx_solidEqn_last) = .FALSE.

    IF ( pore_pressure_flag ) THEN

       implicit_mask(idx_pore) = .TRUE.

    END IF

    CALL equation_partition%initialize( implicit_mask )

    WRITE(*,*) 'Implicit equations =',equation_partition%n_implicit

  RETURN

  END SUBROUTINE init_problem_param


  !******************************************************************************
  !> \brief Local Characteristic speeds x direction
  !
  !> This subroutine computes from the physical variable evaluates the largest
  !> positive and negative characteristic speed in the x-direction.
  !> \param[in]     qpj           array of local physical variables
  !> \param[out]    vel_min       minimum x-velocity
  !> \param[out]    vel_max       maximum x-velocity
  !> @author
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds_x(qpj,grav_coeff,vel_min,vel_max)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: grav_coeff

    REAL(wp), INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL(wp) :: r_h          !< real-value flow thickness [m]
    REAL(wp) :: r_u          !< real-value x-velocity [m s-1]
    REAL(wp) :: r_v          !< real-value y-velocity [m s-1]
    REAL(wp) :: r_Ri         !< real-value Richardson number
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg m-3]
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg m-3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity [m s-2]
    REAL(wp) :: r_celerity
    REAL(wp) :: r_sp_heat_c
    REAL(wp) :: r_sp_heat_mix

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,r_sp_heat_c,             &
         r_sp_heat_mix)

    r_h = qpj(1)
    r_u = qpj(idx_u)
    r_v = qpj(idx_v)

    IF ( r_red_grav * r_h .LT. 0.0_wp ) THEN

       vel_min(1:n_eqns) = r_u
       vel_max(1:n_eqns) = r_u

    ELSE

       r_celerity = SQRT( r_red_grav * r_h * grav_coeff )
       vel_min(1:n_eqns) = r_u - r_celerity
       vel_max(1:n_eqns) = r_u + r_celerity

    END IF

    RETURN

  END SUBROUTINE eval_local_speeds_x

  !******************************************************************************
  !> \brief Local Characteristic speeds y direction
  !
  !> This subroutine computes from the physical variable evaluates the largest
  !> positive and negative characteristic speed in the y-direction.
  !> \param[in]     qpj           array of local physical variables
  !> \param[out]    vel_min       minimum y-velocity
  !> \param[out]    vel_max       maximum y-velocity
  !> @author
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds_y(qpj,grav_coeff,vel_min,vel_max)

    IMPLICIT none

    REAL(wp), INTENT(IN)  :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: grav_coeff
    REAL(wp), INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL(wp) :: r_h          !< real-value flow thickness
    REAL(wp) :: r_u          !< real-value x-velocity
    REAL(wp) :: r_v          !< real-value y-velocity
    REAL(wp) :: r_Ri         !< real-value Richardson number
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity
    REAL(wp) :: r_celerity
    REAL(wp) :: r_sp_heat_c
    REAL(wp) :: r_sp_heat_mix

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,r_sp_heat_c,             &
         r_sp_heat_mix)

    r_h = qpj(1)
    r_u = qpj(idx_u)
    r_v = qpj(idx_v)

    IF ( r_red_grav * r_h .LT. 0.0_wp ) THEN

       vel_min(1:n_eqns) = r_v
       vel_max(1:n_eqns) = r_v

    ELSE

       r_celerity = SQRT( grav_coeff * r_red_grav * r_h )
       vel_min(1:n_eqns) = r_v - r_celerity
       vel_max(1:n_eqns) = r_v + r_celerity

    END IF

    RETURN

  END SUBROUTINE eval_local_speeds_y

  !******************************************************************************
  !> \brief Thermodynamic coefficient of the hydrostatic path
  !>
  !> Gamma=rho_m*g' is recomputed from the final reconstructed primitive state;
  !> no derived thermodynamic quantity is reconstructed independently.
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
  !> \brief Pressure-free conservative flux used by HP-PCCU
  !>
  !> Every retained conservative component is transported by the normal face
  !> velocity.  Hydrostatic pressure is represented exclusively by the path
  !> contribution and must not appear in this flux.
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
  !> \brief Hyperbolic Fluxes
  !
  !> This subroutine evaluates the numerical fluxes given the conservative
  !> variables qcj and physical variables qpj.
  !> \date 01/06/2012
  !> \param[in]     qcj      real local conservative variables
  !> \param[in]     qpj      real local physical variables
  !> \param[in]     dir      direction of the flux (1=x,2=y)
  !> \param[out]    flux     real  fluxes
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_fluxes(qcj,qpj,grav_coeff,dir,flux)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qcj(n_vars)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: grav_coeff
    INTEGER, INTENT(IN) :: dir

    REAL(wp), INTENT(OUT) :: flux(n_eqns)

    REAL(wp) :: r_h          !< real-value flow thickness [m]
    REAL(wp) :: r_u          !< real-value x-velocity [m s-1]
    REAL(wp) :: r_v          !< real-value y-velocity [m s-1]
    REAL(wp) :: r_Ri         !< real-value Richardson number
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg m-3]
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg m-3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity [m s-2]
    REAL(wp) :: r_sp_heat_c
    REAL(wp) :: r_sp_heat_mix

    pos_thick:IF ( qpj(1) .GT. EPSILON(1.0_wp) ) THEN

       r_h = qpj(1)
       r_u = qpj(idx_u)
       r_v = qpj(idx_v)

       CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,                      &
            r_sp_heat_c,r_sp_heat_mix)

       IF ( dir .EQ. 1 ) THEN

          ! Mass flux in x-direction: u * ( rhom * h )
          flux(1) = r_u * qcj(1)

          ! x-momentum flux in x-direction + hydrostatic pressure term
          flux(2) = r_u * qcj(2) + 0.5_wp * r_rho_m *                           &
               grav_coeff * r_red_grav * r_h**2

          ! y-momentum flux in x-direction: u * ( rho * h * v )
          flux(3) = r_u * qcj(3)

          ! Thermal-energy flux in x-direction: u * (rhom * Cp * h * T).
          ! Hydrostatic pressure work is not part of the retained equation.
          flux(4) = r_u * qcj(4)

          ! Solid-component mass flux in x-direction.
          flux(idx_solidEqn_first:idx_solidEqn_last) = r_u *                    &
               qcj(idx_solid_first:idx_solid_last)

          ! Additional-gas component mass flux in x-direction.
          flux(idx_addGasEqn_first:idx_addGasEqn_last) = r_u *                  &
               qcj(idx_add_gas_first:idx_add_gas_last)

          IF ( stoch_transport_flag) THEN

             ! Flux of stochastic variables
             flux(idx_stochEqn) = r_u * qcj(idx_stoch)

          END IF

          IF ( pore_pressure_flag) THEN

             ! Flux of pore pressure
             flux(idx_poreEqn) = r_u * qcj(idx_Pore)

          END IF

          ! Mass flux of liquid in x-direction: u * ( h * alphal * rhol )
          IF ( gas_flag .AND. liquid_flag ) flux(n_vars) = r_u * qcj(n_vars)

       ELSEIF ( dir .EQ. 2 ) THEN

          ! flux G (derivated wrt y in the equations)
          flux(1) = r_v * qcj(1)

          flux(2) = r_v * qcj(2)

          flux(3) = r_v * qcj(3) + 0.5_wp * r_rho_m *                           &
               grav_coeff * r_red_grav * r_h**2

          ! Thermal-energy flux in y-direction: v * (rhom * Cp * h * T).
          ! Hydrostatic pressure work is not part of the retained equation.
          flux(4) = r_v * qcj(4)

          ! Solid-component mass flux in y-direction.
          flux(idx_solidEqn_first:idx_solidEqn_last) = r_v *                    &
               qcj(idx_solid_first:idx_solid_last)

          ! Additional-gas component mass flux in y-direction.
          flux(idx_addGasEqn_first:idx_addGasEqn_last) = r_v *                  &
               qcj(idx_add_gas_first:idx_add_gas_last)

          IF ( stoch_transport_flag) THEN

             ! Flux of stochastic variables
             flux(idx_stochEqn) = r_v * qcj(idx_stoch)

          END IF

          IF ( pore_pressure_flag ) THEN

             ! Flux of pore pressure
             flux(idx_poreEqn) = r_v * qcj(idx_pore)

          END IF

          ! Mass flux of liquid in x-direction: u * ( h * alphal * rhol )
          IF ( gas_flag .AND. liquid_flag ) flux(n_vars) = r_v * qcj(n_vars)

       END IF

    ELSE

       flux(1:n_eqns) = 0.0_wp

    ENDIF pos_thick

    RETURN

  END SUBROUTINE eval_fluxes


  !******************************************************************************
  !> \brief Explicit source term
  !
  !> This subroutine evaluates the non-hyperbolic terms to be treated explicitely
  !> in the DIRK numerical scheme (e.g. gravity,source of mass). The sign of the
  !> terms is taken with the terms on the right-hand side of the equations.
  !> \date 2019/12/13
  !> \param[in]     B_primej_x         local x-slope
  !> \param[in]     B_primej_y         local y_slope
  !> \param[in]     B_secondj_xx       local 2nd derivative in x-direction
  !> \param[in]     B_secondj_xy       local 2nd derivative in xy-direction
  !> \param[in]     B_secondj_yy       local 2nd derivative in y-direction
  !> \param[in]     grav_coeff         correction factor for topography slope
  !> \param[in]     d_grav_coeff_dx    x-derivative of grav_coeff
  !> \param[in]     d_grav_coeff_dy    y-derivative of grav_coeff
  !> \param[in]     qpj                physical variables
  !> \param[in]     time               simlation time (needed for source)
  !> \param[in]     cell_fract_jk      fraction of cell contributing to source
  !> \param[out]    expl_term          explicit term
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_expl_terms( Bprimej_x, Bprimej_y, Bsecondj_xx , Bsecondj_xy , &
       Bsecondj_yy, grav_coeff, d_grav_coeff_dx , d_grav_coeff_dy ,             &
       qpj, expl_term, time, cell_fract_jk,                                    &
       lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, cell_area_jk )

    USE parameters_2d, ONLY : vel_source , T_source , xs_source , xg_source,    &
         xl_source , time_param , bottom_radial_source_flag,                    &
         pore_pressure_flag , pore_pres_fract ,                                &
         n_intervals , t_intervals , vel_intervals ,                            &
         radial_source_flag , h_source

    USE geometry_2d, ONLY : pi_g


    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y
    REAL(wp), INTENT(IN) :: Bsecondj_xx
    REAL(wp), INTENT(IN) :: Bsecondj_xy
    REAL(wp), INTENT(IN) :: Bsecondj_yy
    REAL(wp), INTENT(IN) :: grav_coeff
    REAL(wp), INTENT(IN) :: d_grav_coeff_dx
    REAL(wp), INTENT(IN) :: d_grav_coeff_dy

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)      !< local physical variables
    REAL(wp), INTENT(OUT) :: expl_term(n_eqns) !< local explicit forces

    REAL(wp), INTENT(IN) :: time
    REAL(wp), INTENT(IN) :: cell_fract_jk

    !> Lateral radial-source per-cell geometry. Pass 0 to disable the lateral
    !> injection.
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

    REAL(wp) :: r_tilde_grav
    REAL(wp) :: centr_force_term

    REAL(wp) :: qp_source(n_vars+2)

    REAL(wp) :: exc_pore_pres

    REAL(wp) :: q1

    REAL(wp) :: vel_local
    INTEGER :: i_int


    expl_term(1:n_eqns) = 0.0_wp
    q1 = 0.0_wp

    IF ( ( qpj(1) .LE. EPSILON(1.0_wp) ) .AND. ( cell_fract_jk .EQ. 0.0_wp )    &
         .AND. ( lat_arc_perim_jk .EQ. 0.0_wp ) ) RETURN

    ! Gravity terms - only meaningful when the cell has mass.
    IF ( qpj(1) .GT. EPSILON(1.0_wp) ) THEN

       r_h = qpj(1)
       r_u = qpj(idx_u)
       r_v = qpj(idx_v)

       CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,                      &
            r_sp_heat_c,r_sp_heat_mix)

       q1 = r_h * r_rho_m

       IF ( curvature_term_flag ) THEN

          centr_force_term = Bsecondj_xx * r_u**2 + 2.0_wp*Bsecondj_xy*r_u*r_v &
               + Bsecondj_yy * r_v**2

       ELSE

          centr_force_term = 0.0_wp

       END IF

       ! Hydrostatic pressure and the ordinary bed-slope source are represented
       ! together by the HP path in hyperbolic_2d.  Only the velocity-dependent
       ! curvature force remains here.  Gate E uses G=1; the validated G-weight
       ! is introduced in the following slope-correction gate.
       r_tilde_grav = centr_force_term

       ! units of dqc(2)/dt [kg m-1 s-2]
       expl_term(2) = - r_rho_m * r_tilde_grav * r_h * Bprimej_x

       ! units of dqc(3)/dt [kg m-1 s-2]
       expl_term(3) = - r_rho_m * r_tilde_grav * r_h * Bprimej_y

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

       ! vel_local = 0 means source is off
       IF ( vel_local .GT. 0.0_wp ) THEN

          t_coeff = 1.0_wp

       ELSE

          t_coeff = 0.0_wp

       END IF

       h_dot = cell_fract_jk * vel_local

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

       h_dot = cell_fract_jk * vel_source

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

    expl_term(4) = expl_term(4) + t_coeff * h_dot * r_rho_m * r_sp_heat_mix    &
         * t_source

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
            t_coeff * ( q1 * pore_pres_fract * h_dot * r_rho_m * r_red_grav +  &
            exc_pore_pres * h_dot * r_rho_m )

    END IF

    END IF   ! bottom_radial_source_flag

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
          expl_term(4) = expl_term(4) + t_coeff * h_dot * r_rho_m              &
               * r_sp_heat_mix * T_source

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
             expl_term(idx_poreEqn) = expl_term(idx_poreEqn) +                 &
                  t_coeff * ( q1 * pore_pres_fract * h_dot * r_rho_m           &
                  * r_red_grav + exc_pore_pres * h_dot * r_rho_m )
          END IF

       END IF

    END IF

    RETURN

  END SUBROUTINE eval_expl_terms


  !******************************************************************************
  !> \brief Analytic integration of friction terms
  !
  !> This subroutine integrate analytically the friction term for turbulent
  !> friction only.
  !> \date 2021/12/10
  !> \param[inout]     r_qj            real conservative variables
  !> \param[in]        dt              time step
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE integrate_friction_term( r_qj , dt )

    IMPLICIT NONE

    REAL(wp), INTENT(INOUT) :: r_qj(n_vars)
    REAL(wp), INTENT(IN) :: dt


    REAL(wp) :: r_qp(n_vars+2)
    REAL(wp) :: p_dyn
    REAL(wp) :: r_h
    REAL(wp) :: r_u
    REAL(wp) :: r_v

    REAL(wp) :: mod_vel_hor
    REAL(wp) :: mod_vel

    IF ( r_qj(1) .EQ. 0.0_wp ) RETURN

    CALL qc_to_qp(r_qj , r_qp , p_dyn )

    r_h = r_qp(1)

    r_u = r_qp(idx_u)
    r_v = r_qp(idx_v)

    mod_vel_hor = SQRT( r_u**2 + r_v**2 )

    ! A cell at rest has no velocity for the friction to damp, and the
    ! reciprocal below would divide by zero.
    IF ( mod_vel_hor .LE. EPSILON(1.0_wp) ) RETURN

    mod_vel = 1.0_wp / ( 1.0_wp  / mod_vel_hor + friction_factor / r_h * dt )

    r_u = r_u * ( mod_vel / mod_vel_hor )
    r_v = r_v * ( mod_vel / mod_vel_hor )

    r_qp(2) = r_h*r_u
    r_qp(3) = r_h*r_v

    r_qp(idx_u) = r_u
    r_qp(idx_v) = r_v

    CALL qp_to_qc(r_qp,r_qj)

    RETURN

  END SUBROUTINE integrate_friction_term


  !******************************************************************************
  !> \brief Implicit source terms
  !
  !> This subroutine evaluates the source terms  of the system of equations,
  !> both for real or complex inputs, that are treated implicitely in the DIRK
  !> numerical scheme.
  !> \date 01/06/2012
  !> \param[in]     Bprimej_x       topography slope in x-direction
  !> \param[in]     Bprimej_y       topography slope in y-direction
  !> \param[in]     c_qj            complex conservative variables
  !> \param[in]     r_qj            real conservative variables
  !> \param[out]    c_nh_term_impl  complex non-hyperbolic terms
  !> \param[out]    r_nh_term_impl  real non-hyperbolic terms
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_implicit_terms( Bprimej_x, Bprimej_y, Zij, c_qj,              &
       c_nh_term_impl, r_qj , r_nh_term_impl )

    USE COMPLEXIFY

    USE geometry_2d, ONLY : pi_g

    USE parameters_2d, ONLY : pore_pressure_flag
    USE parameters_2d, ONLY : four_thirds , neg_four_thirds

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y
    REAL(wp), INTENT(IN):: Zij
    COMPLEX(wp), INTENT(IN), OPTIONAL :: c_qj(n_vars)
    COMPLEX(wp), INTENT(OUT), OPTIONAL :: c_nh_term_impl(n_eqns)
    REAL(wp), INTENT(IN), OPTIONAL :: r_qj(n_vars)
    REAL(wp), INTENT(OUT), OPTIONAL :: r_nh_term_impl(n_eqns)

    COMPLEX(wp) :: h                       !< height [m]
    COMPLEX(wp) :: inv_h                   !< 1/height [m-1]
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
    COMPLEX(wp) :: gamma
    COMPLEX(wp) :: rho_g(n_add_gas)
    REAL(wp) :: h_threshold
    COMPLEX(wp) :: rho_c

    INTEGER :: i

    !--- Lahars rheology model variables

    !> Temperature in C
    COMPLEX(wp) :: Tc

    COMPLEX(wp) :: expA , expB

    !> 1st param for fluid viscosity empirical relationship (O'Brian et al, 1993)
    COMPLEX(wp) :: alpha1    ! (units: kg m-1 s-1 )

    !> Fluid dynamic viscosity (units: kg m-1 s-1 )
    COMPLEX(wp) :: fluid_visc

    !> Total friction slope (dimensionless): s_f = s_v+s_td+s_y
    COMPLEX(wp) :: s_f

    !> Viscous slope component of total Friction (dimensionless)
    COMPLEX(wp) :: s_v

    !> Turbulent dispersive slope component of total friction (dimensionless)
    COMPLEX(wp) :: s_td

    COMPLEX(wp) :: temp_term

    COMPLEX(wp) :: c_tau

    COMPLEX(wp) :: Zs

    COMPLEX(wp) :: exc_pore_pres

    COMPLEX(wp) :: D_coeff
    REAL(wp) :: gamma_gas
    COMPLEX(wp) :: gas_compressibility

    COMPLEX(wp) :: rho_gas

    COMPLEX(wp) :: coeff_porosity

    COMPLEX(wp) :: porosity
    COMPLEX(wp) :: f_inhibit
    COMPLEX(wp) :: dyn_visc_c
    COMPLEX(wp) :: red_grav
    COMPLEX(wp) :: hydraulic_permeability_local
    COMPLEX(wp) :: kin_visc_c_local


    !> calculation of permeability
    REAL(wp) :: diam_sauter ! Sauter mean diameter [m]
    REAL(wp) :: sphericity_mean ! < Surface-area-weighted mean sphericity
   !> Inhibit factor for friction term
   COMPLEX(wp) :: x  ! local complex variable for f_inhibit calculation
   COMPLEX(wp) :: s ! local complex variable for f_inhibit calculation
   INTEGER :: n ! summation index for f_inhibit calculation
   REAL(wp) :: tt, term, width
   REAL(wp) :: alpha_trans_dynamic ! changes with height


    COMPLEX(wp) :: turb_stress


    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       qj = c_qj

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       DO i = 1,n_vars

          qj(i) = CMPLX( r_qj(i),0.0_wp,wp )

       END DO

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'
       STOP

    END IF

    ! initialize the source terms
    source_term(1:n_eqns) = CMPLX(0.0_wp,0.0_wp,wp)
    f_inhibit = CMPLX(1.0_wp,0.0_wp,wp)

    IF (rheology_flag) THEN

       CALL c_phys_var(qj,h,u,v,T,rho_m,alphas,alphag,inv_rho_m,Zs,             &
            exc_pore_pres)

       red_grav = ( rho_m - rho_a_amb ) / rho_m * grav

       IF ( slope_correction_flag ) THEN

          w = u * Bprimej_x + v * Bprimej_y

       ELSE

          w = CMPLX( 0.0_wp , 0.0_wp , wp )

       END IF

       mod_vel2 = u**2 + v**2 + w**2
       mod_vel = SQRT( mod_vel2 )
       mod_vel_hor = SQRT( u**2 + v**2 )

       IF ( rheology_model .EQ. 1 .OR. rheology_model .EQ. 12 ) THEN
          ! Voellmy Salm rheology

          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN

             ! Modify xi if using a stochastic model
             IF (stochastic_flag) THEN
               xi_temp = xi + Zij
               ! Limit the boundaies of stochastic friction
               IF (xi_temp .LT. 1._wp) xi_temp = 1._wp
               IF (xi_temp .GT. 1e5_wp) xi_temp = 1e5_wp
             ELSE
               xi_temp = xi
             END IF

             turb_stress = rho_m * red_grav / xi_temp * mod_vel2

             source_term(2) = source_term(2) - turb_stress * ( u / mod_vel )
             source_term(3) = source_term(3) - turb_stress * ( v / mod_vel )

          ENDIF

       ELSEIF ( rheology_model .EQ. 2 ) THEN

          ! Plastic rheology
          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN

             source_term(2) = source_term(2) - rho_m * tau * ( u / mod_vel )

             source_term(3) = source_term(3) - rho_m * tau * ( v / mod_vel )

          ENDIF

       ELSEIF ( rheology_model .EQ. 3 ) THEN

          h_threshold = 1.0E-10_wp

          T_env = 300.0_wp

          ! Temperature dependent rheology
          IF ( REAL(h) .GT. h_threshold ) THEN

             ! Equation 6 from Costa & Macedonio, 2005
             gamma = 3.0_wp * nu_ref / h * EXP( - visc_par * ( T - T_ref ) )

          ELSE

             ! Equation 6 from Costa & Macedonio, 2005
             gamma = 3.0_wp * nu_ref / h_threshold * EXP( - visc_par            &
                  * ( T - T_ref ) )

          END IF

          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN

             ! Last R.H.S. term in equation 2 from Costa & Macedonio, 2005
             source_term(2) = source_term(2) - rho_m * gamma * u

             ! Last R.H.S. term in equation 3 from Costa & Macedonio, 2005
             source_term(3) = source_term(3) - rho_m * gamma * v

          ENDIF

          source_term(4) = source_term(4) - radiative_term_coeff * ( T**4 -     &
               T_env**4 ) - convective_term_coeff * ( T - T_env )


       ELSEIF ( rheology_model .EQ. 4 ) THEN

          ! Lahars rheology (O'Brien 1993, FLO2D)

          ! alpha1 here has units: kg m-1 s-1
          ! in Table 2 from O'Brien 1988, the values reported have different
          ! units ( poises). 1poises = 0.1 kg m-1 s-1

          IF ( wp .EQ. sp ) THEN

             h_threshold = 1.0E-10_wp

          ELSE

             h_threshold = 1.0E-20_wp

          END IF

          ! convert from Kelvin to Celsius
          Tc = T - 273.15_wp

          ! the dependance of viscosity on temperature is modeled with the
          ! equation presented at:
          ! https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781118131473.app3
          !
          ! In addition, we use a reference value provided in input at a
          ! reference temperature. This value is used to scale the equation
          IF ( REAL(Tc) .LT. 20.0_wp ) THEN

             expA = 1301.0_wp / ( 998.333_wp + 8.1855_wp * ( Tc - 20.0_wp )     &
                  + 0.00585_wp * ( Tc - 20.0_wp )**2 ) - 1.30223_wp

             alpha1 = alpha1_coeff * 1.0E-3_wp * 10.0_wp**expA

          ELSE

             expB = ( 1.3272_wp * ( 20.0_wp - Tc ) - 0.001053_wp *              &
                  ( Tc - 20.0_wp )**2 ) / ( Tc + 105.0_wp )

             alpha1 = alpha1_coeff * 1.002E-3_wp * 10.0_wp**expB

          END IF

          ! Fluid dynamic viscosity [kg m-1 s-1]
          fluid_visc = alpha1 * EXP( beta1 * SUM(alphas) )

          IF ( REAL(h) .GT. h_threshold ) THEN

             inv_h = 1.0_wp / h

             ! Viscous slope component (dimensionless)
             s_v = Kappa * fluid_visc * mod_vel * 0.125_wp * inv_rho_m *        &
                  inv_grav * inv_h**2

             ! Turbulent dispersive component (dimensionless)
             s_td = n_td2 * mod_vel2 * inv_h**four_thirds

          ELSE

             ! Viscous slope component (dimensionless)
             s_v = Kappa * fluid_visc * mod_vel / ( 8.0_wp * rho_m * grav *     &
                  h_threshold**2 )

             ! Turbulent dispersive components (dimensionless)
             s_td = n_td2 * mod_vel2 * h_threshold**neg_four_thirds

          END IF

          ! Total implicit friction slope (dimensionless)
          s_f = s_v + s_td

          IF ( REAL(mod_vel) .GT. 0.0_wp ) THEN

             temp_term = grav * rho_m * h * s_f / mod_vel_hor

             ! same units of dqc(2)/dt: kg m-1 s-2
             source_term(2) = source_term(2) - u * temp_term

             ! same units of dqc(3)/dt: kg m-1 s-2
             source_term(3) = source_term(3) - v * temp_term

          END IF

       ELSEIF ( rheology_model .EQ. 5 ) THEN

          c_tau = 1.0E-3_wp / ( 1.0_wp + 10.0_wp * h ) * mod_vel

          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN

             source_term(2) = source_term(2) - rho_m * c_tau * ( u/mod_vel_hor )
             source_term(3) = source_term(3) - rho_m * c_tau * ( v/mod_vel_hor )

          END IF


       ELSEIF ( rheology_model .EQ. 6 ) THEN

          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN

             source_term(2) = source_term(2) - rho_m * ( u / mod_vel_hor ) *    &
                  friction_factor * mod_vel2

             source_term(3) = source_term(3) - rho_m * ( v / mod_vel_hor ) *    &
                  friction_factor * mod_vel2

          ENDIF

       ! Coulomb function rheology: mu(Fr)
       ELSEIF ( rheology_model .EQ. 9 ) THEN

       ENDIF

    ENDIF

    IF ( pore_pressure_flag ) THEN

       h_threshold = 1.0E-10_wp
       hydraulic_permeability_local = CMPLX(hydraulic_permeability,0.0_wp,wp)
       kin_visc_c_local = CMPLX(kin_visc_c,0.0_wp,wp)

       gamma_gas = sp_heat_a / ( sp_heat_a - sp_gas_const_a )
       gas_compressibility = 1.0_wp / ( gamma_gas * pres )
       rho_gas = pres / ( sp_gas_const_a * T )

       porosity = 1.0_wp - SUM(alphas)

 ! ---------------------------------------------------------
       IF (dynamic_permeability_flag) THEN

         ! diameter and density of particles
         IF ( n_solid .GT. 1) THEN ! if more than one solid phase
            ! Sauter diameter
            diam_sauter = sauter_diameter( real(alphas) )
            ! Surface-area-weighted mean sphericity (consistent with Sauter diameter)
            ! Sphericity is weighted by alpha_i * d_i^2 (proportional to surface area)
            sphericity_mean = DOT_PRODUCT( real(alphas) * diam_s**2 , sphericity_s ) / &
                              DOT_PRODUCT( real(alphas) , diam_s**2 )
         ELSE ! if only one solid phase
            diam_sauter = diam_s(1)
            sphericity_mean = sphericity_s(1)
         END IF
         ! calculating permeability using Carman-Kozeny equation with sphericity
         hydraulic_permeability_local = ( porosity**3.0_wp * &
                                   ( diam_sauter * sphericity_mean )**2.0_wp ) / &
                                   ( 150.0_wp * (1.0_wp - porosity)**2.0_wp )

       END IF

       ! ---------------------------------------------------------

       IF ( gas_flag .AND. sutherland_flag ) THEN

          dyn_visc_c = muRef_Suth * ( T / Tref_Suth )**1.5_wp *                 &
               ( Tref_Suth + S_mu ) / ( T + S_mu )

          kin_visc_c_local = dyn_visc_c / rho_gas

       END IF

       ! Eq. 7 from Gueugneau et al, 2017
       D_coeff = hydraulic_permeability_local /                                &
            ( porosity * kin_visc_c_local * rho_gas *                         &
            gas_compressibility )

       ! Equation 12 from Gueugneau et al, 2017
       ! At zero excess pressure the source value is still zero, but the
       ! positive-side derivative is needed by the complex-step Jacobian.
       IF ( ( REAL(exc_pore_pres) .GE. 0.0_wp )                                 &
            .AND. ( REAL(h) .GT. 0.0_wp) ) THEN


         select case ( trim(adjustl(f_inhibit_mode)) )

         case ( 'OFF' )
            f_inhibit = 1.0_wp

         case ('STEP')
            !-------------------------------------------------------------------!
            ! Calculate f_inhibit based on SOLID FRACTION
            ! for STEP function case
            !-------------------------------------------------------------------!
            if ( SUM(alphas) .LT. maximum_solid_packing ) then
               f_inhibit = 1.0_wp
            else
               f_inhibit = 0.0_wp
            end if

         case ('STATIC')

            !f_inhibit = MAX(0.0_wp , 1.0_wp - (SUM(alphas) /                     &
            !CMPLX(maximum_solid_packing,0.0_wp,wp))**alpha_trans )

            !-------------------------------------------------------------------!
            ! Calculate f_inhibit based on SOLID FRACTION
            ! for COMPLEX input
            ! Generalized smoothstep S_N(x) for compx input
            !-------------------------------------------------------------------!

            ! Validate N_inh read from namelist (should be non-negative integer)
            if (N_inh < 0) then
               WRITE(*,*) 'ERROR: N_inh must be >= 0. N_inh=', N_inh
               STOP
            end if
            ! use solid volume fraction source for x
            x = SUM(alphas)
            ! Normalize x to [0,1] using real part (or use abs(x) for magnitude)
            width = maximum_solid_packing - alpha_trans ! right_limit - left_limit
            tt = (real(x) - alpha_trans) / width
            if (tt < 0.0_wp) tt = 0.0_wp
            if (tt > 1.0_wp) tt = 1.0_wp
            ! Initialize smoothstep sum
            s = CMPLX(0.0_wp,0.0_wp,wp)
            ! call the precomputed pascal triangle coefficients
            do n = 0, N_inh
               term = pascal_coeff(n)  * tt**(N_inh + n + 1.0_wp)
               s = s + term
            end do
            ! flip result: [0→1] becomes [1→0]
            s = 1.0_wp - s
            f_inhibit = s

         case ('DYNAMIC')
            !-------------------------------------------------------------------!
            ! Calculate f_inhibit based on SOLID FRACTION and HEIGHT of the flow
            ! for COMPLEX input
            ! Generalized smoothstep S_N(x) for compx input,
            !-------------------------------------------------------------------!

            ! Calcualte dynamic alpha_trans based on height (from empirical fit)
            alpha_trans_dynamic = 10.0_wp ** (-0.28_wp) *                       &
                 REAL(h,wp) ** 0.12_wp

            ! Correct for alpha_trans_dynamic exceeding maximum_solid_packing
            if ( real(alpha_trans_dynamic) .GE. maximum_solid_packing ) then
                alpha_trans_dynamic = 0.99_wp * maximum_solid_packing
            end if

            ! Validate N_inh read from namelist (should be non-negative integer)
            if (N_inh < 0) then
               WRITE(*,*) 'ERROR: N_inh must be >= 0. N_inh=', N_inh
               STOP
            end if
            ! use solid volume fraction source for x
            x = SUM(alphas)
            ! Normalize x to [0,1] using real part (or use abs(x) for magnitude)
            width = maximum_solid_packing - alpha_trans_dynamic ! right_limit - left_limit
            tt = (real(x) - alpha_trans_dynamic) / width
            if (tt < 0.0_wp) tt = 0.0_wp
            if (tt > 1.0_wp) tt = 1.0_wp
            ! Initialize smoothstep sum
            s = CMPLX(0.0_wp,0.0_wp,wp)
            ! call the precomputed pascal triangle coefficients
            do n = 0, N_inh
               term = pascal_coeff(n)  * tt**(N_inh + n + 1.0_wp)
               s = s + term
            end do
            ! flip result: [0→1] becomes [1→0]
            s = 1.0_wp - s
            f_inhibit = s

         end select
         ! -------------------------------------------------------------------!

	 ! calculate source term for pore pressure equatio
	 source_term(idx_poreEqn) = - rho_m * ( pi_g / 2.0_wp )**2 *           &
               D_coeff / MAX(h_threshold,h) * exc_pore_pres * f_inhibit

       END IF

    END IF

    nh_term = source_term

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       c_nh_term_impl = nh_term

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       r_nh_term_impl = REAL( nh_term )

    END IF

    RETURN

  END SUBROUTINE eval_implicit_terms

  !******************************************************************************
  !> \brief Non-Hyperbolic semi-implicit terms
  !
  !> This subroutine evaluates the non-hyperbolic terms that are solved
  !> semi-implicitely by the solver. For example, any discontinuous term that
  !> appears in the friction terms.
  !> \date 20/01/2018
  !> \param[in]     grav3_surf         gravity correction
  !> \param[in]     qcj                real conservative variables
  !> \param[out]    nh_semi_impl_term  real non-hyperbolic terms
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_nh_semi_impl_terms( Bprimej_x , Bprimej_y , Bsecondj_xx ,     &
       Bsecondj_xy , Bsecondj_yy , grav_coeff , qcj , qpj , nh_semi_impl_term , &
       Zj )

    USE parameters_2D, ONLY: pore_pressure_flag

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y
    REAL(wp), INTENT(IN) :: Bsecondj_xx
    REAL(wp), INTENT(IN) :: Bsecondj_xy
    REAL(wp), INTENT(IN) :: Bsecondj_yy
    REAL(wp), INTENT(IN) :: grav_coeff

    REAL(wp), INTENT(IN) :: qcj(n_vars)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: Zj ! value stochastic process

    REAL(wp), INTENT(OUT) :: nh_semi_impl_term(n_eqns)

    REAL(wp) :: source_term(n_eqns)

    REAL(wp) :: mod_vel
    REAL(wp) :: mod_hor_vel

    REAL(wp) :: h_threshold

    !--- Lahars rheology model variables

    !> Yield strenght (units: kg m-1 s-2)
    REAL(wp) :: tau_y

    !> Yield slope component of total friction (dimensionless)
    REAL(wp) :: s_y

    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity
    REAL(wp) :: r_w               !< real_value z-velocity
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(wp) :: r_T               !< real-value temperature [K]
    REAL(wp) :: r_alphal          !< real-value liquid volume fraction
    REAL(wp) :: r_alphag(n_add_gas) !< real-value add. gas volume fractions
    REAL(wp) :: r_red_grav
    REAL(wp) :: r_rho_c
    REAL(wp) :: r_Ri
    REAL(wp) :: Fr                !< Froude number
    REAL(wp) :: U_f               !< norm velocity for frirction
    REAL(wp) :: muF               !< mu(fr) o rmu(U)
    !REAL(wp) :: muU               !< mu(U)
    ! REAL(wp) :: Fr_x                !< Froude number
    ! REAL(wp) :: Fr_y                !< Froude number
    ! REAL(wp) :: mu_Fr_x             !< mu(fr)_x
    ! REAL(wp) :: mu_Fr_y             !< mu(fr)_y

    ! mu(I) rheology variables
    REAL(wp) :: diam_characteristic !< area weighted mean diameter of particles
    REAL(wp) :: rho_particle !< volume weighted mean density of particles
    REAL(wp) :: I !< inertial number
    REAL(wp) :: mu_I !< mu(I)
    REAL(wp) :: shear_rate !< shear rate using horizontal velocity only
    REAL(wp) :: tau_cosA !< shear stress projected to horizontal plane
    REAL(wp) :: vert_stress_eff !< effective vertical stress
    REAL(wp) :: eff_normal_stress !< effective normal stress

    REAL(wp) :: temp_term
    REAL(wp) :: centr_force_term

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

       CALL primitive_to_volume_fractions(qpj, r_alphas, r_alphag, r_alphal)

       r_T = qpj(4)

       CALL mixt_var(qpj, r_Ri, r_rho_m, r_rho_c, r_red_grav,                  &
            r_sp_heat_c, r_sp_heat_mix)


       IF ( slope_correction_flag ) THEN

          r_w = r_u * Bprimej_x + r_v * Bprimej_y

       ELSE

          r_w = 0.0_wp

       END IF

       mod_vel = SQRT( r_u**2 + r_v**2 + r_w**2 )

       mod_hor_vel = SQRT( r_u**2 + r_v**2 )

       ! Voellmy Salm rheology
       rheology_model_if:IF ( rheology_model .EQ. 1 ) THEN

          IF ( mod_vel .GT. 0.0_wp ) THEN

             IF ( curvature_term_flag ) THEN

                ! See Eq. (3) Xia & Liang, 2018 Eng.Geol.
                ! centrifugal force term: (u,v)^T*Hessian*(u,v)
                centr_force_term = Bsecondj_xx * r_u**2 + 2.0_wp * Bsecondj_xy  &
                     * r_u * r_v + Bsecondj_yy * r_v**2

             ELSE

                centr_force_term = 0.0_wp

             END IF

             ! See Eq. (3,4) Xia & Liang, 2018 Eng.Geol.
             ! add the contribution on mu (with coeff for large slope)
             ! and the contribution of centr. force (with coeff for slope)
             ! a = grav_coeff * ( r_red_grav + centr_force_term )
             ! 1/phi = SQRT(grav_coeff)
             temp_term = mu * r_rho_m * ( r_red_grav + centr_force_term ) * r_h &
                  * SQRT(grav_coeff)

             temp_term = MAX(0.0_wp, temp_term)

             IF ( pore_pressure_flag ) THEN

                exc_pore_pres = qpj(idx_pore)

                ! See Eq. (2) Gueugneau et al. 2017, GRL
                ! add the contribution of pore pressure ( with coeff for slope)
                temp_term = MAX(0.0_wp, temp_term - mu * SQRT(grav_coeff)       &
                     * exc_pore_pres )

             END IF

             ! Friction terms cannot accelerate the flow
             ! this term is parallel to the full vel vector (u,v,w)
             ! tangential to the topography
             temp_term = MAX(0.0_wp,temp_term)

             ! horizontal terms
             ! units of dqc(2)/dt=d(rho h u)/dt (kg m-1 s-2)
             source_term(2) = source_term(2) - temp_term * r_u / mod_vel
             ! units of dqc(3)/dt=d(rho h v)/dt (kg m-1 s-2)
             source_term(3) = source_term(3) - temp_term * r_v / mod_vel

          END IF

          ! Plastic rheology
       ELSEIF ( rheology_model .EQ. 2 ) THEN


          ! Temperature dependent rheology
       ELSEIF ( rheology_model .EQ. 3 ) THEN

          IF ( mod_vel .GT. 0.0_wp ) THEN

             ! units of dqc(2)/dt [kg m-1 s-2]
             source_term(2) = source_term(2) - tau0 * r_u / mod_vel

             ! units of dqc(3)/dt [kg m-1 s-2]
             source_term(3) = source_term(3) - tau0 * r_v / mod_vel

          END IF

          ! Lahars rheology (O'Brien 1993, FLO2D)
       ELSEIF ( rheology_model .EQ. 4 ) THEN

          h_threshold = 1.0E-20_wp

          ! Yield strength (units: kg m-1 s-2)
          tau_y = alpha2 * ( EXP( beta2 * SUM(r_alphas) ) - 1.0_wp )

          IF ( r_h .GT. h_threshold ) THEN

             ! Yield slope component (dimensionless)
             s_y = tau_y / ( grav * r_rho_m * r_h )

          ELSE

             ! Yield slope component (dimensionless)
             s_y = tau_y / ( grav * r_rho_m * h_threshold )

          END IF

          IF ( mod_vel .GT. 0.0_wp ) THEN

             temp_term = grav * r_rho_m * r_h * s_y

             ! units of dqc(2)/dt [kg m-1 s-2]
             source_term(2) = source_term(2) - temp_term * r_u / mod_hor_vel

             ! units of dqc(3)/dt [kg m-1 s-2]
             source_term(3) = source_term(3) - temp_term * r_v / mod_hor_vel

          END IF

       ELSEIF ( rheology_model .EQ. 9 ) THEN

          ! From Zhu et al. 2020
          ! (DOI: https://doi.org/10.1007/s10035-020-01053-7)
          ! mu_0: Coulomb friction coefficient at Fr=+inf
          ! mu_inf: Coulomb friction coefficient at Fr=0
          ! Fr : froude number (!!!computed here using the total velocity and
          ! the thickness instead of the particle holdup!!!)
          ! Fr_0 : Renormalization factor controlling the gradient of the
          ! function Zj : Value of the OU process at given cell (they can be
          ! tranformed)

          ! should also use rho or alpha?
          ! must add something for curvature or temperature?

          ! Compute friction only if mass is flowing (this implies that there
          ! is mass)

          IF ( mod_vel .GT. 0.0_wp ) THEN
             ! Computing froude number (The definition in Zhu 2020 et Roche 2021
             ! is sligtlhy different!)
             Fr = mod_vel / SQRT(r_red_grav * r_h)

             ! Mofidy deterministic Fr if needed
             IF ( stochastic_flag ) THEN
               ! Add mean field correction or stochastic noise
               IF (mean_field_flag) THEN
                  Fr = Fr !+ MeanFieldCorrection(grav, r_h, Fr) ! TEMP MODEL, DO NOT USE IT
               ELSE
                  IF (r_h .GT. 0.01_wp) THEN
                     Fr = Fr + Zj
                  END IF
               END IF
               ! Set Foude number to zero if negative
               IF (Fr .LT. 0._wp) Fr = 0._wp
             END IF

             ! Evaluating mu(fr), (mu_0 < mu_inf), it is increasing with Fr
             muF = mu_inf + (mu_0 - mu_inf) * exp(-Fr/Fr_0)

             ! Compute temporal friction term (temp_term = 0 if h=0)
             temp_term = r_rho_m *  ( muF * r_h * grav_coeff * ( r_red_grav +   &
                  centr_force_term ) )

             ! Update the source term
             ! units of dqc(2)/dt=d(rho h u)/dt (kg m-1 s-2)
             source_term(2) = source_term(2) - temp_term * r_u / mod_hor_vel

             ! units of dqc(3)/dt=d(rho h v)/dt (kg m-1 s-2)
             source_term(3) = source_term(3) - temp_term * r_v / mod_hor_vel

          ELSE ! If ||u|| = 0 then there will be no friction in this code

            muF = 0._wp

          END IF

      ELSEIF ( rheology_model .EQ. 10 ) THEN

          ! From Lucas et al. 2014 (DOI: 10.1038/ncomms4417)
          ! mu_0: Coulomb friction coefficient at U=0
          ! mu_inf: Coulomb friction coefficient at U=+inf
          ! U_w : Renormalization factor controlling the gradient of the function
          ! Zj : Value of the OU process at given cell (they can be tranformed)

          ! Compute friction only if mass is flowing (this implies that there is mass)
          IF ( mod_vel .GT. 0.0_wp ) THEN

             ! Mofidy deterministic U if needed
             IF ( stochastic_flag ) THEN
               ! Add mean field correction or stochastic noise
               IF (mean_field_flag) THEN
                  U_f = mod_vel !TEMP MODEL, DO NOT USE IT
               ELSE
                  U_f = mod_vel + Zj
               END IF
             ELSE
               U_f = mod_vel
             END IF

             ! Evaluating mu(U), (mu_0 > mu_inf), it is decreasing with U
             IF ( U_f .LT. U_w) THEN
               muF = mu_0
             ELSE
               ! (for small h ; u->inf but if not true)
               ! Max friction if h -> 0 : else weakening friction
               IF (r_h .LT. 0.1_wp) THEN
                  muF = mu_0
               ELSE
                  muF = mu_inf + (mu_0-mu_inf)/(U_f/U_w)
               END IF
               ! Alternative :
               ! Use exp and limit with sigmoid : (a=15, b=0.5)
               ! muF = mu_0 - ((mu_inf - mu_0) * exp(-(U_f-U_w)/U_w) *             &
               !        (1._wp / (1._wp + exp(-15_wp * (r_h - 0.5_wp))))
             END IF

             ! Compute temporal friction term (temp_term = 0 if h=0)
             temp_term = r_rho_m *  ( muF * r_h * grav_coeff * ( r_red_grav +      &
                  centr_force_term ) )

             ! Update the source term
             ! units of dqc(2)/dt=d(rho h u)/dt (kg m-1 s-2)
             source_term(2) = source_term(2) - temp_term * r_u / mod_hor_vel

             ! units of dqc(3)/dt=d(rho h v)/dt (kg m-1 s-2)
             source_term(3) = source_term(3) - temp_term * r_v / mod_hor_vel

          ELSE ! If ||u|| = 0 then there will be no friction in this code

            muF = 0._wp

          END IF

      ELSEIF ( rheology_model .EQ. 11 ) THEN
         !    mu(I) rheology for dense granular flows   !
         ! -------------------------------------------- !
         ! from https://doi.org/10.1016/j.jnnfm.2015.02.006
            !mu_s = 0.48_wp         ! static friction coefficient
            !mu_2 = 0.73_wp         ! friction at high inertial number above which flow accelerates
            !muI_inf = 1.2_wp        ! friction at high I to avoid plateau (Barker et al. 2017) doi:10.1017/jfm.2017.428
            !I_0 = 0.279_wp         ! reference inertial

         IF ( mod_vel .GT. 0.0_wp ) THEN ! v>0 (avoid div by 0)

            ! diameter and density of particles
            IF ( n_solid .GT. 1) THEN ! if more than one solid phase
               ! Sauter diameter
               diam_characteristic = sauter_diameter( r_alphas )
               ! Particle density
               rho_particle = average_density_solids( r_alphas )
            ELSE ! if only one solid phase
               diam_characteristic = diam_s(1)
               rho_particle = rho_s(1)
            END IF

            ! flow thickness (avoid div by 0)
            IF (r_h .LT. EPSILON(1.0_wp)) THEN
               r_h = EPSILON(1.0_wp)
            END IF
            ! IF (r_h .LT. diam_characteristic) THEN
            !    r_h = diam_characteristic
            ! END IF

            ! Centrifugal force contribution (computed once, used consistently below).
            ! See Eq. (3) Xia & Liang, 2018 Eng.Geol.
            ! centrifugal force term: (u,v)^T*Hessian*(u,v)
            IF ( curvature_term_flag ) THEN
               centr_force_term = Bsecondj_xx * r_u**2 + 2.0_wp * Bsecondj_xy * r_u * r_v + &
                                 Bsecondj_yy * r_v**2
            ELSE
               centr_force_term = 0.0_wp
            END IF

            ! Effective vertical stress: gravity + centrifugal − pore pressure
            ! See Eq. (2) Gueugneau et al. 2017, GRL for pore pressure contribution
            IF ( pore_pressure_flag ) THEN
               ! pore pressure at the base (See Eq. (2) Gueugneau et al. 2017, GRL)
               exc_pore_pres = qpj(idx_pore)
               vert_stress_eff = r_rho_m * ( r_red_grav + centr_force_term ) * r_h &
                           - exc_pore_pres
            ELSE
               vert_stress_eff = r_rho_m * ( r_red_grav + centr_force_term ) * r_h
            END IF
            vert_stress_eff = MAX(vert_stress_eff, 0.0_wp)

           ! Effective normal stress — computed unconditionally; independent of shear rate
            eff_normal_stress = MAX(0.0_wp, vert_stress_eff * SQRT(grav_coeff))
            IF (eff_normal_stress .LT. EPSILON(1.0_wp)) THEN
               eff_normal_stress = EPSILON(1.0_wp)
            END IF

            ! Shear rate at the base (Bouchut et al. 2021, Eqn. 2.15)
            shear_rate = 5.0_wp/2.0_wp * mod_vel / (r_h * SQRT(grav_coeff))

            ! Inertial number
            I = diam_characteristic * shear_rate / SQRT(eff_normal_stress / rho_particle)

            ! mu(I) with regularisation
            mu_I = (mu_s * I_0 + mu_2 * I + muI_inf * I**2) / (I_0 + I)

            ! Friction force per unit horizontal area, parallel to full velocity
            ! vector (u,v,w), tangential to the topography.
            temp_term = mu_I * eff_normal_stress

            ! Apply friction force to momentum equations (only if there's velocity)
            ! Friction force projected onto x direction, opposite to motion
            source_term(2) = source_term(2) - temp_term * (r_u / mod_vel)
            ! Friction force projected onto y direction, opposite to motion
            source_term(3) = source_term(3) - temp_term * (r_v / mod_vel)

         END IF


       ELSEIF ( rheology_model .EQ. 12 ) THEN

         !    mu(I) rheology with Voellmy-Salm style projection        !
         ! ----------------------------------------------------------- !
         ! mu(I) rheology for dense granular flows:
         ! from https://doi.org/10.1016/j.jnnfm.2015.02.006
         !mu_s = 0.7_wp          ! static friction coefficient
         !mu_2 = 1.4_wp          ! friction at high inertial number
         !muI_inf = 0.05_wp      ! friction at high I (Barker et al. 2017)
         !I_0 = 0.3_wp           ! reference inertial number

         IF ( mod_vel .GT. 0.0_wp ) THEN

            ! diameter and density of particles
            IF ( n_solid .GT. 1) THEN ! if more than one solid phase
               ! Sauter diameter
               diam_characteristic = sauter_diameter( r_alphas )
               ! Particle density
               rho_particle = average_density_solids( r_alphas )
            ELSE ! if only one solid phase
               diam_characteristic = diam_s(1)
               rho_particle = rho_s(1)
            END IF

            !Flow thickness (avoid div by 0)
            IF (r_h .LT. EPSILON(1.0_wp)) THEN
               r_h = EPSILON(1.0_wp)
            END IF

            ! Centrifugal force contribution (computed once, used consistently below).
            ! See Eq. (3) Xia & Liang, 2018 Eng.Geol.
            ! centrifugal force term: (u,v)^T * Hessian * (u,v)
            IF ( curvature_term_flag ) THEN
               centr_force_term = Bsecondj_xx * r_u**2 + 2.0_wp * Bsecondj_xy  &
                     * r_u * r_v + Bsecondj_yy * r_v**2
            ELSE
                centr_force_term = 0.0_wp
            END IF

            ! Effective vertical stress: gravity + centrifugal − pore pressure.
            ! See Eq. (2) Gueugneau et al. 2017, GRL for pore pressure contribution.
            IF ( pore_pressure_flag ) THEN
               exc_pore_pres = qpj(idx_pore)
               vert_stress_eff = r_rho_m * ( r_red_grav + centr_force_term ) * r_h &
                                 - exc_pore_pres
            ELSE
               vert_stress_eff = r_rho_m * ( r_red_grav + centr_force_term ) * r_h
            END IF

            vert_stress_eff = MAX(vert_stress_eff, 0.0_wp)

            !Effective normal stress at the base (clamped at 0, then floored above 0)
            eff_normal_stress = MAX(0.0_wp, vert_stress_eff * SQRT(grav_coeff))
            IF (eff_normal_stress .LT. EPSILON(1.0_wp)) THEN
               eff_normal_stress = EPSILON(1.0_wp)
            END IF


            !Shear rate at the base (Eqn. 2.15 from Bouchut et al. 2021)
            shear_rate = 5.0_wp/2.0_wp * mod_vel / (r_h * SQRT(grav_coeff))

            !Inertial number
            I = diam_characteristic*shear_rate/SQRT(eff_normal_stress/rho_particle)

            !Coefficient of friction -> accounting for regularisation
            mu_I = (mu_s * I_0 + mu_2 * I + muI_inf * I**2) / ( I_0 + I )


            ! Friction force per unit horizontal area, parallel to full velocity
            ! vector (u,v,w), tangential to the topography.
            temp_term = mu_I * eff_normal_stress

            ! horizontal terms
            ! units of dqc(2)/dt=d(rho h u)/dt (kg m-1 s-2)
            source_term(2) = source_term(2) - temp_term * (r_u / mod_vel)
            ! units of dqc(3)/dt=d(rho h v)/dt (kg m-1 s-2)
            source_term(3) = source_term(3) - temp_term * (r_v / mod_vel)

         END IF

       ENDIF rheology_model_if

    ENDIF rheology_if

    nh_semi_impl_term = source_term

    RETURN

  END SUBROUTINE eval_nh_semi_impl_terms

  !******************************************************************************
  !> \brief Mass exchange terms
  !
  !> This subroutine evaluates all the terms related to mass exhange between the
  !> flow and the environment: air entrainment, deposition, erosion, entrainment
  !> of water vapour.
  !> \date 20/01/2018
  !> \param[in]     qpj                real physical variables
  !> \param[in]     B_zone             integer defining topo type (ground,water)
  !> \param[in]     B_prime_x          local topography x-slope
  !> \param[in]     B_prime_y          local topography y-slope
  !> \param[in]     erodible           erodible thickness of solid phases
  !> \param[in]     dt                 time step
  !> \param[out]    erosion_term       erosion rate of solid phases
  !> \param[out]    deposition_term    deposition rate of solid phases
  !> \param[out]    continuous_phase_erosion_term
  !> \param[out]    continuous_phase_loss_term
  !> \param[out]    eqns_term          terms associated with mass exchange
  !> \param[out]    topo_term          rate of change of topography
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_mass_exchange_terms( qpj , B_zone , B_prime_x , B_prime_y ,   &
       erodible , dt , erosion_term , deposition_term ,                         &
       continuous_phase_erosion_term , continuous_phase_loss_term , eqns_term , &
       topo_term  )

    USE geometry_2d, ONLY : pi_g

    USE parameters_2d, ONLY : erodible_deposit_flag , liquid_vaporization_flag ,&
         pore_pressure_flag , gas_loss_flag

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)              !< physical variables
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(IN) :: B_prime_x
    REAL(wp), INTENT(IN) :: B_prime_y
    REAL(wp), INTENT(IN) :: erodible(n_solid)          !< erodible thickness
    REAL(wp), INTENT(IN) :: dt

    REAL(wp), INTENT(OUT) :: erosion_term(n_solid)     !< erosion term
    REAL(wp), INTENT(OUT) :: deposition_term(n_solid)  !< deposition term
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

    REAL(wp) :: Tc           !< temperature of carrier pphase [K]

    REAL(wp) :: alpha1       !< viscosity of continuous phase [kg m-1 s-1]
    REAL(wp) :: fluid_visc
    REAL(wp) :: inv_kin_visc
    REAL(wp) :: rhoc
    REAL(wp) :: expA , expB
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

    REAL(wp) :: T_liquid
    REAL(wp) :: T_boiling
    REAL(wp) :: sp_latent_heat
    REAL(wp) :: sp_heat_liq_water
    REAL(wp) :: mass_vap_rate

    REAL(wp) :: pore_pressure_term
    REAL(wp) :: f_inhibit
    REAL(wp) :: vel_loss_gas

   !> permeability calculationn
    REAL(wp) :: diam_sauter ! < Sauter mean diameter
    REAL(wp) :: sphericity_mean ! < Surface-area-weighted mean sphericity
   !> Inhibit factor for friction term
    REAL(wp) :: x  !< local variable for f_inhibit calculation
    REAL(wp) :: s !< local variable for f_inhibit calculation
    INTEGER :: n !< local variable for f_inhibit summation
    REAL(wp) :: tt, term, width !< local variables for f_inhibit
    REAL(wp) :: alpha_trans_dynamic !< dynamic transition solid fraction


    REAL(wp) :: dyn_visc_c
    REAL(wp) :: rho_c
    REAL(wp) :: hydraulic_permeability_local
    REAL(wp) :: kin_visc_c_local

    erosion_term(1:n_solid) = 0.0_wp
    deposition_term(1:n_solid) = 0.0_wp
    continuous_phase_erosion_term = 0.0_wp
    continuous_phase_loss_term = 0.0_wp
    eqns_term(1:n_eqns) = 0.0_wp
    topo_term = 0.0_wp
    hydraulic_permeability_local = hydraulic_permeability
    kin_visc_c_local = kin_visc_c
    r_Zs = 0.0_wp
    r_exc_pore_pres = 0.0_wp
    f_inhibit = 1.0_wp
    pore_pressure_term = 0.0_wp

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

    IF ( alphas_tot .LE. alphastot_min ) RETURN

    r_T = qpj(4)

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,r_sp_heat_c,             &
         r_sp_heat_mix)

    IF ( rheology_model .EQ. 4 ) THEN

       ! alpha1 here has units: kg m-1 s-1
       ! in Table 2 from O'Brien 1988, the values reported have different
       ! units ( poises). 1poises = 0.1 kg m-1 s-1

       ! convert from Kelvin to Celsius
       Tc = r_T - 273.15_wp

       ! the dependance of viscosity on temperature is modeled with the
       ! equation presented at:
       ! https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781118131473.app3
       !
       ! In addition, we use a reference value provided in input at a
       ! reference temperature. This value is used to scale the equation
       IF ( REAL(Tc) .LT. 20.0_wp ) THEN

          expA = 1301.0_wp / ( 998.333_wp + 8.1855_wp * ( Tc - 20.0_wp )        &
               + 0.00585_wp * ( Tc - 20.0_wp )**2 ) - 1.30223_wp

          alpha1 = alpha1_coeff * 1.0E-3_wp * 10.0_wp**expA

       ELSE

          expB = ( 1.3272_wp * ( 20.0_wp - Tc ) - 0.001053_wp *                 &
               ( Tc - 20.0_wp )**2 ) / ( Tc + 105.0_wp )

          alpha1 = alpha1_coeff * 1.002E-3_wp * 10.0_wp**expB

       END IF

       ! Fluid dynamic viscosity [kg m-1 s-1]
       fluid_visc = alpha1 * EXP( beta1 * alphas_tot )
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

          deposition_term(i_solid) = r_alphas(i_solid) * settling_vel

          IF ( rheology_model .NE. 4 ) THEN

             ! Michaels and Bolger (1962) sedimentation correction accounting
             ! for hindered settling due to the presence of particles
             deposition_term(i_solid) = deposition_term(i_solid) *              &
                  ( 1.0_wp - MIN( 1.0_wp , alphas_tot / alpha_max ) )**hind_exp

          END IF

          ! limit the deposition (cannot remove more than particles present
          ! in the flow)
          deposition_term(i_solid) = MIN( deposition_term(i_solid) ,            &
               r_h * r_alphas(i_solid) / dt )

          IF ( deposition_term(i_solid) .LT. 0.0_wp ) THEN

             WRITE(*,*) 'eval_erosion_dep_term'
             WRITE(*,*) 'deposition_term(i_solid)',deposition_term(i_solid)
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

         ! compute f_inhibit based on the selected model
         select case (trim(adjustl(f_inhibit_mode)))
         case ('OFF')
            f_inhibit = 1.0_wp ! no inhibition (pure diffusion)

         case ('STEP')
            !-------------------------------------------------------------------!
            ! Calculate f_inhibit based on SOLID FRACTION
            ! for STEP function case
            !-------------------------------------------------------------------!
            if ( alphas_tot .LT. maximum_solid_packing ) then
               f_inhibit = 1.0_wp
            else
               f_inhibit = 0.0_wp
            end if

         case ('STATIC')
            !-------------------------------------------------------------------!
            ! Calculate f_inhibit based on SOLID FRACTION
            ! for REAL case
            ! Generalized smoothstep S_N(x) for real input, result flipped [1→0]
            !-------------------------------------------------------------------!
            ! Validate N_inh read from namelist (should be non-negative integer)
            if (N_inh < 0) then
               WRITE(*,*) 'ERROR: N_inh must be >= 0. N_inh=', N_inh
               STOP
            end if
            ! use solid volume fraction source for x
            x = alphas_tot
            ! Normalize x to [0,1]
            width = maximum_solid_packing - alpha_trans ! right limit - left limit
            tt = (x - alpha_trans) / width
            if (tt < 0.0_wp) tt = 0.0_wp
            if (tt > 1.0_wp) tt = 1.0_wp
            ! Initialize smoothstep sum
            s = 0.0_wp
            ! call the precomputed pascal triangle coefficients
            do n = 0, N_inh
               term = pascal_coeff(n)  * tt**(N_inh + n + 1.0_wp)
               s = s + term
            end do
            ! Flip the smoothstep to start at 1 → 0
            s = 1.0_wp - s
            ! Assign to f_inhibit (convert real-like complex value)
            f_inhibit = s

         case ('DYNAMIC')
            !-------------------------------------------------------------------!
            ! Calculate f_inhibit based on SOLID FRACTION and FLOW HEIGHT
            ! For a REAL case
            ! Generalized smoothstep S_N(x) for real input, result flipped [1→0]
            !-------------------------------------------------------------------!
            ! Calculate the alpha_trans_dynamic based on flow height
            IF (r_h .LT. EPSILON(1.0_wp)) THEN
               r_h = EPSILON(1.0_wp) ! avoid multiplying by something lower than working precision
            END IF

            alpha_trans_dynamic = 10.0_wp ** (-0.28_wp) * r_h ** 0.12_wp


            ! Correct for alpha_trans_dynamic exceeding maximum_solid_packing
            if ( alpha_trans_dynamic .GE. maximum_solid_packing ) then
                alpha_trans_dynamic = 0.99_wp * maximum_solid_packing
            end if
            ! Validate N_inh read from namelist (should be non-negative integer)
            if (N_inh < 0) then
               WRITE(*,*) 'ERROR: N_inh must be >= 0. N_inh=', N_inh
               STOP
            end if
            ! use solid volume fraction source for x
            x = alphas_tot
            ! Normalize x to [0,1]
            width = maximum_solid_packing - alpha_trans_dynamic ! right limit - left limit
            tt = (x - alpha_trans_dynamic) / width
            if (tt < 0.0_wp) tt = 0.0_wp
            if (tt > 1.0_wp) tt = 1.0_wp
            ! Initialize smoothstep sum
            s = 0.0_wp
            ! call the precomputed pascal triangle coefficients
            do n = 0, N_inh
               term = pascal_coeff(n)  * tt**(N_inh + n + 1.0_wp)
               s = s + term
            end do

           ! Flip the smoothstep to start at 1 → 0
            s = 1.0_wp - s
            ! Assign to f_inhibit (convert real-like complex value)
            f_inhibit = s

         end select
         ! ------------------------------------------------------------------- !

   ! -------------------------------------------------------------------- !
          ! Guard on alphas_tot: Carman-Kozeny below divides by
          ! alphas_tot**2, and the surface-area weights divide by a sum
          ! that vanishes with it, so a gas-only or numerically depleted
          ! cell would divide by zero.
          IF (dynamic_permeability_flag .AND.                                  &
              ( alphas_tot .GT. EPSILON(1.0_wp) )) THEN

            ! diameter and density of particles
            IF ( n_solid .GT. 1) THEN ! if more than one solid phase
               diam_sauter = sauter_diameter( r_alphas )
               ! Surface-area-weighted mean sphericity, consistent with the
               ! Sauter diameter used beside it: the surface area of class i
               ! per unit bulk volume scales as alpha_i/d_i, not
               ! alpha_i*d_i^2, which over-weights the coarse classes.
               sphericity_mean = DOT_PRODUCT( r_alphas / diam_s , sphericity_s ) / &
                                 SUM( r_alphas / diam_s )
            ELSE ! if only one solid phase
               diam_sauter = diam_s(1)
               sphericity_mean = sphericity_s(1)
            END IF
            ! calculating permeability using Carman-Kozeny equation with sphericity
            hydraulic_permeability_local = ( (1.0_wp - alphas_tot)**3.0_wp * &
                                      ( diam_sauter * sphericity_mean )**2.0_wp ) / &
                                      ( 150.0_wp * (alphas_tot)**2.0_wp )


          END IF
	! -------------------------------------------------------------------- !

	! -------------------------------------------------------------------- !
	 ! kinetic viscosity
          IF ( gas_flag .AND. sutherland_flag ) THEN

             dyn_visc_c = muRef_Suth * ( r_T / Tref_Suth )**1.5_wp *            &
                  ( Tref_Suth + S_mu ) / ( r_T + S_mu )

             rho_c = pres / ( sp_gas_const_a * r_T )
             kin_visc_c_local = dyn_visc_c / rho_c

          END IF
	! -------------------------------------------------------------------- !

        ! velocity of gas loss due to pore pressure gradient
        vel_loss_gas = hydraulic_permeability_local /                         &
               ( kin_visc_c_local * r_rho_c ) / MAX(1.0e-5_wp,r_h) *          &
               0.5_wp * pi_g *                                                &
               r_exc_pore_pres

	! add pore pressure driven gas loss term, inhibited by f_inhibit
         pore_pressure_term = vel_loss_gas * f_inhibit

         continuous_phase_loss_term = continuous_phase_loss_term +             &
               pore_pressure_term

       END IF

    END IF

    ! limit the loss accountaing for maximum solid packing
    continuous_phase_loss_term = MIN( continuous_phase_loss_term ,              &
         r_h * MAX( 0.0_wp , maximum_solid_packing - SUM(r_alphas) ) / dt )

    ! loss of continuous phase cannot
    continuous_phase_loss_term = MIN( continuous_phase_loss_term ,              &
         ( r_h *  ( 1.0_wp - SUM(r_alphas) ) - coeff_porosity * ( r_h *         &
         SUM(r_alphas) - dt * SUM( deposition_term(1:n_solid) ) ) ) / dt )


    IF ( erodible_deposit_flag ) T_erodible = r_T

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
    eqns_term(4) = - r_T * ( SUM( rho_s * sp_heat_s * deposition_term )        &
         + r_rho_c * r_sp_heat_c * continuous_phase_loss_term )                &
         + T_erodible * ( SUM( rho_s * sp_heat_s * erosion_term )              &
         + rho_c_sub * r_sp_heat_c * continuous_phase_erosion_term )           &
         + T_ambient * sp_heat_a * rho_a_amb * air_entr

    ! solid phase mass equation source term [kg m-2 s-1]:
    ! due to solid erosion and deposition
    eqns_term(idx_solidEqn_first:idx_solidEqn_last) = rho_s(1:n_solid) *        &
         ( erosion_term(1:n_solid) - deposition_term(1:n_solid) )

    IF ( gas_flag .AND. liquid_vaporization_flag .AND. ( B_zone .NE. 0 ) ) THEN

       ! gamma_steam = 0.114_wp
       T_liquid = 290.0_wp
       T_boiling = 373.15_wp
       sp_latent_heat = 2264705.0_wp
       sp_heat_liq_water = 4184.0_wp

       mass_vap_rate = gamma_steam * MAX(0.0_wp , r_T-T_boiling) * SUM( rho_s * &
            sp_heat_s * deposition_term ) / ( sp_heat_liq_water * ( T_boiling - &
            T_liquid ) + sp_latent_heat )

       eqns_term(1) = eqns_term(1) + mass_vap_rate
       eqns_term(4) = eqns_term(4) + mass_vap_rate * sp_heat_g(1) * T_boiling
       eqns_term(5+n_solid) = eqns_term(5+n_solid) + mass_vap_rate

    END IF

    IF ( stoch_transport_flag) THEN

       ! Evaluate the equation term related to the noise transport equation
       ! (if transport flag is false nothing is done)
       eqns_term(idx_stochEqn) = eqns_term(1) * r_Zs

    END IF

    IF ( pore_pressure_flag ) THEN

       ! Equation for q1*exc_pore_press
       eqns_term(idx_poreEqn) = eqns_term(1) * r_exc_pore_pres

    END IF

    ! erodible layer thickness source terms [m s-1]:
    ! due to erosion and deposition of solid+continuous phase
    topo_term = ( dep_tot - ers_tot ) / ( 1.0_wp - erodible_porosity )

    RETURN

  END SUBROUTINE eval_mass_exchange_terms


  !******************************************************************************
  !> \brief Internal boundary source fluxes
  !
  !> This subroutine evaluates the source terms at the interfaces when an
  !> internal radial source is present, as for a base surge. The terms are
  !> applied as boundary conditions, and thus they have the units of the
  !> physical variable qp
  !> \date 2019/12/01
  !> \param[in]     time         time
  !> \param[in]     vect_x       unit vector velocity x-component
  !> \param[in]     vect_y       unit vector velocity y-component
  !> \param[out]    source_bdry  source terms
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
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
