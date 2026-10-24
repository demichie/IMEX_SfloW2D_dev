!> \brief Common physical source assembly entry points.
!> \details Owns the explicit STEP source-rate bridge independently of the
!> spatial flux and hydrostatic interfaces in equation_terms_2d.
!> The bounded mass-exchange transaction and all inventory updates remain
!> the responsibility of step_source_executor_2d.
MODULE source_terms_2d
  USE mass_boundary_injection_2d
  USE mass_bottom_circular_injection_2d
  USE mass_bottom_fissural_injection_2d
  USE mass_radial_injection_2d
  USE mass_injection_coefficients_2d
  USE parameters_2d
  USE constitutive_parameters_2d
  USE state_conversion_2d, ONLY : mixt_var
  USE thermal_laws_2d, ONLY : eval_inlet_thermal_energy_source
  USE pore_pressure_laws_2d, ONLY : eval_pore_pressure_inlet_source
  USE geometry_2d, ONLY : pi_g, bottom_source_area_ratio
  USE mass_exchange_terms_2d, ONLY : eval_mass_exchange_terms
  USE rheology_2d, ONLY : curvature_acceleration
  USE diagnostics_2d, ONLY : fatal_error
  USE source_catalog_2d, ONLY : SOURCE_METHOD_EXPLICIT, SOURCE_CADENCE_STAGE
  USE source_registry_2d, ONLY : source_group_type, group_contains_name
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_step_source_explicit, eval_stage_source_explicit
  PUBLIC :: eval_geometry_curvature_explicit
  PUBLIC :: eval_geometry_curvature_semi_implicit
  PUBLIC :: eval_geometry_curvature_implicit
  PUBLIC :: solve_geometry_curvature_analytic
  PUBLIC :: bottom_coeff_step, bottom_vel_step, lateral_coeff_step, fissure_coeff_step
  PUBLIC :: eval_boundary_injection_explicit
  PUBLIC :: eval_boundary_injection_semi_implicit
  PUBLIC :: eval_boundary_injection_implicit
  PUBLIC :: solve_boundary_injection_analytic
  PUBLIC :: eval_bottom_circular_injection_explicit
  PUBLIC :: eval_bottom_circular_injection_semi_implicit
  PUBLIC :: eval_bottom_circular_injection_implicit
  PUBLIC :: solve_bottom_circular_injection_analytic
  PUBLIC :: eval_bottom_fissural_injection_explicit
  PUBLIC :: eval_bottom_fissural_injection_semi_implicit
  PUBLIC :: eval_bottom_fissural_injection_implicit
  PUBLIC :: solve_bottom_fissural_injection_analytic
  PUBLIC :: eval_radial_injection_explicit
  PUBLIC :: eval_radial_injection_semi_implicit
  PUBLIC :: eval_radial_injection_implicit
  PUBLIC :: solve_radial_injection_analytic

  ! Shared source-window coefficients updated once before the RK stages.
  ! They remain private to the physics source layer, but are imported by
  ! the legacy window initializer until its migration.
CONTAINS

  SUBROUTINE eval_stage_source_explicit( Bprimej_x, Bprimej_y, Bsecondj_xx , Bsecondj_xy , &
       Bsecondj_yy, grav_coeff, qpj, expl_term, time, cell_fract_jk,             &
       cell_fissure_fract,                                                       &
       lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, cell_area_jk, group )

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
    TYPE(source_group_type), INTENT(IN) :: group

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

    IF (group%method/=SOURCE_METHOD_EXPLICIT .OR. &
         group%cadence/=SOURCE_CADENCE_STAGE .OR. &
         .NOT.ALLOCATED(group%source_names)) THEN
       CALL fatal_error('D17: eval_stage_source_explicit requires a validated explicit stage group')
    END IF

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

       IF ( curvature_term_flag .AND.                                      &
            group_contains_name(group,'geometry_curvature') .AND.     &
            r_h .GT. dry_thickness_tolerance ) THEN

          ! Hydrostatic pressure and the ordinary bed-slope source are already
          ! represented by the HP path. Curvature remains a separate,
          ! velocity-dependent cell source using the fitted Hessian and G_c.
          ! Below the shared dry threshold the auxiliary velocity is not a
          ! resolved physical velocity: do not square it in this source.
          ! This guard neither discards mass nor disables external injection.
          CALL eval_geometry_curvature_explicit( Bprimej_x, Bprimej_y,         &
               Bsecondj_xx, Bsecondj_xy, Bsecondj_yy, grav_coeff, q1,        &
               r_u, r_v, expl_term(2), expl_term(3) )

       END IF

       ! The hydrostatic path/source has no thermal-energy component.
       expl_term(4) = 0.0_wp

    END IF

    ! Each configured source owns exactly one explicit rate calculation.
    IF (group_contains_name(group,'bottom_circular_injection')) &
         CALL eval_bottom_circular_injection_explicit(qpj,time,cell_fract_jk,q1,expl_term)
    IF (group_contains_name(group,'bottom_fissural_injection')) &
         CALL eval_bottom_fissural_injection_explicit(qpj,time,cell_fissure_fract,q1,expl_term)
    IF (group_contains_name(group,'radial_injection')) &
         CALL eval_radial_injection_explicit(qpj,time,lat_arc_perim_jk, &
         lat_n_x_jk,lat_n_y_jk,cell_area_jk,q1,expl_term)

    RETURN

  END SUBROUTINE eval_stage_source_explicit











  PURE SUBROUTINE eval_geometry_curvature_explicit( Bprime_x, Bprime_y,         &
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

  END SUBROUTINE eval_geometry_curvature_explicit

  SUBROUTINE eval_geometry_curvature_semi_implicit(Bprime_x, Bprime_y, &
       Bsecond_xx, Bsecond_xy, Bsecond_yy, grav_coefficient, &
       mixture_mass, velocity_x, velocity_y, source_x, source_y)
    REAL(wp), INTENT(IN) :: Bprime_x, Bprime_y
    REAL(wp), INTENT(IN) :: Bsecond_xx, Bsecond_xy, Bsecond_yy
    REAL(wp), INTENT(IN) :: grav_coefficient
    REAL(wp), INTENT(IN) :: mixture_mass, velocity_x, velocity_y
    REAL(wp), INTENT(OUT) :: source_x, source_y
    CALL fatal_error('geometry_curvature: semi_implicit treatment is unsupported')
  END SUBROUTINE eval_geometry_curvature_semi_implicit

  SUBROUTINE eval_geometry_curvature_implicit(Bprime_x, Bprime_y, &
       Bsecond_xx, Bsecond_xy, Bsecond_yy, grav_coefficient, &
       mixture_mass, velocity_x, velocity_y, source_x, source_y)
    REAL(wp), INTENT(IN) :: Bprime_x, Bprime_y
    REAL(wp), INTENT(IN) :: Bsecond_xx, Bsecond_xy, Bsecond_yy
    REAL(wp), INTENT(IN) :: grav_coefficient
    REAL(wp), INTENT(IN) :: mixture_mass, velocity_x, velocity_y
    REAL(wp), INTENT(OUT) :: source_x, source_y
    CALL fatal_error('geometry_curvature: implicit treatment is unsupported')
  END SUBROUTINE eval_geometry_curvature_implicit

  SUBROUTINE solve_geometry_curvature_analytic(Bprime_x, Bprime_y, &
       Bsecond_xx, Bsecond_xy, Bsecond_yy, grav_coefficient, &
       mixture_mass, velocity_x, velocity_y, source_x, source_y)
    REAL(wp), INTENT(IN) :: Bprime_x, Bprime_y
    REAL(wp), INTENT(IN) :: Bsecond_xx, Bsecond_xy, Bsecond_yy
    REAL(wp), INTENT(IN) :: grav_coefficient
    REAL(wp), INTENT(IN) :: mixture_mass, velocity_x, velocity_y
    REAL(wp), INTENT(OUT) :: source_x, source_y
    CALL fatal_error('geometry_curvature: analytic treatment is unsupported')
  END SUBROUTINE solve_geometry_curvature_analytic

  !> \brief Evaluate the atomic explicit STEP mass-exchange rate.
  !> \details Returns rates only; it never mutates conservative state, bed or inventory.
  !> \param[in] qpj Cell-local primitive state.
  !> \param[in] B_zone Connected-water-region indicator.
  !> \param[in] B_prime_x Filtered bed x derivative.
  !> \param[in] B_prime_y Filtered bed y derivative.
  !> \param[in] erodible Per-class available substrate inventory.
  !> \param[in] dt Time-step length.
  !> \param[out] erosion_term Per-class substrate erosion rate.
  !> \param[out] deposition_term Per-class deposition rate.
  !> \param[out] continuous_phase_erosion_term Carrier entrainment rate.
  !> \param[out] continuous_phase_loss_term Limited carrier loss rate.
  !> \param[out] eqns_term Conservative equation source-rate vector.
  !> \param[out] topo_term Proposed bed-elevation rate.
  SUBROUTINE eval_step_source_explicit(qpj, B_zone, B_prime_x, B_prime_y, &
       erodible, dt, erosion_term, deposition_term, &
       continuous_phase_erosion_term, continuous_phase_loss_term, eqns_term, topo_term)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2), B_prime_x, B_prime_y
    REAL(wp), INTENT(IN) :: erodible(n_solid), dt
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(OUT) :: erosion_term(n_solid), deposition_term(n_solid)
    REAL(wp), INTENT(OUT) :: continuous_phase_erosion_term, continuous_phase_loss_term
    REAL(wp), INTENT(OUT) :: eqns_term(n_eqns), topo_term

    CALL eval_mass_exchange_terms(qpj, B_zone, B_prime_x, B_prime_y, &
         erodible, dt, erosion_term, deposition_term, &
         continuous_phase_erosion_term, continuous_phase_loss_term, eqns_term, topo_term)
  END SUBROUTINE eval_step_source_explicit
END MODULE source_terms_2d
