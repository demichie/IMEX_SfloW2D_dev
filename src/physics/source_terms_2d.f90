!> \brief Common physical source assembly entry points.
!> \details Owns STAGE/STEP explicit, semi-implicit, implicit and analytic
!> source assembly independently of spatial fluxes in equation_terms_2d.
!> The bounded mass-exchange transaction and all inventory updates remain
!> the responsibility of step_source_executor_2d.
MODULE source_terms_2d
    USE parameters_2d, ONLY : idx_h
  USE mass_boundary_injection_2d
  USE mass_injection_2d, ONLY : eval_mass_injection_explicit
  USE mass_injection_coefficients_2d
  USE parameters_2d
  USE constitutive_parameters_2d
  USE state_conversion_2d
  USE thermal_laws_2d, ONLY : eval_thermal_explicit, eval_thermal_implicit, &
       solve_thermal_analytic
  USE pore_pressure_laws_2d, ONLY : &
       eval_pore_pressure_dissipation_implicit, &
       eval_pore_pressure_dissipation_explicit, solve_pore_pressure_dissipation_analytic
  USE mass_exchange_2d, ONLY : eval_mass_exchange_terms
  USE rheology_2d, ONLY : curvature_acceleration, solve_rheology_analytic, &
       eval_rheology_implicit, eval_rheology_semi_implicit, &
       eval_rheology_explicit
  USE diagnostics_2d, ONLY : fatal_error
  USE source_catalog_2d, ONLY : SOURCE_METHOD_NEWTON, SOURCE_ROLE_IMPLICIT, &
       SOURCE_ROLE_SEMI_IMPLICIT, canonical_rheology_source_name, rheology_name_has_role
  USE source_catalog_2d, ONLY : SOURCE_METHOD_EXPLICIT, SOURCE_METHOD_ANALYTIC, SOURCE_CADENCE_STAGE
  USE source_registry_2d, ONLY : source_group_type, group_contains_name
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_source_explicit, eval_cell_source_explicit
  PUBLIC :: cell_analytic_source_context
  ! Source evaluation is dispatched through one public generic: STAGE
  ! carries injection geometry; STEP carries atomic exchange rates/inventories.
  INTERFACE eval_source_explicit
     MODULE PROCEDURE eval_stage_source_explicit, eval_step_source_explicit
  END INTERFACE
  PUBLIC :: eval_source_implicit, eval_source_semi_implicit, eval_source_analytic
  PUBLIC :: analytic_term_context_type
  PUBLIC :: eval_geometry_curvature_explicit
  PUBLIC :: eval_geometry_curvature_semi_implicit
  PUBLIC :: eval_geometry_curvature_implicit
  PUBLIC :: solve_geometry_curvature_analytic
  PUBLIC :: bottom_coeff_step, bottom_vel_step, lateral_coeff_step, fissure_coeff_step
  PUBLIC :: eval_boundary_injection_explicit
  PUBLIC :: eval_boundary_injection_semi_implicit
  PUBLIC :: eval_boundary_injection_implicit
  PUBLIC :: solve_boundary_injection_analytic
  ! Cell-local geometry and stochastic context for analytic source solves.
  TYPE :: analytic_term_context_type
    REAL(wp) :: Bprime_x=0.0_wp, Bprime_y=0.0_wp
    REAL(wp) :: Bsecond_xx=0.0_wp, Bsecond_xy=0.0_wp, Bsecond_yy=0.0_wp
    REAL(wp) :: grav_coeff=0.0_wp, Zij=0.0_wp
  END TYPE analytic_term_context_type
CONTAINS
  ! Shared cell adapters: STAGE and STEP use identical physical geometry.
  SUBROUTINE eval_cell_source_explicit(j,k,qp_cell,time,group,rate)
    USE geometry_2d, ONLY : B_prime_x_geom,B_prime_y_geom, &
         B_second_xx_geom,B_second_xy_geom,B_second_yy_geom,grav_coeff, &
         cell_source_fractions,cell_fissure_fractions,cell_arc_perim, &
         cell_arc_n_x,cell_arc_n_y,dx,dy
    INTEGER, INTENT(IN) :: j,k
    REAL(wp), INTENT(IN) :: qp_cell(n_vars+2),time
    TYPE(source_group_type), INTENT(IN) :: group
    REAL(wp), INTENT(OUT) :: rate(n_eqns)
    CALL eval_source_explicit(B_prime_x_geom(j,k),B_prime_y_geom(j,k), &
         B_second_xx_geom(j,k),B_second_xy_geom(j,k), &
         B_second_yy_geom(j,k),grav_coeff(j,k),qp_cell,rate,time, &
         cell_source_fractions(j,k),cell_fissure_fractions(j,k,:), &
         cell_arc_perim(j,k),cell_arc_n_x(j,k),cell_arc_n_y(j,k), &
         dx*dy,group)
  END SUBROUTINE eval_cell_source_explicit

  FUNCTION cell_analytic_source_context(j,k,Zij) RESULT(context)
    USE geometry_2d, ONLY : B_prime_x_geom,B_prime_y_geom, &
         B_second_xx_geom,B_second_xy_geom,B_second_yy_geom,grav_coeff
    INTEGER, INTENT(IN) :: j,k
    REAL(wp), INTENT(IN) :: Zij
    TYPE(analytic_term_context_type) :: context
    context%Bprime_x=B_prime_x_geom(j,k)
    context%Bprime_y=B_prime_y_geom(j,k)
    context%Bsecond_xx=B_second_xx_geom(j,k)
    context%Bsecond_xy=B_second_xy_geom(j,k)
    context%Bsecond_yy=B_second_yy_geom(j,k)
    context%grav_coeff=grav_coeff(j,k)
    context%Zij=Zij
  END FUNCTION cell_analytic_source_context


  RECURSIVE SUBROUTINE eval_source_analytic(group,q_before,context,interval,q_after,status)
    TYPE(source_group_type), INTENT(IN) :: group
    REAL(wp), INTENT(IN) :: q_before(n_vars),interval
    TYPE(analytic_term_context_type), INTENT(IN) :: context
    REAL(wp), INTENT(OUT) :: q_after(n_vars)
    INTEGER, INTENT(OUT) :: status
    REAL(wp) :: flow_h,flow_u,flow_v,flow_temperature,flow_density
    REAL(wp) :: flow_alphas(n_solid),flow_alphag(n_add_gas)
    REAL(wp) :: flow_red_grav,flow_p_dyn,flow_Zs,flow_pore_pressure
    TYPE(physical_state_type) :: physical_state
    REAL(wp) :: horizontal_speed,speed_ratio,curvature_accel,velocity_after(2)
    REAL(wp) :: heat_capacity,temperature_after,pressure_known,pressure_after
    LOGICAL :: supported
    CHARACTER(LEN=64) :: selected_model_name
    TYPE(source_group_type) :: component
    REAL(wp) :: component_state(n_vars)
    INTEGER :: part

    q_after=q_before
    status=-5
    IF (group%method/=SOURCE_METHOD_ANALYTIC) RETURN
    IF (.NOT.ALLOCATED(group%equation_ids)) RETURN
    ! A mixed group is one ordered sweep of complete physical operators.
    IF (SIZE(group%equation_ids)>1 .AND. &
         .NOT.(SIZE(group%equation_ids)==2 .AND. &
         ANY(group%equation_ids==2) .AND. ANY(group%equation_ids==3))) THEN
      component=group
      DO part=1,3
        SELECT CASE(part)
        CASE(1)
          IF (.NOT.ANY(group%equation_ids==2)) CYCLE
          component%equation_ids=[2,3]
        CASE(2)
          IF (.NOT.ANY(group%equation_ids==4)) CYCLE
          component%equation_ids=[4]
        CASE(3)
          IF (.NOT.pore_pressure_flag) CYCLE
          IF (.NOT.ANY(group%equation_ids==idx_pore)) CYCLE
          component%equation_ids=[idx_pore]
        END SELECT
        CALL eval_source_analytic(component,q_after,context,interval, &
             component_state,status)
        IF (status/=0) RETURN
        q_after=component_state
      END DO
      status=0
      RETURN
    END IF
    ! The group owns the analytic operator and its complete equation target.
    SELECT CASE(SIZE(group%equation_ids))
    CASE(2)
      IF (.NOT.(ANY(group%equation_ids==2) .AND. &
           ANY(group%equation_ids==3))) RETURN
      IF (rheology_flag) THEN

        CALL real_physical_state_from_q(q_before, physical_state)
        flow_h = physical_state%h
        flow_u = physical_state%u
        flow_v = physical_state%v
        flow_temperature = physical_state%T
        flow_density = physical_state%rho_m
        flow_alphas = physical_state%alphas
        flow_alphag = physical_state%alphag
        flow_red_grav = physical_state%red_grav
        flow_Zs = physical_state%Zs
        flow_pore_pressure = physical_state%exc_pore_pres
        horizontal_speed=SQRT(flow_u**2+flow_v**2)
        speed_ratio=1.0_wp
        IF (horizontal_speed>0.0_wp) THEN
          speed_ratio=SQRT(horizontal_speed**2+                           &
               MERGE((flow_u*context%Bprime_x+flow_v*context%Bprime_y)**2,&
               0.0_wp,slope_correction_flag))/horizontal_speed
        END IF
        curvature_accel=0.0_wp
        IF (curvature_term_flag) curvature_accel=curvature_acceleration(  &
             context%Bsecond_xx,context%Bsecond_xy,context%Bsecond_yy,    &
             flow_u,flow_v)
        velocity_after=[flow_u,flow_v]
        selected_model_name=group_rheology_name(group)
        IF (LEN_TRIM(selected_model_name)==0) THEN
          status=-7
          RETURN
        END IF
        CALL solve_rheology_analytic(selected_model_name, interval, physical_state, &
             speed_ratio, context%grav_coeff, curvature_accel, pore_pressure_flag, &
             stochastic_flag, context%Zij, velocity_after, supported)
        IF (.NOT.supported) THEN
          status=-7
          RETURN
        END IF
        q_after(2:3)=q_before(1)*velocity_after
      END IF
      status=0

    CASE(1)
      IF (group%equation_ids(1)==4) THEN
      IF (.NOT.thermal_cooling_flag) THEN
        status=-2
        RETURN
      END IF

      CALL real_physical_state_from_q(q_before, physical_state)
      flow_h = physical_state%h
      flow_u = physical_state%u
      flow_v = physical_state%v
      flow_temperature = physical_state%T
      flow_density = physical_state%rho_m
      flow_alphas = physical_state%alphas
      flow_alphag = physical_state%alphag
      flow_red_grav = physical_state%red_grav
      flow_Zs = physical_state%Zs
      flow_pore_pressure = physical_state%exc_pore_pres
      IF (flow_temperature<=0.0_wp) THEN
        status=-3
        RETURN
      END IF
      heat_capacity=q_before(4)/flow_temperature
      CALL solve_thermal_analytic('thermal_cooling',physical_state,interval, &
           heat_capacity,radiative_term_coeff,convective_term_coeff, &
           T_env,temperature_after)
      q_after(4)=heat_capacity*temperature_after
      status=0

      ELSE IF (group%equation_ids(1)==idx_pore .AND. pore_pressure_flag) THEN
      IF (q_before(1)>0.0_wp) THEN

        CALL real_physical_state_from_q(q_before, physical_state)
        flow_h = physical_state%h
        flow_u = physical_state%u
        flow_v = physical_state%v
        flow_temperature = physical_state%T
        flow_density = physical_state%rho_m
        flow_alphas = physical_state%alphas
        flow_alphag = physical_state%alphag
        flow_red_grav = physical_state%red_grav
        flow_Zs = physical_state%Zs
        flow_pore_pressure = physical_state%exc_pore_pres
        pressure_known=q_before(idx_pore)/q_before(1)
        CALL solve_pore_pressure_dissipation_analytic(pressure_known,interval,flow_h,  &
             flow_temperature,flow_alphas,pressure_after)
        q_after(idx_pore)=q_before(1)*pressure_after
      END IF
      status=0
      ELSE
        RETURN
      END IF
    END SELECT
  END SUBROUTINE eval_source_analytic


  !> Recover the rheology directly from a validated group's canonical sources.
  !> Numerical equation identifiers are not physical model identifiers.
  PURE FUNCTION group_rheology_name(group) RESULT(name)
    TYPE(source_group_type), INTENT(IN) :: group
    CHARACTER(LEN=64) :: name
    INTEGER :: i
    name=''
    IF (.NOT.ALLOCATED(group%source_names)) RETURN
    DO i=1,SIZE(group%source_names)
       IF (INDEX(TRIM(group%source_names(i)),'rheology_')/=1) CYCLE
       name=group%source_names(i)
       RETURN
    END DO
  END FUNCTION group_rheology_name

  PURE LOGICAL FUNCTION group_has_role(group,role) RESULT(found)
    TYPE(source_group_type), INTENT(IN) :: group
    INTEGER, INTENT(IN) :: role
    found=rheology_name_has_role(group_rheology_name(group),role)
  END FUNCTION group_has_role

  SUBROUTINE eval_source_implicit( Bprimej_x, Bprimej_y, Zij, c_qj,              &
       c_nh_term_impl, r_qj , r_nh_term_impl, group )

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
    TYPE(source_group_type), INTENT(IN), OPTIONAL :: group

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
    TYPE(complex_physical_state_type) :: physical_state

    COMPLEX(wp) :: red_grav
    LOGICAL :: selected_rheology, selected_thermal, selected_pore
    CHARACTER(LEN=64) :: selected_model_name

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       qj = c_qj

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       DO i = 1,n_vars

          qj(i) = CMPLX( r_qj(i),0.0_wp,wp )

       END DO

    ELSE

      WRITE(*,*) 'Constitutive, eval_source_implicit: problem with arguments'
       STOP

    END IF

    selected_model_name=canonical_rheology_source_name(rheology_model)
    IF (PRESENT(group)) THEN
       selected_model_name=group_rheology_name(group)
       selected_rheology=group_has_role(group,SOURCE_ROLE_IMPLICIT)
       selected_thermal=group_contains_name(group,'thermal_cooling')
       selected_pore=group_contains_name(group, &
            'pore_pressure_dissipation')
    ELSE
       selected_rheology=rheology_flag
       selected_thermal=thermal_cooling_flag
       selected_pore=pore_pressure_flag
    END IF

    ! initialize the source terms
    source_term(1:n_eqns) = CMPLX(0.0_wp,0.0_wp,wp)
    IF (selected_rheology .OR. selected_thermal .OR. selected_pore) THEN
       CALL physical_state_from_q(qj, physical_state)
       h = physical_state%h
       u = physical_state%u
       v = physical_state%v
       T = physical_state%T
       rho_m = physical_state%rho_m
       alphas = physical_state%alphas
       alphag = physical_state%alphag
       inv_rho_m = physical_state%inv_rhom
       Zs = physical_state%Zs
       exc_pore_pres = physical_state%exc_pore_pres
    END IF

    IF (rheology_flag .AND. selected_rheology) THEN

       red_grav = physical_state%red_grav

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

       CALL eval_rheology_implicit(selected_model_name, physical_state, &
            mod_vel, mod_vel_hor, mod_vel2, moving, stochastic_flag, Zij, source_term(2:3))
    END IF

    IF (selected_thermal) THEN
       CALL eval_thermal_implicit('thermal_cooling',physical_state, &
            radiative_term_coeff,convective_term_coeff,T_env,source_term(4))
    ENDIF

    IF ( pore_pressure_flag .AND. selected_pore ) &
         CALL eval_pore_pressure_dissipation_implicit(h, T,                       &
         rho_m, alphas, exc_pore_pres, source_term(idx_poreEqn))

    nh_term = source_term

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       c_nh_term_impl = nh_term

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       r_nh_term_impl = REAL( nh_term )

    END IF

    RETURN

  END SUBROUTINE eval_source_implicit

  SUBROUTINE eval_source_semi_implicit( Bprimej_x , Bprimej_y , Bsecondj_xx ,    &
       Bsecondj_xy , Bsecondj_yy , grav_coeff , qpj , semi_implicit_term ,     &
       Zj, group )

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

    REAL(wp), INTENT(OUT) :: semi_implicit_term(n_eqns)
    TYPE(source_group_type), INTENT(IN), OPTIONAL :: group

    REAL(wp) :: source_term(n_eqns)

    REAL(wp) :: mod_vel
    REAL(wp) :: mod_hor_vel

    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity
    REAL(wp) :: r_w               !< real_value z-velocity
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(wp) :: r_alphag(n_add_gas) !< real-value add. gas volume fractions
    REAL(wp) :: r_red_grav
    REAL(wp) :: r_rho_c
    REAL(wp) :: r_Ri
    REAL(wp) :: r_sp_heat_c
    REAL(wp) :: r_sp_heat_mix

    REAL(wp) :: exc_pore_pres       !< excess pore pressure
    TYPE(physical_state_type) :: physical_state
    LOGICAL :: selected_rheology
    CHARACTER(LEN=64) :: selected_model_name

    semi_implicit_term(1:n_eqns) = 0.0_wp

    selected_model_name=canonical_rheology_source_name(rheology_model)
    IF (PRESENT(group)) THEN
       selected_model_name=group_rheology_name(group)
       selected_rheology=group_has_role(group,SOURCE_ROLE_SEMI_IMPLICIT)
    ELSE
       selected_rheology=rheology_flag
    END IF

    ! A dry cell carries no friction source.
    IF ( qpj(idx_h) .LE. EPSILON(1.0_wp) ) RETURN

    ! initialize and evaluate the forces terms
    source_term(1:n_eqns) = 0.0_wp

    rheology_if:IF (rheology_flag .AND. selected_rheology) THEN

       r_h = qpj(idx_h)
       r_u = qpj(idx_u)
       r_v = qpj(idx_v)

       ! Every SI law below is zero at exact rest. Return before norms and
       ! divisions so optimized builds with arithmetic traps cannot speculate
       ! a reciprocal of zero from a later moving-state branch.
       IF ((r_u .EQ. 0.0_wp) .AND. (r_v .EQ. 0.0_wp)) RETURN

       CALL physical_state_from_qp(qpj, physical_state)
       r_alphas = physical_state%alphas
       r_alphag = physical_state%alphag
       r_rho_m = physical_state%rho_m
       r_red_grav = physical_state%red_grav


       IF ( slope_correction_flag ) THEN

          r_w = r_u * Bprimej_x + r_v * Bprimej_y

       ELSE

          r_w = 0.0_wp

       END IF

       mod_vel = SQRT( r_u**2 + r_v**2 + r_w**2 )

       mod_hor_vel = SQRT( r_u**2 + r_v**2 )

       exc_pore_pres = 0.0_wp
       IF (pore_pressure_flag) exc_pore_pres = qpj(idx_pore)
       CALL eval_rheology_semi_implicit(selected_model_name, physical_state, &
            mod_vel, mod_hor_vel, grav_coeff, curvature_term_flag, &
            Bsecondj_xx, Bsecondj_xy, Bsecondj_yy, pore_pressure_flag, &
            stochastic_flag, Zj, source_term(2:3))

    ENDIF rheology_if

    semi_implicit_term = source_term

    RETURN

  END SUBROUTINE eval_source_semi_implicit


  SUBROUTINE eval_stage_source_explicit( Bprimej_x, Bprimej_y, Bsecondj_xx , Bsecondj_xy , &
       Bsecondj_yy, grav_coeff, qpj, expl_term, time, cell_fract_jk,             &
       cell_fissure_fract,                                                       &
       lat_arc_perim_jk, lat_n_x_jk, lat_n_y_jk, cell_area_jk, group )

    USE parameters_2d, ONLY : vel_source , T_source , xs_source , xg_source,    &
         time_param , bottom_radial_source_flag,                    &
         pore_pressure_flag , pore_pres_fract ,                                &
         n_intervals , t_intervals , vel_intervals ,                            &
         radial_source_flag , h_source, bottom_fissural_source_flag,            &
         n_fissures, linear_vel_fissures, T_fissures, time_param_fissures,      &
         xs_fissures, xg_fissures



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


    TYPE(physical_state_type) :: physical_state

    REAL(wp) :: exc_pore_pres

    REAL(wp) :: q1

    REAL(wp) :: rheology_rate(2), thermal_rate, pressure_rate
    REAL(wp) :: injection_rate(n_eqns)
    REAL(wp) :: solid_volumes(n_solid), gas_volumes(n_add_gas)
    REAL(wp) :: speed, horizontal_speed, vertical_velocity, speed_squared
    REAL(wp) :: inverse_density, pressure_excess
    CHARACTER(LEN=64) :: active_rheology


    expl_term(1:n_eqns) = 0.0_wp
    physical_state=physical_state_type()
    q1 = 0.0_wp


      IF ( ( qpj(idx_h) .LE. EPSILON(1.0_wp) ) .AND. ( cell_fract_jk .EQ. 0.0_wp )    &
             .AND. ( lat_arc_perim_jk .EQ. 0.0_wp ) .AND.                          &
             ALL(cell_fissure_fract .EQ. 0.0_wp) ) RETURN

    ! Gravity terms - only meaningful when the cell has mass.
    IF ( qpj(idx_h) .GT. EPSILON(1.0_wp) ) THEN

       r_h = qpj(idx_h)
       r_u = qpj(idx_u)
       r_v = qpj(idx_v)

       CALL physical_state_from_qp(qpj, physical_state)
       r_rho_m = physical_state%rho_m
       r_red_grav = physical_state%red_grav

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

       ! Source rates are evaluated at the same unmodified stage state.
       ! The construction of EXPLICIT stays inside each physical module.
       active_rheology=group_rheology_name(group)
       IF (rheology_flag .AND. LEN_TRIM(active_rheology)>0 .AND. &
            r_h>dry_thickness_tolerance) THEN
          solid_volumes = physical_state%alphas
          gas_volumes = physical_state%alphag
          horizontal_speed=SQRT(r_u**2+r_v**2)
          vertical_velocity=0.0_wp
          IF (slope_correction_flag) vertical_velocity=r_u*Bprimej_x+r_v*Bprimej_y
          speed_squared=horizontal_speed**2+vertical_velocity**2
          speed=SQRT(speed_squared)
          inverse_density=1.0_wp/r_rho_m
          pressure_excess=0.0_wp
          IF (pore_pressure_flag) pressure_excess=qpj(idx_pore)
          CALL eval_rheology_explicit(active_rheology, physical_state, &
               speed, horizontal_speed, speed_squared, grav_coeff, curvature_term_flag, &
               Bsecondj_xx, Bsecondj_xy, Bsecondj_yy, pore_pressure_flag, &
               speed>0.0_wp, stochastic_flag, 0.0_wp, rheology_rate)
          expl_term(2:3)=expl_term(2:3)+rheology_rate
        END IF
        IF (thermal_cooling_flag .AND. group_contains_name(group,'thermal_cooling')) THEN
          CALL eval_thermal_explicit('thermal_cooling',physical_state, &
               radiative_term_coeff,convective_term_coeff,T_env,thermal_rate)
          expl_term(4)=expl_term(4)+thermal_rate
       END IF
       IF (pore_pressure_flag .AND. &
            group_contains_name(group,'pore_pressure_dissipation')) THEN
          solid_volumes = physical_state%alphas
          gas_volumes = physical_state%alphag
          CALL eval_pore_pressure_dissipation_explicit(r_h,qpj(idx_T), &
               r_rho_m,solid_volumes,qpj(idx_pore),pressure_rate)
          expl_term(idx_poreEqn)=expl_term(idx_poreEqn)+pressure_rate
       END IF

    END IF

    CALL eval_mass_injection_explicit(physical_state,time,cell_fract_jk, &
         cell_fissure_fract,lat_arc_perim_jk,lat_n_x_jk,lat_n_y_jk,cell_area_jk, &
         group_contains_name(group,'bottom_circular_injection'), &
         group_contains_name(group,'bottom_fissural_injection'), &
         group_contains_name(group,'radial_injection'),injection_rate)
    expl_term=expl_term+injection_rate

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
  !> \param[out] carrier_erosion_term Carrier entrainment rate.
  !> \param[out] carrier_loss_term Limited carrier loss rate.
  !> \param[out] eqns_term Conservative equation source-rate vector.
  !> \param[out] topo_term Proposed bed-elevation rate.
  SUBROUTINE eval_step_source_explicit(qpj,B_zone,B_prime_x,B_prime_y, &
       erodible,dt,rates,eqns_term)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2),B_prime_x,B_prime_y
    REAL(wp), INTENT(IN) :: erodible(n_solid),dt
    INTEGER, INTENT(IN) :: B_zone
    TYPE(mass_exchange_rates_type), INTENT(INOUT) :: rates
    REAL(wp), INTENT(OUT) :: eqns_term(n_eqns)
    CALL eval_mass_exchange_terms(qpj,B_zone,B_prime_x,B_prime_y,erodible,dt,rates,eqns_term)
  END SUBROUTINE eval_step_source_explicit
END MODULE source_terms_2d
