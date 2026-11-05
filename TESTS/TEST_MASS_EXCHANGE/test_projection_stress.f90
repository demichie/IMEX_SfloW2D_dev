!> \brief N8-B prescribed-rate refinement and seeded repeated production geometry checks.
PROGRAM test_projection_stress
  USE parameters_2d
  USE constitutive_parameters_2d
  USE geometry_2d
  USE state_fields, ONLY : state_type,substrate_snapshot_type,conservative_mass_admissible
  USE domain_2d, ONLY : domain_type
  USE source_registry_2d, ONLY : source_registry_type,source_group_type
  USE source_catalog_2d, ONLY : SOURCE_METHOD_EXPLICIT
  USE lie_source_executor_2d, ONLY : execute_lie_source_groups,initialize_lie_source_executor, &
       finalize_lie_source_executor,lie_source_executor_type
  USE state_conversions, ONLY : conservative_to_primitive
  USE source_terms_2d, ONLY : limit_explicit_source_rates
  USE omp_lib
  USE, INTRINSIC :: iso_fortran_env, ONLY : int64, error_unit
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite,ieee_value,ieee_quiet_nan
  IMPLICIT NONE
  INTEGER, PARAMETER :: ns=2, nq=6, rough_steps=32
  REAL(wp), PARAMETER :: tolerance=2048.0_wp*EPSILON(1.0_wp)
  INTEGER :: nx,ny,dimension,level,step,unit,requested,actual,j,k
  INTEGER(int64) :: seed
  CHARACTER(LEN=16) :: arg
  REAL(wp), ALLOCATABLE :: q(:,:,:),qp(:,:,:),rate(:,:),nodes(:,:),weights(:,:)
  TYPE(domain_type) :: domain
  TYPE(state_type) :: state

  CALL get_command_argument(1,arg)
  READ(arg,*) requested
  actual=0
  !$OMP PARALLEL
  !$OMP SINGLE
  actual=omp_get_num_threads()
  !$OMP END SINGLE
  !$OMP END PARALLEL
  IF (actual/=requested) ERROR STOP 'incorrect actual team'
  WRITE(*,*) 'Actual OpenMP team:',actual
  OPEN(NEWUNIT=unit,FILE='projection_snapshot.bin',ACCESS='STREAM', &
       FORM='UNFORMATTED',STATUS='REPLACE')
  CALL initialize_physics
  CALL check_limiter_representability
  CALL check_taal_saturated_species_budget
  CALL check_late_group_recipient_budget
  CALL check_injection_before_settling_budget
  CALL check_preexisting_settling_budget
  CALL check_settling_drag_transition
  CALL check_settling_newton_roots
  CALL check_complete_physical_state
  DO dimension=1,2
     DO level=0,2
        nx=16*2**level; ny=1
        IF (dimension==2) ny=12*2**level
        CALL allocate_grid
        CALL reset_state(.FALSE.)
        IF (level==0) THEN
           CALL check_lie_exchange_transaction
           CALL check_lie_rejected_transaction
        END IF
        DO k=1,ny
           DO j=1,nx
              rate(j,k)=-0.001_wp*(2.0_wp+COS(2.0_wp*ACOS(-1.0_wp)*(REAL(j,wp)-0.5_wp)*dx))
              IF (dimension==2) rate(j,k)=-0.001_wp*(2.0_wp+ &
                   COS(2.0_wp*ACOS(-1.0_wp)*(REAL(j,wp)-0.5_wp)*dx)* &
                   COS(2.0_wp*ACOS(-1.0_wp)*(REAL(k,wp)-0.5_wp)*dy))
           END DO
        END DO
        CALL projection_update('smooth',1)
        CALL reset_state(.TRUE.)
        rate=0.0_wp
        WHERE(q(1,:,:)>0.0_wp) rate=-0.001_wp
        CALL projection_update('wet_dry',1)

        ! Signed, cell-scale random rates exercise the public production
        ! projector without pretending that a rough field must converge.
        CALL reset_state(.FALSE.)
        seed=104729_int64
        DO step=1,rough_steps
           CALL random_rates(.TRUE.)
           CALL projection_update('rough',step)
        END DO

        ! Independently seeded repeated geometry updates for the N8-B oracle.
        CALL reset_state(.FALSE.)
        seed=130363_int64
        DO step=1,64
           CALL random_rates(.TRUE.)
           CALL projection_update('repeated',step)
        END DO
        CALL release_grid
     END DO
  END DO
  CLOSE(unit)
  WRITE(*,*) 'PASS: Lie-compatible prescribed projection refinement and geometry checks'
CONTAINS

  !> The limiter must preserve strict positivity after the exact floating-point
  !> multiplication order used by the Lie executor (including rounded ratios).
  SUBROUTINE check_limiter_representability
    REAL(wp) :: reference(nq),rates(nq),candidate(nq),scale
    REAL(wp) :: qsolid,dt
    LOGICAL :: valid
    INTEGER :: sample
    dt=0.01_wp
    DO sample=1,1000
       reference=0.0_wp
       rates=0.0_wp
       reference(idx_totMassEqn)=1.4739E-4_wp
       qsolid=5.0E-5_wp*(1.0_wp+REAL(sample,wp)*EPSILON(1.0_wp))
       reference(idx_solidEqn_first)=qsolid
       rates(idx_solidEqn_first)=-NEAREST(qsolid/dt,1.0_wp)
       rates(idx_totMassEqn)=rates(idx_solidEqn_first)
       CALL limit_explicit_source_rates(reference,[REAL(wp)::],dt,rates, &
            [REAL(wp)::],scale,valid)
       IF (.NOT.valid) ERROR STOP 'Representable deposition limiter rejected valid state'
       candidate=reference+dt*(scale*rates)
       IF (.NOT.conservative_mass_admissible(candidate)) &
            ERROR STOP 'Deposition limiter returned an inadmissible candidate'
       IF (ANY(candidate(idx_solidEqn_first:idx_solidEqn_last)<0.0_wp)) &
            ERROR STOP 'Deposition limiter left negative solid mass'
    END DO
  END SUBROUTINE check_limiter_representability

  !> Exact conservative witness from EXAMPLE_TAAL group 4, cell (224,379).
  !> The unbuffered active species bound rounded to a negative 1.43E-20.
  !> Treat this as a cross-routine floating-point reproducibility regression.
  SUBROUTINE check_taal_saturated_species_budget
    REAL(wp) :: qref(nq),rates(nq),candidate(nq),scale,step_dt
    LOGICAL :: valid

    qref=0.0_wp
    rates=0.0_wp
    qref(idx_totMassEqn)=1.1179177442718284E-3_wp
    qref(idx_uEqn)=3.8556231966220479E-3_wp
    qref(idx_vEqn)=2.4512212438858446E-3_wp
    qref(idx_engyEqn)=4.0886928508675646E+2_wp
    qref(idx_solidEqn_first)=2.6224851927896702E-4_wp
    step_dt=2.6630953223459757E-2_wp
    rates(idx_totMassEqn)=3.5072004284224628E-1_wp
    rates(idx_uEqn)=-3.3963392791815233E-2_wp
    rates(idx_vEqn)=-2.1592304455133129E-2_wp
    rates(idx_engyEqn)=1.0285961528261482E+5_wp
    rates(idx_solidEqn_first)=-9.8475077883410835E-3_wp
    ! Force a saturated negative species rate to exercise the limiter.
    rates(idx_solidEqn_first)=NEAREST(rates(idx_solidEqn_first),-1.0_wp)
    CALL limit_explicit_source_rates(qref,[REAL(wp)::],step_dt,rates, &
         [REAL(wp)::],scale,valid)
    IF (.NOT.valid) ERROR STOP 'TAAL D19 rate rejected by limiter'
    candidate=qref+step_dt*(scale*rates)
    IF (.NOT.conservative_mass_admissible(candidate)) &
         ERROR STOP 'TAAL D19 limiter returned negative component'
    IF (candidate(idx_solidEqn_first)<=0.0_wp) &
         ERROR STOP 'TAAL D19 depleted species at floating-point boundary'
  END SUBROUTINE check_taal_saturated_species_budget

  !> A later Lie group consumes the state left by earlier groups, not
  !> the array held on entry to the entire source sweep.
  SUBROUTINE check_late_group_recipient_budget
    REAL(wp) :: stale(nq),recipient(nq),rates(nq),candidate(nq)
    REAL(wp) :: stale_scale,current_scale
    LOGICAL :: old_valid,current_valid

    stale=0.0_wp
    recipient=0.0_wp
    rates=0.0_wp
    stale(idx_totMassEqn)=1.0_wp
    stale(idx_solidEqn_first)=0.20_wp
    recipient(idx_totMassEqn)=0.50_wp
    recipient(idx_solidEqn_first)=0.05_wp
    rates(idx_totMassEqn)=-0.30_wp
    rates(idx_solidEqn_first)=-0.15_wp
    CALL limit_explicit_source_rates(stale,[REAL(wp)::],1.0_wp,rates, &
         [REAL(wp)::],stale_scale,old_valid)
    IF (.NOT.old_valid) ERROR STOP 'Late-group stale fixture failed'
    candidate=recipient+stale_scale*rates
    IF (conservative_mass_admissible(candidate)) &
         ERROR STOP 'Late-group stale budget did not reproduce rejection'

    CALL limit_explicit_source_rates(recipient,[REAL(wp)::],1.0_wp,rates, &
         [REAL(wp)::],current_scale,current_valid)
    IF (.NOT.current_valid) ERROR STOP 'Late-group recipient limiter failed'
    candidate=recipient+current_scale*rates
    IF (.NOT.conservative_mass_admissible(candidate)) &
         ERROR STOP 'Late-group recipient yielded inadmissible state'
    IF (current_scale>=stale_scale) &
         ERROR STOP 'Late-group limiter used a stale conservative budget'
  END SUBROUTINE check_late_group_recipient_budget

  !> An initially dry cell receives solids in group 0. With the stale
  !> pre-group reference, group 4 is artificially disabled; the actual
  !> recipient budget permits settling. Characterize this behavioral change
  !> independently of the state physical-state migration.
  SUBROUTINE check_injection_before_settling_budget
    REAL(wp) :: dry(nq),injected(nq),settling_rate(nq),candidate(nq)
    REAL(wp) :: stale_scale,recipient_scale
    LOGICAL :: stale_valid,recipient_valid

    dry=0.0_wp
    injected=0.0_wp
    settling_rate=0.0_wp
    injected(idx_totMassEqn)=1.0_wp
    injected(idx_solidEqn_first)=0.25_wp
    settling_rate(idx_totMassEqn)=-0.20_wp
    settling_rate(idx_solidEqn_first)=-0.15_wp

    CALL limit_explicit_source_rates(dry,[REAL(wp)::],1.0_wp,settling_rate, &
        [REAL(wp)::],stale_scale,stale_valid)
    IF (.NOT.stale_valid .OR. stale_scale/=0.0_wp) &
        ERROR STOP 'Dry pre-injection budget must suppress settling'

    CALL limit_explicit_source_rates(injected,[REAL(wp)::],1.0_wp,settling_rate, &
        [REAL(wp)::],recipient_scale,recipient_valid)
    IF (.NOT.recipient_valid .OR. recipient_scale<=0.0_wp) &
        ERROR STOP 'Current post-injection budget must allow settling'
    candidate=injected+recipient_scale*settling_rate
    IF (.NOT.conservative_mass_admissible(candidate)) &
        ERROR STOP 'Post-injection settling candidate inadmissible'
  END SUBROUTINE check_injection_before_settling_budget

  !> Preserve legacy Lie source splitting without an unsafe stale recipient:
  !> newly injected solids are not available for same-step settling; on a
  !> previously wet cell the existing and current budgets must BOTH apply.
  SUBROUTINE check_preexisting_settling_budget
    REAL(wp) :: before(nq),recipient(nq),rates(nq),candidate(nq),scale
    LOGICAL :: valid

    before=0.0_wp
    recipient=0.0_wp
    rates=0.0_wp
    recipient(idx_totMassEqn)=1.0_wp
    recipient(idx_solidEqn_first)=0.25_wp
    rates(idx_totMassEqn)=-0.20_wp
    rates(idx_solidEqn_first)=-0.15_wp
    CALL limit_explicit_source_rates(recipient,[REAL(wp)::],1.0_wp,rates, &
        [REAL(wp)::],scale,valid,before)
    IF (.NOT.valid .OR. scale/=0.0_wp) &
        ERROR STOP 'Preexisting dry state permitted instant injection settling'
    candidate=recipient+scale*rates
    IF (.NOT.conservative_mass_admissible(candidate)) &
        ERROR STOP 'Dry preexisting budget yielded invalid recipient'

    before(idx_totMassEqn)=1.0_wp
    before(idx_solidEqn_first)=0.20_wp
    recipient(idx_totMassEqn)=0.6_wp
    recipient(idx_solidEqn_first)=0.05_wp
    CALL limit_explicit_source_rates(recipient,[REAL(wp)::],1.0_wp,rates, &
        [REAL(wp)::],scale,valid,before)
    IF (.NOT.valid .OR. scale<=0.0_wp .OR. scale>=1.0_wp) &
        ERROR STOP 'Preexisting/current two-budget limiter failed'
    candidate=recipient+scale*rates
    IF (.NOT.conservative_mass_admissible(candidate)) &
        ERROR STOP 'Two-budget limiter accepted negative solid mass'
  END SUBROUTINE check_preexisting_settling_budget

  !> \brief Configure the existing two-solid thermal-energy liquid exchange law.
  SUBROUTINE initialize_physics
    n_solid=ns; n_add_gas=0; n_stoch_vars=0; n_pore_vars=0
    kt_trace_reconstruction = 'froude_blend'
    n_vars=nq; n_eqns=nq
    CALL configure_canonical_state_indices()
    gas_flag=.FALSE.; liquid_flag=.TRUE.; rheology_flag=.FALSE.; rheology_model='none'
    stochastic_flag=.FALSE.; stoch_transport_flag=.FALSE.; pore_pressure_flag=.FALSE.
    gas_loss_flag=.FALSE.; liquid_vaporization_flag=.FALSE.; entrainment_flag=.FALSE.
    sutherland_flag=.FALSE.; slope_correction_flag=.TRUE.; topo_change_flag=.TRUE.
    n_fissures=0
    erodible_deposit_flag=.FALSE.; bottom_radial_source_flag=.FALSE.
    bottom_fissural_source_flag=.FALSE.; radial_source_flag=.FALSE.
    erosion_coeff=1000.0_wp; settling_flag=.FALSE.; alphastot_min=0.0_wp
    erodible_porosity=0.25_wp; coeff_porosity=erodible_porosity/(1.0_wp-erodible_porosity)
    grav=9.81_wp; rho_a_amb=1.2_wp; T_ambient=300.0_wp; T_erodible=317.0_wp
    rho_l=1000.0_wp; inv_rho_l=1.0_wp/rho_l; rho_c_sub=rho_l
    sp_heat_l=4180.0_wp; sp_heat_a=998.0_wp; sp_gas_const_a=287.05_wp
    pres=101300.0_wp; inv_pres=1.0_wp/pres
    kin_visc_c=1.0E-6_wp; hydraulic_permeability=0.0_wp; maximum_solid_packing=0.6_wp
    eps_sing=1.0E-8_wp; eps_sing4=eps_sing**4; verbose_level=0
    ALLOCATE(rho_s(ns),inv_rho_s(ns),sp_heat_s(ns),diam_s(ns),erodible_fract(ns))
    ALLOCATE(sp_heat_g(0),sp_gas_const_g(0),loss_rate)
    rho_s=[2000.0_wp,3000.0_wp]; inv_rho_s=1.0_wp/rho_s
    sp_heat_s=[900.0_wp,1100.0_wp]; diam_s=[1.0E-4_wp,2.0E-4_wp]
    erodible_fract=[0.6_wp,0.4_wp]; loss_rate=0.0_wp
  END SUBROUTINE initialize_physics

  !> \brief Check continuity of the settling drag closure at physical Re=1000.
  !> \details The threshold belongs to the converged Reynolds number rather
  !>          than the initial C_D=1 estimate. Exercise both supported drag
  !>          laws with safeguarded Newton on the low-Re branch.
  SUBROUTINE check_settling_drag_transition
    USE mass_exchange_deposition_2d, ONLY : settling_velocity
    REAL(wp), PARAMETER :: re_transition=1000.0_wp, delta=1.0E-6_wp
    REAL(wp), PARAMETER :: diam=1.0E-4_wp, solid_density=2000.0_wp, carrier_density=1.2_wp
    REAL(wp) :: cd_transition, const_part, initial_re_transition
    REAL(wp) :: inv_visc, w_below, w_at, w_above, re_below, w_expected
    INTEGER :: mode

    const_part=SQRT((4.0_wp/3.0_wp)*(solid_density/carrier_density-1.0_wp)*diam*grav)
    DO mode=1,2
       collective_settling_flag=(mode==2)
       IF (collective_settling_flag) THEN
          ! Use a positive, smooth collective correction at the transition.
          A_drag=0.01_wp
          B_drag=1.0_wp
       END IF
       cd_transition=24.0_wp*(1.0_wp+0.15_wp*re_transition**0.687_wp)/re_transition
       IF (collective_settling_flag) cd_transition=cd_transition+ &
            A_drag*(LOG(re_transition)-B_drag)
       initial_re_transition=re_transition*SQRT(cd_transition)
       inv_visc=initial_re_transition/(diam*const_part)

       w_below=settling_velocity(diam,solid_density,carrier_density, &
            inv_visc*(1.0_wp-delta))
       w_at=settling_velocity(diam,solid_density,carrier_density,inv_visc)
       w_above=settling_velocity(diam,solid_density,carrier_density, &
            inv_visc*(1.0_wp+delta))
       w_expected=const_part/SQRT(cd_transition)

       ! The limiting velocities must agree to within the 1e-6 fixed-point
       ! stopping criterion plus the O(delta) parameter perturbation.
       IF (ABS(w_below-w_at)>1.0E-5_wp*w_expected .OR. &
           ABS(w_above-w_at)>1.0E-5_wp*w_expected) &
           ERROR STOP 'discontinuous settling velocity at physical Re=1000'
       IF (ABS(w_at-w_expected)>64.0_wp*EPSILON(1.0_wp)*w_expected) &
           ERROR STOP 'incorrect settling drag at high-Re transition'

       re_below=diam*w_below*inv_visc*(1.0_wp-delta)
       IF (re_below>re_transition*(1.0_wp+1.0E-5_wp)) &
           ERROR STOP 'low-Re solver crossed settling transition'
       ! The high-Re branch must have a constant, positive drag.
       IF (ABS(w_above-w_expected)>64.0_wp*EPSILON(1.0_wp)*w_expected) &
           ERROR STOP 'incorrect constant high-Re settling drag'
    END DO
    collective_settling_flag=.FALSE.
    IF (settling_velocity(diam,carrier_density,carrier_density,inv_visc)/=0.0_wp) &
         ERROR STOP 'no buoyant settling at equal density'
  END SUBROUTINE check_settling_drag_transition

  !> \brief Recover prescribed Reynolds roots of both piecewise drag laws.
  !> \details Construct inverse viscosity from Re^2*C_D(Re)=Re0^2.
  !>          Thus the expected solution is independent of the implementation
  !>          of Newton, including the crossover and the constant-drag branch.
  SUBROUTINE check_settling_newton_roots
    USE mass_exchange_deposition_2d, ONLY : settling_velocity
    REAL(wp), PARAMETER :: diam=1.0E-4_wp, solid_density=2000.0_wp, carrier_density=1.2_wp
    REAL(wp), PARAMETER :: reference_re(*)=[1.0E-5_wp,1.0E-2_wp,1.0_wp,10.0_wp, &
         100.0_wp,500.0_wp,999.0_wp,1000.0_wp,1100.0_wp,3000.0_wp]
    REAL(wp) :: cd_transition, cd, reference, inv_visc, speed, actual_re, k
    INTEGER :: mode,i

    k=SQRT((4.0_wp/3.0_wp)*(solid_density/carrier_density-1.0_wp)*diam*grav)
    DO mode=1,2
       collective_settling_flag=(mode==2)
       A_drag=0.01_wp
       B_drag=1.0_wp
       cd_transition=24.0_wp*(1.0_wp+0.15_wp*1000.0_wp**0.687_wp)/1000.0_wp
       IF (collective_settling_flag) cd_transition=cd_transition+ &
            A_drag*(LOG(1000.0_wp)-B_drag)
       DO i=1,SIZE(reference_re)
          reference=reference_re(i)
          IF (reference<=1000.0_wp) THEN
             cd=24.0_wp*(1.0_wp+0.15_wp*reference**0.687_wp)/reference
             IF (collective_settling_flag) cd=cd+A_drag*(LOG(reference)-B_drag)
          ELSE
             cd=cd_transition
          END IF
          inv_visc=reference*SQRT(cd)/(diam*k)
          speed=settling_velocity(diam,solid_density,carrier_density,inv_visc)
          actual_re=diam*speed*inv_visc
          IF (ABS(actual_re-reference)>1.0E-8_wp*reference) THEN
             WRITE(*,*) 'FAIL: settling root',mode,reference,actual_re
             ERROR STOP 'incorrect safeguarded settling root'
          END IF
       END DO
    END DO
    collective_settling_flag=.FALSE.
  END SUBROUTINE check_settling_newton_roots

  !> \brief Verify that cached real and complex physical properties are complete.
  SUBROUTINE check_complete_physical_state
    USE state_mixture, ONLY : evaluate_mixture_phase_properties, evaluate_mixture_heat_capacity
    USE state_types, ONLY : physical_state_real, physical_state_complex, &
        initialize_physical_state_real, promote_physical_state
    USE state_physical, ONLY : build_physical_state_from_qp
    TYPE(physical_state_real) :: wet,dry,roundtrip
    TYPE(physical_state_complex) :: complex_state
    REAL(wp) :: primitive(n_vars),ri,rho_m,rho_c,red_grav,cp_c,cp_m
    REAL(wp) :: expected_cp,inv_rho_m,inv_rho_c,xc
    REAL(wp) :: xs(n_solid),xg(n_add_gas),alphas(n_solid),alphag(n_add_gas)
    REAL(wp) :: speed_squared,buoyancy_depth

    primitive=0.0_wp
    primitive(idx_h)=2.0_wp
    primitive(idx_u)=3.0_wp
    primitive(idx_v)=4.0_wp
    primitive(idx_T)=320.0_wp
    primitive(idx_solid_first)=0.15_wp
    primitive(idx_solid_first+1)=0.05_wp
    CALL initialize_physical_state_real(wet)
    CALL initialize_physical_state_real(dry)
    CALL build_physical_state_from_qp(primitive,wet)
    xs=primitive(idx_solid_first:idx_solid_last)
    xg=primitive(idx_add_gas_first:idx_add_gas_last)
    CALL evaluate_mixture_phase_properties(primitive(idx_T),xs,xg, &
        rho_m,inv_rho_m,rho_c,inv_rho_c,alphas,alphag)
    CALL evaluate_mixture_heat_capacity(xs,xg,cp_m,cp_c)
    red_grav=(rho_m-rho_a_amb)*inv_rho_m*grav
    speed_squared=primitive(idx_u)**2+primitive(idx_v)**2
    buoyancy_depth=red_grav*primitive(idx_h)
    ri=0.0_wp
    IF (speed_squared>0.0_wp) THEN
        IF (buoyancy_depth/1.0E15_wp>speed_squared) THEN
            ri=1.0E15_wp
        ELSE
            ri=MIN(1.0E15_wp,buoyancy_depth/speed_squared)
        END IF
    END IF
    IF (ABS(wet%rho_c-rho_c)>tolerance*MAX(1.0_wp,ABS(rho_c))) &
        ERROR STOP 'physical carrier density mismatch'
    IF (ABS(wet%inv_rho_c-1.0_wp/rho_c)>tolerance) &
        ERROR STOP 'physical inverse carrier density mismatch'
    IF (ABS(wet%cp_c-cp_c)>tolerance*cp_c .OR. &
        ABS(wet%cp_m-cp_m)>tolerance*cp_m) &
        ERROR STOP 'physical heat capacity mismatch'
    IF (ABS(wet%Ri-ri)>tolerance*MAX(1.0_wp,ABS(ri))) &
        ERROR STOP 'physical Richardson mismatch'
    IF (ABS(wet%rho_m-rho_m)>tolerance*rho_m) &
        ERROR STOP 'physical mixture density mismatch'

    CALL promote_physical_state(wet,complex_state)
    roundtrip%rho_c=REAL(complex_state%rho_c,wp)
    roundtrip%inv_rho_c=REAL(complex_state%inv_rho_c,wp)
    roundtrip%cp_c=REAL(complex_state%cp_c,wp)
    roundtrip%cp_m=REAL(complex_state%cp_m,wp)
    roundtrip%Ri=REAL(complex_state%Ri,wp)
    IF (roundtrip%rho_c/=wet%rho_c .OR. roundtrip%inv_rho_c/=wet%inv_rho_c .OR. &
        roundtrip%cp_c/=wet%cp_c .OR. roundtrip%cp_m/=wet%cp_m .OR. roundtrip%Ri/=wet%Ri) &
        ERROR STOP 'physical real-complex roundtrip mismatch'

    primitive=0.0_wp
    CALL build_physical_state_from_qp(primitive,dry)
    expected_cp=sp_heat_l
    IF (.NOT.liquid_flag) expected_cp=sp_heat_a
    IF (dry%rho_c/=rho_a_amb .OR. dry%inv_rho_c/=1.0_wp/rho_a_amb .OR. &
        dry%cp_c/=expected_cp .OR. dry%cp_m/=expected_cp .OR. dry%Ri/=0.0_wp) &
        ERROR STOP 'physical dry state properties mismatch'
  END SUBROUTINE check_complete_physical_state

  !> \brief Allocate the rectangular production grid and independent oracle buffers.
  SUBROUTINE allocate_grid
    comp_cells_x=nx; comp_cells_y=ny; comp_interfaces_x=nx+1; comp_interfaces_y=ny+1
    dx=1.0_wp/REAL(nx,wp); dy=1.0_wp/REAL(ny,wp); cell_size=dx
    ALLOCATE(q(nq,nx,ny),qp(nq+2,nx,ny),rate(nx,ny),nodes(nx+1,ny+1),weights(nx+1,ny+1))
    ALLOCATE(B_vertex(nx+1,ny+1),B_cent(nx,ny),B_face_x(nx+1,ny),B_face_y(nx,ny+1))
    ALLOCATE(B_prime_x_geom(nx,ny),B_prime_y_geom(nx,ny))
    ALLOCATE(B_second_xx_geom(nx,ny),B_second_xy_geom(nx,ny),B_second_yy_geom(nx,ny))
    ALLOCATE(grav_coeff(nx,ny),grav_coeff_stag_x(nx+1,ny),grav_coeff_stag_y(nx,ny+1))
    ALLOCATE(state%deposit(nx,ny,ns),state%erosion(nx,ny,ns), &
         state%erodible(ns,nx,ny),state%lava_deposit(nx,ny))
    ALLOCATE(state%deposit_tot(nx,ny),state%erosion_tot(nx,ny))
    state%deposit_tot=0.0_wp;state%erosion_tot=0.0_wp
    ALLOCATE(B_zone(nx,ny),cell_source_fractions(nx,ny))
    ALLOCATE(cell_fissure_fractions(nx,ny,n_fissures))
    ALLOCATE(cell_arc_perim(nx,ny),cell_arc_n_x(nx,ny),cell_arc_n_y(nx,ny))
    cell_fissure_fractions=0.0_wp
    cell_arc_perim=0.0_wp;cell_arc_n_x=0.0_wp;cell_arc_n_y=0.0_wp
    CALL domain%initialize
  END SUBROUTINE allocate_grid

  !> \brief Reset a known bed/flow and select a full or wet-plus-halo workset.
  !> \param[in] wet_dry Restrict flow to the left half or southwest quadrant.
  SUBROUTINE reset_state(wet_dry)
    LOGICAL, INTENT(IN) :: wet_dry
    INTEGER :: a,b,l
    REAL(wp) :: x,y,carrier_mass
    state%deposit=0.0_wp; state%erosion=0.0_wp; state%erodible=0.0_wp; state%lava_deposit=0.0_wp
    B_zone=0; cell_source_fractions=0.0_wp
    qp=0.0_wp; qp(4,:,:)=T_ambient
    l=0
    DO b=1,ny
       DO a=1,nx
          carrier_mass=rho_l*0.86_wp
          q(5:6,a,b)=rho_s*[0.10_wp,0.04_wp]
          q(1,a,b)=carrier_mass+SUM(q(5:6,a,b))
          q(2,a,b)=q(1,a,b); q(3,a,b)=-0.3_wp*q(1,a,b)
          q(4,a,b)=350.0_wp*(carrier_mass*sp_heat_l+DOT_PRODUCT(q(5:6,a,b),sp_heat_s))
          IF (wet_dry) THEN
             IF (a>nx/2 .OR. (dimension==2 .AND. b>ny/2)) q(:,a,b)=0.0_wp
             IF (a>nx/2+1 .OR. (dimension==2 .AND. b>ny/2+1)) CYCLE
          END IF
          l=l+1; domain%j_cent(l)=a; domain%k_cent(l)=b
       END DO
    END DO
    domain%solve_cells=l
    DO b=1,ny+1
       y=REAL(b-1,wp)*dy
       DO a=1,nx+1
          x=REAL(a-1,wp)*dx
          B_vertex(a,b)=2.0_wp+0.03_wp*x+0.002_wp*x*x
          IF (dimension==2) B_vertex(a,b)=B_vertex(a,b)-0.02_wp*y+0.001_wp*x*y+0.003_wp*y*y
       END DO
    END DO
    CALL refresh_topography_geometry
    CALL check_geometry
  END SUBROUTINE reset_state

  !> Lie group-0 exchange: verify complete conservative/inventory accounting.
  !> Use the production rate evaluator, executor and geometry update.
  SUBROUTINE check_lie_exchange_transaction
    TYPE(source_registry_type) :: registry
    TYPE(lie_source_executor_type) :: executor
    TYPE(source_group_type) :: groups(1)
    REAL(wp) :: q_before(nq,nx,ny),qp_before(nq+2,nx,ny)
    REAL(wp) :: deposited_before(nx,ny,ns),eroded_before(nx,ny,ns)
    REAL(wp) :: bed_before(nx+1,ny+1),available_before(ns,nx,ny)
    REAL(wp) :: pressure,fluid_delta,solid_delta,inventory_delta
    INTEGER :: a,b,l
    LOGICAL :: failed,accepted
    CHARACTER(LEN=256) :: message,error
    groups(1)%id=0
    groups(1)%method=SOURCE_METHOD_EXPLICIT
    ! A mixed explicit group exercises reuse of the same worker physical state.
    groups(1)%source_names=[CHARACTER(LEN=64)::'mass_exchange_erosion','geometry_curvature']
    groups(1)%equation_ids=[(l,l=1,nq)]
    CALL registry%initialize(nq,0,groups,error)
    IF (LEN_TRIM(error)>0) ERROR STOP 'Lie exchange registry initialization'
    DO b=1,ny
       DO a=1,nx
          CALL conservative_to_primitive(q(:,a,b),qp(:,a,b),accepted,pressure, &
               numerical_reconstruction=.TRUE.)
          IF (.NOT.accepted) ERROR STOP 'Lie exchange initial state conversion'
          state%erodible(:,a,b)=[0.002_wp,0.001_wp]
       END DO
    END DO
    q_before=q;qp_before=qp
    deposited_before=state%deposit;eroded_before=state%erosion
    bed_before=B_vertex;available_before=state%erodible
    CALL initialize_lie_source_executor(executor,registry)
    CALL execute_lie_source_groups(executor,state,q,qp,1.0_wp,domain, &
         q_original=q_before,qp_original=qp_before, &
         step_failed=failed,failure_message=message)
    IF (failed) THEN
       WRITE(error_unit,'(A)') 'Lie exchange transaction rejected: '//TRIM(message)
       WRITE(error_unit,'(A,I0,A,I0)') 'Lie exchange workset=',domain%solve_cells,', group=',groups(1)%id
       FLUSH(error_unit)
       ERROR STOP 'Lie exchange transaction rejected'
    END IF
    fluid_delta=SUM(q(1,:,:)-q_before(1,:,:))
    solid_delta=SUM(q(5:6,:,:)-q_before(5:6,:,:))
    inventory_delta=SUM(state%erosion-eroded_before)
    CALL require_small('Lie solid mass gain / state%erosion inventory', &
         ABS(solid_delta-DOT_PRODUCT(rho_s, &
         [SUM(state%erosion(:,:,1)-eroded_before(:,:,1)), &
          SUM(state%erosion(:,:,2)-eroded_before(:,:,2))])), &
         MAX(1.0_wp,ABS(solid_delta)))
    IF (ANY(state%erodible<0.0_wp)) ERROR STOP 'Lie erodible inventory negative'
    IF (ANY(state%erodible>available_before)) ERROR STOP 'Lie erodible inventory increased'
    DO l=1,ns
       CALL require_small('Lie erodible transfer balance', &
            MAXVAL(ABS(state%erodible(l,:,:)+state%erosion(:,:,l)-eroded_before(:,:,l) &
            -available_before(l,:,:))),1.0_wp)
    END DO
    IF (.NOT.ALL(ieee_is_finite(q)) .OR. .NOT.ALL(ieee_is_finite(B_vertex))) &
         ERROR STOP 'Lie exchange non-finite state'
    IF (fluid_delta<0.0_wp .OR. solid_delta<0.0_wp .OR. inventory_delta<0.0_wp) &
         ERROR STOP 'Lie erosion transfer changed sign'
    CALL check_geometry
    CALL finalize_lie_source_executor(executor)
    CALL registry%finalize()
    CALL reset_state(.FALSE.)
  END SUBROUTINE check_lie_exchange_transaction

  !> A later-cell failure invalidates the whole Lie time-step.
  !> The caller must restore substrate and geometry snapshots, while the
  !> partially updated fluid trial is never published.
  SUBROUTINE check_lie_rejected_transaction
    TYPE(source_registry_type) :: registry
    TYPE(lie_source_executor_type) :: executor
    TYPE(source_group_type) :: groups(1)
    TYPE(substrate_snapshot_type) :: substrate_before
    REAL(wp) :: q_before(nq,nx,ny),qp_before(nq+2,nx,ny)
    REAL(wp) :: q_trial(nq,nx,ny),qp_trial(nq+2,nx,ny)
    REAL(wp) :: qp_invalid(nq+2,nx,ny)
    REAL(wp) :: bed_before(nx+1,ny+1)
    REAL(wp) :: deposited_before(nx,ny,ns),eroded_before(nx,ny,ns)
    REAL(wp) :: erodible_before(ns,nx,ny),lava_before(nx,ny)
    REAL(wp) :: centers_before(nx,ny),slope_before(nx,ny)
    REAL(wp) :: pressure
    INTEGER :: a,b,l
    LOGICAL :: failed,accepted
    CHARACTER(LEN=256) :: error,message
    CALL reset_state(.FALSE.)
    groups(1)%id=0
    groups(1)%method=SOURCE_METHOD_EXPLICIT
    groups(1)%source_names=[CHARACTER(LEN=64)::'mass_exchange_erosion']
    groups(1)%equation_ids=[(l,l=1,nq)]
    CALL registry%initialize(nq,0,groups,error)
    IF (LEN_TRIM(error)>0) ERROR STOP 'Lie rejection registry initialization'
    DO b=1,ny
       DO a=1,nx
          CALL conservative_to_primitive(q(:,a,b),qp(:,a,b),accepted,pressure, &
               numerical_reconstruction=.TRUE.)
          IF (.NOT.accepted) ERROR STOP 'Lie rejection primitive conversion'
          state%erodible(:,a,b)=[0.002_wp,0.001_wp]
       END DO
    END DO
    q_before=q;qp_before=qp
    bed_before=B_vertex
    deposited_before=state%deposit;eroded_before=state%erosion
    erodible_before=state%erodible;lava_before=state%lava_deposit
    centers_before=B_cent;slope_before=B_prime_x_geom
    q_trial=q_before;qp_trial=qp_before
    ! Cell 2 is rejected after cell 1 has been considered for transfer.
    ! An independent immutable pre-spatial input avoids aliasing an INTENT(INOUT)
    ! candidate with the executor's INTENT(IN) GROUP-0 reference.
    qp_invalid=qp_before
    qp_invalid(:,2,1)=ieee_value(0.0_wp,ieee_quiet_nan)
    CALL state%snapshot_substrate(substrate_before)
    CALL begin_geometry_transaction
    CALL initialize_lie_source_executor(executor,registry)
    CALL execute_lie_source_groups(executor,state,q_trial,qp_trial,1.0_wp,domain, &
         q_original=q_before,qp_original=qp_invalid, &
         step_failed=failed,failure_message=message)
    IF (.NOT.failed) ERROR STOP 'Lie exchange rejection not reported'
    IF (LEN_TRIM(message)==0) ERROR STOP 'Lie rejection missing diagnostic'
    CALL state%restore_substrate(substrate_before)
    CALL rollback_geometry_transaction
    ! Fluid workspaces are not published by the Lie driver on rejection.
    q_trial=q_before;qp_trial=qp_before
    IF (ANY(q_trial/=q_before) .OR. ANY(qp_trial/=qp_before)) &
         ERROR STOP 'Rejected Lie fluid state changed'
    IF (ANY(B_vertex/=bed_before) .OR. ANY(B_cent/=centers_before) .OR. &
         ANY(B_prime_x_geom/=slope_before)) ERROR STOP 'Lie bed rollback mismatch'
    IF (ANY(state%deposit/=deposited_before) .OR. ANY(state%erosion/=eroded_before) .OR. &
         ANY(state%erodible/=erodible_before) .OR. ANY(state%lava_deposit/=lava_before)) &
         ERROR STOP 'Lie inventory rollback mismatch'
    CALL finalize_lie_source_executor(executor)
    CALL registry%finalize()
    CALL reset_state(.FALSE.)
  END SUBROUTINE check_lie_rejected_transaction

  !> \brief Generate reproducible rates with a portable nonoverflowing integer generator.
  !> \param[in] signed_rate Allow either deposition-like or erosion-like bed proposals.
  SUBROUTINE random_rates(signed_rate)
    LOGICAL, INTENT(IN) :: signed_rate
    INTEGER :: a,b
    REAL(wp) :: value
    DO b=1,ny
       DO a=1,nx
          seed=MODULO(16807_int64*seed,2147483647_int64)
          value=REAL(seed,wp)/2147483647.0_wp
          rate(a,b)=-0.001_wp*(0.5_wp+value)
          IF (signed_rate) rate(a,b)=0.002_wp*(value-0.5_wp)
       END DO
    END DO
  END SUBROUTINE random_rates

  !> \brief Independently scatter AC/4 contributions, including all dry cells and boundaries.
  SUBROUTINE reference_nodes
    INTEGER :: a,b
    nodes=0.0_wp; weights=0.0_wp
    DO b=1,ny
       DO a=1,nx
          nodes(a:a+1,b:b+1)=nodes(a:a+1,b:b+1)+0.25_wp*dx*dy*rate(a,b)
          weights(a:a+1,b:b+1)=weights(a:a+1,b:b+1)+0.25_wp*dx*dy
       END DO
    END DO
    nodes=nodes/weights
  END SUBROUTINE reference_nodes

  !> \brief Advance the public signed-rate projector and refresh actual production geometry.
  !> \param[in] name Scenario label.
  !> \param[in] iteration Update number.
  SUBROUTINE projection_update(name,iteration)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: iteration
    REAL(wp) :: actual_nodes(nx+1,ny+1),bed0(nx+1,ny+1),q0(nq,nx,ny)
    CALL reference_nodes
    bed0=B_vertex; q0=q
    CALL project_cell_field_to_vertices(rate,actual_nodes)
    CALL require_small('signed production projection',MAXVAL(ABS(actual_nodes-nodes)),1.0_wp)
    B_vertex=B_vertex+actual_nodes
    CALL refresh_topography_geometry
    CALL require_small('projector does not change flow',MAXVAL(ABS(q-q0)),1.0_wp)
    CALL check_update(name,iteration,bed0)
  END SUBROUTINE projection_update

  !> \brief Apply saturated erosion with an independent complete conservative transaction oracle.
  !> \param[in] name Scenario label.
  !> \param[in] iteration Update number.


  !> \brief Verify bed volume, front localization and all refreshed caches after each update.
  !> \param[in] name Scenario label.
  !> \param[in] iteration Update number.
  !> \param[in] bed0 Authoritative nodal bed before the transaction.
  SUBROUTINE check_update(name,iteration,bed0)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: iteration
    REAL(wp), INTENT(IN) :: bed0(nx+1,ny+1)
    REAL(wp) :: delta(nx+1,ny+1),geom(nx,ny),mismatch(nx,ny),metrics(6),x,y
    LOGICAL :: front_band
    INTEGER :: a,b
    delta=B_vertex-bed0
    CALL require_small('limited nodal proposal',MAXVAL(ABS(delta-nodes)),MAXVAL(ABS(bed0)))
    CALL require_small('weighted global nodal volume',ABS(SUM(weights*delta)-dx*dy*SUM(rate)),1.0_wp)
    DO b=1,ny
       DO a=1,nx
          geom(a,b)=0.25_wp*SUM(delta(a:a+1,b:b+1))
       END DO
    END DO
    mismatch=geom-rate
    metrics=[dx*dy*SUM(rate),dx*dy*SUM(geom),dx*dy*SUM(mismatch), &
         dx*dy*SUM(ABS(mismatch)),MAXVAL(ABS(mismatch)),0.0_wp]
    CALL require_small('signed global mismatch',ABS(metrics(3)),1.0_wp)
    IF (name=='wet_dry') THEN
       DO b=1,ny
          y=(REAL(b,wp)-0.5_wp)*dy
          DO a=1,nx
             x=(REAL(a,wp)-0.5_wp)*dx
             front_band=ABS(x-0.5_wp)<=dx
             IF (dimension==2) front_band=(front_band .AND. y<=0.5_wp+dy) .OR. &
                  (ABS(y-0.5_wp)<=dy .AND. x<=0.5_wp+dx)
             IF (.NOT.front_band) metrics(6)=MAX(metrics(6),ABS(mismatch(a,b)))
          END DO
       END DO
       CALL require_small('front-localized mismatch',metrics(6),1.0_wp)
       IF (ABS(mismatch(nx/2+1,1))<=tolerance) ERROR STOP 'missing dry geometric redistribution'
    END IF
    IF (metrics(4)<=tolerance) ERROR STOP 'mismatch stress not exercised'
    CALL check_geometry
    WRITE(*,'(A,1X,A,4(1X,I0),6(1X,ES24.16))') 'N8B_RECORD',name,dimension,nx,ny,iteration,metrics
    ! Capture intermediate fields too: final-only equality could hide a transient race.
    WRITE(unit) q,qp,B_vertex,B_cent,B_face_x,B_face_y,B_prime_x_geom,B_prime_y_geom, &
         B_second_xx_geom,B_second_xy_geom,B_second_yy_geom,grav_coeff, &
         grav_coeff_stag_x,grav_coeff_stag_y,state%deposit,state%erosion,state%erodible
  END SUBROUTINE check_update

  !> \brief Independently differentiate Q1 cell means with retained five-point LS kernels.
  SUBROUTINE check_geometry
    REAL(wp), PARAMETER :: c1(5)=[-2.0_wp,-1.0_wp,0.0_wp,1.0_wp,2.0_wp]
    REAL(wp), PARAMETER :: c2(5)=[2.0_wp,-1.0_wp,-2.0_wp,-1.0_wp,2.0_wp]
    REAL(wp) :: centers(nx,ny),bx,by,bxx,byy,bxy,G,bed_scale
    INTEGER :: a,b,ac,bc,s,t
    centers=0.25_wp*(B_vertex(1:nx,1:ny)+B_vertex(2:nx+1,1:ny)+ &
         B_vertex(1:nx,2:ny+1)+B_vertex(2:nx+1,2:ny+1))
    ! Freeze cancellation-aware operation scales before refinement runs:
    ! the epsilon multiplier is unchanged, while derivatives divide bed
    ! elevations by dx/dx**2. Do not scale by a nearly vanishing derivative.
    bed_scale=MAXVAL(ABS(centers))
    CALL require_small('Q1 centers',MAXVAL(ABS(B_cent-centers)),MAXVAL(ABS(centers)))
    CALL require_small('Q1 x faces',MAXVAL(ABS(B_face_x- &
         0.5_wp*(B_vertex(:,1:ny)+B_vertex(:,2:ny+1)))),MAXVAL(ABS(B_vertex)))
    CALL require_small('Q1 y faces',MAXVAL(ABS(B_face_y- &
         0.5_wp*(B_vertex(1:nx,:)+B_vertex(2:nx+1,:)))),MAXVAL(ABS(B_vertex)))
    DO b=1,ny
       bc=1
       IF (dimension==2) bc=MAX(3,MIN(ny-2,b))
       DO a=1,nx
          ac=MAX(3,MIN(nx-2,a))
          bx=DOT_PRODUCT(c1,centers(ac-2:ac+2,bc))/(10.0_wp*dx)
          bxx=DOT_PRODUCT(c2,centers(ac-2:ac+2,bc))/(7.0_wp*dx*dx)
          by=0.0_wp; byy=0.0_wp; bxy=0.0_wp
          IF (dimension==2) THEN
             by=DOT_PRODUCT(c1,centers(ac,bc-2:bc+2))/(10.0_wp*dy)
             byy=DOT_PRODUCT(c2,centers(ac,bc-2:bc+2))/(7.0_wp*dy*dy)
             DO t=1,5
                DO s=1,5
                   bxy=bxy+c1(s)*c1(t)*centers(ac+s-3,bc+t-3)
                END DO
             END DO
             bxy=bxy/(100.0_wp*dx*dy)
          END IF
          G=1.0_wp/(1.0_wp+bx*bx+by*by)
          CALL require_small('filtered x slope',ABS(B_prime_x_geom(a,b)-bx),bed_scale*SUM(ABS(c1))/(10.0_wp*dx))
          CALL require_small('filtered y slope',ABS(B_prime_y_geom(a,b)-by),bed_scale*SUM(ABS(c1))/(10.0_wp*dy))
          CALL require_small('filtered xx Hessian',ABS(B_second_xx_geom(a,b)-bxx),bed_scale*SUM(ABS(c2))/(7.0_wp*dx**2))
          CALL require_small('filtered xy Hessian',ABS(B_second_xy_geom(a,b)-bxy), &
               bed_scale*SUM(ABS(c1))**2/(100.0_wp*dx*dy))
          CALL require_small('filtered yy Hessian',ABS(B_second_yy_geom(a,b)-byy),bed_scale*SUM(ABS(c2))/(7.0_wp*dy**2))
          CALL require_small('center G',ABS(grav_coeff(a,b)-G),1.0_wp)
       END DO
    END DO
    CALL require_small('interior x-face G',MAXVAL(ABS(grav_coeff_stag_x(2:nx,:)- &
         0.5_wp*(grav_coeff(1:nx-1,:)+grav_coeff(2:nx,:)))),1.0_wp)
    IF (ny>1) CALL require_small('interior y-face G',MAXVAL(ABS(grav_coeff_stag_y(:,2:ny)- &
         0.5_wp*(grav_coeff(:,1:ny-1)+grav_coeff(:,2:ny)))),1.0_wp)
    CALL require_small('west G',MAXVAL(ABS(grav_coeff_stag_x(1,:)-grav_coeff(1,:))),1.0_wp)
    CALL require_small('east G',MAXVAL(ABS(grav_coeff_stag_x(nx+1,:)-grav_coeff(nx,:))),1.0_wp)
    CALL require_small('south G',MAXVAL(ABS(grav_coeff_stag_y(:,1)-grav_coeff(:,1))),1.0_wp)
    CALL require_small('north G',MAXVAL(ABS(grav_coeff_stag_y(:,ny+1)-grav_coeff(:,ny))),1.0_wp)
  END SUBROUTINE check_geometry

  !> \brief Reject nonfinite or excessive errors with the unchanged N8 roundoff threshold.
  !> \param[in] label Assertion label.
  !> \param[in] error Absolute error.
  !> \param[in] scale Reference magnitude.
  SUBROUTINE require_small(label,error,scale)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: error,scale
    IF (.NOT.ieee_is_finite(error) .OR. error>tolerance*MAX(1.0_wp,scale)) THEN
       WRITE(*,*) 'FAIL: ',label,error,tolerance*MAX(1.0_wp,scale)
       ERROR STOP 1
    END IF
  END SUBROUTINE require_small

  !> \brief Release test-grid allocations before changing resolution or dimension.
  SUBROUTINE release_grid
    CALL domain%finalize
    DEALLOCATE(q,qp,rate,nodes,weights,B_vertex,B_cent,B_face_x,B_face_y)
    DEALLOCATE(B_prime_x_geom,B_prime_y_geom,B_second_xx_geom,B_second_xy_geom,B_second_yy_geom)
    DEALLOCATE(grav_coeff,grav_coeff_stag_x,grav_coeff_stag_y,B_zone,cell_source_fractions)
    DEALLOCATE(cell_fissure_fractions,cell_arc_perim,cell_arc_n_x,cell_arc_n_y)
    DEALLOCATE(state%deposit,state%erosion,state%erodible,state%lava_deposit, &
         state%deposit_tot,state%erosion_tot)
  END SUBROUTINE release_grid
END PROGRAM test_projection_stress
