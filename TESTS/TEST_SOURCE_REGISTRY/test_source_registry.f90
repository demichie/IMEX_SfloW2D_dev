PROGRAM test_source_registry
  USE source_catalog_2d
  USE layer_contacts_2d, ONLY : layer_contacts_type, evaluate_layer_contacts
  USE source_registry_2d
  USE source_configuration_2d, ONLY : initialize_configured_sources
  USE rheology_parameters_2d, ONLY : rheology_l1,rheology_l2
  USE thermal_parameters_2d, ONLY : thermal_l1, thermal_l2,                                              &
          thermal_radiative_coefficient, thermal_convective_coefficient
  USE thermal_laws_2d, ONLY : eval_thermal_implicit_temperature
  USE source_input_2d
  USE parameters_2d
  USE mass_exchange_parameters_2d, ONLY : mass_exchange_l1,mass_exchange_l2
  USE interlayer_mass_exchange_parameters_2d, ONLY : interlayer_mass_exchange,interlayer_mass_exchange_config_type
  USE interlayer_momentum_exchange_2d, ONLY : interlayer_momentum_exchange, &
       interlayer_momentum_exchange_config_type
  USE pore_pressure_parameters_2d, ONLY : pore_pressure_l1,pore_pressure_l2
  IMPLICIT NONE
  TYPE(layer_contacts_type) :: contact_result
  TYPE(source_registry_type) :: registry
  TYPE(source_group_type), ALLOCATABLE :: groups(:)
  CHARACTER(LEN=256) :: error
  LOGICAL :: contact_valid
  REAL(wp) :: thermal_rate_l1,thermal_rate_l2
  COMPLEX(wp) :: thermal_rate_complex
  n_eqns=4; n_vars=4
  idx_totMassEqn=1;idx_uEqn=2;idx_vEqn=3;idx_engyEqn=4
  idx_pore=0;idx_poreEqn=0
  rheology_flag=.TRUE.;rheology_l1%rheology_model='voellmy'
  thermal_l1%thermal_cooling_flag=.TRUE.;pore_pressure_l1%pore_pressure_flag=.FALSE.
  curvature_term_flag=.FALSE.;bottom_radial_source_flag=.FALSE.
  bottom_fissural_source_flag=.FALSE.;radial_source_flag=.FALSE.
  lateral_source_flag=.FALSE.;mass_exchange_l1%vaporization%liquid_vaporization_flag=.FALSE.
  lava_flag=.FALSE.;gas_flag=.FALSE.;liquid_flag=.TRUE.
  mass_exchange_l1%deposition%settling_flag=.FALSE.
  mass_exchange_l1%entrainment%flag=.FALSE.
  mass_exchange_l1%erosion%coefficient=-1.0_wp
  rheology_l1%rheology_method='NEWTON';thermal_l1%thermal_method='NEWTON'
  rheology_l1%rheology_group=1;thermal_l1%thermal_group=1
  ! The layer-specific rheology stochastic switches are independent.
  rheology_l1%stochastic_rheology_flag=.FALSE.
  rheology_l2%stochastic_rheology_flag=.TRUE.
  CALL require(.NOT.rheology_l1%stochastic_rheology_flag .AND. rheology_l2%stochastic_rheology_flag, &
       'independent stochastic rheology activation per layer')
  ! Keep the automatically active liquid-loss source in explicit GROUP 0.
  ! This fixture tests only coupling/separation of the Newton source groups.
  mass_exchange_l1%method='EXPLICIT';mass_exchange_l1%group=0
  ! Layer-specific coefficients must remain independent and derived rates must
  ! match the established single-layer algebra without duplicate global caches.
  thermal_l1%emissivity=0.8_wp
  thermal_l1%atm_heat_transf_coeff=4.0_wp
  thermal_l1%exp_area_fract=0.5_wp
  thermal_l2%emissivity=0.3_wp
  thermal_l2%atm_heat_transf_coeff=2.0_wp
  thermal_l2%exp_area_fract=0.25_wp
  CALL require(thermal_radiative_coefficient(thermal_l1)/=thermal_radiative_coefficient(thermal_l2), &
       'independent thermal radiation coefficients')
  CALL require(thermal_convective_coefficient(thermal_l1)==2.0_wp, &
       'L1 convective coefficient matches prior algebra')
  CALL require(thermal_convective_coefficient(thermal_l2)==0.5_wp, &
       'L2 convective coefficient is independent')
  CALL eval_thermal_implicit_temperature(SOURCE_ID_THERMAL_COOLING,600.0_wp,thermal_l1,300.0_wp,thermal_rate_l1)
  CALL eval_thermal_implicit_temperature(SOURCE_ID_THERMAL_COOLING,600.0_wp,thermal_l2,300.0_wp,thermal_rate_l2)
  CALL require(thermal_rate_l1/=thermal_rate_l2,'thermal dispatcher uses the supplied layer config')
  CALL eval_thermal_implicit_temperature(SOURCE_ID_THERMAL_COOLING,CMPLX(600.0_wp,0.0_wp,wp), &
       thermal_l1,300.0_wp,thermal_rate_complex)
  CALL require(ABS(REAL(thermal_rate_complex,wp)-thermal_rate_l1)<= &
       64.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,ABS(thermal_rate_l1)), &
       'REAL/COMPLEX thermal dispatchers use the same coefficients')
  ! TL04-C: independently test all depth-averaged column contacts.
  CALL check_contacts([1.0_wp,2.0_wp],[.FALSE.,.TRUE.],[.TRUE.,.FALSE.],.TRUE., &
       'wet L1 and wet L2')
  CALL check_contacts([1.0_wp,0.0_wp],[.TRUE.,.FALSE.],[.TRUE.,.FALSE.],.FALSE., &
       'wet L1 and dry L2')
  CALL check_contacts([0.0_wp,2.0_wp],[.FALSE.,.TRUE.],[.FALSE.,.TRUE.],.FALSE., &
       'dry L1 and wet L2')
  CALL check_contacts([0.0_wp,0.0_wp],[.FALSE.,.FALSE.],[.FALSE.,.FALSE.],.FALSE., &
       'both layers dry')
  CALL check_contacts([1.0_wp],[.TRUE.,.FALSE.],[.TRUE.,.FALSE.],.FALSE., &
       'single layer wet')
  CALL check_contacts([0.0_wp],[.FALSE.,.FALSE.],[.FALSE.,.FALSE.],.FALSE., &
       'single layer dry')
  CALL check_contacts([dry_thickness_tolerance,2.0_wp],[.FALSE.,.TRUE.], &
       [.FALSE.,.TRUE.],.FALSE.,'canonical dry tolerance boundary')
  CALL check_contacts([2.0_wp*dry_thickness_tolerance],[.TRUE.,.FALSE.], &
       [.TRUE.,.FALSE.],.FALSE.,'canonical wet tolerance boundary')
  CALL evaluate_layer_contacts([-1.0_wp,1.0_wp],contact_result,contact_valid)
  CALL require(.NOT.contact_valid,'negative thickness rejected')
  CALL require(.NOT.ANY(contact_result%terrain_exposed) .AND. &
       .NOT.ANY(contact_result%atmosphere_exposed) .AND. &
       .NOT.contact_result%interlayer_contact,'invalid column has no contacts')
  CALL evaluate_layer_contacts([1.0_wp,2.0_wp,3.0_wp],contact_result,contact_valid)
  CALL require(.NOT.contact_valid,'unsupported layer count rejected')
  CALL require(source_catalog_names_unique(),'unique catalog names')
  CALL require(ALL(source_newton_affected_equations('rheology_voellmy')==[2,3]), &
       'rheology momentum footprint')
  CALL require(ALL(source_newton_affected_equations('thermal_cooling')==[4]), &
       'thermal footprint')
  CALL initialize_configured_sources(registry,error, &
       rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1)
  CALL require(LEN_TRIM(error)==0,'configure coupled Newton group: '//TRIM(error))
  groups=registry%groups()
  CALL require(SIZE(groups)==2,'group zero plus coupled group')
  CALL require(groups(1)%id==0 .AND. groups(2)%id==1,'Lie group order')
  CALL require(groups(2)%method==SOURCE_METHOD_NEWTON,'Newton group method')
  CALL require(SIZE(groups(2)%equation_ids)==3,'union has three unique equations')
  CALL require(ALL(groups(2)%equation_ids==[2,3,4]),'complete conservative footprint')
  CALL require(group_contains_name(groups(2),'thermal_cooling'),'cooling registered')
  CALL require(BTEST(groups(2)%source_mask,RHEOLOGY_ID_VOELLMY-1), &
       'numeric source mask contains rheology')
  CALL require(BTEST(groups(2)%source_mask,SOURCE_ID_THERMAL_COOLING-1), &
       'numeric source mask contains cooling')
  CALL require(groups(2)%source_mask==IOR(source_id_mask(RHEOLOGY_ID_VOELLMY), &
       source_id_mask(SOURCE_ID_THERMAL_COOLING)), 'Newton group source mask is exact')
  CALL require(.NOT.BTEST(groups(2)%source_mask,SOURCE_ID_GEOMETRY_CURVATURE-1), &
       'numeric source mask excludes inactive curvature')
   CALL require(groups(2)%rheology_name=='rheology_voellmy','cached rheology identity')
   CALL require(groups(2)%rheology_id==RHEOLOGY_ID_VOELLMY,'cached numeric rheology ID')
   CALL require(groups(2)%rheology_implicit,'cached implicit role')
   CALL require(groups(2)%thermal_active,'cached Newton thermal term')
  CALL registry%finalize()
  thermal_l1%thermal_group=2
  CALL initialize_configured_sources(registry,error, &
       rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1)
  CALL require(LEN_TRIM(error)==0,'separated Newton groups')
  groups=registry%groups()
  CALL require(SIZE(groups)==3,'three ordered Lie groups')
  CALL require(ALL(groups(2)%equation_ids==[2,3]),'momentum group footprint')
  CALL require(ALL(groups(3)%equation_ids==[4]),'energy group footprint')
   CALL require(groups(2)%rheology_implicit .AND. .NOT.groups(2)%thermal_active, &
        'separated Newton rheology metadata')
   CALL require(groups(3)%thermal_active .AND. .NOT.groups(3)%rheology_implicit, &
        'separated Newton cooling metadata')
  CALL registry%finalize()
  thermal_l1%thermal_method='ANALYTIC'
  CALL initialize_configured_sources(registry,error, &
       rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1)
  CALL require(LEN_TRIM(error)==0,'mixed Newton/analytic groups')
  groups=registry%groups()
  CALL require(groups(3)%method==SOURCE_METHOD_ANALYTIC,'analytic method preserved')
   CALL require(groups(3)%thermal_active,'analytic cooling metadata')
   CALL require(groups(3)%rheology_id==RHEOLOGY_ID_NONE,'analytic cooling no rheology')
  CALL registry%finalize()
  ! Historical Etna lava uses a liquid carrier with the lava rheology.
  ! Cooling and bottom fissural injection must remain valid in this phase.
  rheology_l1%rheology_model='lava'
  thermal_l1%thermal_method='NEWTON'
  thermal_l1%thermal_group=1
  bottom_fissural_source_flag=.TRUE.
  CALL initialize_configured_sources(registry,error, &
       rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1)
  CALL require(LEN_TRIM(error)==0,'Etna liquid-carrier lava configuration: '//TRIM(error))
  groups=registry%groups()
  CALL require(group_contains_name(groups(1),'bottom_fissural_injection'), &
       'Etna fissural injection in explicit group zero')
  CALL require(group_contains_name(groups(2),'rheology_lava'), &
       'Etna lava rheology in Newton group')
  CALL require(group_contains_name(groups(2),'thermal_cooling'), &
       'Etna thermal cooling in Newton group')
  CALL registry%finalize()
  ! The selected L1 configuration alone determines Darcy registration,
  ! even if L2 enables Darcy in the common structural layout.
  rheology_flag=.FALSE.
  thermal_l1%thermal_cooling_flag=.FALSE.
  bottom_fissural_source_flag=.FALSE.
  n_eqns=5; n_vars=5; idx_pore=5; idx_poreEqn=5
  pore_pressure_l1%pore_pressure_flag=.TRUE.
  pore_pressure_l2%pore_pressure_flag=.FALSE.
  pore_pressure_l1%pore_pressure_method='ANALYTIC'
  pore_pressure_l1%pore_pressure_group=3
  CALL initialize_configured_sources(registry,error, &
       rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1)
  CALL require(LEN_TRIM(error)==0,'selected L1 Darcy source: '//TRIM(error))
  groups=registry%groups()
  CALL require(SIZE(groups)==2,'Darcy source added to shared layout')
  CALL require(group_contains_name(groups(2),'pore_pressure_dissipation'), &
       'local Darcy source registered from L1 configuration')
  CALL registry%finalize()
  pore_pressure_l1%pore_pressure_flag=.FALSE.
  pore_pressure_l2%pore_pressure_flag=.TRUE.
  CALL require(.NOT.pore_pressure_l1%pore_pressure_flag .AND. pore_pressure_l2%pore_pressure_flag, &
       'independent Darcy activation per layer')
  CALL initialize_configured_sources(registry,error, &
       rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1)
  CALL require(LEN_TRIM(error)==0,'unselected Darcy in shared layout: '//TRIM(error))
  groups=registry%groups()
  CALL require(SIZE(groups)==1,'unselected Darcy has no source group')
  CALL require(.NOT.group_contains_name(groups(1),'pore_pressure_dissipation'), &
       'unselected Darcy does not register an active source')
  CALL registry%finalize()
  ! TL07-B: equal Lie group IDs can carry one source for each layer.
  DEALLOCATE(groups)
  ALLOCATE(groups(3))
  groups(1)%id=0
  groups(1)%method=SOURCE_METHOD_EXPLICIT
  groups(1)%source_names=[CHARACTER(LEN=64)::]
  groups(1)%equation_ids=[1,2,3,4,5]
  groups(2)%id=1
  groups(2)%method=SOURCE_METHOD_EXPLICIT
  groups(2)%scope=SOURCE_SCOPE_L1
  groups(2)%source_names=[CHARACTER(LEN=64)::'thermal_cooling']
  groups(2)%equation_ids=[1,2,3,4,5]
  groups(3)=groups(2)
  groups(3)%scope=SOURCE_SCOPE_L2
  CALL registry%initialize(5,0,groups,error)
  CALL require(LEN_TRIM(error)==0,'independent L1/L2 group identity')
  groups=registry%groups()
  CALL require(ALL(groups%scope==[SOURCE_SCOPE_L1,SOURCE_SCOPE_L1,SOURCE_SCOPE_L2]), &
       'ordered source scopes')
  CALL registry%finalize()
  groups(3)%scope=SOURCE_SCOPE_L1
  CALL registry%initialize(5,0,groups,error)
  CALL require(LEN_TRIM(error)>0,'duplicate group and source within L1 rejected')
  groups(3)%scope=SOURCE_SCOPE_INTERLAYER
  CALL registry%initialize(5,0,groups,error)
  CALL require(INDEX(error,'only paired mass exchange or drag')>0, &
       'layer-local source rejected in INTERLAYER scope')
  CALL require(.NOT.registry%is_initialized(), &
       'invalid interlayer group must not leave an active registry')
  ! An empty group or non-explicit method cannot bypass the paired contract.
  groups(3)%source_names=[CHARACTER(LEN=64)::]
  CALL registry%initialize(5,0,groups,error)
  CALL require(LEN_TRIM(error)>0, &
       'empty interlayer group rejected without fallback')
  CALL require(.NOT.registry%is_initialized(),'empty interlayer group leaves no registry')
  groups(3)%source_names=[CHARACTER(LEN=64)::'thermal_cooling']
  groups(3)%method=SOURCE_METHOD_NEWTON
  CALL registry%initialize(5,0,groups,error)
  CALL require(LEN_TRIM(error)>0, &
       'Newton method rejected for non-Newton interlayer group')
  CALL require(.NOT.registry%is_initialized(),'interlayer Newton leaves no registry')
  groups(3)%method=SOURCE_METHOD_EXPLICIT
  groups(3)%source_names=[CHARACTER(LEN=64)::'interlayer_mass_exchange']
  CALL registry%initialize(5,0,groups,error)
  CALL require(LEN_TRIM(error)==0,'registered paired interlayer mass exchange: '//TRIM(error))
  CALL require(registry%is_initialized(),'paired source preserved by registry')
  groups=registry%groups()
  CALL require(group_contains_name(groups(3),'interlayer_mass_exchange'), &
       'paired source identity is preserved')
  CALL registry%finalize()
  groups(3)%scope=SOURCE_SCOPE_L2
  CALL registry%initialize(5,0,groups,error)
  CALL require(INDEX(error,'requires INTERLAYER scope')>0, &
       'paired source rejected in layer-local scope')
  CALL require(.NOT.registry%is_initialized(),'wrong paired source scope leaves no registry')
  groups(3)%scope=99
  CALL registry%initialize(5,0,groups,error)
  CALL require(LEN_TRIM(error)>0,'invalid source scope rejected')
  CALL registry%finalize()
  ! Registry construction selects both layer configurations without cloning source logic.
  n_eqns=5; n_vars=5; idx_pore=5; idx_poreEqn=5
  rheology_flag=.FALSE.
  thermal_l1%thermal_cooling_flag=.FALSE.
  thermal_l2%thermal_cooling_flag=.FALSE.
  pore_pressure_l1%pore_pressure_flag=.FALSE.
  pore_pressure_l2%pore_pressure_flag=.TRUE.
  pore_pressure_l2%pore_pressure_method='ANALYTIC'
  pore_pressure_l2%pore_pressure_group=3
  CALL initialize_configured_sources(registry,error, &
       rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1, &
       rheology_l2,thermal_l2,pore_pressure_l2,mass_exchange_l2)
  CALL require(LEN_TRIM(error)==0,'two-layer source configuration: '//TRIM(error))
  groups=registry%groups()
  CALL require(COUNT(groups%id==0)==2,'one explicit group zero for each layer')
  CALL require(ANY(groups%scope==SOURCE_SCOPE_L2 .AND. groups%id==3), &
       'L2 Darcy source configured at its own group')
  CALL registry%finalize()
  interlayer_mass_exchange%elutriation_flag=.TRUE.
  interlayer_mass_exchange%method='EXPLICIT'
  interlayer_mass_exchange%group=4
  CALL initialize_configured_sources(registry,error, &
       rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1, &
       rheology_l2,thermal_l2,pore_pressure_l2,mass_exchange_l2)
  CALL require(LEN_TRIM(error)==0,'two-layer coupled group configuration: '//TRIM(error))
  groups=registry%groups()
  CALL require(COUNT(groups%scope==SOURCE_SCOPE_INTERLAYER)==1, &
       'one paired interlayer group for the two transfer laws')
  CALL require(ANY(groups%scope==SOURCE_SCOPE_INTERLAYER .AND. groups%id==4), &
       'interlayer mass exchange preserves configured Lie group')
  CALL registry%finalize()
  ! The same Lie registry supports a second paired analytic source,
  ! independently scheduled from the explicitly integrated mass exchange.
  interlayer_momentum_exchange%chezy_drag_coeff=0.03_wp
  interlayer_momentum_exchange%method='ANALYTIC'
  interlayer_momentum_exchange%group=5
  CALL initialize_configured_sources(registry,error, &
       rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1, &
       rheology_l2,thermal_l2,pore_pressure_l2,mass_exchange_l2)
  CALL require(LEN_TRIM(error)==0,'independent mass and drag groups: '//TRIM(error))
  groups=registry%groups()
  CALL require(COUNT(groups%scope==SOURCE_SCOPE_INTERLAYER)==2,'two paired Lie groups')
  CALL require(ANY(groups%scope==SOURCE_SCOPE_INTERLAYER .AND. groups%id==5 .AND. &
       groups%method==SOURCE_METHOD_ANALYTIC),'analytic drag uses its own Lie group')
  CALL require(group_contains_name(groups(SIZE(groups)),'interlayer_drag'), &
       'drag is the last configured source')
  CALL registry%finalize()
  interlayer_momentum_exchange%group=4
  CALL initialize_configured_sources(registry,error, &
       rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1, &
       rheology_l2,thermal_l2,pore_pressure_l2,mass_exchange_l2)
  CALL require(INDEX(error,'distinct Lie groups')>0,'same-group mass and drag are rejected')
  CALL require(.NOT.registry%is_initialized(),'conflicting groups leave no registry')
  interlayer_momentum_exchange%group=5
  interlayer_momentum_exchange%chezy_drag_coeff=0.0_wp
  interlayer_mass_exchange%elutriation_flag=.FALSE.
  CALL initialize_configured_sources(registry,error, &
       rheology_l1,thermal_l1,pore_pressure_l1,mass_exchange_l1, &
       rheology_l2=rheology_l2)
  CALL require(LEN_TRIM(error)>0,'incomplete L2 configurations rejected')
  CALL require(.NOT.registry%is_initialized(),'failed L2 setup leaves registry uninitialized')
  CALL check_interlayer_input
  PRINT *, 'PASS Lie registry complete Newton footprints, ordering and methods'
CONTAINS
  !> @brief Verify one canonical pair-wide NAMELIST with independent mass and drag settings.
  SUBROUTINE check_interlayer_input
    TYPE(interlayer_mass_exchange_config_type) :: saved_mass
    TYPE(interlayer_momentum_exchange_config_type) :: saved_drag
    INTEGER :: unit,saved_layers
    CHARACTER(LEN=256) :: input_error
    CHARACTER(LEN=256) :: line

    saved_layers=n_layers
    saved_mass=interlayer_mass_exchange
    saved_drag=interlayer_momentum_exchange
    n_layers=2
    interlayer_mass_exchange=interlayer_mass_exchange_config_type()
    interlayer_momentum_exchange=interlayer_momentum_exchange_config_type()

    OPEN(NEWUNIT=unit,STATUS='SCRATCH',FORM='FORMATTED')
    WRITE(unit,'(A)') '&INTERLAYER_PARAMETERS'
    WRITE(unit,'(A)') ' interlayer_mass_exchange%group=3,'
    WRITE(unit,'(A)') ' interlayer_momentum_exchange%group=7,'
    WRITE(unit,'(A)') ' interlayer_momentum_exchange%chezy_drag_coeff=0.125,'
    WRITE(unit,'(A)') ' /'
    REWIND(unit)
    CALL read_interlayer_parameters(unit,input_error)
    CALL require(LEN_TRIM(input_error)==0,'canonical interlayer input: '//TRIM(input_error))
    CALL require(interlayer_mass_exchange%group==3,'interlayer mass group read from canonical pair')
    CALL require(interlayer_momentum_exchange%group==7,'interlayer drag group read independently')
    CALL require(interlayer_momentum_exchange%chezy_drag_coeff==0.125_wp, &
         'interlayer Chezy coefficient read from canonical pair')
    CLOSE(unit)

    OPEN(NEWUNIT=unit,STATUS='SCRATCH',FORM='FORMATTED')
    CALL write_interlayer_parameters(unit)
    REWIND(unit)
    READ(unit,'(A)') line
    CALL require(INDEX(uppercase_source_text(ADJUSTL(line)),'&INTERLAYER_PARAMETERS')==1, &
         'interlayer backup uses canonical NAMELIST')
    CLOSE(unit)

    OPEN(NEWUNIT=unit,STATUS='SCRATCH',FORM='FORMATTED')
    WRITE(unit,'(A)') '&INTERLAYER_PARAMETERS /'
    WRITE(unit,'(A)') '&INTERLAYER_PARAMETERS /'
    REWIND(unit)
    CALL read_interlayer_parameters(unit,input_error)
    CALL require(INDEX(input_error,'Duplicate INTERLAYER_PARAMETERS')>0, &
         'duplicate pair-wide NAMELIST rejected')
    CLOSE(unit)

    n_layers=1
    OPEN(NEWUNIT=unit,STATUS='SCRATCH',FORM='FORMATTED')
    WRITE(unit,'(A)') '&INTERLAYER_PARAMETERS /'
    REWIND(unit)
    CALL read_interlayer_parameters(unit,input_error)
    CALL require(INDEX(input_error,'requires N_LAYERS=2')>0, &
         'interlayer NAMELIST rejected for single-layer input')
    CLOSE(unit)

    n_layers=saved_layers
    interlayer_mass_exchange=saved_mass
    interlayer_momentum_exchange=saved_drag
  END SUBROUTINE check_interlayer_input

  !> @brief Check the complete physical column-contact contract.
  !> @param[in] thickness Candidate layer thicknesses in metres.
  !> @param[in] atmosphere Expected atmosphere exposure for L1/L2.
  !> @param[in] terrain Expected ground exposure for L1/L2.
  !> @param[in] interface_contact Expected wet L1/L2 interface contact.
  !> @param[in] description Regression case name.
  SUBROUTINE check_contacts(thickness,atmosphere,terrain,interface_contact,description)
    REAL(wp), INTENT(IN) :: thickness(:)
    LOGICAL, INTENT(IN) :: atmosphere(2),terrain(2),interface_contact
    CHARACTER(LEN=*), INTENT(IN) :: description
    TYPE(layer_contacts_type) :: contacts
    LOGICAL :: valid
    CALL evaluate_layer_contacts(thickness,contacts,valid)
    CALL require(valid,'valid contact column: '//description)
    CALL require(ALL(contacts%atmosphere_exposed .EQV. atmosphere), &
         'atmosphere contact: '//description)
    CALL require(ALL(contacts%terrain_exposed .EQV. terrain), &
         'terrain contact: '//description)
    CALL require(contacts%interlayer_contact .EQV. interface_contact, &
         'interlayer contact: '//description)
  END SUBROUTINE check_contacts

  SUBROUTINE require(ok,message)
    LOGICAL,INTENT(IN) :: ok
    CHARACTER(*),INTENT(IN) :: message
    IF (.NOT.ok) THEN
       PRINT *, 'FAIL: ',TRIM(message)
       ERROR STOP 1
    END IF
  END SUBROUTINE require
END PROGRAM test_source_registry
