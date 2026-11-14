PROGRAM test_source_registry
  USE source_catalog_2d
  USE layer_contacts_2d, ONLY : layer_contacts_type, evaluate_layer_contacts
  USE source_registry_2d
  USE source_configuration_2d, ONLY : initialize_configured_sources
  USE rheology_parameters_2d, ONLY : rheology_l1
  USE thermal_parameters_2d, ONLY : thermal_l1, thermal_l2,                                              &
          thermal_radiative_coefficient, thermal_convective_coefficient
  USE thermal_laws_2d, ONLY : eval_thermal_implicit_temperature
  USE source_input_2d
  USE parameters_2d
  USE mass_exchange_parameters_2d, ONLY : mass_exchange_l1
  IMPLICIT NONE
  TYPE(layer_contacts_type) :: contact_result
  TYPE(source_registry_type) :: registry
  TYPE(source_group_type), ALLOCATABLE :: groups(:)
  CHARACTER(LEN=256) :: error
  INTEGER :: i,idx
  LOGICAL :: contact_valid
  REAL(wp) :: thermal_rate_l1,thermal_rate_l2
  COMPLEX(wp) :: thermal_rate_complex
  n_eqns=4; n_vars=4
  idx_totMassEqn=1;idx_uEqn=2;idx_vEqn=3;idx_engyEqn=4
  idx_pore=0;idx_poreEqn=0
  rheology_flag=.TRUE.;rheology_l1%rheology_model='voellmy'
  thermal_l1%thermal_cooling_flag=.TRUE.;pore_pressure_flag=.FALSE.
  curvature_term_flag=.FALSE.;bottom_radial_source_flag=.FALSE.
  bottom_fissural_source_flag=.FALSE.;radial_source_flag=.FALSE.
  lateral_source_flag=.FALSE.;mass_exchange_l1%vaporization%liquid_vaporization_flag=.FALSE.
  lava_flag=.FALSE.;gas_flag=.FALSE.;liquid_flag=.TRUE.
  mass_exchange_l1%deposition%settling_flag=.FALSE.
  mass_exchange_l1%entrainment%flag=.FALSE.
  mass_exchange_l1%erosion%coefficient=-1.0_wp
  rheology_l1%rheology_method='NEWTON';thermal_l1%thermal_method='NEWTON'
  rheology_l1%rheology_group=1;thermal_l1%thermal_group=1
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
  CALL require(canonical_rheology_catalog_unique(),'unique rheology names')
  CALL require(ALL(source_newton_affected_equations('rheology_voellmy')==[2,3]), &
       'rheology momentum footprint')
  CALL require(ALL(source_newton_affected_equations('thermal_cooling')==[4]), &
       'thermal footprint')
  CALL initialize_configured_sources(registry,error)
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
  CALL require(.NOT.BTEST(groups(2)%source_mask,SOURCE_ID_GEOMETRY_CURVATURE-1), &
       'numeric source mask excludes inactive curvature')
   CALL require(groups(2)%rheology_name=='rheology_voellmy','cached rheology identity')
   CALL require(groups(2)%rheology_id==RHEOLOGY_ID_VOELLMY,'cached numeric rheology ID')
   CALL require(groups(2)%rheology_implicit,'cached implicit role')
   CALL require(groups(2)%thermal_active,'cached Newton thermal term')
  CALL registry%finalize()
  thermal_l1%thermal_group=2
  CALL initialize_configured_sources(registry,error)
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
  CALL initialize_configured_sources(registry,error)
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
  CALL initialize_configured_sources(registry,error)
  CALL require(LEN_TRIM(error)==0,'Etna liquid-carrier lava configuration: '//TRIM(error))
  groups=registry%groups()
  CALL require(group_contains_name(groups(1),'bottom_fissural_injection'), &
       'Etna fissural injection in explicit group zero')
  CALL require(group_contains_name(groups(2),'rheology_lava'), &
       'Etna lava rheology in Newton group')
  CALL require(group_contains_name(groups(2),'thermal_cooling'), &
       'Etna thermal cooling in Newton group')
  CALL registry%finalize()
  PRINT *, 'PASS Lie registry complete Newton footprints, ordering and methods'
CONTAINS
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
