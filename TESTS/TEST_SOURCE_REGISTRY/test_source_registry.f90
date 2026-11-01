PROGRAM test_source_registry
  USE source_catalog_2d
  USE source_registry_2d
  USE source_configuration_2d, ONLY : initialize_configured_sources
  USE source_input_2d
  USE parameters_2d
  USE constitutive_parameters_2d, ONLY : erosion_coeff,settling_flag,entrainment_flag
  IMPLICIT NONE
  TYPE(source_registry_type) :: registry
  TYPE(source_group_type), ALLOCATABLE :: groups(:)
  CHARACTER(LEN=256) :: error
  INTEGER :: i,idx
  n_eqns=4; n_vars=4
  idx_totMassEqn=1;idx_uEqn=2;idx_vEqn=3;idx_engyEqn=4
  idx_pore=0;idx_poreEqn=0
  rheology_flag=.TRUE.;rheology_model='voellmy'
  thermal_cooling_flag=.TRUE.;pore_pressure_flag=.FALSE.
  curvature_term_flag=.FALSE.;bottom_radial_source_flag=.FALSE.
  bottom_fissural_source_flag=.FALSE.;radial_source_flag=.FALSE.
  lateral_source_flag=.FALSE.;liquid_vaporization_flag=.FALSE.
  lava_flag=.FALSE.;gas_flag=.FALSE.;liquid_flag=.TRUE.
  settling_flag=.FALSE.;entrainment_flag=.FALSE.
  IF (ALLOCATED(erosion_coeff)) DEALLOCATE(erosion_coeff)
  rheology_method='NEWTON';thermal_method='NEWTON'
  rheology_group=1;thermal_group=1
  ! Keep the automatically active liquid-loss source in explicit GROUP 0.
  ! This fixture tests only coupling/separation of the Newton source groups.
  mass_exchange_method='EXPLICIT';mass_exchange_group=0
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
  CALL registry%finalize()
  thermal_group=2
  CALL initialize_configured_sources(registry,error)
  CALL require(LEN_TRIM(error)==0,'separated Newton groups')
  groups=registry%groups()
  CALL require(SIZE(groups)==3,'three ordered Lie groups')
  CALL require(ALL(groups(2)%equation_ids==[2,3]),'momentum group footprint')
  CALL require(ALL(groups(3)%equation_ids==[4]),'energy group footprint')
  CALL registry%finalize()
  thermal_method='ANALYTIC'
  CALL initialize_configured_sources(registry,error)
  CALL require(LEN_TRIM(error)==0,'mixed Newton/analytic groups')
  groups=registry%groups()
  CALL require(groups(3)%method==SOURCE_METHOD_ANALYTIC,'analytic method preserved')
  CALL registry%finalize()
  PRINT *, 'PASS Lie registry complete Newton footprints, ordering and methods'
CONTAINS
  SUBROUTINE require(ok,message)
    LOGICAL,INTENT(IN) :: ok
    CHARACTER(*),INTENT(IN) :: message
    IF (.NOT.ok) THEN
       PRINT *, 'FAIL: ',TRIM(message)
       ERROR STOP 1
    END IF
  END SUBROUTINE require
END PROGRAM test_source_registry
