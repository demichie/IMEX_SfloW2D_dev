!> \brief D17 metadata regressions: legacy schedules, coupling and invalid native groups.
PROGRAM test_source_registry
  USE source_catalog_2d
  USE source_registry_2d
  USE source_configuration_2d, ONLY : initialize_configured_sources
  USE equation_metadata_2d, ONLY : equation_partition_type, &
       IM_GROUP_METHOD_ANALYTIC, IM_GROUP_METHOD_NEWTON
  USE parameters_2d, ONLY : rheology_id_lava
  IMPLICIT NONE
  TYPE(source_registry_type) :: registry,trial
  TYPE(equation_partition_type) :: partition
  TYPE(source_group_type), ALLOCATABLE :: groups(:),steps(:),badgroups(:),badsteps(:),snapshot(:)
  TYPE(source_solver_type), ALLOCATABLE :: solves(:),badsolves(:)
  TYPE(source_term_type) :: term
  INTEGER, ALLOCATABLE :: terms(:),family(:),bundle(:),expected(:),expected_targets(:)
  INTEGER :: a,b,c,bits,mode,g,t,k,j,position,checks,methods(3),original(3),pore,nphysical,negative,model_id
  LOGICAL :: enabled(3),mask(7)
  CHARACTER(LEN=256) :: error
  checks=0; negative=0
  nphysical=0
  DO k=1,SOURCE_TERM_COUNT
    term=source_term_descriptor(k)
    CALL require(term%id==k .AND. LEN_TRIM(term%name)>0,'catalog ID/name')
    IF (term%schedulable) nphysical=nphysical+1
  END DO
  CALL require(nphysical==SOURCE_PHYSICAL_TERM_COUNT,'physical catalog count')
  term=source_term_descriptor(0); CALL require(term%id==0,'invalid catalog lookup')
  ! All 864 legacy assignments/method combinations, with optional thermal/pore targets.
  DO mode=0,3
    enabled=[.TRUE.,BTEST(mode,0),BTEST(mode,1)]
    pore=0
    IF (enabled(3)) pore=7
    family=rheology_source_terms(MERGE(3,1,enabled(2)))
    terms=[family,SOURCE_TERM_GEOMETRY_CURVATURE,SOURCE_TERM_BOTTOM_CIRCULAR_INJECTION, &
         SOURCE_TERM_BOTTOM_FISSURAL_INJECTION,SOURCE_TERM_RADIAL_INJECTION, &
         SOURCE_TERM_BOUNDARY_INJECTION,SOURCE_TERM_MASS_EXCHANGE_EROSION, &
         SOURCE_TERM_MASS_EXCHANGE_DEPOSITION,SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS, &
         SOURCE_TERM_MASS_EXCHANGE_CARRIER_DEPOSITION_LOSS,SOURCE_TERM_MASS_EXCHANGE_CARRIER_DARCY_LOSS, &
         SOURCE_TERM_MASS_EXCHANGE_ENTRAINMENT,SOURCE_TERM_MASS_EXCHANGE_VAPORIZATION]
    IF (enabled(2)) terms=[terms,SOURCE_TERM_THERMAL_RADIATIVE_COOLING,SOURCE_TERM_THERMAL_CONVECTIVE_COOLING]
    IF (enabled(3)) terms=[terms,SOURCE_TERM_PORE_PRESSURE_DISSIPATION]
    DO a=1,3
    DO b=1,3
    DO c=1,3
    DO bits=0,7
      original=[a,b,c]
      DO g=1,3
        methods(g)=MERGE(IM_GROUP_METHOD_ANALYTIC,IM_GROUP_METHOD_NEWTON,BTEST(bits,g-1))
      END DO
      CALL make_partition()
      CALL registry%initialize_legacy(partition,original,methods,pore,terms,error)
      CALL require(error=='' .AND. registry%is_initialized(),'legacy normalization: '//TRIM(error))
      groups=registry%stage_groups(); steps=registry%step_groups(); solves=registry%solvers()
      position=0
      DO g=1,3
        expected_targets=PACK([1,2,3],enabled .AND. original==g)
        IF (SIZE(expected_targets)==0) CYCLE
        IF (methods(g)==IM_GROUP_METHOD_NEWTON) THEN
          position=position+1
          CALL require(groups(position)%method==SOURCE_METHOD_NEWTON,'legacy Newton method')
          CALL require(groups(position)%id==g .AND. groups(position)%legacy_id==g,'retained Newton ID')
          CALL require(SIZE(groups(position)%implicit_solver_ids)==1,'joint Newton solve retained')
          k=groups(position)%implicit_solver_ids(1)
          expected=[INTEGER::]
          DO j=1,SIZE(expected_targets)
            expected=[expected,target_eqs(expected_targets(j))]
          END DO
          CALL require(same_members(solves(k)%equation_ids,expected),'joint equation/derivative block')
          CALL require(same_members(solves(k)%target_ids,expected_targets),'joint target coverage')
        ELSE
          DO t=1,3
            IF (.NOT.enabled(t) .OR. original(t)/=g) CYCLE
            position=position+1
            CALL require(groups(position)%legacy_id==g,'analytic original ID')
            CALL require(groups(position)%method==SOURCE_METHOD_ANALYTIC,'analytic method')
            CALL require(SIZE(groups(position)%analytic_solver_ids)==1,'one analytic solve per group')
            k=groups(position)%analytic_solver_ids(1)
            CALL require(ALL(solves(k)%target_ids==t),'legacy velocity-temperature-pressure order')
            IF (SIZE(expected_targets)>1) THEN
              CALL require(groups(position)%id>4,'split identity avoids EX/legacy collisions')
            ELSE
              CALL require(groups(position)%id==g,'single analytic identity retained')
            END IF
          END DO
        END IF
      END DO
      CALL require(SIZE(groups)==position+1,'only one additional EX group')
      CALL require(groups(position+1)%id==4,'reserved EX ID')
      CALL require(groups(position+1)%method==SOURCE_METHOD_EXPLICIT,'EX homogeneous method')
      CALL require(SIZE(steps)==1 .AND. steps(1)%id==1,'independent step namespace')
      CALL require(group_contains(steps(1),SOURCE_TERM_MASS_EXCHANGE_DEPOSITION),'mass membership')
      CALL registry%validate(error); CALL require(error=='','repeat validation')
      snapshot=registry%stage_groups(); snapshot(1)%id=999
      snapshot=registry%stage_groups(); CALL require(snapshot(1)%id/=999,'detached immutable storage')
    END DO
    END DO
    END DO
    END DO
  END DO
  ! One configured bridge smoke test covers all injection and mass-exchange gates.
  CALL configure_bridge()
  enabled=.TRUE.; original=[3,1,2]; methods=[IM_GROUP_METHOD_ANALYTIC, &
       IM_GROUP_METHOD_ANALYTIC,IM_GROUP_METHOD_NEWTON]; pore=7
  CALL make_partition()
  CALL initialize_configured_sources(registry,partition,error)
  CALL require(error=='','configured source bridge: '//TRIM(error))
  groups=registry%stage_groups();steps=registry%step_groups();solves=registry%solvers()
  CALL require(SIZE(groups)==4 .AND. SIZE(solves)==3,'configured default groups')
  CALL require(SIZE(groups(4)%term_ids)==5,'curvature/circular/fissural/radial/boundary inventory')
  CALL require(SIZE(steps(1)%term_ids)==7,'all bounded mass proposals registered')
  CALL require(group_contains(groups(4),SOURCE_TERM_BOUNDARY_INJECTION),'boundary hook retained')
  DO mode=1,21
    badgroups=groups;badsteps=steps;badsolves=solves
    SELECT CASE(mode)
    CASE(1);badgroups(1)%id=4 ! Collision with EX group.
    CASE(2);badgroups(4)%term_ids=[badgroups(4)%term_ids,SOURCE_TERM_GEOMETRY_CURVATURE]
    CASE(3);badsteps(1)%term_ids=[badsteps(1)%term_ids,SOURCE_TERM_GEOMETRY_CURVATURE]
    CASE(4);badgroups(1)%analytic_solver_ids=[1,2]
    CASE(5);badgroups(3)%implicit_solver_ids=[INTEGER::]
    CASE(6);badgroups(4)%term_ids=[SOURCE_TERM_INLET_THERMAL_COUPLING]
    CASE(7);badgroups(4)%term_ids=[999]
    CASE(8);badsolves(1)%target_ids=[4]
    CASE(9);badsolves(1)%equation_ids=[99]
    CASE(10);DEALLOCATE(badgroups(1)%equation_ids)
    CASE(11);badgroups(1)%id=0
    CASE(12);badgroups(4)%cadence=SOURCE_CADENCE_STEP
    CASE(13)
      badgroups(1)%term_ids=[SOURCE_TERM_THERMAL_RADIATIVE_COOLING]
      badsolves(1)%term_ids=badgroups(1)%term_ids
    CASE(14)
      badsolves=[badsolves,badsolves(3)];badsolves(4)%id=4
    CASE(15)
      badsolves=[badsolves,badsolves(3)];badsolves(4)%id=4;badsolves(4)%execution_order=2
      badgroups(3)%implicit_solver_ids=[3,4]
    CASE(16)
      badgroups(3)%term_ids=[badgroups(3)%term_ids,SOURCE_TERM_RHEOLOGY_PLASTIC_YIELD]
      badsolves(3)%term_ids=badgroups(3)%term_ids
    CASE(17)
      badsteps=[badsteps,badsteps];badsteps(2)%id=2;badsteps(2)%execution_order=2
      badsteps(1)%term_ids=steps(1)%term_ids(1:2);badsteps(2)%term_ids=steps(1)%term_ids(3:)
    CASE(18);badgroups(2)%execution_order=1
    CASE(19);badgroups(1)%method=99
    CASE(20);badgroups(1)%analytic_solver_ids=[99]
    CASE(21);badsolves(3)%term_ids=[INTEGER::]
    END SELECT
    CALL trial%initialize(7,7,badgroups,badsteps,badsolves,error)
    CALL require(error/='' .AND. .NOT.trial%is_initialized(),'invalid native registry rejected')
    snapshot=trial%stage_groups();CALL require(SIZE(snapshot)==0,'failed initialization releases storage')
    negative=negative+1
  END DO
  ! Native metadata can represent independent Newton blocks without merging equations.
  badgroups=[groups(1),groups(3),groups(4)];badsolves=solves
  badgroups(1)%method=SOURCE_METHOD_NEWTON
  badgroups(1)%equation_ids=[4,7]
  badgroups(1)%term_ids=[groups(1)%term_ids,groups(2)%term_ids]
  badgroups(1)%analytic_solver_ids=[INTEGER::];badgroups(1)%implicit_solver_ids=[1,2]
  badgroups(2)%execution_order=2
  badsolves(1)%kind=SOURCE_SOLVE_NEWTON_BLOCK;badsolves(2)%kind=SOURCE_SOLVE_NEWTON_BLOCK
  badsolves(2)%execution_order=2
  CALL trial%initialize(7,7,badgroups,steps,badsolves,error)
  CALL require(error=='','multiple independent Newton blocks: '//TRIM(error))
  ! A partial analytic velocity bundle must never silently evaluate the full law.
  badgroups=groups;badsolves=solves
  badgroups(3)%method=SOURCE_METHOD_ANALYTIC
  badgroups(3)%implicit_solver_ids=[INTEGER::];badgroups(3)%analytic_solver_ids=[3]
  badgroups(3)%term_ids=[SOURCE_TERM_RHEOLOGY_LAVA_VISCOUS]
  badsolves(3)%kind=SOURCE_SOLVE_ANALYTIC_VELOCITY;badsolves(3)%term_ids=badgroups(3)%term_ids
  CALL trial%initialize(7,7,badgroups,steps,badsolves,error)
  CALL require(INDEX(error,'complete law bundle')>0,'partial analytic velocity bundle rejected')
  negative=negative+1
  ! Model membership and atomic batches are derived from the catalog.
  DO k=1,SOURCE_TERM_COUNT
    term=source_term_descriptor(k)
    IF (term%model_id==0) CYCLE
    model_id=term%model_id
    family=rheology_source_terms(model_id)
    CALL require(SIZE(family)>0,'catalogued rheology family present')
    DO j=1,SIZE(family)
      term=source_term_descriptor(family(j))
      CALL require(term%model_id==model_id, &
           'rheology lookup returns only its registered family')
    END DO
  END DO
  DO k=1,SOURCE_TERM_COUNT
    term=source_term_descriptor(k)
    IF (term%bundle_id==0) CYCLE
    bundle=source_bundle_terms(term%bundle_id)
    CALL require(SIZE(bundle)>=2,'catalogued atomic batch has multiple terms')
    IF (k/=bundle(1)) CYCLE
    enabled=[.TRUE.,.FALSE.,.FALSE.]; original=[1,1,1]; methods=IM_GROUP_METHOD_NEWTON; pore=0
    DO j=1,SIZE(bundle)
      term=source_term_descriptor(bundle(j))
      IF (term%target>0) enabled(term%target)=.TRUE.
    END DO
    IF (enabled(3)) pore=7
    CALL make_partition()
    terms=bundle(:SIZE(bundle)-1)
    CALL trial%initialize_legacy(partition,original,methods,pore,terms,error)
    CALL require(INDEX(error,'batch')>0,'incomplete catalogued batch rejected')
    negative=negative+1
  END DO
  enabled=[.TRUE.,.FALSE.,.FALSE.]; pore=0
  original=[3,1,2]; methods=[IM_GROUP_METHOD_ANALYTIC,IM_GROUP_METHOD_ANALYTIC,IM_GROUP_METHOD_NEWTON]
  CALL make_partition()
  terms=[SOURCE_TERM_THERMAL_RADIATIVE_COOLING,SOURCE_TERM_THERMAL_CONVECTIVE_COOLING]
  CALL trial%initialize_legacy(partition,original,methods,pore,terms,error)
  CALL require(INDEX(error,'enabled legacy target')>0,'active source cannot disappear on an inactive target')
  negative=negative+1
  ! Input/source activation may be empty while momentum remains an implicit block.
  CALL disable_sources()
  enabled=[.TRUE.,.FALSE.,.FALSE.];pore=0
  CALL make_partition()
  CALL initialize_configured_sources(registry,partition,error)
  CALL require(error=='','disabled physics still has valid metadata')
  groups=registry%stage_groups();solves=registry%solvers()
  CALL require(SIZE(groups)==2 .AND. SIZE(solves)==1,'zero-rate momentum block retained')
  CALL require(SIZE(solves(1)%term_ids)==0,'disabled rheology has no physical rates')
  CALL require(ALL(solves(1)%equation_ids==[2,3]),'zero-rate momentum target preserved')
  mask=.FALSE.
  CALL partition%initialize(mask)
  CALL initialize_configured_sources(registry,partition,error)
  CALL require(error=='','empty implicit partition accepted')
  solves=registry%solvers();CALL require(SIZE(solves)==0,'no orphan solve in explicit-only partition')
  mask=.FALSE.;mask(2)=.TRUE.
  CALL partition%initialize(mask)
  CALL trial%initialize_legacy(partition,original,methods,0,[INTEGER::],error)
  CALL require(INDEX(error,'two-equation block')>0,'split legacy momentum rejected');negative=negative+1
  enabled=[.TRUE.,.FALSE.,.FALSE.];CALL make_partition()
  CALL trial%initialize_legacy(partition,[0,1,2],methods,0,[INTEGER::],error)
  CALL require(INDEX(error,'group IDs')>0,'invalid legacy ID rejected');negative=negative+1
  CALL trial%initialize_legacy(partition,original,[99,99,99],0,[INTEGER::],error)
  CALL require(INDEX(error,'legacy method')>0,'invalid legacy method rejected');negative=negative+1
  CALL registry%finalize();CALL registry%finalize()
  CALL require(.NOT.registry%is_initialized(),'idempotent finalization')
  snapshot=registry%stage_groups();CALL require(SIZE(snapshot)==0,'released stage storage')
  WRITE(*,'(A,I0,A,I0,A)') 'PASS: D17 registry, 864 legacy configurations; ',checks,' assertions; ',negative, &
       ' rejected invalid configurations'
CONTAINS
  SUBROUTINE require(condition,label)
    LOGICAL, INTENT(IN) :: condition
    CHARACTER(*), INTENT(IN) :: label
    IF (.NOT.condition) THEN
      WRITE(*,*) 'FAIL: ',label
      ERROR STOP 1
    END IF
    checks=checks+1
  END SUBROUTINE require
  FUNCTION target_eqs(target) RESULT(ids)
    INTEGER, INTENT(IN) :: target
    INTEGER, ALLOCATABLE :: ids(:)
    SELECT CASE(target)
    CASE(1);ids=[2,3]
    CASE(2);ids=[4]
    CASE(3);ids=[7]
    END SELECT
  END FUNCTION target_eqs
  LOGICAL FUNCTION same_members(left,right) RESULT(equal)
    INTEGER, INTENT(IN) :: left(:),right(:)
    INTEGER :: index
    equal=SIZE(left)==SIZE(right)
    DO index=1,SIZE(left)
      equal=equal .AND. COUNT(right==left(index))==1
    END DO
  END FUNCTION same_members
  SUBROUTINE make_partition()
    INTEGER :: membership(7),compact_methods(3),count_groups,group,target
    mask=.FALSE.;mask(2:3)=.TRUE.;mask(4)=enabled(2);mask(7)=enabled(3)
    membership=0;count_groups=0
    DO group=1,3
      IF (.NOT.ANY(enabled .AND. original==group)) CYCLE
      count_groups=count_groups+1;compact_methods(count_groups)=methods(group)
      DO target=1,3
        IF (enabled(target) .AND. original(target)==group) membership(target_eqs(target))=count_groups
      END DO
    END DO
    CALL partition%initialize(mask,membership,compact_methods(1:count_groups))
  END SUBROUTINE make_partition
  SUBROUTINE configure_bridge()
    USE parameters_2d, ONLY : rheology_flag,rheology_model_id,pore_pressure_flag,idx_pore, &
         curvature_term_flag,bottom_radial_source_flag,bottom_fissural_source_flag, &
         radial_source_flag,lateral_source_flag,liquid_flag,gas_flag,gas_loss_flag, &
         liquid_vaporization_flag,imex_velocity_group,imex_energy_group,imex_pore_pressure_group,imex_group_method
    USE constitutive_parameters_2d, ONLY : erosion_coeff,settling_flag,coeff_porosity,entrainment_flag
    rheology_flag=.TRUE.;rheology_model_id=rheology_id_lava;pore_pressure_flag=.TRUE.;idx_pore=7
    curvature_term_flag=.TRUE.;bottom_radial_source_flag=.TRUE.;bottom_fissural_source_flag=.TRUE.
    radial_source_flag=.TRUE.;lateral_source_flag=.TRUE.;liquid_flag=.TRUE.;gas_flag=.TRUE.
    gas_loss_flag=.TRUE.;liquid_vaporization_flag=.TRUE.;erosion_coeff=0.1;settling_flag=.TRUE.
    coeff_porosity=0.3;entrainment_flag=.TRUE.
    imex_velocity_group=3;imex_energy_group=1;imex_pore_pressure_group=2
    imex_group_method=[CHARACTER(LEN=16)::'ANALYTIC','ANALYTIC','NEWTON']
  END SUBROUTINE configure_bridge
  SUBROUTINE disable_sources()
    USE parameters_2d, ONLY : rheology_flag,rheology_model_id,rheology_id_none,pore_pressure_flag, &
         curvature_term_flag,bottom_radial_source_flag,bottom_fissural_source_flag, &
         radial_source_flag,lateral_source_flag,liquid_flag,gas_flag,gas_loss_flag,liquid_vaporization_flag
    USE constitutive_parameters_2d, ONLY : erosion_coeff,settling_flag,coeff_porosity,entrainment_flag
    rheology_flag=.FALSE.;rheology_model_id=rheology_id_none;pore_pressure_flag=.FALSE.
    curvature_term_flag=.FALSE.;bottom_radial_source_flag=.FALSE.;bottom_fissural_source_flag=.FALSE.
    radial_source_flag=.FALSE.;lateral_source_flag=.FALSE.;liquid_flag=.FALSE.;gas_flag=.FALSE.
    gas_loss_flag=.FALSE.;liquid_vaporization_flag=.FALSE.;erosion_coeff=0.0;settling_flag=.FALSE.
    coeff_porosity=0.0;entrainment_flag=.FALSE.
  END SUBROUTINE disable_sources
END PROGRAM test_source_registry
