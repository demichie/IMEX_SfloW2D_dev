!> \brief D17 metadata regressions: legacy schedules, coupling and invalid native groups.
PROGRAM test_source_registry
  USE source_catalog_2d
  USE source_registry_2d
  USE source_configuration_2d, ONLY : initialize_configured_sources
  USE source_configuration_2d, ONLY : initialize_native_source_configuration
  USE source_input_2d
  USE equation_metadata_2d, ONLY : equation_partition_type, &
       IM_GROUP_METHOD_ANALYTIC, IM_GROUP_METHOD_NEWTON
  USE parameters_2d, ONLY : rheology_id_lava,n_eqns,idx_pore
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
    ! Every catalog entry must resolve to the same descriptor by canonical name.
    ! This also covers non-schedulable dependent contributions.
    BLOCK
      TYPE(source_term_type) :: named_term
      named_term=source_term_descriptor_by_name(term%name)
      CALL require(named_term%id==k .AND. named_term%target==term%target, &
           'catalog name round trip')
      CALL require(named_term%has_explicit .EQV. term%has_explicit, &
           'catalog explicit capability round trip')
      CALL require(named_term%has_newton .EQV. term%has_newton, &
           'catalog Newton capability round trip')
      CALL require(named_term%has_analytic .EQV. term%has_analytic, &
           'catalog analytic capability round trip')
    END BLOCK
    IF (term%schedulable) nphysical=nphysical+1
  END DO
  CALL require(nphysical==SOURCE_PHYSICAL_TERM_COUNT,'physical catalog count')
  ! Transitional regression: preserve the legacy method matrix while the
  ! numerical capabilities are moved from roles into the catalog.
  DO k=1,SOURCE_TERM_COUNT
    term=source_term_descriptor(k)
    SELECT CASE(term%role)
    CASE(SOURCE_ROLE_IMPLICIT,SOURCE_ROLE_THERMAL,SOURCE_ROLE_PORE)
      CALL require(.NOT.term%has_explicit,'legacy implicit source is not explicit')
      CALL require(term%has_newton .EQV. term%complex_rate, &
           'legacy complex Newton capability')
      CALL require(term%has_analytic,'legacy analytic eligibility')
    CASE(SOURCE_ROLE_SEMI_IMPLICIT)
      CALL require(.NOT.term%has_explicit .AND. term%has_newton &
           .AND. term%has_analytic,'semi-implicit contribution in Newton contract')
    CASE(SOURCE_ROLE_CELL,SOURCE_ROLE_BOUNDARY,SOURCE_ROLE_MASS)
      CALL require(term%has_explicit .AND. .NOT.term%has_newton &
           .AND. .NOT.term%has_analytic,'explicit source capabilities')
    CASE(SOURCE_ROLE_DERIVED)
      CALL require(.NOT.term%has_explicit .AND. .NOT.term%has_newton &
           .AND. .NOT.term%has_analytic,'dependent contributions not schedulable')
    END SELECT
  END DO
  CALL require(source_catalog_names_unique(),'catalog canonical names are unique')
  term=source_term_descriptor_by_name('SOURCE_TERM_MASS_EXCHANGE_DEPOSITION')
  CALL require(term%id==SOURCE_TERM_MASS_EXCHANGE_DEPOSITION, &
       'catalog lookup by canonical source name')
  term=source_term_descriptor_by_name('UNKNOWN_SOURCE_NAME')
  CALL require(term%id==0,'unknown source name rejected')
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
      CALL require(ALLOCATED(steps(1)%term_names),'step group materializes source names')
      CALL require(group_contains_name(steps(1),'SOURCE_TERM_MASS_EXCHANGE_DEPOSITION'), &
           'name-based mass membership')
      CALL require(.NOT.group_contains_name(steps(1),'UNKNOWN_SOURCE_NAME'), &
           'unknown name is not a group member')
      CALL registry%validate(error); CALL require(error=='','repeat validation')
      DO j=1,SIZE(groups)
        CALL require(ALLOCATED(groups(j)%term_names),'stage group has canonical source names')
        CALL require(SIZE(groups(j)%term_names)==SIZE(groups(j)%term_ids), &
             'stage group name and legacy counts agree')
      END DO
      DO j=1,SIZE(solves)
        CALL require(ALLOCATED(solves(j)%term_names),'solver has canonical source names')
        CALL require(SIZE(solves(j)%term_names)==SIZE(solves(j)%term_ids), &
             'solver name and legacy counts agree')
      END DO
      snapshot=registry%stage_groups(); snapshot(1)%id=999
      snapshot=registry%stage_groups(); CALL require(snapshot(1)%id/=999,'detached immutable storage')
    END DO
    END DO
    END DO
    END DO
  END DO
  CALL test_configured_default_schedule()
  CALL test_native_source_input_parser()
  ! One configured bridge smoke test covers all injection and mass-exchange gates.
  CALL configure_bridge()
  CALL test_native_source_configuration()
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
  SUBROUTINE test_configured_default_schedule()
    USE parameters_2d, ONLY : wp,rheology_flag,rheology_model_id,pore_pressure_flag, &
         idx_pore,imex_velocity_group,imex_energy_group,imex_pore_pressure_group, &
         imex_group_method,curvature_term_flag,bottom_radial_source_flag, &
         bottom_fissural_source_flag,radial_source_flag,lateral_source_flag, &
         liquid_flag,gas_flag,gas_loss_flag,liquid_vaporization_flag
    USE constitutive_parameters_2d, ONLY : erosion_coeff,settling_flag, &
         coeff_porosity,entrainment_flag
    INTEGER :: default_mask(7)
    LOGICAL :: active_mask(7)
    INTEGER :: default_methods(3)
    CHARACTER(LEN=16) :: default_method_names(3)
    INTEGER :: expected_terms(2)

    idx_pore=7
    CALL require(imex_velocity_group==3 .AND. imex_energy_group==1 .AND. &
         imex_pore_pressure_group==2,'compiled legacy IMEX group defaults')
    default_method_names=[CHARACTER(LEN=16)::'ANALYTIC','ANALYTIC','NEWTON']
    CALL require(ALL(imex_group_method==default_method_names), &
         'compiled legacy IMEX method defaults')
    default_methods=[IM_GROUP_METHOD_ANALYTIC,IM_GROUP_METHOD_ANALYTIC, &
         IM_GROUP_METHOD_NEWTON]

    rheology_flag=.FALSE.;pore_pressure_flag=.FALSE.
    curvature_term_flag=.FALSE.;bottom_radial_source_flag=.FALSE.
    bottom_fissural_source_flag=.FALSE.;radial_source_flag=.FALSE.
    lateral_source_flag=.FALSE.;liquid_flag=.FALSE.;gas_flag=.TRUE.
    gas_loss_flag=.FALSE.;liquid_vaporization_flag=.FALSE.
    erosion_coeff=0.0_wp;settling_flag=.FALSE.;coeff_porosity=0.0_wp
    entrainment_flag=.FALSE.
    active_mask=.FALSE.;active_mask(2:3)=.TRUE.
    default_mask=0;default_mask(2:3)=1
    CALL partition%initialize(active_mask,default_mask,[IM_GROUP_METHOD_NEWTON])
    CALL initialize_configured_sources(registry,partition,error)
    CALL require(error=='','inactive compiled default schedule: '//TRIM(error))
    groups=registry%stage_groups();steps=registry%step_groups();solves=registry%solvers()
    CALL require(SIZE(groups)==2 .AND. SIZE(steps)==1 .AND. SIZE(solves)==1, &
         'inactive compiled default group counts')
    CALL require(groups(1)%id==3 .AND. groups(1)%legacy_id==3 .AND. &
         groups(1)%method==SOURCE_METHOD_NEWTON .AND. &
         ALL(groups(1)%equation_ids==[2,3]),'default velocity group mapping')
    CALL require(SIZE(groups(1)%term_ids)==0 .AND. &
         SIZE(steps(1)%term_ids)==0,'inactive defaults register no source terms')
    CALL require(groups(2)%id==4 .AND. &
         groups(2)%method==SOURCE_METHOD_EXPLICIT,'default EX group reservation')

    rheology_flag=.TRUE.;rheology_model_id=rheology_id_lava
    pore_pressure_flag=.TRUE.
    active_mask(4)=.TRUE.;active_mask(idx_pore)=.TRUE.
    default_mask(2:3)=3;default_mask(4)=1;default_mask(idx_pore)=2
    CALL partition%initialize(active_mask,default_mask,default_methods)
    CALL initialize_configured_sources(registry,partition,error)
    CALL require(error=='','active compiled default schedule: '//TRIM(error))
    groups=registry%stage_groups();steps=registry%step_groups();solves=registry%solvers()
    CALL require(SIZE(groups)==4 .AND. SIZE(steps)==1 .AND. SIZE(solves)==3, &
         'active compiled default group counts')
    CALL require(groups(1)%id==1 .AND. groups(1)%method==SOURCE_METHOD_ANALYTIC .AND. &
         ALL(groups(1)%equation_ids==[4]),'default analytic energy group')
    CALL require(groups(2)%id==2 .AND. groups(2)%method==SOURCE_METHOD_ANALYTIC .AND. &
         ALL(groups(2)%equation_ids==[idx_pore]),'default analytic pore group')
    CALL require(groups(3)%id==3 .AND. groups(3)%method==SOURCE_METHOD_NEWTON .AND. &
         ALL(groups(3)%equation_ids==[2,3]),'default Newton velocity group')
    expected_terms=[SOURCE_TERM_THERMAL_RADIATIVE_COOLING, &
         SOURCE_TERM_THERMAL_CONVECTIVE_COOLING]
    CALL require(same_members(groups(1)%term_ids,expected_terms), &
         'default energy cooling bundle')
    expected_terms=rheology_source_terms(rheology_id_lava)
    CALL require(same_members(groups(3)%term_ids,expected_terms), &
         'default lava rheology terms')
    CALL require(SIZE(groups(2)%term_ids)==1 .AND. &
         groups(2)%term_ids(1)==SOURCE_TERM_PORE_PRESSURE_DISSIPATION, &
         'default pore-pressure term')
    CALL require(SIZE(groups(4)%term_ids)==0 .AND. SIZE(steps(1)%term_ids)==0, &
         'default explicit and step groups remain empty when disabled')
  END SUBROUTINE test_configured_default_schedule
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
  SUBROUTINE test_native_source_configuration()
    INTEGER :: native_terms(7),term_index
    CHARACTER(LEN=64) :: term_names(7)

    n_eqns=7
    idx_pore=7
    native_source_configuration_present=.TRUE.
    CALL clear_source_input()
    native_source_configuration_present=.TRUE.
    source_group_count=5
    CALL set_group(1,10,'STAGE','ANALYTIC',1, &
         [SOURCE_TERM_THERMAL_RADIATIVE_COOLING,SOURCE_TERM_THERMAL_CONVECTIVE_COOLING], &
         ['TEMPERATURE'])
    CALL set_group(2,11,'STAGE','ANALYTIC',2, &
         [SOURCE_TERM_PORE_PRESSURE_DISSIPATION],['PORE_PRESSURE'])
    CALL set_group(3,12,'STAGE','NEWTON',3, &
         rheology_source_terms(rheology_id_lava),['VELOCITY'])
    CALL set_group(4,4,'STAGE','EXPLICIT',1, &
         [SOURCE_TERM_GEOMETRY_CURVATURE,SOURCE_TERM_BOTTOM_CIRCULAR_INJECTION, &
         SOURCE_TERM_BOTTOM_FISSURAL_INJECTION,SOURCE_TERM_RADIAL_INJECTION, &
         SOURCE_TERM_BOUNDARY_INJECTION],[CHARACTER(LEN=16)::])
    CALL set_group(5,1,'STEP','EXPLICIT',1, &
         [SOURCE_TERM_MASS_EXCHANGE_EROSION,SOURCE_TERM_MASS_EXCHANGE_DEPOSITION, &
         SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS, &
         SOURCE_TERM_MASS_EXCHANGE_CARRIER_DEPOSITION_LOSS, &
         SOURCE_TERM_MASS_EXCHANGE_CARRIER_DARCY_LOSS, &
         SOURCE_TERM_MASS_EXCHANGE_ENTRAINMENT,SOURCE_TERM_MASS_EXCHANGE_VAPORIZATION], &
         [CHARACTER(LEN=16)::])
    idx_pore=0
    CALL initialize_native_source_configuration(trial,partition,error)
    CALL require(INDEX(error,'pore-pressure equation index')>0 .AND. &
         .NOT.trial%is_initialized() .AND. .NOT.partition%is_initialized(), &
         'out-of-range pore equation index rejected safely')
    idx_pore=7
    n_eqns=0
    CALL initialize_native_source_configuration(trial,partition,error)
    CALL require(INDEX(error,'valid equation layout')>0 .AND. &
         .NOT.trial%is_initialized() .AND. .NOT.partition%is_initialized(), &
         'uninitialized equation layout rejected safely')
    n_eqns=7
    CALL initialize_native_source_configuration(registry,partition,error)
    CALL require(error=='' .AND. registry%is_initialized(), &
         'native named source groups accepted: '//TRIM(error))
    groups=registry%stage_groups();steps=registry%step_groups();solves=registry%solvers()
    CALL require(SIZE(groups)==4 .AND. SIZE(steps)==1 .AND. SIZE(solves)==3, &
         'native source group and solver counts')
    CALL require(groups(1)%method==SOURCE_METHOD_ANALYTIC .AND. &
         groups(2)%method==SOURCE_METHOD_ANALYTIC .AND. &
         groups(3)%method==SOURCE_METHOD_NEWTON .AND. &
         groups(4)%method==SOURCE_METHOD_EXPLICIT, &
         'native phase and execution order')
    CALL require(ALL(groups(3)%equation_ids==[2,3]) .AND. &
         SIZE(steps(1)%term_ids)==7,'native targets and complete step batch')

    source_group_term_count(5)=6
    CALL initialize_native_source_configuration(trial,partition,error)
    CALL require(error/='' .AND. .NOT.trial%is_initialized() .AND. &
         .NOT.partition%is_initialized(),'native schedule omission rejected transactionally')
    source_group_term_count(5)=7
    native_terms=[SOURCE_TERM_MASS_EXCHANGE_EROSION,SOURCE_TERM_MASS_EXCHANGE_DEPOSITION, &
         SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS, &
         SOURCE_TERM_MASS_EXCHANGE_CARRIER_DEPOSITION_LOSS, &
         SOURCE_TERM_MASS_EXCHANGE_CARRIER_DARCY_LOSS, &
         SOURCE_TERM_MASS_EXCHANGE_ENTRAINMENT,SOURCE_TERM_MASS_EXCHANGE_VAPORIZATION]
    DO term_index=1,SIZE(native_terms)
      term=source_term_descriptor(native_terms(term_index))
      term_names(term_index)=term%name
    END DO
    source_group_terms(1,5)=source_group_terms(2,5)
    CALL initialize_native_source_configuration(trial,partition,error)
    CALL require(error/='' .AND. .NOT.trial%is_initialized(), &
         'duplicate native source term rejected')
    source_group_terms(1,5)=term_names(1)
    CALL clear_source_input()
  END SUBROUTINE test_native_source_configuration

  SUBROUTINE test_native_source_input_parser()
    CHARACTER(LEN=128) :: valid_lines(13)

    valid_lines=[CHARACTER(LEN=128) :: &
         '&source_groups', &
         ' source_group_count=1,', &
         ' source_group_id(1)=5,', &
         ' source_group_order(1)=1,', &
         " source_group_cadence(1)='STAGE',", &
         " source_group_method(1)='NEWTON',", &
         ' source_group_term_count(1)=1,', &
         " source_group_terms(1,1)='SOURCE_TERM_RHEOLOGY_VOELLMY_COULOMB',", &
         ' source_group_solver_count(1)=1,', &
         ' source_solver_order(1,1)=1,', &
         ' source_solver_target_count(1,1)=1,', &
         " source_solver_targets(1,1,1)='VELOCITY',", &
         '/']
    CALL check_native_input(valid_lines,.FALSE.,'',.TRUE.)
    CALL require(source_group_count==1 .AND. source_group_id(1)==5 .AND. &
         source_group_terms(1,1)=='SOURCE_TERM_RHEOLOGY_VOELLMY_COULOMB' .AND. &
         source_solver_targets(1,1,1)=='VELOCITY','native namelist values parsed')
    CALL check_native_input([CHARACTER(LEN=128) :: '! no native source groups'], &
         .FALSE.,'',.FALSE.)
    CALL require(source_group_count==0,'absent native namelist preserves legacy path')
    CALL check_native_input([CHARACTER(LEN=128) :: &
         '&source_groups', ' source_group_count=1', '/', &
         '&SOURCE_GROUPS', ' source_group_count=1', '/'], &
         .FALSE.,'may appear only once',.FALSE.)
    CALL check_native_input([CHARACTER(LEN=128) :: &
         '&source_groups', ' unknown_source_field=1', '/'], &
         .FALSE.,'Invalid SOURCE_GROUPS',.FALSE.)
    CALL check_native_input(valid_lines,.TRUE.,'cannot be combined',.TRUE.)
    CALL clear_source_input()
  END SUBROUTINE test_native_source_input_parser

  SUBROUTINE check_native_input(lines,legacy_imex_explicit,expected_error,expected_present)
    CHARACTER(*), INTENT(IN) :: lines(:),expected_error
    LOGICAL, INTENT(IN) :: legacy_imex_explicit,expected_present
    INTEGER :: input_unit,line_index
    OPEN(NEWUNIT=input_unit,STATUS='scratch',ACTION='readwrite',FORM='formatted')
    DO line_index=1,SIZE(lines)
      WRITE(input_unit,'(A)') TRIM(lines(line_index))
    END DO
    CALL clear_source_input()
    CALL read_native_source_groups(input_unit,legacy_imex_explicit,error)
    IF (LEN_TRIM(expected_error)==0) THEN
      CALL require(error=='','native input reader: '//TRIM(error))
    ELSE
      CALL require(INDEX(error,TRIM(expected_error))>0, &
           'native input error diagnostic: '//TRIM(error))
    END IF
    CALL require(native_source_configuration_present .EQV. expected_present, &
         'native input presence flag')
    CLOSE(input_unit)
  END SUBROUTINE check_native_input

  SUBROUTINE set_group(slot,id,cadence,method,order,term_ids,targets)
    INTEGER, INTENT(IN) :: slot,id,order,term_ids(:)
    CHARACTER(*), INTENT(IN) :: cadence,method,targets(:)
    TYPE(source_term_type) :: descriptor
    INTEGER :: i
    source_group_id(slot)=id
    source_group_cadence(slot)=cadence
    source_group_method(slot)=method
    source_group_order(slot)=order
    source_group_term_count(slot)=SIZE(term_ids)
    DO i=1,SIZE(term_ids)
      descriptor=source_term_descriptor(term_ids(i))
      source_group_terms(i,slot)=descriptor%name
    END DO
    source_group_solver_count(slot)=SIZE(targets)
    IF (SIZE(targets)>0) THEN
      source_solver_order(1,slot)=1
      source_solver_target_count(1,slot)=SIZE(targets)
      source_solver_targets(1:SIZE(targets),1,slot)=targets
    END IF
  END SUBROUTINE set_group
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
