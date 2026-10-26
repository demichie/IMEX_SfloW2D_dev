PROGRAM test_source_registry
  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit
  USE source_catalog_2d
  USE source_registry_2d
  USE equation_metadata_2d, ONLY : equation_partition_type, &
       IM_GROUP_METHOD_NEWTON, IM_GROUP_METHOD_ANALYTIC
  USE source_configuration_2d, ONLY : initialize_native_source_configuration, initialize_configured_sources
  USE source_input_2d
  USE parameters_2d, ONLY : wp, n_eqns, idx_pore, &
       rheology_flag, rheology_model, pore_pressure_flag, thermal_cooling_flag, &
       curvature_term_flag, bottom_radial_source_flag, bottom_fissural_source_flag, &
       radial_source_flag, lateral_source_flag, liquid_flag, gas_flag, &
       gas_loss_flag, liquid_vaporization_flag, imex_velocity_group, &
       imex_energy_group, imex_pore_pressure_group, imex_group_method
  USE constitutive_parameters_2d, ONLY : erosion_coeff, settling_flag, &
       coeff_porosity, entrainment_flag
  IMPLICIT NONE
  TYPE(equation_partition_type) :: partition
  TYPE(source_registry_type) :: registry,invalid
  TYPE(source_group_type), ALLOCATABLE :: stages(:),steps(:),wrong(:)
  TYPE(source_term_type) :: term
  CHARACTER(LEN=64), ALLOCATABLE :: active(:)
  CHARACTER(LEN=256) :: error
  INTEGER :: i,mode,a,b,c,bits,ids(3),methods(3),membership(7),nm
  LOGICAL :: mask(7),enabled(3)
  CALL require(SOURCE_TERM_COUNT==24,'physical catalog contains 24 indivisible sources including extension fixture')
  CALL require(source_catalog_names_unique(),'unique physical source names')
  CALL require(canonical_rheology_catalog_unique(),'unique complete rheology names')
  term=source_term_descriptor_by_name('rheology_newmodel')
  CALL require(term%has_newton .AND. .NOT.term%has_explicit .AND. &
       .NOT.term%has_analytic,'Newton-only extension capabilities')
  DO i=1,SOURCE_TERM_COUNT
    term=source_term_descriptor_by_name(SOURCE_CATALOG(i)%name)
    CALL require(term%name==SOURCE_CATALOG(i)%name,'physical catalog lookup')
    CALL require(LEN_TRIM(term%target)>0,'physical target declared')
  END DO
  term=source_term_descriptor_by_name('SOURCE_TERM_RHEOLOGY_LAVA_VISCOUS')
  CALL require(term%name=='','legacy mathematical addend removed')
  term=source_term_descriptor_by_name('mass_exchange_carrier_deposition_loss')
  CALL require(term%name=='','dependent deposition carrier loss is not independently schedulable')
  term=source_term_descriptor_by_name('SOURCE_TERM_THERMAL_RADIATIVE_COOLING')
  CALL require(term%name=='','thermal formula component is not schedulable')
  DO mode=0,3
    enabled=[.TRUE.,BTEST(mode,0),BTEST(mode,1)]
    mask=.FALSE.;mask(2:3)=.TRUE.;mask(4)=enabled(2);mask(7)=enabled(3)
    active=[CHARACTER(LEN=64)::'rheology_lava','geometry_curvature', &
         'bottom_circular_injection','boundary_injection','mass_exchange_deposition']
    IF (enabled(2)) active=[active,'thermal_cooling']
    IF (enabled(3)) active=[active,'pore_pressure_dissipation']
    DO a=1,3
    DO b=1,3
    DO c=1,3
    DO bits=0,7
      ids=[a,b,c]
      DO i=1,3
        methods(i)=MERGE(IM_GROUP_METHOD_ANALYTIC,IM_GROUP_METHOD_NEWTON,BTEST(bits,i-1))
      END DO
      membership=0;nm=0
      BLOCK
        INTEGER :: selected(3)
        DO i=1,3
          IF (.NOT.ANY(enabled .AND. ids==i)) CYCLE
          nm=nm+1;selected(nm)=methods(i)
          IF (enabled(1) .AND. ids(1)==i) membership(2:3)=nm
          IF (enabled(2) .AND. ids(2)==i) membership(4)=nm
          IF (enabled(3) .AND. ids(3)==i) membership(7)=nm
        END DO
        CALL partition%initialize(mask,membership,selected(1:nm))
      END BLOCK
      CALL registry%initialize_legacy(partition,ids,methods,MERGE(7,0,enabled(3)),active,error)
      CALL require(error=='' .AND. registry%is_initialized(), &
           'legacy combination normalizes: '//TRIM(error))
      stages=registry%stage_groups();steps=registry%step_groups()
      CALL require(SIZE(steps)==1,'one post-RK source group')
      CALL require(group_contains_name(steps(1),'mass_exchange_deposition'),'mass source ownership')
      CALL require(.NOT.group_contains_name(steps(1),'rheology_lava'),'no stage source in step')
      CALL registry%validate(error)
      CALL require(error=='','repeat physical registry validation')
      wrong=stages
      DO i=1,SIZE(wrong)
        IF (wrong(i)%method==SOURCE_METHOD_EXPLICIT) CYCLE
        IF (.NOT.ANY(wrong(i)%equation_ids==2)) CYCLE
        wrong(i)%source_names=[CHARACTER(LEN=64)::'SOURCE_TERM_RHEOLOGY_LAVA_VISCOUS']
        CALL invalid%initialize(7,MERGE(7,0,enabled(3)),wrong,steps,error)
        CALL require(error/='' .AND. .NOT.invalid%is_initialized(), &
             'mathematical addend rejected in velocity group')
        EXIT
      END DO
    END DO
    END DO
    END DO
    END DO
  END DO
  CALL verify_default_schedule()
  CALL verify_thermal_activation_matrix()
  CALL verify_native_parser()
  CALL verify_native_schedule()
  WRITE(*,'(A)') 'D17 physical catalog, native configuration and 864 legacy mappings verified'
CONTAINS
  SUBROUTINE configure_lava_case()
    n_eqns=7;idx_pore=7
    rheology_flag=.TRUE.;rheology_model='lava'
    thermal_cooling_flag=.TRUE.
    pore_pressure_flag=.TRUE.
    curvature_term_flag=.FALSE.
    bottom_radial_source_flag=.FALSE.
    bottom_fissural_source_flag=.FALSE.
    radial_source_flag=.FALSE.
    lateral_source_flag=.TRUE.
    liquid_flag=.FALSE.;gas_flag=.FALSE.
    gas_loss_flag=.FALSE.;liquid_vaporization_flag=.FALSE.
    erosion_coeff=0.0_wp;settling_flag=.TRUE.
    coeff_porosity=0.0_wp;entrainment_flag=.FALSE.
    imex_velocity_group=3;imex_energy_group=1;imex_pore_pressure_group=2
    imex_group_method=[CHARACTER(LEN=16)::'ANALYTIC','ANALYTIC','NEWTON']
  END SUBROUTINE configure_lava_case

  SUBROUTINE verify_default_schedule()
    TYPE(source_group_type), ALLOCATABLE :: g(:),s(:)
    LOGICAL :: enabled_eq(7)
    INTEGER :: assigned(7)
    CALL configure_lava_case()
    enabled_eq=.FALSE.;enabled_eq(2:4)=.TRUE.;enabled_eq(7)=.TRUE.
    assigned=0;assigned(2:3)=3;assigned(4)=1;assigned(7)=2
    CALL partition%initialize(enabled_eq,assigned, &
         [IM_GROUP_METHOD_ANALYTIC,IM_GROUP_METHOD_ANALYTIC,IM_GROUP_METHOD_NEWTON])
    CALL initialize_configured_sources(registry,partition,error)
    CALL require(error=='','configured default physical schedule: '//TRIM(error))
    g=registry%stage_groups();s=registry%step_groups()
    CALL require(SIZE(g)==4 .AND. SIZE(s)==1,'default schedule group count')
    CALL require(SIZE(g(1)%source_names)==1 .AND. &
         group_contains_name(g(1),'thermal_cooling'),'one complete cooling source')
    CALL require(SIZE(g(2)%source_names)==1 .AND. &
         group_contains_name(g(2),'pore_pressure_dissipation'),'pore dissipation scheduling')
    CALL require(SIZE(g(3)%source_names)==1 .AND. &
         group_contains_name(g(3),'rheology_lava'),'complete lava law scheduling')
    CALL require(group_contains_name(g(4),'boundary_injection'), &
         'explicit boundary scheduling')
    CALL require(group_contains_name(s(1),'mass_exchange_deposition'), &
         'post-RK deposition scheduling')
  END SUBROUTINE verify_default_schedule

  SUBROUTINE verify_thermal_activation_matrix()
    TYPE(source_group_type), ALLOCATABLE :: groups(:)
    LOGICAL :: enabled_eq(7),has_cooling
    INTEGER :: assigned(7),model,choice,method(2),i
    DO model=1,2
      DO choice=1,2
        CALL configure_lava_case()
        rheology_model=MERGE('lava   ','voellmy',model==2)
        thermal_cooling_flag=(choice==2)
        pore_pressure_flag=.FALSE.
        lateral_source_flag=.FALSE.
        settling_flag=.FALSE.
        enabled_eq=.FALSE.; enabled_eq(2:3)=.TRUE.
        enabled_eq(4)=thermal_cooling_flag
        assigned=0
        IF (thermal_cooling_flag) THEN
          assigned(4)=1
          assigned(2:3)=2
          method=[IM_GROUP_METHOD_ANALYTIC,IM_GROUP_METHOD_NEWTON]
          CALL partition%initialize(enabled_eq,assigned,method)
        ELSE
          assigned(2:3)=1
          CALL partition%initialize(enabled_eq,assigned,[IM_GROUP_METHOD_NEWTON])
        END IF
        CALL initialize_configured_sources(registry,partition,error)
        CALL require(error=='','cooling on/off independent of rheology: '//TRIM(error))
        groups=registry%stage_groups()
        has_cooling=.FALSE.
        DO i=1,SIZE(groups)
          has_cooling=has_cooling .OR. group_contains_name(groups(i),'thermal_cooling')
        END DO
        CALL require(has_cooling .EQV. thermal_cooling_flag, &
             'physical cooling membership follows dedicated flag')
        CALL require(COUNT(enabled_eq)==MERGE(3,2,thermal_cooling_flag), &
             'legacy equation layout changes only when cooling active')
      END DO
    END DO
  END SUBROUTINE verify_thermal_activation_matrix

  SUBROUTINE verify_native_parser()
    INTEGER :: unit,i
    CHARACTER(LEN=160), PARAMETER :: lines(10)=[CHARACTER(LEN=160):: &
         '&source_groups', &
         ' source_group_count=1,', &
         ' source_group_id(1)=5,', &
         " source_group_cadence(1)='STAGE',", &
         " source_group_method(1)='NEWTON',", &
         ' source_group_term_count(1)=1,', &
         " source_group_terms(1,1)='RHEOLOGY_LAVA',", &
         ' source_group_target_count(1)=1,', &
         " source_group_targets(1,1)='VELOCITY',", &
         '/']
    OPEN(NEWUNIT=unit,STATUS='scratch',ACTION='readwrite',FORM='formatted')
    DO i=1,SIZE(lines)
      WRITE(unit,'(A)') TRIM(lines(i))
    END DO
    CALL read_native_source_groups(unit,.FALSE.,error)
    CALL require(error=='' .AND. native_source_configuration_present, &
         'native SOURCE_GROUPS namelist parsed')
    CALL require(source_group_terms(1,1)=='RHEOLOGY_LAVA', &
         'native parser preserves user source spelling')
    REWIND(unit)
    CALL read_native_source_groups(unit,.TRUE.,error)
    CALL require(INDEX(error,'cannot be combined')>0, &
         'native and explicit legacy IMEX input conflict')
    CLOSE(unit)
    CALL clear_source_input()
    OPEN(NEWUNIT=unit,STATUS='scratch',ACTION='readwrite',FORM='formatted')
    WRITE(unit,'(A)') '&source_groups'
    WRITE(unit,'(A)') ' source_group_count=1'
    WRITE(unit,'(A)') '/'
    WRITE(unit,'(A)') '&source_groups'
    WRITE(unit,'(A)') ' source_group_count=1'
    WRITE(unit,'(A)') '/'
    CALL read_native_source_groups(unit,.FALSE.,error)
    CALL require(INDEX(error,'may appear only once')>0, &
         'duplicate native namelist rejected')
    CLOSE(unit)
    CALL clear_source_input()
  END SUBROUTINE verify_native_parser

  SUBROUTINE set_native_group(slot,id,cadence,method,name,target)
    INTEGER, INTENT(IN) :: slot,id
    CHARACTER(*), INTENT(IN) :: cadence,method,name,target
    source_group_id(slot)=id
    source_group_cadence(slot)=cadence
    source_group_method(slot)=method
    source_group_term_count(slot)=1
    source_group_terms(1,slot)=name
    source_group_target_count(slot)=0
    IF (LEN_TRIM(target)>0) THEN
      source_group_target_count(slot)=1
      source_group_targets(1,slot)=target
    END IF
  END SUBROUTINE set_native_group

  SUBROUTINE compare_groups(left,right,label)
    TYPE(source_group_type), INTENT(IN) :: left,right
    CHARACTER(*), INTENT(IN) :: label
    INTEGER :: i
    CALL require(left%method==right%method .AND. &
         left%cadence==right%cadence, &
         label//' method/cadence')
    CALL require(ALLOCATED(left%equation_ids) .EQV. &
         ALLOCATED(right%equation_ids),label//' equations allocated')
    IF (ALLOCATED(left%equation_ids)) THEN
      CALL require(SIZE(left%equation_ids)==SIZE(right%equation_ids), &
           label//' equation count')
      CALL require(ALL(left%equation_ids==right%equation_ids), &
           label//' equation indices')
    END IF
    CALL require(ALLOCATED(left%source_names) .AND. &
         ALLOCATED(right%source_names),label//' source names allocated')
    CALL require(SIZE(left%source_names)==SIZE(right%source_names), &
         label//' source count')
    DO i=1,SIZE(left%source_names)
      CALL require(ANY(right%source_names==left%source_names(i)), &
           label//' canonical source membership')
    END DO
  END SUBROUTINE compare_groups

  SUBROUTINE verify_native_schedule()
    TYPE(source_group_type), ALLOCATABLE :: g(:),s(:),legacy_g(:),legacy_s(:)
    INTEGER :: i
    CALL configure_lava_case()
    ! Obtain an independently normalized legacy schedule before building
    ! the native equivalent. Compare full equation/source footprints below.
    BLOCK
      LOGICAL :: active_eq(7)
      INTEGER :: assigned(7)
      active_eq=.FALSE.;active_eq(2:4)=.TRUE.;active_eq(7)=.TRUE.
      assigned=0;assigned(2:3)=3;assigned(4)=1;assigned(7)=2
      CALL partition%initialize(active_eq,assigned, &
           [IM_GROUP_METHOD_ANALYTIC,IM_GROUP_METHOD_ANALYTIC,IM_GROUP_METHOD_NEWTON])
    END BLOCK
    CALL initialize_configured_sources(registry,partition,error)
    CALL require(error=='','legacy schedule for native equivalence: '//TRIM(error))
    legacy_g=registry%stage_groups()
    legacy_s=registry%step_groups()
    CALL clear_source_input()
    native_source_configuration_present=.TRUE.
    source_group_count=5
    CALL set_native_group(1,10,'STAGE','ANALYTIC','THERMAL_COOLING','TEMPERATURE')
    CALL set_native_group(2,11,'STAGE','ANALYTIC','PORE_PRESSURE_DISSIPATION','PORE_PRESSURE')
    CALL set_native_group(3,12,'STAGE','NEWTON','RHEOLOGY_LAVA','VELOCITY')
    CALL set_native_group(4,13,'STAGE','EXPLICIT','BOUNDARY_INJECTION','')
    CALL set_native_group(5,14,'STEP','EXPLICIT','MASS_EXCHANGE_DEPOSITION','')
    CALL initialize_native_source_configuration(registry,partition,error)
    CALL require(error=='' .AND. registry%is_initialized(), &
         'native names case-insensitive and complete: '//TRIM(error))
    g=registry%stage_groups();s=registry%step_groups()
    CALL require(SIZE(g)==4 .AND. SIZE(s)==1,'native schedule retains five groups')
    CALL require(group_contains_name(g(1),'thermal_cooling') .AND. &
         group_contains_name(g(3),'rheology_lava'),'native canonical physical names')
    CALL require(group_contains_name(s(1),'mass_exchange_deposition'), &
         'native post-RK physical group')
    CALL require(SIZE(g)==SIZE(legacy_g) .AND. SIZE(s)==SIZE(legacy_s), &
         'legacy/native same group counts')
    DO i=1,SIZE(g)
      CALL compare_groups(legacy_g(i),g(i),'legacy/native STAGE')
    END DO
    DO i=1,SIZE(s)
      CALL compare_groups(legacy_s(i),s(i),'legacy/native STEP')
    END DO
    source_group_terms(1,3)='SOURCE_TERM_RHEOLOGY_LAVA_VISCOUS'
    CALL initialize_native_source_configuration(invalid,partition,error)
    CALL require(error/='' .AND. .NOT.invalid%is_initialized(), &
         'native mathematical addend rejected')
    source_group_terms(1,3)='RHEOLOGY_LAVA'
    source_group_targets(1,3)='UNKNOWN'
    CALL initialize_native_source_configuration(invalid,partition,error)
    CALL require(INDEX(error,'unknown or duplicate group target')>0, &
         'native unsupported target rejected')
    source_group_targets(1,3)='VELOCITY'
    source_group_terms(1,5)='BOUNDARY_INJECTION'
    CALL initialize_native_source_configuration(invalid,partition,error)
    CALL require(error/='' .AND. .NOT.invalid%is_initialized(), &
         'native missing and duplicated physical source rejected')
    CALL clear_source_input()
  END SUBROUTINE verify_native_schedule

  SUBROUTINE require(valid,message)
    LOGICAL, INTENT(IN) :: valid
    CHARACTER(*), INTENT(IN) :: message
    IF (valid) RETURN
    WRITE(error_unit,'(A)') 'FAIL: '//message
    ERROR STOP 1
  END SUBROUTINE require
END PROGRAM test_source_registry
